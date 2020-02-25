import re
from bisect import bisect_left
import itertools as it
from collections import defaultdict, Counter
from functools import lru_cache, partial
from multiprocessing import Pool

import numpy as np
from scipy import stats

import pysam
import click

from .utils import write_model

RC = str.maketrans('ACGT', 'TGCA')


@lru_cache(maxsize=128)
def revcomp(seq):
    return seq.translate(RC)[::-1]


CS_SPLITTER = '([-+*~=:])'


def parse_cs_tag_to_alignment(cs_tag, strand):
    '''
    generalisable function for parsing minimap2 cs tag (long form only) into pw alignment
    '''
    cs_tag = re.split(CS_SPLITTER, cs_tag)[1:]
    cs_ops = cs_tag[::2]
    cs_info = cs_tag[1::2]
    if strand == '-':
        cs_ops = cs_ops[::-1]
        cs_info = cs_info[::-1]
    qur_seq = []
    ref_seq = []
    i = 0
    for op, info in zip(cs_ops, cs_info):
        if op == '=':
            # long frm match
            if strand == '-':
                info = revcomp(info)
            qur_seq.append(info)
            ref_seq.append(info)
            i += len(info)
        elif op == ':':
            # short form match
            raise ValueError('Need long form CS')
        elif op == '*':
            # mismatch
            info = info.upper()
            if strand == '-':
                info = info.translate(RC)
            ref = info[0]
            alt = info[1]
            qur_seq.append(alt)
            ref_seq.append(ref)
            i += 1
        elif op == '+':
            if strand == '-':
                info = revcomp(info)
            qur_seq.append(info.upper())
            ref_seq.append('-' * len(info))
        elif op == '-':
            if strand == '-':
                info = revcomp(info)
            qur_seq.append('-' * len(info))
            ref_seq.append(info.upper())
            i += len(info)
        elif op == '~':
            # ignore
            pass
    return ''.join(qur_seq), ''.join(ref_seq)


class ref_ignore_ins:

    def __init__(self):
        self.prev = None


    def __call__(self, aln_col):
        qur_base, ref_base = aln_col
        if self.prev is None:
            # the first base shouldn't be an insertion
            assert ref_base != '-'
            self.prev = ref_base

        if ref_base == '-':
            return self.prev
        else:
            self.prev = ref_base
            return ref_base


@lru_cache(maxsize=128)
def identify_state(ref_base, basecall):
    if len(basecall) > 1:
        return '+'
    elif basecall == '-':
        return '-'
    elif ref_base != basecall:
        return '*'
    else:
        return '='


class PairwiseAlignment:

    def __init__(self, cs_tag=None, strand=None, qur_seq=None, ref_seq=None):
        if cs_tag is not None:
            if strand is None:
                strand = '+'
            self.qur_seq, self.ref_seq = parse_cs_tag_to_alignment(cs_tag, strand)
        else:
            if qur_seq is None or ref_seq is None:
                raise ValueError()
            self.validate_seqs(qur_seq, ref_seq)
            # RNA is sequenced 3' -> 5' so reverse it to build the model
            self.qur_seq = qur_seq[::-1]
            self.ref_seq = ref_seq[::-1]
            

    def validate_seqs(self, qur_seq, ref_seq):
        assert len(qur_seq) == len(ref_seq), 'Seq lengths not equal'
        assert len(qur_seq) == len(re.findall('[ACGTN-]', qur_seq)), 'Query contains non ACGT- chars'
        assert len(ref_seq) == len(re.findall('[ACGTN-]', ref_seq)), 'Ref contains non ACGT- chars'

    def ref_rle(self):
        for base, hp_aln in it.groupby(zip(self.qur_seq, self.ref_seq), key=ref_ignore_ins()):
            qur_hp, ref_hp = zip(*hp_aln)
            qur_hp = ''.join(qur_hp)
            ref_hp = ''.join(ref_hp)
            yield base, qur_hp, ref_hp

    def iter_ref_cols(self, pad_size=0, pad_val='A'):
        qur_padded = pad_val * pad_size + self.qur_seq
        ref_padded = pad_val * pad_size + self.ref_seq
        q_with_ins = []
        for q, r in zip(qur_padded, ref_padded):
            q_with_ins.append(q)
            if r != '-':
                q = ''.join(q_with_ins)
                state = identify_state(r, q)
                yield q, r, state
                q_with_ins = []

    def iter_ref_kmers_and_bc(self, k=5):
        qur_cols, ref_cols, states = zip(*self.iter_ref_cols(pad_size=k - 1))
        for i in range(k, len(ref_cols) + 1):
            ref_kmer = ''.join(ref_cols[i - k: i])
            prev_states = ''.join(states[i - k: i - 1])
            next_state = states[i - 1]
            bc = qur_cols[i - 1]
            yield ref_kmer, prev_states, bc, next_state


def count_and_compress_homopolymers(p_aln, min_hp_length=3):
    qur_comp = []
    ref_comp = []
    homopolymer_counts = defaultdict(Counter)
    for n, qur_hp, ref_hp in p_aln.ref_rle():
        ref_hp_no_ins = ref_hp.replace('-', '')
        hp_ln = len(ref_hp_no_ins)
        if hp_ln >= min_hp_length:
            qur_hp_no_ins = qur_hp.replace('-', '')
            homopolymer_counts[ref_hp_no_ins][qur_hp_no_ins] += 1
            qur_comp.append(qur_hp_no_ins)
            ref_comp.append(qur_hp_no_ins)
        else:
            qur_comp.append(qur_hp)
            ref_comp.append(ref_hp)
    p_aln_comp = PairwiseAlignment(
        qur_seq=''.join(qur_comp),
        ref_seq=''.join(ref_comp)
    )
    return p_aln_comp, homopolymer_counts


def build_prob_tree(tag_generator, *, kmer_size, max_ins, min_hp_ln):
    basecall_counts = defaultdict(partial(defaultdict, Counter))
    homopolymer_counts = defaultdict(Counter)
    aln_lengths = defaultdict(list)
    for i, (ref_name, aln_len, cs_tag, strand) in enumerate(tag_generator):
        aln_lengths[ref_name].append(aln_len)

        p_aln = PairwiseAlignment(cs_tag, strand)
        p_aln, hp = count_and_compress_homopolymers(p_aln, min_hp_ln)
        nested_dd_of_counters_update(homopolymer_counts, hp)

        for ref_kmer, prev_states, bc, next_state in p_aln.iter_ref_kmers_and_bc(k=kmer_size):
            if len(bc) > (max_ins + 1):
                # if there is a huge insertion we make the assumption (perhaps wrongly)
                # that it is an alignment failure not a basecall failure
                bc = bc[-1]
                next_state = identify_state(ref_kmer[-1], bc)

            basecall_counts[ref_kmer][prev_states][(bc, next_state)] += 1

    return basecall_counts, homopolymer_counts, aln_lengths


def chunk_bam_tags(bam_fn, chunk_size=10_000):
    with pysam.AlignmentFile(bam_fn) as bam:
        tags = []
        i = 0
        mapped = bam.mapped
        with click.progressbar(bam.fetch(), length=mapped) as fetch_iter:
            for aln in fetch_iter:
                if aln.is_secondary:
                    continue
                try:
                    cs_tag = aln.get_tag('cs')
                except KeyError:
                    continue
                ref_name = aln.reference_name
                aln_len = aln.query_alignment_length
                strand = '+-'[aln.is_reverse]
                tags.append((ref_name, aln_len, cs_tag, strand))
                i += 1
                if i == chunk_size:
                    yield tags
                    i = 0
                    tags = []
            else:
                if tags:
                    yield tags


def nested_dd_of_counters_update(d1, d2):
    for k in d2:
        if isinstance(d2[k], Counter):
            assert isinstance(d1[k], Counter)
            d1[k] += d2[k]
        elif isinstance(d2[k], list):
            assert isinstance(d1[k], list)
            d1[k] += d2[k]
        elif isinstance(d2[k], defaultdict):
            assert isinstance(d1[k], defaultdict)
            nested_dd_of_counters_update(d1[k], d2[k])
        else:
            raise ValueError()


def cumsum_counts(dd):
    cum_counts = defaultdict(dict)
    for k in dd:
        if isinstance(dd[k], Counter):
            ns, cs = zip(*dd[k].most_common())
            cs = np.cumsum(cs) / sum(cs)
            cs = np.insert(cs, 0, 0)
            cum_counts[k] = ns, cs
        elif isinstance(dd[k], defaultdict):
            cum_counts[k] = cumsum_counts(dd[k])
        else:
            raise ValueError()
    return cum_counts


def estimate_p_fragmentation(aln_lens, perc=95, tol=500, min_reads=100):
    max_lengths = []
    frag_frac = []
    for ref, a in aln_lens.items():
        if len(a) >= min_reads:
            a = np.array(a)
            mlen = np.percentile(a, perc)
            n_frag = sum(a < (mlen - tol))
            max_lengths.append(mlen)
            frag_frac.append(n_frag / len(a))
    model = stats.linregress(max_lengths, frag_frac)._asdict()
    return model


def parallel_build_prob_tree(bam_fn, processes=12, chunk_size=1000,
                             kmer_size=5, max_ins=15, min_hp_len=5,
                             frag_full_length_perc=95,
                             frag_full_length_tol=500,
                             frag_min_reads=100):
    basecall_counts = defaultdict(lambda: defaultdict(Counter))
    homopolymer_counts = defaultdict(Counter)
    aln_lengths = defaultdict(list)
    with Pool(processes=processes) as pool:
        chunk_iter = chunk_bam_tags(bam_fn, chunk_size)
        build_prob_tree_ = partial(
            build_prob_tree,
            kmer_size=kmer_size,
            max_ins=max_ins,
            min_hp_ln=min_hp_len
        )
        for i, (b, h, a) in enumerate(pool.imap_unordered(build_prob_tree_, chunk_iter)):
            nested_dd_of_counters_update(basecall_counts, b)
            nested_dd_of_counters_update(homopolymer_counts, h)
            nested_dd_of_counters_update(aln_lengths, a)

    frag_model = estimate_p_fragmentation(
        aln_lengths,
        frag_full_length_perc,
        frag_min_reads
    )
    return cumsum_counts(basecall_counts), cumsum_counts(homopolymer_counts), frag_model