import itertools as it
from functools import partial
from multiprocessing import Pool

import numpy as np
import click

from .utils import read_fasta, load_quantification


def random_choice(bases, probs):
    p = np.random.random()
    idx = np.searchsorted(probs, p) - 1
    return bases[idx]


def mutate_basecalls(ref_seq, basecall_ptree, k=5):
    state = '=' * (k - 1)
    m = []
    for i in range(k, len(ref_seq) + 1):
        ref_kmer = ref_seq[i - k: i]
        try:
            possible_basecalls, probs = basecall_ptree[ref_kmer][state]
            bc, new_state = random_choice(*basecall_ptree[ref_kmer][state])
        except KeyError:
            # just use ref seq as basecall
            bc = ref_kmer[-1]
            new_state = '='
        m.append(bc)
        state = state[1:] + new_state
    return ''.join(m).replace('-', '')


def mutate_homopolymers(ref_seq, hp_ptree, min_hp_len=5):
    m = []
    for a, g in it.groupby(ref_seq):
        hp = ''.join(list(g))
        if len(hp) >= min_hp_len:
            while True:
                if hp in hp_ptree:
                    sim_hp = random_choice(*hp_ptree[hp])
                    break
                else:
                    # if homopolymer is so long that it hasn't been seen before,
                    # shorten by one and try again
                    hp = hp[:-1]
        else:
            sim_hp = hp
        m.append(sim_hp)
    return ''.join(m)


def random_fragment(ref_seq, frag_model):
    slen = len(ref_seq)
    p_frag = max(frag_model['slope'] * slen + frag_model['intercept'], 0)
    is_frag = np.random.random() <= p_frag
    if is_frag:
        frag_point = np.random.randint(0, slen)
        ref_seq = ref_seq[frag_point:]
    return ref_seq


def simulate_read(ref_seq, *, basecall_ptree, hp_ptree, frag_model,
                  model_frag=True, polya_len=10, five_prime_loss=False):
    if polya_len:
        ref_seq = ref_seq + 'A' * polya_len

    # sometimes the read is fragmented at a random point
    if model_frag:
        ref_seq = random_fragment(ref_seq, frag_model)
    # reverse sequence for simulation as RNA is sequenced 3' -> 5'
    ref_seq = ref_seq[::-1]
    # mutate the read to match the errors in the model
    sim_read = mutate_basecalls(ref_seq, basecall_ptree)
    # mutate the read to match the homopolymer errors in the model
    sim_read = mutate_homopolymers(sim_read, hp_ptree)
    # re-reverse sequence to 5' -> 3' direction
    sim_read = sim_read[::-1]
    # Direct RNA reads are known to generally be missing 11nt from the 5' end
    if five_prime_loss:
        sim_read = sim_read[11:]
    return sim_read


def get_seqs_to_sim(fasta_fn, quantification):
    for ref_name, seq in read_fasta(fasta_fn):
        for i in range(1, quantification[ref_name] + 1):
            read_id = f'{ref_name}_sim{i}'
            yield read_id, seq


def get_reads(fasta_fn, quantification_fn, chunk_size=10_000):
    quantification = load_quantification(quantification_fn)
    nsim = sum(list(quantification.values()))
    chunk = []
    seqs = get_seqs_to_sim(fasta_fn, quantification)
    with click.progressbar(seqs, length=nsim) as seqs:
        for read_id, seq in seqs:
            chunk.append((read_id, seq))
            if len(chunk) == chunk_size:
                yield chunk
                chunk = []
        else:
            if len(chunk):
                yield chunk


def _imap_simulate(chunk, **sim_kwargs):
    sim = []
    for read_id, seq in chunk:
        sim.append((read_id, simulate_read(seq, **sim_kwargs)))
    return sim


def parallel_simulate_read(fasta_fn, quantification_fn, processes,
                           chunk_size=10_000, **sim_kwargs):
    with Pool(processes=processes) as pool:
        for res in pool.imap_unordered(
                partial(_imap_simulate, **sim_kwargs),
                get_reads(fasta_fn, quantification_fn, chunk_size)):
            yield from res