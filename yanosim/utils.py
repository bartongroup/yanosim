import numpy as np
import gzip
import json
import itertools as it
import re


def serialise(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj


def write_model(output_fn, basecall_ptree, homopolymer_ptree, frag_model):
    model = {
        'basecalls': basecall_ptree,
        'homopolymers': homopolymer_ptree,
        'fragmentation': frag_model,
    }
    with gzip.open(output_fn, 'wt', encoding="ascii") as f:
        json.dump(model, f, default=serialise)


def load_model(input_fn):
    with gzip.open(input_fn, 'rt', encoding='ascii') as f:
        model = json.load(f)
    return model['basecalls'], model['homopolymers'], model['fragmentation']


def load_quantification(quant_fn):
    quantification = {}
    with open(quant_fn) as q:
        for record in q:
            ref_name, n_reads = record.split()
            n_reads = int(n_reads)
            quantification[ref_name] = n_reads
    return quantification


def read_fasta(fasta_fn):
    with open(fasta_fn) as f:
        header_grouped = it.groupby(f, key=lambda line: line.startswith('>'))
        for _, read_id in header_grouped:
            read_id = list(read_id)
            assert len(read_id) == 1
            read_id = read_id[0][1:].strip()
            read_id = re.split('[\s|]', read_id)[0]
            _, seq = next(header_grouped)
            seq = ''.join([line.strip() for line in seq])
            yield read_id, seq