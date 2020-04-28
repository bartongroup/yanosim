import re
import pysam
import click
import numpy as np

from .model import parallel_build_prob_tree
from .simulate import parallel_simulate_read
from .utils import write_model, load_model

@click.group()
def yanosim():
    pass


@yanosim.command()
@click.option('-b', '--bam-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-p', '--processes', default=4)
def model(bam_fn, output_fn, processes):
    '''
    Creates an model of mismatches, insertions and deletions based on
    an alignment of nanopore DRS reads to a reference. Reads should
    be aligned to a transcriptome i.e. without spliced alignment, using
    minimap2. They should have the cs tag.
    '''
    basecall_model, homopolymer_model, frag_model = parallel_build_prob_tree(
        bam_fn, processes=processes
    )
    write_model(output_fn, basecall_model, homopolymer_model, frag_model)

    
@yanosim.command()
@click.option('-b', '--bam-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-g', '--gtf-fn', required=False)
@click.option('-f', '--filtered-gtf-output-fn', required=False)
@click.option('-r', '--remove-ensembl-version', default=False, is_flag=True)
def quantify(bam_fn, output_fn, gtf_fn, filtered_gtf_output_fn, remove_ensembl_version):
    '''
    Quantify the number of reads mapping to each transcript in a reference, so that
    the right number of reads can be simulated.
    '''
    if gtf_fn and not filtered_gtf_output_fn:
        raise click.BadParameter('If -g is specified, must also provide -f')
    expressed = set()
    with open(output_fn, 'w') as f, pysam.AlignmentFile(bam_fn) as bam:
        for ref in bam.references:
            i = 0
            for aln in bam.fetch(ref):
                if not aln.is_secondary:
                    i += 1
            f.write(f'{ref}\t{i}\n')
            if i:
                if remove_ensembl_version:
                    ref = ref.split('.')[0]
                expressed.add(ref)
    if gtf_fn:
        with open(gtf_fn) as g, open(filtered_gtf_output_fn, 'w') as f:
            for record in g:
                try:
                    transcript_id = re.search('transcript_id "(.*?)";', record).group(1)
                except AttributeError:
                    continue
                if transcript_id in expressed:
                    f.write(record)


@yanosim.command()
@click.option('-f', '--fasta-fn', required=True)
@click.option('-m', '--model-fn', required=True)
@click.option('-q', '--quantification-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('--model-frag/--no-model-frag', default=True)
@click.option('-p', '--processes', default=4)
@click.option('-s', '--seed', default=None, type=int)
def simulate(fasta_fn, model_fn, quantification_fn, output_fn, model_frag, processes, seed):
    '''
    Given a model created using yanosim model, and per-transcript read counts created using
    yanosim simulate, simulate error-prone long-reads from the given fasta file.
    '''
    basecall_model, homopolymer_model, fragmentation_model = load_model(model_fn)
    if seed is not None:
        np.random.seed(seed)
    with open(output_fn, 'w') as f:
        sim_reads = parallel_simulate_read(
            fasta_fn, quantification_fn, processes,
            basecall_ptree=basecall_model,
            hp_ptree=homopolymer_model,
            frag_model=fragmentation_model,
            model_frag=model_frag,
            polya_len=10,
            five_prime_loss=False,
            chunk_size=1000,
        )
        for read_id, seq in sim_reads:
            f.write(f'>{read_id}\n{seq}\n')
            f.flush()