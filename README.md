# Yanosim (yet another nanopore simulator)

Read simulator for nanopore DRS datasets.

## Installation:

`pip install git+git://github.com/bartongroup/yanosim`

## Usage:

### `model`:

```
$ yanosim model --help
Usage: yanosim model [OPTIONS]

  Creates an model of mismatches, insertions and deletions based on an
  alignment of nanopore DRS reads to a reference. Reads should be aligned to
  a transcriptome i.e. without spliced alignment, using minimap2. They
  should have the cs tag.

Options:
  -b, --bam-fn TEXT        [required]
  -o, --output-fn TEXT     [required]
  -p, --processes INTEGER
  --help                   Show this message and exit.
```

### `quantify`:

```
$ yanosim quantify --help
Usage: yanosim quantify [OPTIONS]

  Quantify the number of reads mapping to each transcript in a reference, so
  that the right number of reads can be simulated.

Options:
  -b, --bam-fn TEXT               [required]
  -o, --output-fn TEXT            [required]
  -g, --gtf-fn TEXT
  -f, --filtered-gtf-output-fn TEXT
  -r, --remove-ensembl-version
  --help                          Show this message and exit.
```

### `simulate`:

```
$ yanosim simulate --help
Usage: yanosim simulate [OPTIONS]

  Given a model created using yanosim model, and per-transcript read counts
  created using yanosim simulate, simulate error-prone long-reads from the
  given fasta file.

Options:
  -f, --fasta-fn TEXT             [required]
  -m, --model-fn TEXT             [required]
  -q, --quantification-fn TEXT    [required]
  -o, --output-fn TEXT            [required]
  --model-frag / --no-model-frag
  -p, --processes INTEGER
  -s, --seed INTEGER
  --help                          Show this message and exit.
```