# rna_secstruct_design

[![CI](https://github.com/jyesselm/rna_secstruct_design/actions/workflows/ci.yml/badge.svg)](https://github.com/jyesselm/rna_secstruct_design/actions/workflows/ci.yml)
[![PYPI status]( https://badge.fury.io/py/rna_secstruct_design.png)](http://badge.fury.io/py/rna_secstruct_design)[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Common secondary structure design algorithms for RNA: pick out the parts of a
construct you care about, then mutate, resize or redesign everything else.

## Install

```shell
python -m pip install git+https://github.com/jyesselm/rna_secstruct_design
```

Folding goes through the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) python
bindings, which come in as a dependency, so there is no separate install step.

## Quick start

Sequences and structures are handled by `SecStruct` from
[rna_secstruct](https://github.com/jyesselm/rna_secstruct). Positions are always
0-based.

```python
from rna_secstruct import SecStruct
from rna_secstruct_design.helix_randomizer import HelixRandomizer
from rna_secstruct_design.selection import get_selection

secstruct = SecStruct("AAGGGGGAAACCCC", "..((((....))))")

# hold the GAAA tetraloop fixed
exclude = get_selection(secstruct, {"motif": {"name": "gaaa_tetraloop"}})

# redesign every helix around it, keeping the same fold
ens_defect, sequence = HelixRandomizer().run(secstruct, exclude)
```

## Features

**Selections** (`selection.py`) turn a description of "which positions" into a
list of indices. Select by motif type, by a named motif (`gaaa_tetraloop`,
`tlr`, `tlr_extended`, `ref_hp`), by an exact sequence + structure, by a range,
or by the flanking pairs of every loop. Combine several, and `invert` to get
everything else. The same dict can be written as yaml and passed to the CLI.

**Mutations** (`mutations.py`) build variant libraries: exhaustive or sampled
single nucleotide scans, and basepair mutations that change both sides of a
helix so it stays paired.

**Structure edits** (`mutations.py`) change the length of the RNA: grow and
shrink helices, scan a helix across a range of lengths, insert bulges, delete
nucleotides.

**Helix design** (`helix_randomizer.py`) randomizes the helices of a construct,
leaving every other motif untouched, and keeps the variant that folds into the
target structure with the best ensemble defect.

**Replacement** (`replace.py`) finds a sequence + structure pattern, including
multi-strand ones, and swaps in a different one.

**Quality control** (`dataframe.py`, `constraints.py`) scores a library for GC
stretches and repeated nucleotides.

## Examples

Runnable, self-checking scripts live in [`examples/`](examples/) — see
[examples/README.md](examples/README.md) for what each one covers.

```shell
python examples/01_selections.py
```

## Command line

Installing the package provides `rna-struct-design`:

```shell
# every single mutant of a sequence, folded
rna-struct-design mut-scan -s GGGGAAAACCCC -ss "((((....))))" -o muts.csv

# 10 helix designs, protecting whatever select.yml matches
rna-struct-design helix-rand -s GGGGAAAACCCC -pf select.yml -n 10 -o designs.csv

# redesign a whole csv of constructs, on 4 cores
rna-struct-design helix-rand -csv constructs.csv -n 10 -p 4 -o designs.csv

# swap a motif across a csv
rna-struct-design replace constructs.csv swap.yml -o swapped.csv
```

`mut-scan` and `helix-rand` take a `-pf/--param-file` yaml selecting the
positions to **protect**:

```yaml
motif_1:
  name: gaaa_tetraloop
  extend_flank: 1
motif_2:
  name: tlr
  extend_flank: 2
```

`examples/07_cli.py` runs all of these end to end.

## Development

```shell
python -m pip install -e .
python -m pip install pytest black
python -m pytest test/
```
