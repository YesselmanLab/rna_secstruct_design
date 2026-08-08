# Examples

Each script is standalone and self-checking: it prints what it is doing and
asserts the result, so if one runs to `ok` the feature works. They run in CI for
exactly that reason.

Run one:

```shell
python examples/01_selections.py
```

Run all of them:

```shell
for f in examples/*.py; do python "$f" || exit 1; done
```

| Script | What it covers |
| --- | --- |
| `01_selections.py` | Picking positions by motif type, name, range, sequence/structure, and from yaml |
| `02_mutations.py` | Single nucleotide and basepair mutant libraries |
| `03_structure_edits.py` | Changing helix lengths, adding bulges, deleting nucleotides |
| `04_helix_design.py` | Redesigning helices with `HelixRandomizer`, and the sequence constraints |
| `05_replace_motifs.py` | Finding and swapping motifs, including across two strands |
| `06_full_workflow.py` | End to end: design a library, check it folds, QC it, write a csv |
| `07_cli.py` | The same operations from the `rna-struct-design` command line |

Scripts 04, 06 and 07 fold sequences with ViennaRNA, which comes in as a
dependency, so there is nothing extra to install. They are also the slow ones,
though all seven finish in a few seconds.
