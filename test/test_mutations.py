from rna_secstruct_design.mutations import (
    possible_nucleotide_mutations,
    find_mutations,
    find_multiple_mutations,
)

import pytest


def test_possible_nucleotide_mutations():
    """
    test possible_nucleotide_mutations
    """
    assert possible_nucleotide_mutations("A") == ["U", "C", "G"]
    assert possible_nucleotide_mutations("U") == ["A", "C", "G"]
    assert possible_nucleotide_mutations("C") == ["A", "U", "G"]
    assert possible_nucleotide_mutations("G") == ["A", "U", "C"]

    with pytest.raises(ValueError):
        possible_nucleotide_mutations("X")


def test_find_mutations():
    seq = "AU"
    exclude = [0]
    seqs = find_mutations(seq, exclude)
    assert len(seqs) == 3

def test_find_multiple_mutations():
    seq = "AUC"
    seqs = find_multiple_mutations(seq, 1, [])
    print(seqs)