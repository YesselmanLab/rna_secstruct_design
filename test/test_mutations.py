from rna_secstruct.secstruct import SecStruct
from rna_secstruct_design.mutations import (
    possible_nucleotide_mutations,
    find_mutations,
    find_multiple_mutations,
    change_helix_length,
)
from rna_secstruct_design.selection import get_selection

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
    results = find_multiple_mutations(seq, 1, [])
    assert len(results) == 9
    assert results[0].name == "A1U"
    assert results[0].sequence == "UUC"


def test_find_multiple_mutations_2():
    seq = "AUC"
    results = find_multiple_mutations(seq, 2, [])
    assert len(results) == 27
    assert results[0].name == "A1U_U2A"
    assert results[0].sequence == "UAC"


def test_find_mutations_mttr():
    seq = (
        "GUUGAUAUGGAUUUACUCCGAGGAGACGAACUACCACGAACAGGGGAAACUCUACCCGUGGCGUCUCCGUU"
        "UGACGAGUAAGUCCUAAGUCAACAAAGUCCGCGAGUAGCGGACAC"
    )
    ss = (
        "((((((..((((((((((((((((((((.....(((((...((((....))))...))))))))))))..)"
        "))..))))))))))...))))))...((((((.....)))))).."
    )
    secstruct = SecStruct(seq, ss)
    params = {"motif": {"name": "tlr", "extend_flank": 2}, "invert": True}
    exclude = get_selection(secstruct, params)
    results = find_multiple_mutations(seq, 1, exclude)
    assert len(results) == 51
    assert results[0].name == "G4A"


def test_change_helix_length():
    seqstruct = SecStruct("GGGGAAAACCCC", "((((....))))")
    new_secstruct = change_helix_length(seqstruct, 0, 5)
    assert new_secstruct.structure == "(((((....)))))"
    seqstruct = SecStruct("GGGAAAAAUCCC", "((((....))))")
    new_secstruct = change_helix_length(seqstruct, 0, 1)
    assert new_secstruct.structure == "(....)"
    assert new_secstruct.sequence == "AAAAAU"
    seqstruct = SecStruct("GGGG&CCCC", "((((&))))")
    new_secstruct = change_helix_length(seqstruct, 0, 1)
    assert new_secstruct.structure == "(&)"
    assert new_secstruct.sequence == "G&C"
    seqstruct = SecStruct("CGGAAAAAUCCG", "((((....))))")
    new_secstruct = change_helix_length(seqstruct, 0, 2)
    assert new_secstruct.structure == "((....))"
    assert new_secstruct.sequence == "CAAAAAUG"
    seqstruct = SecStruct("CGGAAAAAUCCG", "((((....))))")
    new_secstruct = change_helix_length(seqstruct, 0, 3)
    assert new_secstruct.structure == "(((....)))"
    assert new_secstruct.sequence == "CGAAAAAUCG"
