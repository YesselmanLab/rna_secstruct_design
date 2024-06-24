from rna_secstruct.secstruct import SecStruct
from seq_tools.structure import SequenceStructure
from rna_secstruct_design.mutations import (
    possible_nucleotide_mutations,
    find_mutations,
    find_multiple_mutations,
    get_basepair_mutation,
    get_basepair_mutations,
    get_basepair_mutuations_random,
    change_helix_length,
    scan_helix_lengths,
    scan_all_helix_lengths,
    add_unpaired,
    add_unpaired_sweep,
    remove_nucleotides,
    remove_unpaired_nucleotide_sweep,
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


def test_scan_helix_lengths():
    seqstruct = SecStruct("GGGGAAAACCCC", "((((....))))")
    results = scan_helix_lengths(seqstruct, 0, 1, 10)
    assert len(results) == 10
    assert results[0].structure == "(....)"
    assert results[-1].structure == "((((((((((....))))))))))"

    seqstruct = SecStruct("ACCAUCGGAAACGAUGGU", "(((((((....)))))))")
    new_secstruct = change_helix_length(seqstruct, 0, 5)


def test_scan_all_helix_lengths():
    seqstruct = SecStruct("GGAGGAAAACCCC", "((.((....))))")
    h_ranges = {
        0: [2, 10],
        2: [2, 10],
    }
    results = scan_all_helix_lengths(seqstruct, h_ranges)
    assert len(results) == 81


def test_add_unpaired():
    seqstruct = SequenceStructure("GGGGAAAACCCC", "((((....))))")
    new_secstructs = add_unpaired(seqstruct, 1, 2, all_nucleotides=True)
    assert len(new_secstructs) == 16


def test_add_unpaired_sweep():
    seqstruct = SequenceStructure("GGGGAAAACCCC", "((((....))))")
    exclude = [2, 3, 4, 5, 6, 7, 8]
    new_secstructs = add_unpaired_sweep(seqstruct, 2, exclude, all_nucleotides=True)
    assert len(new_secstructs) == 32


def test_remove_nucleotides():
    seqstruct = SequenceStructure("GGGGAAAACCCC", "((((....))))")
    new_secstruct = remove_nucleotides(seqstruct, [4, 5])
    assert new_secstruct.sequence == "GGGGAACCCC"
    assert new_secstruct.structure == "((((..))))"
    seqstruct = SequenceStructure("GGAAGGAAAACCCC", "((..((....))))")
    new_secstruct = remove_nucleotides(seqstruct, [2, 3, 6])
    assert new_secstruct.sequence == "GGGGAAACCCC"
    assert new_secstruct.structure == "((((...))))"


def test_remove_nucleotide_sweep():
    seqstruct = SequenceStructure("GGGGAAAACCCC", "((((....))))")
    new_secstructs = remove_unpaired_nucleotide_sweep(seqstruct, 1)
    assert len(new_secstructs) == 4
    seq = "CGACAUGGAGUUUCGCCGAGCCUGCGAACUACAGCGAACACUCUUCGGAGUACCCGCUGCGUAGGCGUUUGACGCGAGGCUCCUAAAUCG"
    ss = "(((...((((((((((((((((((((.....(((((...((((....))))...))))))))))))..)))..))))))))))....)))"
    seqstruct = SequenceStructure(seq, ss)
    new_secstructs = remove_unpaired_nucleotide_sweep(seqstruct, 2)


def test_mutate_basepair():
    seq = "GGGAAACCC"
    ss = "(((...)))"
    struct = SecStruct(seq, ss)
    new_struct = get_basepair_mutation(struct, 0)
    assert new_struct.sequence[-1] != "C"


def test_mutate_basepairs():
    seq = "GGGGAAACCCC"
    ss = "((((...))))"
    struct = SecStruct(seq, ss)
    new_seqs = get_basepair_mutations(struct, 2)


def test_mutate_basepair_random():
    seq = "GGGGAAACCCC"
    ss = "((((...))))"
    struct = SecStruct(seq, ss)
    new_seqs = get_basepair_mutuations_random(struct, 3)
    assert new_seqs[-1] != "C"
