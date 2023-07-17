from rna_secstruct_design.helix_randomizer import (
    generate_helix_sequence,
    HelixRandomizer,
)
from rna_secstruct_design.util import can_form_helix
from rna_secstruct_design.selection import get_selection
from rna_secstruct.secstruct import SecStruct


def test_generate_helix_sequence():
    ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
    helix = ss.motifs[1]
    seqs = generate_helix_sequence(helix, None).split("&")
    assert can_form_helix(seqs[0], seqs[1])
    seqs = generate_helix_sequence(helix, [2, 3]).split("&")
    assert can_form_helix(seqs[0], seqs[1])
    assert seqs[0][0:2] == "GG"
    assert seqs[1][-2:] == "CC"


class TestHelixRandomizer:
    def test_simple(self):
        hr = HelixRandomizer()
        secstruct = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        ens_defect, seq = hr.run(secstruct)
        assert ens_defect < 1
        assert seq[5:11] == "GAAAAC"

    def test_long_helix(self):
        hr = HelixRandomizer()
        secstruct = SecStruct(
            "G" * 20 + "UUCG" + "C" * 20, "(" * 20 + "...." + ")" * 20
        )
        ens_defect, seq = hr.run(secstruct)
        assert ens_defect < 1

    def test_mttr(self):
        seq = (
            "GUUGAUAUGGAUUUACUCCGAGGAGACGAACUACCACGAACAGGGGAAACUCUACCCGUGGCGUCUCCGUU"
            "UGACGAGUAAGUCCUAAGUCAACAAAGUCCGCGAGUAGCGGACAC"
        )
        ss = (
            "((((((..((((((((((((((((((((.....(((((...((((....))))...))))))))))))..)"
            "))..))))))))))...))))))...((((((.....)))))).."
        )
        secstruct = SecStruct(seq, ss)
        params = {"motif": {"max_id": 1, "extend_flank": 2}}
        exclude = get_selection(secstruct, params)
        hr = HelixRandomizer()
        ens_defect, seq = hr.run(secstruct, exclude)
        assert ens_defect < 5
        assert seq.startswith("GUUGAUAUGG")

    def test_tlr_extended(self):
        seq = "AAGAUAUGGCGUGGAGUACCGGAAACGGUGACUUUGAUACUGCCUAAGUCUU"
        ss = "(((((..((((((((((((((....)))).))))...))).)))...)))))"
        secstruct = SecStruct(seq, ss)
        params = {
            "seq_struct_1": {
                "name": "tlr_extended",
            },
            "motif_2": {"name": "gaaa_tetraloop", "extend_flank": 1},
        }
        exclude = get_selection(secstruct, params)
        correct = [
            3,
            4,
            5,
            6,
            7,
            8,
            42,
            43,
            44,
            45,
            46,
            47,
            48,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
        ]
        assert exclude == correct
