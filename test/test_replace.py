import pytest
from rna_secstruct.secstruct import SecStruct

from rna_secstruct_design.replace import replace_seq_structures, replace_gaaa_w_uucg


class TestReplaceSeqStructure:
    def test_hairpin(self):
        org_ss = SecStruct("GGGGAAAACCCC", "((((....))))")
        sub_ss = SecStruct("GAAAAC", "(....)")
        new_sub_struct = SecStruct("GUUUUC", "((..))")
        new_ss = replace_seq_structures(org_ss, sub_ss, new_sub_struct)
        assert new_ss.sequence == "GGGGUUUUCCCC"
        assert new_ss.structure == "(((((..)))))"

    def test_helices(self):
        org_ss = SecStruct("GGGGAAAACCCC", "((((....))))")
        sub_ss = SecStruct("GGGG&CCCC", "((((&))))")
        new_sub_struct = SecStruct("GACG&CCAC", "(..(&)..)")
        new_ss = replace_seq_structures(org_ss, sub_ss, new_sub_struct)
        assert new_ss.sequence == "GACGAAAACCAC"
        assert new_ss.structure == "(..(....)..)"

    def test_returns_secstruct(self):
        org_ss = SecStruct("GGGGAAAACCCC", "((((....))))")
        new_ss = replace_seq_structures(
            org_ss, SecStruct("GAAAAC", "(....)"), SecStruct("GUUUUC", "((..))")
        )
        assert isinstance(new_ss, SecStruct)
        # the original is left untouched
        assert org_ss.sequence == "GGGGAAAACCCC"

    def test_not_found(self):
        org_ss = SecStruct("GGGGAAAACCCC", "((((....))))")
        with pytest.raises(ValueError, match="cannot find substructure"):
            replace_seq_structures(
                org_ss, SecStruct("UUUUUU", "(....)"), SecStruct("GUUUUC", "((..))")
            )

    def test_multiple_matches(self):
        org_ss = SecStruct("GAAACGAAAC", "(...)(...)")
        with pytest.raises(ValueError, match="found multiple substructures"):
            replace_seq_structures(
                org_ss, SecStruct("GAAAC", "(...)"), SecStruct("GUUUC", "(...)")
            )

    def test_strand_count_mismatch(self):
        org_ss = SecStruct("GGGGAAAACCCC", "((((....))))")
        with pytest.raises(ValueError, match="strands"):
            replace_seq_structures(
                org_ss, SecStruct("GGGG&CCCC", "((((&))))"), SecStruct("GACG", "(..(")
            )


def test_replace_gaaa_to_uucg():
    ss = SecStruct("GGGGGAAACCCC", "((((....))))")
    new_ss = replace_gaaa_w_uucg(ss)
    assert new_ss.sequence == "GGGCUUCGGCCC"
    assert new_ss.structure == "((((....))))"
