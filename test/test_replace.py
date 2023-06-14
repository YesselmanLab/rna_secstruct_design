from seq_tools.structure import SequenceStructure
from rna_secstruct_design.replace import replace_seq_structures, replace_gaaa_w_uucg


class TestReplaceSeqStructure:
    def test_hairpin(self):
        org_ss = SequenceStructure("GGGGAAAACCCC", "((((....))))")
        sub_ss = SequenceStructure("GAAAAC", "(....)")
        new_sub_struct = SequenceStructure("GUUUUC", "((..))")
        new_ss = replace_seq_structures(org_ss, sub_ss, new_sub_struct)
        assert new_ss.sequence == "GGGGUUUUCCCC"
        assert new_ss.structure == "(((((..)))))"

    def test_helices(self):
        org_ss = SequenceStructure("GGGGAAAACCCC", "((((....))))")
        sub_ss = SequenceStructure("GGGG&CCCC", "((((&))))")
        new_sub_struct = SequenceStructure("GACG&CCAC", "(..(&)..)")
        new_ss = replace_seq_structures(org_ss, sub_ss, new_sub_struct)
        assert new_ss.sequence == "GACGAAAACCAC"
        assert new_ss.structure == "(..(....)..)"


def test_replace_gaaa_to_uucg():
    ss = SequenceStructure("GGGGGAAACCCC", "((((....))))")
    new_ss = replace_gaaa_w_uucg(ss)
    assert new_ss.sequence == "GGGCUUCGGCCC"
    assert new_ss.structure == "((((....))))"
