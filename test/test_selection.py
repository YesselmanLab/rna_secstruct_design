from rna_secstruct import SecStruct
from rna_secstruct_design.selection import get_selection, get_selection_from_motifs


class TestMotifSelection:
    def test_by_type(self):
        params = {"m_type": "HELIX"}
        ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        selection = get_selection_from_motifs(ss, params)
        assert selection == [2, 3, 4, 5, 10, 11, 12, 13]

    def test_by_pos(self):
        params = {"max_pos": 1}
        ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        selection = get_selection_from_motifs(ss, params)
        assert selection == [0, 1]

    def test_extend_flank(self):
        params = {"m_type": "HELIX", "extend_flank": 2}
        ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        selection = get_selection_from_motifs(ss, params)
        assert selection == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

    def test_named_motif(self):
        params = {"name": "gaaa_tetraloop"}
        ss = SecStruct("AAGGGGGAAACCCC", "..((((....))))")
        selection = get_selection_from_motifs(ss, params)
        assert selection == [5, 6, 7, 8, 9, 10]


class TestGetSelection:
    def test_motif(self):
        params = {"motif": {"m_type": "HELIX"}}
        ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        selection = get_selection(ss.sequence, ss.structure, params)
        assert selection == [2, 3, 4, 5, 10, 11, 12, 13]

    def test_flanks(self):
        params = {"flanks": {}}
        ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        selection = get_selection(ss.sequence, ss.structure, params)
        assert selection == [0, 1, 5, 10]

    def test_range(self):
        params = {"range": "1-5"}
        ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        selection = get_selection(ss.sequence, ss.structure, params)
        assert selection == [0, 1, 2, 3, 4]