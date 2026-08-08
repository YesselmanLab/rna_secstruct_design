from rna_secstruct.secstruct import SecStruct

from rna_secstruct_design.util import (
    can_form_helix,
    find_seq_struct,
    random_helix,
    str_to_range,
    max_repeating_nucleotides,
    max_gc_stretch,
)


class TestStrToRange:
    def test_empty_string(self):
        """
        Test that an empty string returns an empty list.
        """
        assert str_to_range("") == []

    def test_single_number(self):
        """
        Test that a single number returns a list with that number.
        """
        assert str_to_range("1") == [1]

    def test_comma_separated_numbers(self):
        """
        Test that a comma separated list of numbers returns a list with those numbers.
        """
        assert str_to_range("1,2,3") == [1, 2, 3]

    def test_range_of_numbers(self):
        """
        Test that a range of numbers returns a list with those numbers.
        """
        assert str_to_range("1-3") == [1, 2, 3]

    def test_mixed_comma_separated_and_range(self):
        """
        Test that a mixed list of comma separated and range of numbers returns a list with
        those numbers.
        """
        assert str_to_range("1-3,5") == [1, 2, 3, 5]

    def test_mixed_comma_separated_and_range_with_spaces(self):
        """
        Test that a mixed list of comma separated and range of numbers with spaces
        returns a list with those numbers.
        """
        assert str_to_range("1-3, 5") == [1, 2, 3, 5]

    def test_two_ranges(self):
        """
        Test that two ranges of numbers returns a list with those numbers.
        """
        assert str_to_range("1-3,5-7") == [1, 2, 3, 5, 6, 7]


class TestMaxRepeatingNucleotides:
    def test_empty_sequence(self):
        """Test empty sequence"""
        sequence = ""
        expected = {"A": 0, "C": 0, "G": 0, "U": 0}
        result = max_repeating_nucleotides(sequence)
        assert result == expected

    def test_single_nucleotide(self):
        """Test sequence with only one nucleotide"""
        sequence = "A" * 10
        expected = {"A": 10, "C": 0, "G": 0, "U": 0}
        result = max_repeating_nucleotides(sequence)
        assert result == expected

    def test_multiple_nucleotides(self):
        """Test sequence with multiple nucleotides"""
        sequence = "ACGUA" * 10
        expected = {"A": 2, "C": 1, "G": 1, "U": 1}
        result = max_repeating_nucleotides(sequence)
        assert result == expected

    def test_long_repeating_nucleotides(self):
        """Test sequence with long repeating nucleotides"""
        sequence = "AAAGGUUCC" * 10
        expected = {"A": 3, "C": 2, "G": 2, "U": 2}
        result = max_repeating_nucleotides(sequence)
        assert result == expected

    def test_mixed_repeating_nucleotides(self):
        """Test sequence with mixed repeating nucleotides"""
        sequence = "AAAGGUUCCAGGU" * 10
        expected = {"A": 3, "C": 2, "G": 2, "U": 2}
        result = max_repeating_nucleotides(sequence)
        assert result == expected


def test_max_gc_stretch():
    """Test max_gc_stretch"""
    assert max_gc_stretch("GGGGAAAACCCC", "((((....))))") == 4
    assert max_gc_stretch("CAGGAAAACCUG", "((((....))))") == 2


def test_can_form_helix():
    """Test can_form_helix"""
    assert can_form_helix("GAC", "GUC")
    assert can_form_helix("AGAC", "GUCU")
    assert not can_form_helix("GUC", "GUG")


class TestFindSeqStruct:
    def test_single_strand(self):
        """A single strand match returns one bound"""
        full = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        sub = SecStruct("GGGGAAAACCCC", "((((....))))")
        assert find_seq_struct(full, sub) == [[(2, 14)]]

    def test_two_strands(self):
        """A two strand match returns one bound per strand"""
        full = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        sub = SecStruct("GGGG&CCCC", "((((&))))")
        assert find_seq_struct(full, sub) == [[(2, 6), (10, 14)]]

    def test_structure_must_match(self):
        """The sequence matching is not enough, the structure has to match too"""
        full = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        assert find_seq_struct(full, SecStruct("GGGG", "((((")) == [[(2, 6)]]
        assert find_seq_struct(full, SecStruct("GGGG", "....")) == []

    def test_multiple_matches(self):
        """Every placement is returned"""
        full = SecStruct("GAAACGAAAC", "(...)(...)")
        assert find_seq_struct(full, SecStruct("GAAAC", "(...)")) == [
            [(0, 5)],
            [(5, 10)],
        ]

    def test_strands_do_not_overlap(self):
        """Strands are matched in order and cannot reuse the same positions"""
        full = SecStruct("GAAAC", "(...)")
        assert find_seq_struct(full, SecStruct("G&G", "(&(")) == []

    def test_no_match(self):
        """A missing substructure gives an empty list rather than raising"""
        full = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
        assert find_seq_struct(full, SecStruct("UUUU", "....")) == []


def test_random_helix():
    """random_helix builds a valid two strand SecStruct"""
    helix = random_helix(5)
    assert isinstance(helix, SecStruct)
    assert helix.structure == "(((((&)))))"
    strand_1, strand_2 = helix.sequence.split("&")
    assert can_form_helix(strand_1, strand_2)


def test_random_helix_gu():
    """Asking for all GU pairs gives a helix of only GU pairs"""
    helix = random_helix(4, gu=4)
    strand_1, strand_2 = helix.sequence.split("&")
    assert can_form_helix(strand_1, strand_2)
    for i, nt in enumerate(strand_1):
        assert f"{nt}{strand_2[-i - 1]}" in ("GU", "UG")
