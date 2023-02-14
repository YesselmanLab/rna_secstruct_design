from rna_secstruct_design.util import (
    can_form_helix,
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
