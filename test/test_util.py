from rna_secstruct_design.util import str_to_range


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

