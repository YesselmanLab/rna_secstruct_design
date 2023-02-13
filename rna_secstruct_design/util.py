import re


def str_to_range(x):
    """
    Convert a string representation of a range of numbers to a list of integers.

    Given a string representation of a range of numbers, this function returns a
    list of integers corresponding to the numbers in the range. The string can
    contain single numbers separated by commas, and ranges of numbers separated by
    a hyphen.

    :param x: A string representation of a range of numbers.
    :return: A list of integers corresponding to the numbers in the range.
    """
    return sum(
        (
            i if len(i) == 1 else list(range(i[0], i[1] + 1))
            for i in (
                [int(j) for j in i if j] for i in re.findall(r"(\d+),?(?:-(\d+))?", x)
            )
        ),
        [],
    )


def max_repeating_nucleotides(sequence: str) -> dict:
    """
    Return max number of repeating nucleotides in a row for each.

    Given a RNA sequence, returns a dictionary where the keys are the nucleotides
    and the values are the maximum number of repeating nucleotides in a row for
    each nucleotide.

    :param sequence: A RNA sequence as a string.
    :return: A dictionary where the keys are the nucleotides and the values are
        the max number of repeating nucleotides in a row for each nucleotide.
    """
    result = {nucleotide: 0 for nucleotide in set(sequence)}
    for nucleotide in set(sequence):
        count = 0
        max_count = 0
        for char in sequence:
            if char == nucleotide:
                count += 1
            else:
                max_count = max(max_count, count)
                count = 0
        result[nucleotide] = max(max_count, count)
    return result
