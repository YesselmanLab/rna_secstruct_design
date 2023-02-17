import re
import random
from rna_secstruct.secstruct import SecStruct
from seq_tools import SequenceStructure

BASEPAIRS = ["AU", "UA", "GC", "CG", "GU", "UG"]
BASEPAIRS_WC = ["AU", "UA", "GC", "CG"]
BASEPAIRS_GU = ["GU", "UG"]


def random_helix(length, gu=0) -> SequenceStructure:
    """
    generate a random helix
    """
    seq_1 = ""
    seq_2 = ""
    basepairs = ["AU", "UA", "GC", "CG", "GU", "UG"]
    basepairs_wc = ["AU", "UA", "GC", "CG"]
    bps = []
    for _ in range(0, gu):
        bps.append(random.choice(basepairs))
    for _ in range(0, length - gu):
        bps.append(random.choice(basepairs_wc))
    random.shuffle(bps)
    for bp in bps:
        seq_1 += bp[0]
        seq_2 = bp[1] + seq_2
    seq = seq_1 + "&" + seq_2
    ss = "(" * length + "&" + ")" * length
    return SequenceStructure(seq, ss)


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


def random_basepair():
    """
    Return a random base pair.
    :return: A random base pair as a string.
    """
    return random.choice(BASEPAIRS)


def random_gu_basepair():
    return random.choice(BASEPAIRS_GU)


def random_wc_basepair():
    return random.choice(BASEPAIRS_WC)


def random_weighted_basepair(frac_gu=0.3):
    if random.random() > frac_gu:
        return random_wc_basepair()
    else:
        return random_gu_basepair()


def max_repeating_nucleotides(sequence: str) -> dict:
    """
    Return max number of repeating nucleotides in a row for each.

    Given a RNA sequence, returns a dictionary where the keys are the nucleotides
    (A, C, G, U) and the values are the maximum number of repeating nucleotides in a row
    for each nucleotide.

    :param sequence: A RNA sequence as a string.
    :return: A dictionary where the keys are the nucleotides (A, C, G, U) and the values
     are the max number of repeating nucleotides in a row for each nucleotide.
    """
    result = {nucleotide: 0 for nucleotide in "ACGU"}
    for nucleotide in "ACGU":
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


def max_gc_stretch(sequence, structure) -> int:
    """
    Return the maximum number of GC basepairs in a row
    :param sequence: sequence of RNA
    :param structure: structure of RNA
    """

    def helix_max_gc_stretch(s1, s2):
        i = -1
        j = len(s2)
        max_count = 0
        count = 0
        while i < len(s1) - 1:
            i += 1
            j -= 1
            flag = 0
            if s1[i] == "G" and s2[j] == "C":
                flag = 1
            elif s1[i] == "C" and s2[j] == "G":
                flag = 1
            if flag:
                count += 1
                continue
            else:
                if count > max_count:
                    max_count = count
                count = 0
        if count > max_count:
            max_count = count
        return max_count

    ss = SecStruct(sequence, structure)

    longest_gc_stretch = 0
    for h in ss.get_helices():
        seq = h.sequence.split("&")
        current_gc_stretch = helix_max_gc_stretch(seq[0], seq[1])
        if current_gc_stretch > longest_gc_stretch:
            longest_gc_stretch = current_gc_stretch
    return longest_gc_stretch


def hamming(a, b):
    """hamming distance between two strings"""
    dist = 0
    for i, j in zip(a, b):
        if i != j:
            dist += 1
    return dist


def can_form_helix(sequence1: str, sequence2: str) -> bool:
    """
    Check if two RNA sequences can form a helix.

    Given two RNA sequences, returns True if the sequences can form a helix and
    False otherwise. A helix is formed when the base pairs of the two sequences
    complement each other in the manner sequence1[i] and sequence2[-i-1].

    :param sequence1: The first RNA sequence as a string.
    :param sequence2: The second RNA sequence as a string.
    :return: True if the sequences can form a helix, False otherwise.
    """
    if len(sequence1) != len(sequence2):
        return False

    base_pairs = {"AU", "UA", "GC", "CG", "GU", "UG"}
    for i in range(len(sequence1)):
        if f"{sequence1[i]}{sequence2[-i - 1]}" not in base_pairs:
            return False

    return True
