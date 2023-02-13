import click
import random

from seq_tools import SequenceStructure

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


def possible_nucleotide_mutations(nucleotide: str) -> list:
    """
    Given a RNA nucleotide, returns all possible other nucleotides it could be.

    :param nucleotide: A single RNA nucleotide as a string (must be one of
    'A', 'U', 'C', 'G').
    :return: A list of possible other nucleotides as strings.
    :raises ValueError: If the input nucleotide is not one of the four valid RNA
    nucleotides ('A', 'U', 'C', 'G').
    """
    nucleotides = ["A", "U", "C", "G"]
    if nucleotide not in nucleotides:
        raise ValueError("Invalid nucleotide. Must be one of 'A', 'U', 'C', 'G'.")
    return [nt for nt in nucleotides if nt != nucleotide]


def find_mutations(sequence: str, mutation_pos: list) -> list:
    """
    Given a RNA sequence and a list indicating mutation positions, returns all new
    sequences with mutations at allowed positions.

    :param sequence: A RNA sequence as a string.
    :param mutation_pos: A list of 1s and 0s indicating the positions in the sequence
     where mutations are allowed (1) or not allowed (0).
    :return: A list of new sequences with mutations at allowed positions as strings.
    """
    if len(sequence) != len(mutation_pos):
        raise ValueError(
            "The length of the mutation_pos list must be the same as the length of the "
            "sequence."
        )
    result = []
    for i, nucleotide in enumerate(sequence):
        if mutation_pos[i] == 0:
            continue
        for new_nucleotide in possible_nucleotide_mutations(nucleotide):
            new_sequence = sequence[:i] + new_nucleotide + sequence[i + 1 :]
            result.append(new_sequence)
    return result


def find_double_mutations(sequence: str, mutation_pos: list) -> list:
    """
    Given a RNA sequence and a list indicating mutation positions, returns all new
    sequences with two mutations at different allowed positions.

    :param sequence: A RNA sequence as a string.
    :param mutation_pos: A list of 1s and 0s indicating the positions in the sequence
    where mutations are allowed (1) or not allowed (0).
    :return: A list of new sequences with two mutations at different allowed positions
    as strings.
    """
    if len(sequence) != len(mutation_pos):
        raise ValueError(
            "The length of the mutation_pos list must be the same as the length of the "
            "sequence."
        )
    result = []
    for i, nucleotide1 in enumerate(sequence):
        if mutation_pos[i] == 0:
            continue
        for new_nucleotide1 in possible_nucleotide_mutations(nucleotide1):
            new_sequence1 = sequence[:i] + new_nucleotide1 + sequence[i + 1 :]
            for j, nucleotide2 in enumerate(new_sequence1):
                if j == i or mutation_pos[j] == 0:
                    continue
                for new_nucleotide2 in possible_nucleotide_mutations(nucleotide2):
                    new_sequence2 = (
                        new_sequence1[:j] + new_nucleotide2 + new_sequence1[j + 1 :]
                    )
                    result.append(new_sequence2)
    return result


def valid_base_pairs(base_pair: str) -> list:
    """Return other valid base pairs if supplied base pair is valid.

    Given a base pair as a string, returns the other valid base pairs if the
    supplied base pair is valid.

    :param base_pair: A base pair as a string.
    :return: A list of other valid base pairs as strings.
    """
    if len(base_pair) != 2:
        raise ValueError("Base pair must be a 2-letter string.")
    valid_pairs = ["CG", "GC", "AU", "UA"]
    if base_pair not in valid_pairs:
        raise ValueError("Supplied base pair is not valid.")
    return [pair for pair in valid_pairs if pair != base_pair]


def valid_base_pairs_GU(base_pair: str) -> list:
    """Return other valid base pairs if supplied base pair is valid.

    Given a base pair as a string, returns the other valid base pairs if the
    supplied base pair is valid. The valid base pairs are `CG`, `GC`, `AU`,
    `UA`, `GU`, and `UG`.

    :param base_pair: A base pair as a string.
    :return: A list of other valid base pairs as strings.
    """
    if len(base_pair) != 2:
        raise ValueError("Base pair must be a 2-letter string.")
    valid_pairs = ["CG", "GC", "AU", "UA", "GU", "UG"]
    if base_pair not in valid_pairs:
        raise ValueError("Supplied base pair is not valid.")
    return [pair for pair in valid_pairs if pair != base_pair]


def mutate_base_pairs(sequence: str, secondary_structure: str) -> list:
    """Return new sequences with mutated base pairs.

    Given a RNA sequence and its secondary structure as a string of `(` and
    `)` characters, returns all new sequences with mutations in their base
    pairs using the `valid_base_pairs_GU` function.

    :param sequence: A RNA sequence as a string.
    :param secondary_structure: A secondary structure as a string of `(` and
        `)` characters.
    :return: A list of new sequences with mutated base pairs as strings.
    """
    if len(sequence) != len(secondary_structure):
        raise ValueError("Seq and sec str must have same length.")
    stack = []
    result = []
    for i, char in enumerate(secondary_structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            j = stack.pop()
            base_pair = sequence[j] + sequence[i]
            for new_base_pair in valid_base_pairs_GU(base_pair):
                new_sequence = (sequence[:j] + new_base_pair[0] +
                                sequence[j+1:i] + new_base_pair[1] +
                                sequence[i+1:])
                result.append(new_sequence)
    return result


"""
# multi commmand format
@click.group()
def cli():
    pass


@cli.command()
def func1():
    pass
"""


# @click.command()
def cli():
    print(mutate_base_pairs("A&U", "(&)"))


if __name__ == "__main__":
    cli()
