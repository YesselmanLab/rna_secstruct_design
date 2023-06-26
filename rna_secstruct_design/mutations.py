import itertools
from typing import List
from dataclasses import dataclass

from rna_secstruct import SecStruct
from seq_tools.structure import SequenceStructure

from rna_secstruct_design.util import random_helix


@dataclass(frozen=True, order=True)
class Mutation:
    name: str
    sequence: str


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


def find_mutations(sequence: str, exclude: list) -> list:
    """
    Given a RNA sequence and a list indicating mutation positions, returns all new
    sequences with mutations at allowed positions.

    :param sequence: A RNA sequence as a string.
    :param mutation_pos: A list of 1s and 0s indicating the positions in the sequence
     where mutations are allowed (1) or not allowed (0).
    :return: A list of new sequences with mutations at allowed positions as strings.
    """
    result = []
    for i, nucleotide in enumerate(sequence):
        if i in exclude:
            continue
        for new_nucleotide in possible_nucleotide_mutations(nucleotide):
            new_sequence = sequence[:i] + new_nucleotide + sequence[i + 1 :]
            result.append(new_sequence)
    return result


def find_double_mutations(sequence: str, exclude: list) -> list:
    """
    Given a RNA sequence and a list indicating mutation positions, returns all new
    sequences with two mutations at different allowed positions.

    :param sequence: A RNA sequence as a string.
    :param mutation_pos: A list of 1s and 0s indicating the positions in the sequence
    where mutations are allowed (1) or not allowed (0).
    :return: A list of new sequences with two mutations at different allowed positions
    as strings.
    """
    result = []
    for i, nucleotide1 in enumerate(sequence):
        if i in exclude:
            continue
        for new_nucleotide1 in possible_nucleotide_mutations(nucleotide1):
            new_sequence1 = sequence[:i] + new_nucleotide1 + sequence[i + 1 :]
            for j, nucleotide2 in enumerate(new_sequence1):
                if j == i or j in exclude:
                    continue
                for new_nucleotide2 in possible_nucleotide_mutations(nucleotide2):
                    new_sequence2 = (
                        new_sequence1[:j] + new_nucleotide2 + new_sequence1[j + 1 :]
                    )
                    result.append(new_sequence2)
    return result


def find_multiple_mutations(sequence: str, num: int, exclude: list) -> List[Mutation]:
    """
    Given a RNA sequence and a list indicating mutation positions, returns all new
    sequences with multiple mutations at different allowed positions.

    :param sequence: A RNA sequence as a string.
    :param num: The number of mutations to make.
    :param exclude: A list of positions in the sequence where mutations are not allowed.
    """
    allowed_pos = []
    for i in range(0, len(sequence)):
        if i not in exclude:
            allowed_pos.append(i)
    combos = itertools.product(allowed_pos, repeat=num)
    count = 0
    seen = []
    results = []
    for muts in combos:
        if len(set(muts)) != len(muts):
            continue
        muts = sorted(muts)
        key = "-".join([str(x) for x in muts])
        if key in seen:
            continue
        seen.append(key)
        mut_pos = []
        for i in muts:
            muts_at_pos = []
            for new_nucleotide in possible_nucleotide_mutations(sequence[i]):
                muts_at_pos.append([i, new_nucleotide])
            mut_pos.append(muts_at_pos)
        mut_combos = list(itertools.product(*mut_pos))
        for mut_combo in mut_combos:
            names = []
            new_sequence = sequence[:]
            for mut in mut_combo:
                names.append(new_sequence[mut[0]] + str(mut[0] + 1) + mut[1])
                new_sequence = (
                    new_sequence[: mut[0]] + mut[1] + new_sequence[mut[0] + 1 :]
                )
            results.append(Mutation("_".join(names), new_sequence))
        count += 1
    return results


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
        if char == "(":
            stack.append(i)
        elif char == ")":
            j = stack.pop()
            base_pair = sequence[j] + sequence[i]
            for new_base_pair in valid_base_pairs_GU(base_pair):
                new_sequence = (
                    sequence[:j]
                    + new_base_pair[0]
                    + sequence[j + 1 : i]
                    + new_base_pair[1]
                    + sequence[i + 1 :]
                )
                result.append(new_sequence)
    return result


def change_helix_length(struct: SecStruct, pos, new_length) -> SecStruct:
    """
    Change the length of a helix in a secondary structure.
    :param struct: a secondary structure
    :param pos: the position of the helix to change
    :param new_length: the new length of the helix
    :return: a new secondary structure with the helix length changed
    """
    if new_length < 1:
        raise ValueError("new length must be greater than 0")
    m = struct[pos]
    if not m.is_helix():
        raise ValueError("motif must be a helix to change its length")
    length = len(m.sequence.split("&")[0])
    if new_length == length:
        return struct
    struct = SecStruct(struct.sequence, struct.structure)
    if new_length < length:
        # weird special cases I think
        if new_length == 1 and m.has_parent() and m.has_children():
            raise ValueError("cannot shorten helix to 1 if it has parent and children")
        # keep flanking pair if child
        elif new_length == 1 and m.has_children():
            m_struct = SequenceStructure(m.sequence, m.structure).split_strands()
            new_struct = m_struct[0][-1].join(m_struct[1][0])
            struct.change_motif(pos, new_struct.sequence, new_struct.structure)
        # keep flanking pair of parent if there is one
        elif new_length == 1:
            m_struct = SequenceStructure(m.sequence, m.structure).split_strands()
            new_struct = m_struct[0][0].join(m_struct[1][-1])
            struct.change_motif(pos, new_struct.sequence, new_struct.structure)
        # probably a better way to do this but takes the first basepairs on the 5' end
        # and then the last basepair on the 3' end
        else:
            m_struct = SequenceStructure(m.sequence, m.structure).split_strands()
            new_structs = [
                m_struct[0][0 : new_length - 1] + m_struct[0][-1],
                m_struct[1][0 : new_length - 1] + m_struct[1][-1],
            ]
            new_struct = new_structs[0].join(new_structs[1])
            struct.change_motif(pos, new_struct.sequence, new_struct.structure)
    else:
        diff = new_length - length
        m_struct = SequenceStructure(m.sequence, m.structure).split_strands()
        new_helix = random_helix(diff).split_strands()
        mid = len(m_struct[0]) // 2
        new_structs = [
            m_struct[0].insert(mid, new_helix[0]),
            m_struct[1].insert(len(m_struct[1]) - mid, new_helix[1]),
        ]
        new_struct = new_structs[0].join(new_structs[1])
        struct.change_motif(pos, new_struct.sequence, new_struct.structure)
    return struct
