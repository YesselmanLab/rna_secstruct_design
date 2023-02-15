import itertools


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


def find_multiple_mutations(sequence: str, num: int, exclude: list) -> list:
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
    seqs = []
    for muts in combos:
        if(len(set(muts)) != len(muts)):
            continue
        key = "-".join([str(x) for x in sorted(muts)])
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
            new_sequence = sequence[:]
            for mut in mut_combo:
                new_sequence = new_sequence[:mut[0]] + mut[1] + new_sequence[mut[0] + 1 :]
            seqs.append(new_sequence)
        count += 1
    return seqs
