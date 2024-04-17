import itertools
import random
from typing import List, Dict, Tuple
from dataclasses import dataclass

from rna_secstruct import SecStruct
from rna_secstruct.parser import ConnectivityList
from seq_tools.structure import SequenceStructure

from rna_secstruct_design.util import random_helix


# introduce mutations into the sequence at allowed positions #########################


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


# mutate basepairs ###################################################################


def possible_basepair_mutations(bp: str, gu=True) -> list:
    """
    Given a RNA basepair such as GC or CG returns all possible other basepairs it could be.
    """
    bps = ["AU", "UA", "CG", "GC", "GU", "UG"]
    if not gu:
        bps = ["AU", "UA", "CG", "GC"]
    if bp not in bps:
        raise ValueError(f"Invalid basepair: {bp}")
    return [b for b in bps if b != bp]


def get_basepair_mutation(struct: SecStruct, pos, new_bp=None, gu=True) -> SecStruct:
    """
    Given a secondary structure and a basepair position, returns a new secondary structure
    with a new basepair at the position specified if the indentity is not specified pick
    at random
    """
    cl = ConnectivityList(struct.sequence, struct.structure)
    if not cl.is_nucleotide_paired(pos):
        raise ValueError("position must be a basepair")

    max_val = pos
    min_val = cl.get_paired_nucleotide(pos)
    if cl.get_paired_nucleotide(pos) > pos:
        max_val = cl.get_paired_nucleotide(pos)
        min_val = pos
    if new_bp is None:
        new_bp = random.choice(
            possible_basepair_mutations(cl.get_basepair(min_val), gu)
        )
    sequence = struct.sequence
    new_sequence = (
        sequence[:min_val]
        + new_bp[0]
        + sequence[min_val + 1 : max_val]
        + new_bp[1]
        + sequence[max_val + 1 :]
    )
    return SecStruct(new_sequence, struct.structure)


def get_basepair_mutations(
    struct: SecStruct, num: int, exclude=None, gu=True, flank_bp=False
) -> List[str]:
    if exclude is None:
        exclude = []
    allowed_pos = []
    for i in range(0, len(struct.structure)):
        # exclude positions in exclude
        if i in exclude:
            continue
        # has to be a opening pair
        if struct.structure[i] != "(":
            continue
        # if we are allowing flanking basepairs then accept
        if flank_bp:
            allowed_pos.append(i)
            continue
        if i == 0:
            allowed_pos.append(i)
            continue
        # are we flanking a non basepair then dont accept
        if struct.structure[i - 1] == "." or struct.structure[i + 1] == ".":
            continue
        allowed_pos.append(i)
    cl = ConnectivityList(struct.sequence, struct.structure)
    seen = []
    combos = itertools.product(allowed_pos, repeat=num)
    sequence = struct.sequence
    sequences = []
    for muts in combos:
        muts = sorted(muts)
        if len(set(muts)) != len(muts):
            continue
        mut_key = "-".join([str(x) for x in muts])
        if mut_key in seen:
            continue
        for m in muts:
            new_bp = random.choice(possible_basepair_mutations(cl.get_basepair(m), gu))
            sequence = (
                sequence[:m]
                + new_bp[0]
                + sequence[m + 1 : cl.get_paired_nucleotide(m)]
                + new_bp[1]
                + sequence[cl.get_paired_nucleotide(m) + 1 :]
            )
        seen.append(mut_key)
        sequences.append(sequence)
    return sequences


# change helix length ###############################################################


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
            # has to flip the sequence to make the algorithm easier
            seq_2 = m_struct[1].sequence[::-1]
            seq_2 = seq_2[: new_length - 1] + seq_2[-1]
            seq_2 = seq_2[::-1]
            strand_2 = SequenceStructure(seq_2, ")" * len(seq_2))
            new_structs = [
                m_struct[0][0 : new_length - 1] + m_struct[0][-1],
                strand_2,
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


def scan_helix_lengths(
    struct: SecStruct, pos, min_length, max_length
) -> List[SecStruct]:
    """
    Given a secondary structure and a helix position, returns all possible
    secondary structures with helices of lengths between min_length and
    max_length.
    """
    structs = []
    for i in range(min_length, max_length + 1):
        structs.append(change_helix_length(struct, pos, i))
    return structs


def scan_all_helix_lengths(
    struct: SecStruct, h_ranges: Dict[int, Tuple[int, int]]
) -> List[SecStruct]:
    """
    Iterates though all possible helical lengths in a given secstruct
    object. The helical lengths are defined by a dictionary of the position
    of the helix in secstruct.
    :param struct: a secondary structure
    :param h_ranges: a dictionary of the positions of the helices and the
        minimum and maximum lengths of the helices
    :return: a list of secondary structures with all possible helical lengths
    """
    ranges = []
    helix_pos = []
    structs = []
    for k, v in h_ranges.items():
        ranges.append(list(range(v[0], v[1] + 1)))
        helix_pos.append(k)
    for i in itertools.product(*ranges):
        new_struct = struct.get_copy()
        for j, pos in enumerate(helix_pos):
            new_struct = change_helix_length(new_struct, pos, i[j])
        structs.append(new_struct)
    return structs


# unpaired insertions ################################################################


def add_unpaired(
    struct: SequenceStructure, pos, bulge_size, all_nucleotides=False
) -> List[SequenceStructure]:
    """
    Add an unpaired residue to an RNA.
    :param struct: a secondary structure
    :param pos: the position of the helix to change
    :param bulge_size: the size of the bulge to add
    :return: a new secondary structure with the bulge added
    """
    if pos < 1:
        raise ValueError("position must be greater than 0")
    if bulge_size < 1:
        raise ValueError("bulge size must be greater than 0")
    sequence = list(struct.sequence)
    structure = list(struct.structure)
    nucs = itertools.product("AUCG", repeat=bulge_size)
    nucs = ["".join(x) for x in list(nucs)]
    structs = []
    for nuc in nucs:
        new_seq = sequence[:]
        new_str = structure[:]
        new_seq.insert(pos, nuc)
        new_str.insert(pos, "." * len(nuc))
        new_struct = SequenceStructure("".join(new_seq), "".join(new_str))
        structs.append(new_struct)
        if not all_nucleotides:
            break
    return structs


def add_unpaired_sweep(
    struct: SequenceStructure, n_include, exclude=None, all_nucleotides=False
) -> List[SequenceStructure]:
    """
    Create a list of secondary structures each with n unpaired nucleotide added.
    :param struct: a secondary structure
    :param n_include: the number of unpaired nucleotides to add
    :param exclude: a list of positions to exclude
    :param all_nucleotides: if True, add all possible nucleotides at each position
    :return: a list of secondary structures with n unpaired nucleotides added
    """
    if exclude is None:
        exclude = []
    if n_include < 1:
        raise ValueError("n_include must be greater than 0")
    unpaired = []
    for i, char in enumerate(struct.structure):
        if i < 1:
            continue
        if i in exclude:
            continue
        if char != ".":
            unpaired.append(i)
    if len(unpaired) < n_include:
        raise ValueError(
            "n_include must be less than the number of unpaired nucleotides"
        )
    combos = itertools.combinations(unpaired, n_include)
    structs = []
    for combo in combos:
        combo = sorted(combo)
        new_structs = [SequenceStructure(struct.sequence, struct.structure)]
        cur_structs = []
        while len(combo) > 0:
            pos = combo.pop()
            for s in new_structs:
                new_structs = add_unpaired(s, pos, 1, all_nucleotides)
                for ns in new_structs:
                    if ns.structure.find(".(.") != -1 or ns.structure.find(".).") != -1:
                        continue
                    cur_structs.append(ns)
            new_structs = cur_structs
            cur_structs = []
        structs.extend(new_structs)
    return structs


# deletions ##########################################################################


def remove_nucleotides(struct: SequenceStructure, pos) -> SequenceStructure:
    """
    Remove residues from an RNA.
    :param struct: a secondary structure
    :param pos: a list of positions to remove
    :return: a new secondary structure with residues removed
    """
    pos = sorted(pos, reverse=True)
    sequence = list(struct.sequence)
    structure = list(struct.structure)
    for p in pos:
        sequence = sequence[:p] + sequence[p + 1 :]
        structure = structure[:p] + structure[p + 1 :]
    new_struct = SequenceStructure("".join(sequence), "".join(structure))
    return new_struct


def remove_unpaired_nucleotide_sweep(
    struct: SequenceStructure, n_remove: int, exclude=None
) -> List[SequenceStructure]:
    """
    Create a list of secondary structures each with n unpaired nucleotide removed.
    :param struct: a secondary structure
    :param n_remove: the number of unpaired nucleotides to remove
    """
    if exclude is None:
        exclude = []
    if n_remove < 1:
        raise ValueError("n_remove must be greater than 0")
    sequence = list(struct.sequence)
    structure = list(struct.structure)
    unpaired = []
    for i, char in enumerate(structure):
        if i in exclude:
            continue
        if char == ".":
            unpaired.append(i)
    if len(unpaired) < n_remove:
        raise ValueError(
            "n_remove must be less than the number of unpaired nucleotides"
        )
    combos = itertools.combinations(unpaired, n_remove)
    structs = []
    for combo in combos:
        new_seq = sequence[:]
        new_str = structure[:]
        for pos in combo:
            new_seq[pos] = "X"
            new_str[pos] = "X"
        new_seq = [x for x in new_seq if x != "X"]
        new_str = [x for x in new_str if x != "X"]
        new_struct = SequenceStructure("".join(new_seq), "".join(new_str))
        structs.append(new_struct)
    return structs


def remove_nucleotide_sweep(
    struct: SequenceStructure, n_remove: int, exclude: None
) -> List[SequenceStructure]:
    """
    Create a list of secondary structures each with n nucleotides removed.
    :param struct: a secondary structure
    :param n_remove: the number of nucleotides to remove
    """
    if exclude is None:
        exclude = []
    if n_remove < 1:
        raise ValueError("n_remove must be greater than 0")
    sequence = list(struct.sequence)
    structure = list(struct.structure)
    nucleotides = []
    for i, char in enumerate(structure):
        if i in exclude:
            continue
        nucleotides.append(i)
    if len(nucleotides) < n_remove:
        raise ValueError("n_remove must be less than the number of nucleotides")
    combos = itertools.combinations(nucleotides, n_remove)
    structs = []
    for combo in combos:
        new_seq = sequence[:]
        new_str = structure[:]
        for pos in combo:
            new_seq[pos] = "X"
            new_str[pos] = "X"
        new_seq = [x for x in new_seq if x != "X"]
        new_str = [x for x in new_str if x != "X"]
        new_struct = SequenceStructure("".join(new_seq), "".join(new_str))
        structs.append(new_struct)
    return structs
