"""
Making mutant libraries.

Two flavours: single nucleotide mutations, which ignore structure, and basepair
mutations, which keep a helix paired by changing both sides at once.

Run with: python examples/02_mutations.py
"""

import random

from rna_secstruct import SecStruct

from rna_secstruct_design.mutations import (
    find_mutations,
    find_multiple_mutations,
    get_basepair_mutation,
    get_basepair_mutations,
    get_basepair_mutuations_random,
    possible_nucleotide_mutations,
)
from rna_secstruct_design.selection import get_selection
from rna_secstruct_design.util import can_form_helix

SEQ = "GGGGAAAACCCC"
SS = "((((....))))"


def single_nucleotide_mutations():
    """Every alternative nucleotide at every allowed position."""
    print("\n--- single nucleotide mutations ---")
    print(f"at one position: A -> {possible_nucleotide_mutations('A')}")

    # a bare list of sequences
    seqs = find_mutations("AUC", exclude=[])
    print(f"all mutants of AUC: {seqs}")
    assert len(seqs) == 3 * 3

    # exclude keeps positions fixed
    seqs = find_mutations("AUC", exclude=[0, 1])
    print(f"only position 3 mutable: {seqs}")
    assert seqs == ["AUA", "AUU", "AUG"]

    # find_multiple_mutations names each mutant, which is what you want in a csv
    muts = find_multiple_mutations("AUC", 1, [])
    print(f"named: {[(m.name, m.sequence) for m in muts[:3]]}")
    assert muts[0].name == "A1U" and muts[0].sequence == "UUC"
    # names are 1-based: <original><position><new>
    assert all(m.name[0] == "AUC"[int(m.name[1]) - 1] for m in muts)

    # two mutations at a time grows fast, so exclude aggressively
    doubles = find_multiple_mutations("AUC", 2, [])
    print(f"double mutants of AUC: {len(doubles)}")
    assert len(doubles) == 27


def mutations_with_a_selection():
    """Combine a selection with a scan to mutate only the part you care about."""
    print("\n--- mutating only the hairpin loop ---")
    ss = SecStruct(SEQ, SS)
    # protect everything that is not the loop
    exclude = get_selection(ss, {"motif": {"m_type": "HAIRPIN"}, "invert": True})
    muts = find_multiple_mutations(SEQ, 1, exclude)
    print(f"{len(muts)} mutants, e.g. {[m.name for m in muts[:4]]}")
    # the hairpin motif is 6 nt (loop plus closing pair), 6 * 3 = 18
    assert len(muts) == 18
    # the helix arms are untouched in every mutant
    assert all(m.sequence[:3] == "GGG" and m.sequence[-3:] == "CCC" for m in muts)


def basepair_mutations():
    """Mutate both halves of a basepair so the helix stays paired."""
    print("\n--- basepair mutations ---")
    random.seed(0)
    struct = SecStruct(SEQ, SS)
    print(f"start  {struct.sequence}  {struct.structure}")

    # change the pair that position 0 takes part in
    mutated = get_basepair_mutation(struct, 0, new_bp="AU")
    print(f"pos 0 -> AU  {mutated.sequence}")
    assert mutated.sequence == "AGGGAAAACCCU"
    assert mutated.structure == struct.structure

    # leave new_bp out and one is picked at random, but never the original pair
    for _ in range(20):
        mutated = get_basepair_mutation(struct, 0)
        assert (mutated.sequence[0], mutated.sequence[-1]) != ("G", "C")
        assert can_form_helix(mutated.sequence[:4], mutated.sequence[-4:])

    # exhaustively walk every combination of n pairs
    singles = get_basepair_mutations(struct, 1)
    print(f"all single pair mutants: {singles}")
    # pairs 0, 1 and 2 are mutable, the pair closing the loop is skipped because
    # by default a pair sitting next to an unpaired region is left alone
    assert len(singles) == 3

    # flank_bp=True lets that closing pair change too
    with_flanks = get_basepair_mutations(struct, 1, flank_bp=True)
    print(f"including flanking pairs: {len(with_flanks)} mutants")
    assert len(with_flanks) == 4

    # gu=False keeps the library Watson-Crick only
    wc_only = get_basepair_mutations(struct, 1, gu=False, flank_bp=True)
    for seq in wc_only:
        for i in range(4):
            assert f"{seq[i]}{seq[-i - 1]}" in ("AU", "UA", "GC", "CG")
    print(f"Watson-Crick only: {wc_only}")


def random_basepair_mutations():
    """For big helices, sample instead of enumerating."""
    print("\n--- random basepair mutants ---")
    random.seed(0)
    struct = SecStruct("G" * 20 + "GAAAC" + "C" * 20, "(" * 20 + "(...)" + ")" * 20)
    sampled = get_basepair_mutuations_random(struct, 3, max_muts=5)
    print(f"{len(sampled)} sampled mutants of a 20 bp helix")
    assert len(sampled) == 5
    for seq in sampled:
        assert len(seq) == len(struct.sequence)
        # still a valid helix after mutating 3 pairs
        assert can_form_helix(seq[:20], seq[-20:])


def main():
    single_nucleotide_mutations()
    mutations_with_a_selection()
    basepair_mutations()
    random_basepair_mutations()
    print("\nok")


if __name__ == "__main__":
    main()
