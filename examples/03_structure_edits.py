"""
Editing the structure itself: helix lengths, bulges and deletions.

These change how long the RNA is, unlike the mutations in example 02 which keep
the length fixed.

Run with: python examples/03_structure_edits.py
"""

import random

from rna_secstruct import SecStruct

from rna_secstruct_design.mutations import (
    add_unpaired,
    add_unpaired_sweep,
    change_helix_length,
    remove_nucleotides,
    remove_unpaired_nucleotide_sweep,
    scan_all_helix_lengths,
    scan_helix_lengths,
)


def show(struct):
    return f"{struct.sequence}  {struct.structure}"


def helix_lengths():
    """Grow or shrink a helix, keeping the closing pairs it needs."""
    print("\n--- changing helix length ---")
    random.seed(0)
    struct = SecStruct("GGGGAAAACCCC", "((((....))))")
    print(f"start (4 bp)  {show(struct)}")

    # growing inserts random Watson-Crick pairs in the middle of the helix
    longer = change_helix_length(struct, 0, 6)
    print(f"grown to 6 bp {show(longer)}")
    assert longer.structure == "((((((....))))))"
    # the original closing pairs are still there
    assert longer.sequence.startswith("GG") and longer.sequence.endswith("CC")

    # shrinking keeps the outermost pair on each end
    shorter = change_helix_length(struct, 0, 2)
    print(f"cut to 2 bp   {show(shorter)}")
    assert shorter.structure == "((....))"
    assert shorter.sequence == "GGAAAACC"

    # motif 0 here is the helix, ask the SecStruct which motifs exist
    for motif in struct:
        print(f"  motif {motif.m_id}: {motif.m_type} {motif.sequence}")


def scanning_helix_lengths():
    """Build a length ladder, one construct per helix length."""
    print("\n--- scanning helix lengths ---")
    random.seed(0)
    struct = SecStruct("GGGGAAAACCCC", "((((....))))")

    ladder = scan_helix_lengths(struct, 0, 2, 6)
    for s in ladder:
        print(f"  {show(s)}")
    assert len(ladder) == 5
    assert [s.structure.count("(") for s in ladder] == [2, 3, 4, 5, 6]

    # with two helices you get every combination of the two ranges
    two_helix = SecStruct("GGAGGAAAACCCC", "((.((....))))")
    grid = scan_all_helix_lengths(two_helix, {0: [2, 4], 2: [2, 4]})
    print(f"2 helices x 3 lengths each -> {len(grid)} constructs")
    assert len(grid) == 9


def bulges():
    """Insert unpaired nucleotides."""
    print("\n--- adding unpaired nucleotides ---")
    struct = SecStruct("GGGGAAAACCCC", "((((....))))")

    # a single 2 nt insertion after position 1, every possible sequence
    variants = add_unpaired(struct, 1, bulge_size=2, all_nucleotides=True)
    print(f"all 2 nt bulges at position 1: {len(variants)}")
    assert len(variants) == 16  # 4 * 4
    print(f"  first  {show(variants[0])}")
    assert variants[0].sequence == "GAAGGGAAAACCCC"
    assert variants[0].structure == "(..(((....))))"

    # all_nucleotides=False just gives one representative, useful when you only
    # care about the structural change
    one = add_unpaired(struct, 1, bulge_size=1)
    assert len(one) == 1

    # sweep across positions instead of picking one
    exclude = [2, 3, 4, 5, 6, 7, 8]
    swept = add_unpaired_sweep(struct, 2, exclude, all_nucleotides=True)
    print(f"2 bulges swept over the allowed positions: {len(swept)}")
    assert len(swept) == 32
    assert all(len(s.sequence) == len(struct.sequence) + 2 for s in swept)


def deletions():
    """Remove nucleotides."""
    print("\n--- removing nucleotides ---")
    struct = SecStruct("GGGGAAAACCCC", "((((....))))")

    trimmed = remove_nucleotides(struct, [4, 5])
    print(f"drop positions 4 and 5: {show(trimmed)}")
    assert trimmed.sequence == "GGGGAACCCC"
    assert trimmed.structure == "((((..))))"

    # sweep every way of removing n unpaired nucleotides
    swept = remove_unpaired_nucleotide_sweep(struct, 1)
    for s in swept:
        print(f"  {show(s)}")
    assert len(swept) == 4  # the 4 loop nucleotides
    assert all(len(s.sequence) == len(struct.sequence) - 1 for s in swept)


def main():
    helix_lengths()
    scanning_helix_lengths()
    bulges()
    deletions()
    print("\nok")


if __name__ == "__main__":
    main()
