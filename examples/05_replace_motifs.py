"""
Swapping one motif for another.

replace_seq_structures finds a sequence + structure pattern and puts a different
one in its place. The replacement does not have to be the same length or have the
same structure, only the same number of strands.

Run with: python examples/05_replace_motifs.py
"""

from rna_secstruct import SecStruct

from rna_secstruct_design.replace import replace_gaaa_w_uucg, replace_seq_structures
from rna_secstruct_design.util import find_seq_struct


def show(struct):
    return f"{struct.sequence}  {struct.structure}"


def finding_patterns():
    """find_seq_struct is the search underneath replace, useful on its own."""
    print("\n--- finding a pattern ---")
    full = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
    print(f"searching {show(full)}")

    # each match is a list of (start, end) bounds, one per strand of the pattern
    matches = find_seq_struct(full, SecStruct("GGGGAAAACCCC", "((((....))))"))
    print(f"one strand pattern  -> {matches}")
    assert matches == [[(2, 14)]]

    # a '&' in the pattern matches two strands that need not be adjacent
    matches = find_seq_struct(full, SecStruct("GGGG&CCCC", "((((&))))"))
    print(f"two strand pattern  -> {matches}")
    assert matches == [[(2, 6), (10, 14)]]

    # the structure has to match too, not just the sequence
    assert find_seq_struct(full, SecStruct("GGGG", "....")) == []
    print("same sequence but wrong structure -> no match")


def swapping_a_loop():
    """Replace a hairpin loop with a different one."""
    print("\n--- swapping a hairpin loop ---")
    org = SecStruct("GGGGAAAACCCC", "((((....))))")
    print(f"start        {show(org)}")

    # the replacement can change the structure, here two loop nucleotides pair up
    new = replace_seq_structures(
        org,
        SecStruct("GAAAAC", "(....)"),
        SecStruct("GUUUUC", "((..))"),
    )
    print(f"loop swapped {show(new)}")
    assert new.sequence == "GGGGUUUUCCCC"
    assert new.structure == "(((((..)))))"
    # the original is not modified
    assert org.sequence == "GGGGAAAACCCC"

    # GAAA to UUCG is common enough to have its own helper
    gaaa = SecStruct("GGGGGAAACCCC", "((((....))))")
    uucg = replace_gaaa_w_uucg(gaaa)
    print(f"GAAA {gaaa.sequence} -> UUCG {uucg.sequence}")
    assert uucg.sequence == "GGGCUUCGGCCC"


def swapping_both_sides_of_a_helix():
    """A two strand pattern replaces both halves of a helix at once."""
    print("\n--- swapping a helix ---")
    org = SecStruct("GGGGAAAACCCC", "((((....))))")
    new = replace_seq_structures(
        org,
        SecStruct("GGGG&CCCC", "((((&))))"),
        SecStruct("GACG&CCAC", "(..(&)..)"),
    )
    print(f"start  {show(org)}")
    print(f"result {show(new)}")
    assert new.sequence == "GACGAAAACCAC"
    assert new.structure == "(..(....)..)"


def errors():
    """Replacement is deliberately strict, it will not guess for you."""
    print("\n--- error cases ---")
    org = SecStruct("GGGGAAAACCCC", "((((....))))")

    try:
        replace_seq_structures(
            org, SecStruct("UUUUUU", "(....)"), SecStruct("GUUUUC", "((..))")
        )
    except ValueError as e:
        print(f"pattern absent   -> {e}")

    # ambiguity is an error rather than a silent "first match wins"
    two_loops = SecStruct("GAAACGAAAC", "(...)(...)")
    try:
        replace_seq_structures(
            two_loops, SecStruct("GAAAC", "(...)"), SecStruct("GUUUC", "(...)")
        )
    except ValueError as e:
        print(f"pattern ambiguous -> {e}")


def main():
    finding_patterns()
    swapping_a_loop()
    swapping_both_sides_of_a_helix()
    errors()
    print("\nok")


if __name__ == "__main__":
    main()
