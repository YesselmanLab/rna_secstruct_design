"""
Selecting positions in a secondary structure.

A "selection" is just a list of 0-based positions. Every design function in this
package takes one as its ``exclude`` argument, so selections are how you say
"redesign this part, leave that part alone".

Run with: python examples/01_selections.py
"""

import yaml
from rna_secstruct import SecStruct

from rna_secstruct_design.selection import get_selection, get_selection_from_motifs

# a two helix construct with a GAAA tetraloop and a tetraloop receptor
MTTR_SEQ = (
    "GUUGAUAUGGAUUUACUCCGAGGAGACGAACUACCACGAACAGGGGAAACUCUACCCGUGGCGUCUCCGUU"
    "UGACGAGUAAGUCCUAAGUCAACAAAGUCCGCGAGUAGCGGACAC"
)
MTTR_SS = (
    "((((((..((((((((((((((((((((.....(((((...((((....))))...))))))))))))..)"
    "))..))))))))))...))))))...((((((.....)))))).."
)


def show(label, positions):
    print(f"{label:<42} {positions}")


def by_motif_type():
    """Select whole motifs by what kind of motif they are."""
    print("\n--- selecting by motif type ---")
    ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
    print(f"sequence  {ss.sequence}")
    print(f"structure {ss.structure}")

    helices = get_selection_from_motifs(ss, {"m_type": "HELIX"})
    show("HELIX", helices)
    assert helices == [2, 3, 4, 5, 10, 11, 12, 13]

    hairpins = get_selection_from_motifs(ss, {"m_type": "HAIRPIN"})
    show("HAIRPIN", hairpins)
    assert hairpins == [5, 6, 7, 8, 9, 10]

    single = get_selection_from_motifs(ss, {"m_type": "SINGLESTRAND"})
    show("SINGLESTRAND", single)
    assert single == [0, 1]

    # extend_flank grows each motif outwards by n positions
    padded = get_selection_from_motifs(ss, {"m_type": "HAIRPIN", "extend_flank": 2})
    show("HAIRPIN + extend_flank 2", padded)
    assert padded == [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]


def by_name():
    """Common motifs have shorthand names so you do not retype their sequences."""
    print("\n--- selecting by name ---")
    ss = SecStruct("AAGGGGGAAACCCC", "..((((....))))")
    print(f"sequence  {ss.sequence}")
    print(f"structure {ss.structure}")

    gaaa = get_selection_from_motifs(ss, {"name": "gaaa_tetraloop"})
    show("gaaa_tetraloop", gaaa)
    assert gaaa == [5, 6, 7, 8, 9, 10]
    # the selected positions really are the GAAA tetraloop
    assert "".join(ss.sequence[i] for i in gaaa) == "GGAAAC"

    # the other names are: ref_hp, tlr and tlr_extended
    tlr = get_selection_from_motifs(SecStruct(MTTR_SEQ, MTTR_SS), {"name": "tlr"})
    show("tlr (in the mttr construct)", tlr)
    assert "".join(MTTR_SEQ[i] for i in tlr) == "UAUGCUAAG"


def combining_selections():
    """get_selection takes a dict and concatenates every sub selection."""
    print("\n--- combining selections ---")
    ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")

    # keys are matched on their prefix, so you can have several of the same kind
    params = {
        "motif_1": {"m_type": "SINGLESTRAND"},
        "motif_2": {"m_type": "HAIRPIN"},
    }
    combined = get_selection(ss, params)
    show("singlestrand + hairpin", combined)
    assert combined == [0, 1, 5, 6, 7, 8, 9, 10]

    # "range" uses 1-based inclusive ranges because that is how people talk about
    # sequences, but the returned positions are 0-based like everything else
    show("range 1-5", get_selection(ss, {"range": "1-5"}))
    assert get_selection(ss, {"range": "1-5"}) == [0, 1, 2, 3, 4]
    show("range 1-3,7", get_selection(ss, {"range": "1-3,7"}))
    assert get_selection(ss, {"range": "1-3,7"}) == [0, 1, 2, 6]

    # "flanks" grabs the closing pair of every non helix motif, these are the
    # positions you almost always want to leave alone when redesigning helices
    show("flanks", get_selection(ss, {"flanks": {}}))
    assert get_selection(ss, {"flanks": {}}) == [0, 1, 5, 10]

    # "invert" flips the whole selection, useful for "everything except ..."
    inverted = get_selection(ss, {"motif": {"m_type": "HELIX"}, "invert": True})
    show("everything except the helices", inverted)
    assert inverted == [0, 1, 6, 7, 8, 9]


def by_sequence_and_structure():
    """seq_struct matches an exact sequence + structure, including across strands."""
    print("\n--- selecting by sequence and structure ---")
    ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")

    params = {"seq_struct": {"sequence": "GGGGAAAACCCC", "structure": "((((....))))"}}
    whole = get_selection(ss, params)
    show("one strand", whole)
    assert whole == list(range(2, 14))

    # a '&' splits the pattern into strands that are matched in order, so this
    # picks up both sides of the helix without the loop in between
    params = {"seq_struct": {"sequence": "GGGG&CCCC", "structure": "((((&))))"}}
    two_strand = get_selection(ss, params)
    show("two strands", two_strand)
    assert two_strand == [2, 3, 4, 5, 10, 11, 12, 13]


def from_yaml():
    """The same params can live in a yaml file, which is what the cli reads."""
    print("\n--- selections from yaml ---")
    text = """
    motif_1:
        name: gaaa_tetraloop
        extend_flank: 1
    motif_2:
        name: tlr
        extend_flank: 2
    invert: True
    """
    params = yaml.safe_load(text)
    ss = SecStruct(MTTR_SEQ, MTTR_SS)
    positions = get_selection(ss, params)
    # the tetraloop is 6 nt plus 1 either side, the receptor is 9 nt over two
    # strands plus 2 either side of each, so 8 + 17 = 25 positions are protected
    print(
        f"everything outside the tetraloop and its receptor: {len(positions)} of "
        f"{len(MTTR_SEQ)} positions"
    )
    assert len(positions) == len(MTTR_SEQ) - 25


def main():
    by_motif_type()
    by_name()
    combining_selections()
    by_sequence_and_structure()
    from_yaml()
    print("\nok")


if __name__ == "__main__":
    main()
