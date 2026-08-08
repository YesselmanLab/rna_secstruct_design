"""
Redesigning helices so the RNA still folds into the structure you asked for.

HelixRandomizer keeps every non helix motif exactly as it is, randomizes the
helices, folds the result with ViennaRNA and keeps the best scoring sequence.

Run with: python examples/04_helix_design.py
"""

import random

from rna_secstruct import SecStruct
from vienna import fold

from rna_secstruct_design.constraints import (
    MaxGCStretchConstraint,
    MaxRepeatingConstraint,
)
from rna_secstruct_design.helix_randomizer import (
    HelixRandomizer,
    generate_helix_sequence,
)
from rna_secstruct_design.selection import get_selection
from rna_secstruct_design.util import can_form_helix

MTTR_SEQ = (
    "GUUGAUAUGGAUUUACUCCGAGGAGACGAACUACCACGAACAGGGGAAACUCUACCCGUGGCGUCUCCGUU"
    "UGACGAGUAAGUCCUAAGUCAACAAAGUCCGCGAGUAGCGGACAC"
)
MTTR_SS = (
    "((((((..((((((((((((((((((((.....(((((...((((....))))...))))))))))))..)"
    "))..))))))))))...))))))...((((((.....)))))).."
)


def one_helix():
    """The building block: a new sequence for a single helix motif."""
    print("\n--- randomizing one helix ---")
    random.seed(0)
    ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
    helix = ss.motifs[1]
    print(f"original helix {helix.sequence}")

    new_seq = generate_helix_sequence(helix, exclude=None)
    print(f"randomized     {new_seq}")
    strand_1, strand_2 = new_seq.split("&")
    assert can_form_helix(strand_1, strand_2)

    # positions in exclude keep their original nucleotides
    new_seq = generate_helix_sequence(helix, exclude=[2, 3])
    strand_1, strand_2 = new_seq.split("&")
    print(f"first 2 pairs held fixed {new_seq}")
    assert strand_1[:2] == "GG" and strand_2[-2:] == "CC"
    assert can_form_helix(strand_1, strand_2)


def whole_construct():
    """Redesign every helix in a construct and check that it folds correctly."""
    print("\n--- redesigning a hairpin ---")
    random.seed(0)
    secstruct = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
    print(f"target structure {secstruct.structure}")

    hr = HelixRandomizer()
    ens_defect, seq = hr.run(secstruct)
    print(f"designed         {seq}  ens_defect {ens_defect:.3f}")

    # the design really folds into the structure we asked for
    assert fold(seq).dot_bracket == secstruct.structure
    # the loop and its closing pair were left alone
    assert seq[5:11] == "GAAAAC"
    assert ens_defect < 1


def protecting_a_motif():
    """Use a selection to protect the parts of a construct that must not change."""
    print("\n--- redesigning around a tetraloop receptor ---")
    random.seed(0)
    secstruct = SecStruct(MTTR_SEQ, MTTR_SS)

    # keep the tetraloop, its receptor and the first helix as they are
    exclude = get_selection(
        secstruct,
        {
            "motif_1": {"name": "gaaa_tetraloop", "extend_flank": 1},
            "motif_2": {"name": "tlr", "extend_flank": 2},
            "motif_3": {"max_id": 1, "extend_flank": 2},
        },
    )
    print(f"holding {len(exclude)} of {len(MTTR_SEQ)} positions fixed")

    hr = HelixRandomizer()
    ens_defect, seq = hr.run(secstruct, exclude)
    print(f"designed  {seq}")
    print(f"ens_defect {ens_defect:.3f}")

    # an empty sequence means the randomizer gave up before finding a design
    assert seq, "helix randomizer could not satisfy the constraints"
    assert fold(seq).dot_bracket == MTTR_SS
    # every protected position kept its original nucleotide
    for i in exclude:
        assert seq[i] == MTTR_SEQ[i], f"position {i} changed"


def constraints():
    """The constraints the randomizer applies, usable on their own."""
    print("\n--- sequence constraints ---")

    # no more than n of the same nucleotide in a row
    no_long_runs = MaxRepeatingConstraint(3)
    print(f"AAGGCCUU passes max 3 repeats: {no_long_runs.satisifes('AAGGCCUU')}")
    print(f"AAAAGGCC passes max 3 repeats: {no_long_runs.satisifes('AAAAGGCC')}")
    assert no_long_runs.satisifes("AAGGCCUU")
    assert not no_long_runs.satisifes("AAAAGGCC")

    # no more than n GC pairs stacked in a row, these make helices too stable
    no_gc_block = MaxGCStretchConstraint(2)
    print(
        f"CAGGAAAACCUG passes max 2 GC: "
        f"{no_gc_block.satisifes('CAGGAAAACCUG', '((((....))))')}"
    )
    print(
        f"GGGGAAAACCCC passes max 2 GC: "
        f"{no_gc_block.satisifes('GGGGAAAACCCC', '((((....))))')}"
    )
    assert no_gc_block.satisifes("CAGGAAAACCUG", "((((....))))")
    assert not no_gc_block.satisifes("GGGGAAAACCCC", "((((....))))")


def main():
    one_helix()
    whole_construct()
    protecting_a_motif()
    constraints()
    print("\nok")


if __name__ == "__main__":
    main()
