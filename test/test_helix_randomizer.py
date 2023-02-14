from rna_secstruct_design.helix_randomizer import (
    generate_helix_sequence,
    HelixRandomizer,
)
from rna_secstruct_design.util import can_form_helix
from rna_secstruct.secstruct import SecStruct


def test_generate_helix_sequence():
    ss = SecStruct("AAGGGGAAAACCCC", "..((((....))))")
    helix = ss.motifs[1]
    seqs = generate_helix_sequence(helix, None)
    assert can_form_helix(seqs[0], seqs[1])
    seqs = generate_helix_sequence(helix, [2, 3])
    assert can_form_helix(seqs[0], seqs[1])
    assert seqs[0][0:2] == "GG"
    assert seqs[1][-2:] == "CC"
