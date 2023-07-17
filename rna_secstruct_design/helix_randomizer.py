from rna_secstruct import SecStruct
from rna_secstruct.motif import Motif
from vienna import fold, cofold

from rna_secstruct_design.constraints import (
    MaxRepeatingConstraint,
    MaxGCStretchConstraint,
    MaxRepeatingIncreaseConstraint,
    MaxGCStretchIncreaseConstraint,
)
from rna_secstruct_design.logger import get_logger
from rna_secstruct_design.selection import get_selection
from rna_secstruct_design.util import random_weighted_basepair

log = get_logger("HELIX-RANDOMIZER")


def generate_helix_sequence(helix: Motif, exclude, frac_gu=0.3):
    """
    Generates a random sequence for a helix motif
    :param helix: Motif object from rna_secstruct
    :param exclude: list of indices to not change sequence
    :param frac_gu: fraction of gu basepairs
    :return: new sequence for helix
    """
    if exclude is None:
        exclude = []
    if not helix.is_helix():
        raise ValueError("motif is not a helix!")
    strand1, strand2 = helix.strands
    org_seq = helix.sequence.split("&")
    # flip second strand so we can just add like the first strand
    org_seq[1] = org_seq[1][::-1]
    seq1, seq2 = "", ""
    for i, (s1, s2) in enumerate(zip(strand1, strand2[::-1])):
        if s1 in exclude or s2 in exclude:
            seq1 += org_seq[0][i]
            seq2 += org_seq[1][i]
        else:
            bp = random_weighted_basepair(frac_gu)
            seq1 += bp[0]
            seq2 += bp[1]
    return seq1 + "&" + seq2[::-1]


def get_designable_sequence(secstruct: SecStruct, exclude=None):
    """
    Returns a sequence with all designable bases as N's and all other bases
    :param secstruct: SecStruct object contains sequence and strucutre
    :param exclude: list of indices to not change sequence
    :return: sequence with designable bases as N's
    """
    if exclude is None:
        exclude = []
    design_sequence = ""
    for i, (seq, ss) in enumerate(zip(secstruct.sequence, secstruct.structure)):
        if i in exclude:
            design_sequence += seq
        elif ss != ".":
            design_sequence += "N"
        else:
            design_sequence += seq
    return design_sequence


class HelixRandomizer(object):
    def __init__(self):
        # do not want to repeat a base more than 4 times in a helix
        self.h_repeat_constraint = MaxRepeatingConstraint(4)
        # do not want more than 3 gcs in a row
        self.h_gc_constraint = MaxGCStretchConstraint(3)

    def __get_randomized_helix_sequence(self, h, exclude):
        for i in range(100):
            h_seq = generate_helix_sequence(h, exclude)
            if not self.h_repeat_constraint.satisifes(h_seq):
                continue
            if not self.h_gc_constraint.satisifes(h_seq, h.structure):
                continue
            if i == 99:
                log.warning(
                    "Could not find a helix sequence that satisfies constraints"
                )
            return h_seq

    def run(self, secstruct, exclude=None, attempts=10):
        log.debug("running helix randomizer")
        log.debug(f"exclude: {exclude}")
        secstruct = SecStruct(secstruct.sequence, secstruct.structure)
        if exclude is None:
            log.debug("no exclude given, using flanks")
            exclude = get_selection(secstruct, {"flanks": ""})
        else:
            log.debug("using exclude given but ensuring flanks are addded")
            log.debug(f"initial exclude length : {len(exclude)}")
            flank_exclude = get_selection(secstruct, {"flanks": ""})
            exclude = list(set(exclude + flank_exclude))

        designable_sequence = get_designable_sequence(secstruct, exclude)
        repeat_constraint = MaxRepeatingIncreaseConstraint(4, designable_sequence)
        gc_constraint = MaxGCStretchIncreaseConstraint(
            3, designable_sequence, secstruct.structure
        )
        best = 1000
        best_seq = ""
        count = 0
        seq_count = 0
        use_cofold = False
        if secstruct.sequence.count("&") > 0:
            use_cofold = True
            log.debug("using cofold for design")
        while True:
            seq_count += 1
            if seq_count > 1000:
                log.warn("Could not find a sequence that satisfies constraints")
                break
            for h in secstruct.get_helices():
                new_seq = self.__get_randomized_helix_sequence(h, exclude)
                secstruct.change_motif(h.m_id, new_seq, h.structure)
            if not repeat_constraint.satisifes(secstruct.sequence):
                continue
            if not gc_constraint.satisifes(secstruct.sequence, secstruct.structure):
                continue
            if not use_cofold:
                r = fold(secstruct.sequence)
            else:
                r = cofold(secstruct.sequence)
            if r.dot_bracket != secstruct.structure:
                continue
            if r.ens_defect < best:
                best = r.ens_defect
                best_seq = secstruct.sequence
            count += 1
            if count >= attempts:
                break
        log.debug(f"sequence: {best_seq} and ens_defect: {best}")
        return best, best_seq
