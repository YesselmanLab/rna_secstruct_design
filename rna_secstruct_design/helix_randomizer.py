from rna_secstruct import SecStruct
from rna_secstruct.motif import Motif
from vienna import fold

from rna_secstruct_design.constraints import (
    MaxRepeatingIncreaseConstraint,
    MaxGCStretchIncreaseConstraint,
)

from rna_secstruct_design.util import random_weighted_basepair


def generate_helix_sequence(helix: Motif, exclude, frac_gu=0.3):
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
    return seq1, seq2[::-1]


def get_excludes(exclude_seqs, exclude, sequence):
    if exclude_seqs is not None:
        for es in exclude_seqs:
            pos = sequence.find(es)
            if pos != -1:
                exclude.extend(list(range(pos, pos + len(es) + 1)))
    return exclude


class HelixRandomizer(object):
    def __init__(self):
        pass

    def __get_randomized_helix_sequence(self, h, exclude):
        for i in range(100):
            org_seq = h.sequence.split("&")
            org_s = util.compute_stretches(org_seq[0], org_seq[1])
            seq1, seq2 = self.__generate_helix_sequence(h, exclude)
            new_s = util.compute_stretches(seq1, seq2)
            if org_s.max_gc_stretch > 3 and new_s.max_gc_stretch > org_s.max_gc_stretch:
                continue
            elif new_s.max_gc_stretch > 3:
                continue
            if org_s.max_stretch_1 > 3 and new_s.max_stretch_1 > org_s.max_stretch_1:
                continue
            elif new_s.max_stretch_1 > 3:
                continue
            if org_s.max_stretch_2 > 3 and new_s.max_stretch_2 > org_s.max_stretch_2:
                continue
            elif new_s.max_stretch_2 > 3:
                continue
            return seq1 + "&" + seq2

    def run(self, sequence, structure, exclude):
        # TODO move validation to another function
        sequence = sequence.replace("T", "U")
        # log.debug("excluded: " + str(exclude))
        if exclude is None:
            exclude = []
            # log.debug("final excluded: " + str(exclude))
        s = SecStruct(sequence, structure)
        best = 1000
        best_seq = ""
        for _ in range(100):
            for h in s.get_helices():
                new_seq = self.__get_randomized_helix_sequence(h, exclude)
                s.change_motif(h.m_id, new_seq, h.structure)
            r = fold(s.sequence)
            if r.dot_bracket != structure:
                continue
            if r.ens_defect < best:
                best = r.ens_defect
                best_seq = s.sequence
        if len(best_seq) == 0:
            return DesignResults(best, rna_structure("A", "."))
        else:
            return DesignResults(best, rna_structure(best_seq, structure))
