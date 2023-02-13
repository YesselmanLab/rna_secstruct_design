

class HelixRandomizer(object):
    def __init__(self):
        pass

    def __randomize_helix(self, h, exclude):
        for i in range(100):
            org_seq = h.sequence.split("&")
            org_s = util.compute_stretches(org_seq[0], org_seq[1])
            seq1, seq2 = self.__generate_helix_sequence(h, exclude)
            new_s = util.compute_stretches(seq1, seq2)
            if (
                org_s.max_gc_stretch > 3
                and new_s.max_gc_stretch > org_s.max_gc_stretch
            ):
                continue
            elif new_s.max_gc_stretch > 3:
                continue
            if (
                org_s.max_stretch_1 > 3
                and new_s.max_stretch_1 > org_s.max_stretch_1
            ):
                continue
            elif new_s.max_stretch_1 > 3:
                continue
            if (
                org_s.max_stretch_2 > 3
                and new_s.max_stretch_2 > org_s.max_stretch_2
            ):
                continue
            elif new_s.max_stretch_2 > 3:
                continue
            return seq1 + "&" + seq2

    def __generate_helix_sequence(self, h, exclude):
        if exclude is None:
            exclude = []
        strand1, strand2 = h.strands
        seq1, seq2 = "", ""
        for i, (s1, s2) in enumerate(zip(strand1, strand2[::-1])):
            # dont change end pairs
            if i == 0 or i == len(strand1) - 1:
                seq1 += self.sequence[s1]
                seq2 += self.sequence[s2]
            elif s1 in exclude or s2 in exclude:
                seq1 += self.sequence[s1]
                seq2 += self.sequence[s2]
            else:
                bp = util.random_weighted_basepair()
                seq1 += bp[0]
                seq2 += bp[1]
        return seq1, seq2[::-1]

    def run(self, sequence, structure, exclude=None, exclude_seqs=None):
        sequence = sequence.replace("T", "U")
        # log.debug("excluded: " + str(exclude))
        if exclude is None:
            exclude = []
        if exclude_seqs is not None:
            for es in exclude_seqs:
                pos = sequence.find(es)
                if pos != -1:
                    exclude.extend(list(range(pos, pos + len(es) + 1)))
        # log.debug("final excluded: " + str(exclude))
        self.sequence, self.structure = sequence, structure
        s = rl.SecStruct(sequence, structure)
        best = 1000
        best_seq = ""
        for _ in range(100):
            for h in s:
                if not h.is_helix():
                    continue
                new_seq = self.__randomize_helix(h, exclude)
                s.change_motif(h.m_id, new_seq, h.structure)
            r = vienna.fold(s.sequence)
            if r.dot_bracket != structure:
                continue
            if r.ens_defect < best:
                best = r.ens_defect
                best_seq = s.sequence
        if len(best_seq) == 0:
            return DesignResults(best, rna_structure("A", "."))
        else:
            return DesignResults(best, rna_structure(best_seq, structure))
