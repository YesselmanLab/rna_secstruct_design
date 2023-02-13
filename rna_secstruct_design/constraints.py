from rna_secstruct_design.util import max_repeating_nucleotides, max_gc_stretch


class SequenceConstraint:
    def __init__(self):
        pass

    def satisifes(self, sequence):
        raise NotImplementedError("apply method not implemented")


class SequenceStructureConstraint:
    def __init__(self):
        pass

    def satisifes(self, sequence, structure):
        raise NotImplementedError("apply method not implemented")


class MaxRepeatingConstraint(SequenceConstraint):
    def __init__(self, max_value):
        super().__init__()
        self.max_value = max_value

    def satisifes(self, sequence):
        max_repeating = max_repeating_nucleotides(sequence)
        for val in max_repeating.values():
            if val > self.max_value:
                return False
        return True


class MaxRepeatingIncreaseConstraint(SequenceConstraint):
    def __init__(self, max_value, org_sequence):
        super().__init__()
        self.max_value = max_value
        self.org_values = max_repeating_nucleotides(org_sequence)

    def satisifes(self, sequence):
        max_repeating = max_repeating_nucleotides(sequence)
        for key, val in max_repeating.items():
            if val <= self.max_value:
                continue
            elif val > self.org_values[key]:
                return False
        return True


class MaxGCStretchConstraint(SequenceStructureConstraint):
    def __init__(self, max_value):
        super().__init__()
        self.max_value = max_value

    def satisifes(self, sequence, structure):
        gc_stretch = max_gc_stretch(sequence, structure)
        if gc_stretch > self.max_value:
            return False
        return True


class MaxGCStretchIncreaseConstraint(SequenceStructureConstraint):
    def __init__(self, max_value, sequence, structure):
        super().__init__()
        self.max_value = max_value
        self.org_value = max_gc_stretch(sequence, structure)

    def satisifes(self, sequence, structure):
        gc_stretch = max_gc_stretch(sequence, structure)
        if gc_stretch <= self.max_value:
            return True
        elif gc_stretch > self.org_value:
            return False
        return True
