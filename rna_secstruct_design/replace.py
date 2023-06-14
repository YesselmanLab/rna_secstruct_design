from seq_tools.structure import SequenceStructure
from seq_tools.structure import find as find_structure


def replace_seq_structures(org_seq_struct, sub_seq_struct, new_seq_struct):
    def replace_substrings(s, bounds, replacements):
        new_s = ""
        prev_pos = 0
        for i, pos in enumerate(bounds):
            new_s += s[prev_pos : pos[0]] + replacements[i]
            prev_pos = pos[1]
        new_s += s[prev_pos:]
        return new_s

    bounds = find_structure(org_seq_struct, sub_seq_struct)
    if len(bounds) == 0:
        raise ValueError("cannot find substructure in original sequence/structure")
    elif len(bounds) > 1:
        raise ValueError("found multiple substructures in original sequence/structure")
    bounds = bounds[0]
    seqs = [seq_ss.sequence for seq_ss in new_seq_struct.split_strands()]
    sss = [seq_ss.structure for seq_ss in new_seq_struct.split_strands()]
    new_sequence = replace_substrings(org_seq_struct.sequence, bounds, seqs)
    new_structure = replace_substrings(org_seq_struct.structure, bounds, sss)
    return SequenceStructure(new_sequence, new_structure)


def replace_gaaa_w_uucg(seq_struct):
    sub_seq_struct = SequenceStructure("GGAAAC", "(....)")
    new_seq_struct = SequenceStructure("CUUCGG", "(....)")
    return replace_seq_structures(seq_struct, sub_seq_struct, new_seq_struct)
