from rna_secstruct.secstruct import SecStruct

from rna_secstruct_design.util import find_seq_struct


def replace_seq_structures(org_seq_struct, sub_seq_struct, new_seq_struct) -> SecStruct:
    """
    Swap a sub sequence/structure for a new one.

    :param org_seq_struct: the SecStruct to modify
    :param sub_seq_struct: the SecStruct to search for, may be multi strand
    :param new_seq_struct: the SecStruct to put in its place, must have the same
        number of strands as ``sub_seq_struct``
    :return: a new SecStruct with the replacement made
    """

    def replace_substrings(s, bounds, replacements):
        new_s = ""
        prev_pos = 0
        for i, pos in enumerate(bounds):
            new_s += s[prev_pos : pos[0]] + replacements[i]
            prev_pos = pos[1]
        new_s += s[prev_pos:]
        return new_s

    matches = find_seq_struct(org_seq_struct, sub_seq_struct)
    if len(matches) == 0:
        raise ValueError("cannot find substructure in original sequence/structure")
    elif len(matches) > 1:
        raise ValueError("found multiple substructures in original sequence/structure")
    bounds = matches[0]
    new_strands = new_seq_struct.split_strands()
    if len(new_strands) != len(bounds):
        raise ValueError(
            f"replacement has {len(new_strands)} strands but the substructure "
            f"being replaced has {len(bounds)}"
        )
    seqs = [strand.sequence for strand in new_strands]
    sss = [strand.structure for strand in new_strands]
    new_sequence = replace_substrings(org_seq_struct.sequence, bounds, seqs)
    new_structure = replace_substrings(org_seq_struct.structure, bounds, sss)
    return SecStruct(new_sequence, new_structure)


def replace_gaaa_w_uucg(seq_struct) -> SecStruct:
    """
    Replace a GAAA tetraloop with a UUCG tetraloop.

    :param seq_struct: the SecStruct to modify
    :return: a new SecStruct with the tetraloop swapped
    """
    sub_seq_struct = SecStruct("GGAAAC", "(....)")
    new_seq_struct = SecStruct("CUUCGG", "(....)")
    return replace_seq_structures(seq_struct, sub_seq_struct, new_seq_struct)
