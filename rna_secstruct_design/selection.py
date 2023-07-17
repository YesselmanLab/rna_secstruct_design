import yaml
from seq_tools.structure import SequenceStructure, find
from rna_secstruct.secstruct import SecStruct, MotifSearchParams
from rna_secstruct_design.util import str_to_range


def flatten(l):
    """Recursively flatten a list of lists of integers."""
    flattened = []
    for item in l:
        if isinstance(item, int):
            flattened.append(item)
        else:
            flattened.extend(flatten(item))
    return flattened


def selection_from_file(filename):
    with open(filename, "r") as f:
        selection = yaml.safe_load(f)
    return selection


def get_selection(secstruct, params):
    pos = []
    for k, v in params.items():
        if k.startswith("motif"):
            pos.extend(get_selection_from_motifs(secstruct, v))
        elif k.startswith("seq_struct"):
            pos.extend(get_seq_struct(secstruct, v))
        elif k.startswith("flanks"):
            pos.extend(get_all_flanking_pairs(secstruct))
        elif k.startswith("range"):
            pos.extend([x - 1 for x in str_to_range(v)])
    if "invert" in params:
        pos = invert_exclude_list(pos, len(secstruct.sequence))
    return pos


def get_named_motif(params):
    type_name = params.pop("name", None)
    if type_name is None:
        return
    if type_name == "ref_hp":
        params["sequence"] = "CGAGUAG"
        params["structure"] = "(.....)"
    elif type_name == "gaaa_tetraloop":
        params["sequence"] = "GGAAAC"
        params["structure"] = "(....)"
    elif type_name == "tlr":
        params["sequence"] = "UAUG&CUAAG"
        params["structure"] = "(..(&)...)"
    elif type_name == "tlr_extended":
        params["sequence"] = "AUAUGG&CCUAAGU"
        params["structure"] = "((..((&))...))"
    else:
        raise ValueError(f"Unknown motif name: {type_name}")


def get_selection_from_motifs(secstruct: SecStruct, params):
    def extend_strands(strands, extend, seq_len):
        """
        Extend strands by a given number of positions
        :param strands: strands from a motif.strands
        :param extend: number of pos to extend
        :param max_pos: the number of nucleotides in the sequence
        :return: strands with extended flanks
        """
        new_strands = []
        for s in strands:
            min_val, max_val = s[0], s[-1]
            r1 = list(range(min_val - extend, min_val))
            r2 = list(range(max_val + 1, max_val + extend + 1))
            new_strand = r1 + s + r2
            new_strand_filtered = [x for x in new_strand if x < seq_len and x >= 0]
            new_strands.append(new_strand_filtered)
        return new_strands

    pos = []
    extend_flank = params.pop("extend_flank", 0)
    get_named_motif(params)
    msg = MotifSearchParams(**params)
    motifs = secstruct.get_motifs(msg)
    for motif in motifs:
        strands = motif.strands.copy()
        if extend_flank > 0:
            strands = extend_strands(strands, extend_flank, len(secstruct.sequence))
        for s in strands:
            pos += s
    return pos


def get_seq_struct(secstruct: SecStruct, v):
    seq = secstruct.sequence
    struct = secstruct.structure
    full = SequenceStructure(seq, struct)
    if "name" in v:
        get_named_motif(v)
    sub = SequenceStructure(v["sequence"], v["structure"])
    bounds = find(full, sub)[0]
    pos = []
    for r in bounds:
        pos.extend(list(range(r[0], r[1])))
    return pos


# TODO helix after or before single strand count as flank?
def get_all_flanking_pairs(secstruct: SecStruct):
    pos = []
    for motif in secstruct:
        if motif.is_helix():
            continue
        for strand in motif.strands:
            pos.extend([strand[0], strand[-1]])
    return pos


def invert_exclude_list(exclude, seq_len):
    allowed = list(range(seq_len))
    for e in exclude:
        allowed.remove(e)
    return allowed
