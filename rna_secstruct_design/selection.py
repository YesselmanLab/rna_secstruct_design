import yaml
from rna_secstruct.secstruct import SecStruct, MotifSearchParams
from rna_secstruct_design.util import str_to_range


def selection_from_file(filename):
    with open(filename, "r") as f:
        selection = yaml.safe_load(f)
    return selection


def get_selection(sequence, structure, params):
    ss = SecStruct(sequence, structure)
    pos = []
    for k, v in params.items():
        if k.startswith("motif"):
            pos.extend(get_selection_from_motifs(ss, v))
        elif k.startswith("flanks"):
            pos.extend(get_all_flanking_pairs(ss))
        elif k.startswith("range"):
            pos.extend([x-1 for x in str_to_range(v)])
    return pos


def get_selection_from_motifs(ss: SecStruct, params):
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

    def get_named_motif(params):
        type_name = params.pop("name", None)
        if type_name == "ref_hp":
            params["sequence"] = "CGAGUAG"
            params["structure"] = "(.....)"
        elif type_name == "gaaa_tetraloop":
            params["sequence"] = "GGAAAC"
            params["structure"] = "(....)"
        elif type_name == "tlr":
            params["sequence"] = "UAUG&CUAAG"
            params["structure"] = "(..(&)...)"

    pos = []
    extend_flank = params.pop("extend_flank", 0)
    get_named_motif(params)
    msg = MotifSearchParams(**params)
    motifs = ss.get_motifs(msg)
    for motif in motifs:
        strands = motif.strands.copy()
        if extend_flank > 0:
            strands = extend_strands(strands, extend_flank, len(ss.sequence))
        for s in strands:
            pos += s
    return pos


def get_all_flanking_pairs(ss: SecStruct):
    pos = []
    for motif in ss:
        if motif.is_helix():
            continue
        for strand in motif.strands:
            pos.extend([strand[0], strand[-1]])
    return pos
