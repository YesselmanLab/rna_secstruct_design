import pandas as pd
import numpy as np

from rna_secstruct_design.util import max_repeating_nucleotides, max_gc_stretch


def get_max_gc_stretch(df):
    """Get the max gc stretch for each sequence and structure in a dataframe"""
    # apply max_gc_stretch to each row in the dataframe uses sequence and structure columns
    return df.apply(
        lambda row: max_gc_stretch(row["sequence"], row["structure"]), axis=1
    )


def get_max_repeating_nucleotides(df):
    """
    Get the max repeating nucleotides for each sequence and structure in a dataframe
    """
    data = []
    for _, row in df.iterrows():
        r = max_repeating_nucleotides(row["sequence"])
        data.append([r["A"], r["C"], r["G"], r["U"]])
    return pd.DataFrame(
        data,
        columns=["max_A_stretch", "max_C_stretch", "max_G_stretch", "max_U_stretch"],
    )


def has_gc_streches_less_than(df, max_stretch):
    """
    Get a boolean series of whether each sequence in a dataframe has a gc stretch
    less than max_stretch"""
    return np.max(get_max_gc_stretch(df)) < max_stretch


def has_repeating_nucleotides_less_than(df, max_repeats):
    """
    Get a boolean of whether each sequence in a dataframe has a repeating
    nucleotide less than max_repeats
    """
    df_nuc_repeats = get_max_repeating_nucleotides(df)
    for c in ["max_A_stretch", "max_C_stretch", "max_G_stretch", "max_U_stretch"]:
        if df_nuc_repeats[c].max() >= max_repeats:
            return False
    return True
