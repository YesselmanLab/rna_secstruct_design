import pandas as pd

from rna_secstruct_design.dataframe import (
    get_max_gc_stretch,
    get_max_repeating_nucleotides,
    has_gc_streches_less_than,
    has_repeating_nucleotides_less_than,
)

# generate test data ################################################################


def get_test_data_rna() -> pd.DataFrame:
    """
    get test data for dna
    :return: pd.DataFrame
    """
    return pd.DataFrame(
        [
            ["seq_0", "GGGGUUUUCCCC", "((((....))))"],
            ["seq_1", "GAGGUAUUCCUC", "((((....))))"],
        ],
        columns=["name", "sequence", "structure"],
    )


# tests  ##########################################################################


def test_get_max_gc_stretch():
    """
    test get_max_gc_stretch
    """
    df = get_test_data_rna()
    gc_stretches = get_max_gc_stretch(df)
    assert gc_stretches[0] == 4
    assert gc_stretches[1] == 2


def test_get_max_repeating_nucleotides():
    """
    test get_max_repeating_nucleotides
    """
    df = get_test_data_rna()
    df_nuc_repeats = get_max_repeating_nucleotides(df)
    assert len(df_nuc_repeats) == 2
    assert df_nuc_repeats["max_A_stretch"][0] == 0


def test_has_gc_streches_less_than():
    """
    test has_gc_streches_less_than
    """
    df = get_test_data_rna()
    assert has_gc_streches_less_than(df, 5)
    assert not has_gc_streches_less_than(df, 4)


def test_has_repeating_nucleotides_less_than():
    """
    test has_repeating_nucleotides_less_than
    """
    df = get_test_data_rna()
    assert has_repeating_nucleotides_less_than(df, 5)
    assert not has_repeating_nucleotides_less_than(df, 4)
