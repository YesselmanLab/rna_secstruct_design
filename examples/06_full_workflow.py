"""
An end to end library design, the way you would actually use this package.

Take a construct, protect the functional parts, redesign everything else into a
set of sequence variants, quality check them and write a csv.

Run with: python examples/06_full_workflow.py
"""

import random
import tempfile
from pathlib import Path

import pandas as pd
from rna_secstruct import SecStruct
from vienna import fold

from rna_secstruct_design.dataframe import (
    get_max_gc_stretch,
    get_max_repeating_nucleotides,
    has_gc_streches_less_than,
    has_repeating_nucleotides_less_than,
)
from rna_secstruct_design.helix_randomizer import HelixRandomizer
from rna_secstruct_design.selection import get_selection

# a GAAA tetraloop docked into its receptor, the two motifs must stay intact
SEQ = "AAGAUAUGGCGUGGAGUACCGGAAACGGUGACUUUGAUACUGCCUAAGUCUU"
SS = "(((((..((((((((((((((....)))).))))...))).)))...)))))"
N_DESIGNS = 6


def build_library():
    """Redesign the helices N times, keeping the two functional motifs fixed."""
    print("--- designing ---")
    print(f"scaffold  {SEQ}")
    print(f"structure {SS}")

    secstruct = SecStruct(SEQ, SS)
    exclude = get_selection(
        secstruct,
        {
            "seq_struct_1": {"name": "tlr_extended"},
            "motif_2": {"name": "gaaa_tetraloop", "extend_flank": 1},
        },
    )
    print(f"protecting {len(exclude)} positions: {sorted(exclude)}")

    hr = HelixRandomizer()
    rows = []
    for i in range(N_DESIGNS):
        ens_defect, seq = hr.run(secstruct, exclude)
        # an empty sequence means the randomizer gave up before finding a design
        assert seq, f"helix randomizer failed on design {i + 1}"
        rows.append(
            {
                "name": f"design_{i + 1}",
                "sequence": seq,
                "structure": SS,
                "ens_defect": ens_defect,
            }
        )
    df = pd.DataFrame(rows)
    print(f"built {len(df)} designs")
    return df, exclude


def check_library(df, exclude):
    """Every design has to fold correctly and keep the protected positions."""
    print("\n--- checking ---")
    for _, row in df.iterrows():
        assert len(row["sequence"]) == len(SEQ)
        assert fold(row["sequence"]).dot_bracket == SS, f"{row['name']} misfolds"
        for i in exclude:
            assert row["sequence"][i] == SEQ[i], f"{row['name']} changed position {i}"
    print(f"all {len(df)} designs fold into the target structure")
    print(f"all {len(df)} designs kept their protected motifs")

    # the designs should differ from each other, otherwise the library is useless
    assert df["sequence"].nunique() == len(df), "designs are not unique"
    print(f"all {len(df)} designs are distinct")


def quality_control(df):
    """The dataframe helpers score a whole library at once."""
    print("\n--- quality control ---")

    df = df.copy()
    df["max_gc_stretch"] = get_max_gc_stretch(df)
    repeats = get_max_repeating_nucleotides(df)
    df = pd.concat([df.reset_index(drop=True), repeats], axis=1)

    cols = ["name", "ens_defect", "max_gc_stretch", "max_A_stretch", "max_G_stretch"]
    print(df[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # these return a single bool for the whole library, not per row
    print(f"\nno GC stretch of 5 or more:      {has_gc_streches_less_than(df, 5)}")
    print(
        f"no nucleotide repeated 6 times:  {has_repeating_nucleotides_less_than(df, 6)}"
    )
    assert has_gc_streches_less_than(df, 5)
    assert has_repeating_nucleotides_less_than(df, 6)
    return df


def write_csv(df):
    """What you would hand off to an oligo order."""
    print("\n--- writing ---")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "library.csv"
        df[["name", "sequence", "structure", "ens_defect"]].to_csv(out, index=False)
        written = pd.read_csv(out)
        print(f"wrote {len(written)} rows to {out.name}")
        assert len(written) == N_DESIGNS


def main():
    random.seed(0)
    df, exclude = build_library()
    check_library(df, exclude)
    df = quality_control(df)
    write_csv(df)
    print("\nok")


if __name__ == "__main__":
    main()
