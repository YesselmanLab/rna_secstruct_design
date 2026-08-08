"""
Driving the same functionality from the command line.

Installing the package puts a `rna-struct-design` command on your PATH. This
script shells out to it so you can see the exact commands and their output.

Run with: python examples/07_cli.py
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SEQ = "GGGGAAAACCCC"
SS = "((((....))))"


def run(args, cwd):
    """Run the cli the same way the console script does"""
    printable = "rna-struct-design " + " ".join(args)
    print(f"\n$ {printable}")
    result = subprocess.run(
        [sys.executable, "-m", "rna_secstruct_design.cli", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise SystemExit(f"command failed: {printable}")
    return result


def mut_scan(tmp):
    """Every single mutant of a sequence, folded."""
    run(["mut-scan", "-s", SEQ, "-ss", SS, "-o", "muts.csv"], tmp)
    df = pd.read_csv(tmp / "muts.csv")
    print(df.head(4).to_string(index=False))
    print(f"... {len(df)} mutants total")
    assert len(df) == 36


def mut_scan_with_params(tmp):
    """The same, but holding a motif fixed.

    A param file selects the positions to *protect*, so this leaves the hairpin
    alone and mutates only the helix arms.
    """
    (tmp / "select.yml").write_text("motif:\n  m_type: HAIRPIN\n")
    print("\nselect.yml:")
    print((tmp / "select.yml").read_text().rstrip())
    run(
        ["mut-scan", "-s", SEQ, "-ss", SS, "-pf", "select.yml", "-o", "arm_muts.csv"],
        tmp,
    )
    df = pd.read_csv(tmp / "arm_muts.csv")
    print(f"{len(df)} mutants, all keeping the hairpin intact")
    # the hairpin motif spans positions 3-8, leaving 6 mutable positions
    assert len(df) == 18
    assert all(s[3:9] == "GAAAAC" for s in df["sequence"])


def helix_rand(tmp):
    """Redesign the helices n times."""
    run(["helix-rand", "-s", SEQ, "-ss", SS, "-n", "4", "-o", "designs.csv"], tmp)
    df = pd.read_csv(tmp / "designs.csv")
    print(df.to_string(index=False))
    assert len(df) == 4
    assert all(df["structure"] == SS)


def helix_rand_from_csv(tmp):
    """Feed it a csv to redesign many constructs at once."""
    pd.DataFrame(
        [["hp_a", SEQ, SS], ["hp_b", "GGGGGAAACCCC", SS]],
        columns=["name", "sequence", "structure"],
    ).to_csv(tmp / "constructs.csv", index=False)
    run(["helix-rand", "-csv", "constructs.csv", "-n", "2", "-o", "many.csv"], tmp)
    df = pd.read_csv(tmp / "many.csv")
    print(df.to_string(index=False))
    assert list(df["name"]) == ["hp_a_1", "hp_a_2", "hp_b_1", "hp_b_2"]


def replace(tmp):
    """Swap a motif across a whole csv."""
    pd.DataFrame(
        [["hp", "GGGGGAAACCCC", SS]], columns=["name", "sequence", "structure"]
    ).to_csv(tmp / "gaaa.csv", index=False)
    (tmp / "swap.yml").write_text(
        "gaaa_to_uucg:\n"
        "  sequence: GGAAAC\n"
        "  structure: (....)\n"
        "  r_sequence: CUUCGG\n"
        "  r_structure: (....)\n"
    )
    print("\nswap.yml:")
    print((tmp / "swap.yml").read_text().rstrip())
    run(["replace", "gaaa.csv", "swap.yml", "-o", "swapped.csv"], tmp)
    df = pd.read_csv(tmp / "swapped.csv")
    print(df.to_string(index=False))
    assert df["sequence"][0] == "GGGCUUCGGCCC"


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        mut_scan(tmp)
        mut_scan_with_params(tmp)
        helix_rand(tmp)
        helix_rand_from_csv(tmp)
        replace(tmp)
    print("\nok")


if __name__ == "__main__":
    main()
