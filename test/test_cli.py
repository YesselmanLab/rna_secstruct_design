import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from rna_secstruct_design.cli import cli

# a small hairpin that folds to a single stable structure
SEQ = "GGGGAAAACCCC"
STRUCT = "((((....))))"


@pytest.fixture(name="runner")
def fixture_runner():
    return CliRunner()


def run_ok(runner, args):
    """Invoke the cli and fail loudly with the traceback if it errored"""
    result = runner.invoke(cli, args)
    assert result.exit_code == 0, result.output + str(result.exception)
    return result


class TestMutScan:
    def test_single_mutations(self, runner, tmp_path):
        out = tmp_path / "muts.csv"
        run_ok(runner, ["mut-scan", "-s", SEQ, "-ss", STRUCT, "-o", str(out)])
        df = pd.read_csv(out)
        # 3 alternative nucleotides at each of the 12 positions
        assert len(df) == 36
        assert set(df.columns) == {"name", "sequence", "structure", "ens_defect"}
        assert df["name"][0] == "G1A"
        assert all(len(s) == len(SEQ) for s in df["sequence"])

    def test_excludes_from_param_file(self, runner, tmp_path):
        # a param file selects the positions to protect, not the ones to mutate
        param_file = tmp_path / "params.yml"
        param_file.write_text(yaml.dump({"motif": {"m_type": "HAIRPIN"}}))
        out = tmp_path / "muts.csv"
        run_ok(
            runner,
            [
                "mut-scan",
                "-s",
                SEQ,
                "-ss",
                STRUCT,
                "-pf",
                str(param_file),
                "-o",
                str(out),
            ],
        )
        df = pd.read_csv(out)
        # the hairpin covers positions 3-8, so those 6 stay fixed and 6 can change
        assert len(df) == 18
        assert all(s[3:9] == "GAAAAC" for s in df["sequence"])


class TestHelixRand:
    def test_from_sequence(self, runner, tmp_path):
        out = tmp_path / "rand.csv"
        run_ok(
            runner, ["helix-rand", "-s", SEQ, "-ss", STRUCT, "-n", "3", "-o", str(out)]
        )
        df = pd.read_csv(out)
        assert len(df) == 3
        # the structure is preserved and only the helix is redesigned
        assert all(df["structure"] == STRUCT)
        assert all(s[4:8] == "AAAA" for s in df["sequence"])

    def test_from_csv(self, runner, tmp_path):
        csv_in = tmp_path / "in.csv"
        pd.DataFrame(
            [["hp", SEQ, STRUCT]], columns=["name", "sequence", "structure"]
        ).to_csv(csv_in, index=False)
        out = tmp_path / "rand.csv"
        run_ok(runner, ["helix-rand", "-csv", str(csv_in), "-n", "2", "-o", str(out)])
        df = pd.read_csv(out)
        assert len(df) == 2
        assert list(df["name"]) == ["hp_1", "hp_2"]

    def test_seq_and_csv_conflict(self, runner, tmp_path):
        csv_in = tmp_path / "in.csv"
        pd.DataFrame(
            [["hp", SEQ, STRUCT]], columns=["name", "sequence", "structure"]
        ).to_csv(csv_in, index=False)
        result = runner.invoke(
            cli, ["helix-rand", "-s", SEQ, "-csv", str(csv_in), "-o", "out.csv"]
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)


def test_replace(runner, tmp_path):
    csv_in = tmp_path / "in.csv"
    pd.DataFrame(
        [["hp", "GGGGGAAACCCC", "((((....))))"]],
        columns=["name", "sequence", "structure"],
    ).to_csv(csv_in, index=False)
    param_file = tmp_path / "params.yml"
    param_file.write_text(
        yaml.dump(
            {
                "gaaa_to_uucg": {
                    "sequence": "GGAAAC",
                    "structure": "(....)",
                    "r_sequence": "CUUCGG",
                    "r_structure": "(....)",
                }
            }
        )
    )
    out = tmp_path / "replaced.csv"
    run_ok(runner, ["replace", str(csv_in), str(param_file), "-o", str(out)])
    df = pd.read_csv(out)
    assert len(df) == 1
    assert df["sequence"][0] == "GGGCUUCGGCCC"
