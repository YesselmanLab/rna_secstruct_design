import click
import pandas as pd
import numpy as np
from multiprocessing import Pool

from vienna import fold
from rna_secstruct.secstruct import SecStruct
from rna_secstruct_design.selection import selection_from_file, get_selection
from rna_secstruct_design.logger import setup_applevel_logger, get_logger
from rna_secstruct_design.mutations import find_multiple_mutations
from rna_secstruct_design.helix_randomizer import HelixRandomizer

log = get_logger("CLI")


def validate_dataframe(df) -> None:
    """
    validates a dataframe to have a column named `sequence` and `name`
    :param df: dataframe with sequences
    :return: None
    """
    if "sequence" not in df.columns:
        raise ValueError("sequence column not found")
    if "structure" not in df.columns:
        raise ValueError("structure column not found")
    if "name" not in df.columns:
        df["name"] = [f"seq_{i}" for i in range(len(df))]


def get_input_dataframe(
        seq,
        struct,
        csv_file,
) -> pd.DataFrame:
    """
    returns a dataframe from a sequence or a file
    :param data: can be a seqeunce or a file
    :return: pd.DataFrame
    """
    if csv_file is not None and seq is not None:
        raise ValueError("cannot specify both a sequence and a csv file")
    elif csv_file is not None:
        log.info(f"reading file {csv_file}")
        df = pd.read_csv(csv_file)
        log.info(f"csv file contains {len(df)} sequences")
    else:
        log.info(f"reading sequence {seq}")
        if struct is None:
            struct = fold(seq).dot_bracket
        data_df = [["seq", seq, struct]]
        df = pd.DataFrame(data_df, columns=["name", "sequence", "structure"])
    validate_dataframe(df)
    return df


def randomize_helices(df, params, num_seqs):
    data = []
    hr = HelixRandomizer()
    for _, row in df.iterrows():
        secstruct = SecStruct(row["sequence"], row["structure"])
        exclude = get_selection(secstruct, params)
        log.info(row["name"])
        print(row['name'])
        for i in range(num_seqs):
            ens_defect, seq = hr.run(secstruct, exclude)
            data.append(
                    {
                        "name"      : row["name"] + "_" + str(i + 1),
                        "sequence"  : seq,
                        "structure" : secstruct.structure,
                        "ens_defect": ens_defect,
                    }
            )
    return pd.DataFrame(data)


def fold_sequences(df):
    data = []
    count = 0
    for _, row in df.iterrows():
        rf = fold(row["sequence"])
        if count % 10 == 0:
            print(count)
        data.append(
                {
                    "name"      : row["name"],
                    "sequence"  : row["sequence"],
                    "structure" : rf.dot_bracket,
                    "ens_defect": rf.ens_defect,
                }
        )
    return pd.DataFrame(data)


# multi commmand format
@click.group()
def cli():
    pass


@cli.command()
@click.option("-s", "--seq", type=str, required=True)
@click.option("-ss", "--struct", type=str, default=None)
@click.option("-n", "--num-muts", type=int, default=1)
@click.option("-pf", "--param-file", type=click.Path(exists=True), default=None)
@click.option("-o", "--output", type=click.Path(exists=False), default="output.csv")
@click.option("-p", "--num-processes", type=int, default=1)
def mut_scan(seq, struct, num_muts, param_file, num_processes, output):
    setup_applevel_logger()
    if num_processes > 1:
        log.info(f"running with multiprocess! {num_processes} processes")
    if struct is None:
        struct = fold(seq).dot_bracket
    secstruct = SecStruct(seq, struct)
    if param_file is not None:
        params = selection_from_file(param_file)
        exclude = get_selection(secstruct, params)
    else:
        exclude = []
    results = find_multiple_mutations(secstruct.sequence, num_muts, exclude)
    if num_processes > 1:
        with Pool(num_processes) as p:
            results = p.starmap_async(fold_sequences,
                                [results_s for results_s in np.array_split(results, num_processes)]
                                )
            df = pd.concat(results)
    else:
        df = fold_sequences(results)
    df.to_csv(output, index=False)


@cli.command()
@click.option("-s", "--seq", type=str, required=False)
@click.option("-ss", "--struct", type=str, default=False)
@click.option("-csv", "--csv-file", type=click.Path(exists=True), default=None)
@click.option("-pf", "--param-file", type=click.Path(exists=True), default=None)
@click.option("-o", "--output", type=click.Path(exists=False), default="output.csv")
@click.option("-n", "--num-seqs", type=int, default=10)
@click.option("-p", "--num-processes", type=int, default=1)
def helix_rand(seq, struct, csv_file, param_file, num_seqs, num_processes, output):
    setup_applevel_logger()
    df = get_input_dataframe(seq, struct, csv_file)
    if param_file is not None:
        params = selection_from_file(param_file)
    else:
        params = {}
    if num_processes > 1:
        with Pool(num_processes) as p:
            results = p.starmap_async(
                    randomize_helices,
                    [
                        (df_s, params, num_seqs)
                        for df_s in np.array_split(df, num_processes)
                    ],
            )
        df = pd.concat(results)
    else:
        df = randomize_helices(df, params, num_seqs)
    # df = pd.DataFrame(data)
    df.to_csv(output, index=False)


if __name__ == "__main__":
    cli()
