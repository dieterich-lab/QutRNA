import click
import pandas as pd
import re
import yaml
from itertools import repeat

import os


@click.command()
@click.option("--sprinzl", required=True, help="Sequence to Sprinzl.")
@click.option("--output", required=True, help="Output FNAME")
@click.option("--linker5", default=0, help=("Length of 5' linker sequence"))
@click.argument("jacusa2", type=click.Path(exists=True))
def transform(sprinzl, output, linker5, jacusa2):
    """Add Sprinzl coordinates to JACUS2A output"""

    sprinzl = pd.read_csv(sprinzl, sep="\t")
    jacusa = pd.read_csv(jacusa2, sep="\t")

    i = sprinzl["id"].isin(jacusa["Ref"].unique())
    sprinzl = sprinzl.loc[i, ["id", "seq_pos", "sprinzl"]]

    jacusa["n_pos"] = jacusa["Pos3"] - linker5
    jacusa = (jacusa.merge(sprinzl,
                           how="left",
                           left_on=("Ref", "n_pos" ),
                           right_on=("id", "seq_pos"),
                           indicator=True)
              .drop(columns=["n_pos", "seq_pos", "id"]))
    jacusa["sprinzl"] = jacusa["sprinzl"].fillna(".")

    jacusa.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    transform()
