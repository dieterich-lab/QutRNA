import click
import pandas as pd
import re
import yaml
from itertools import repeat

import os


@click.command()
@click.option("--numbering", required=True, help="Sequence to universal numbering.")
@click.option("--output", required=True, help="Output FNAME")
@click.option("--aln_prefix", show_default=True, default="Homo_sapiens_", help="Alignment prefix")
@click.option("--linker5", default=0, help=("Length of 5' linker sequence"))
@click.argument("jacusa2", type=click.Path(exists=True))
def transform(numbering, output, linker5, aln_prefix, jacusa2):
    """Transform JACUS2A output to universal conventional tRNAposition"""

    numbering = pd.read_csv(numbering, sep="\t")

    jacusa = pd.read_csv(jacusa2, sep=",")

    i = numbering["ref"].isin(jacusa["Ref"].unique())
    numbering = numbering.loc[i]

    jacusa["n_pos"] = jacusa["Pos3"] - linker5 - 1
    jacusa["aa"] = jacusa["Ref"].apply(_clean_name, prefix=aln_prefix)
    jacusa = (
            jacusa.merge(numbering,
                          how="left", #
                          left_on=("Ref", "n_pos" ),
                          right_on=("ref", "seq_pos"),
                          indicator=True)
                  .drop(columns=["n_pos", "seq_pos", "ref"]) )
    jacusa["u_pos"] = jacusa["u_pos"].fillna(".")

    jacusa.to_csv(output, sep=",", index=False)


def _clean_name(name, prefix=""):
    if prefix:
        pat = "^" + prefix + "tRNA-"
        name = re.sub(pat, "", name)

    return re.sub(r"-.+$", "", name)


if __name__ == "__main__":
    transform()
