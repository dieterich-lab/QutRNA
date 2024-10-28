from os.path import basename

import pandas as pd
import click
import re
import sys


@click.command()
@click.option("-o", "--output", type=str, required=True)
@click.argument("coverages", nargs=-1)
def merge(output, coverages):
    dfs = list()
    for fname in coverages:
        df = pd.read_csv(fname, sep="\t")
        df = df.rename(columns={"#rname": "rname"})
        df["fname"] = basename(fname)
        df = df[["rname", "numreads", "fname"]]
        df["sample"] = df["fname"].str.replace("\_coverage.txt", "", regex=True)
        dfs.append(df)
    df = pd.concat(dfs)
    df.to_csv(output, sep="\t", index=False, quoting=False)


if __name__ == '__main__':
    merge()
