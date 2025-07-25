from pickle import FALSE

import click
import pandas as pd


@click.command()
@click.option("--output", "-o", type=click.Path(), required=True)
@click.option("--read-types", "-r", type=str, multiple=True, required=True)
@click.option("--conditions", "-c", type=str, multiple=True, required=True)
@click.option("--samples", "-s", type=str, multiple=True, required=True)
@click.option("--subsamples", "-u", type=str, multiple=True, required=True)
@click.option("--base-calls", "-b", type=str, multiple=True, required=True)
@click.argument("FNAMES", type=click.Path(), nargs=-1)
def process(output, read_types, conditions, samples, subsamples, base_calls, fnames):
  dfs = []
  for read_type, sample, subsample, condition, base_call, fname in zip(read_types, samples, subsamples, conditions, base_calls, fnames):
    df = pd.read_csv(fname, sep="\t")

    df["read_type"] = read_type
    df["sample"] = sample
    df["subsample"] = subsample
    df["base_calling"] = base_call
    df["condition"] = condition
    df["fname"] = fname

    dfs.append(df)
  df = pd.concat(dfs, ignore_index=True)
  df.to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
  process()
