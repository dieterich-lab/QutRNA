import click
import pandas as pd


@click.command()
@click.option("--output", "-o", type=click.Path(), required=True)
@click.option("--data", "-d", type=(str, str, str, str, str), multiple=True, required=True)
@click.argument("FNAMES", type=click.Path(), nargs=-1)
def process(fnames, data, output):
  dfs = []
  for (read_type, condition, sample, subsample, base_call), fname in zip(data, fnames):
    df = pd.read_csv(fname, sep="\t", dtype=str)

    df["read_type"] = read_type
    df["sample"] = sample
    df["subsample"] = subsample
    df["base_calling"] = base_call
    df["condition"] = condition
    df["fname"] = fname

    dfs.append(df)
  df = pd.concat(dfs, ignore_index=True)
  df.to_csv(output, sep="\t", index=False, quoting=False)


if __name__ == '__main__':
  process()
