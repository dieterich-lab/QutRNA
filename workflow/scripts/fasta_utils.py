import click
from contextlib import contextmanager
import enum
import gzip
import re
import pandas as pd
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO


class BaseChange(enum.Enum):
    U2T = enum.auto()
    T2U = enum.auto()


@click.group()
def cli():
    pass


def read_records(fname):
    encoding = guess_type(fname)[1]
    open_in = partial(gzip.open, mode="rt") if encoding == "gzip" else partial(open, mode="r")
    with open_in(fname) as fin:
        records = {record.name: record for record in SeqIO.parse(fin, "fasta")}

    return records


def open_out(output):
    f = partial(gzip.open, mode="wt") if output.endswith("gz") else partial(open, mode="w")

    return f(output)


@cli.command()
@click.option("-o", "--output", type=click.Path())
@click.argument("FASTA", type=click.Path(exists=True))
def extract_seqids(fasta, output):
    """Extract tRNAs from fasta"""

    records = read_records(fasta)
    with open_out(output) as fout:
        for record in records:
            fout.write(record)


def do_remove_linker(record, linker5, linker3):
    if linker5 > 0:
        record.seq = record.seq[linker5:]

    if linker3 > 0:
        record.seq = record.seq[:-linker3]

    return record


def do_base_change(record, base1, base2):
    record.seq = record.seq.replace(base1, base2)

    return record


def do_reverse(record):
    record.id = f"{record.id}_rev"
    record.seq = record.seq[::-1]

    return record


@cli.command()
@click.option("--remove-linker5", type=int)
@click.option("--remove-linker3", type=int)
@click.option("-c", "--base-change", type = click.Choice(BaseChange, case_sensitive=False), help="Base change: U2T or T2U.")
@click.option("-i", "--ignore", type=str, multiple=True)
@click.option("-r", "--reverse", is_flag=True, default=False, help="Reverse sequence.")
@click.option("-o", "--output", required=True, type=click.Path())
@click.argument("FASTA", type=click.Path())
def transform(fasta, remove_linker5, remove_linker3, base_change, ignore, reverse, output):
    """Transform fasta"""

    tasks = []
    if remove_linker5 or remove_linker3:
        def helper(record):
            return do_remove_linker(record, remove_linker5, remove_linker3)
        tasks.append(helper)

    if base_change:
        bc = str(base_change).split("2")
        def helper(record):
            return do_base_change(record, bc[0], bc[1])
        tasks.append(helper)

    if reverse:
        tasks.append(do_reverse)

    records = read_records(fasta)
    if ignore:
        records = ignore_trnas(records, ignore)

    with open_out(output) as fout:
        for record in records.values():
            for task in tasks:
                record = task(record)
            SeqIO.write(record, fout, "fasta")


def ignore_trnas(records, ignore):
    """ Remove trnas from fasta"""

    filtered = dict(records)

    if len(ignore) == 1 and ignore[0].startswith("/"):
        pattern = re.compile(rf"{ignore[0]}")
        trnas = [trna for trna in filtered if pattern.search(trna)]
    else:
        trnas = list(ignore)

    # remove trnas such as defined by ignore
    for trna in trnas:
        del filtered[trna]

    if not filtered:
        raise Exception("No tRNAs left!")

    if ignore and len(filtered) == len(records):
        raise Exception("No tRNAs were filtered. Check pep file!")

    return filtered

@cli.command()
@click.option("-o", "--output", type=click.Path())
@click.argument("FASTA", type=click.Path(exists=True))
def infer_annotation(fasta, output):
    """Infer tRNA annotation from fasta"""

    trnas = read_records(fasta)

    df = pd.DataFrame.from_dict({"trna": trnas.keys()})
    df["is_mt"] = df["trna"].str.startswith("MT")
    df["is_nuclear"] = ~df["is_mt"]
    info = df["trna"].str.extract(r'.*tRNA-(?P<amino_acid>[^-]+)-(?P<anti_codon>[A-Z]+)')
    df = pd.concat([df, info], axis=1)
    df["isotype"] = df["amino_acid"]
    df["isoacceptor"] = df["amino_acid"].str.cat(df["anti_codon"], sep="-")
    # TODO df["isodecoder"]

    df.to_csv(output, index=False, sep="\t")

if __name__ == "__main__":
    cli()