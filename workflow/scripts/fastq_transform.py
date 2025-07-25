import click
import enum
import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO


class BaseChange(enum.Enum):
    U2T = enum.auto()
    T2U = enum.auto()


@click.command()
@click.option("-c", "--base-change", type = click.Choice(BaseChange, case_sensitive=False), help="Base change: U2T or T2U.")
@click.option("-r", "--reverse",is_flag=True, default=False, help="Reverse sequence.")
@click.option("-o", "--output", required=True, type=click.Path())
@click.argument("fastq", type=click.Path())
def process(fastq, base_change, reverse, output):
    """Generate FASTQ for parasail"""

    encoding = guess_type(fastq)[1]
    open_in = partial(gzip.open, mode="rt") if encoding == "gzip" else partial(open, mode="r")
    open_out = partial(gzip.open, mode="wt") if output.endswith("gz") else partial(open, mode="w")

    def do_reverse(record):
        record.id = f"{record.id}_rev"
        record.seq = record.seq[::-1]

        return record

    def do_base_change(base1, base2):
        def helper(record):
            record.seq = record.seq.replace(base1, base2)

            return record

        return helper

    tasks = []
    if base_change:
        bc = base_change.name.split("2")
        tasks.append(do_base_change(bc[0], bc[1]))
    if reverse:
        tasks.append(do_reverse)

    with open_in(fastq) as fin:
        with open_out(output) as fout:
            for record in SeqIO.parse(fin, 'fastq'):
                for task in tasks:
                    record = task(record)
                SeqIO.write(record, fout, "fastq")


if __name__ == '__main__':
    process()
