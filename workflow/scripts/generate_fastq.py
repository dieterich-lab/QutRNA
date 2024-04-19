import click
import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO


@click.command()
@click.option("-t", "--u2t",is_flag=True, default=False, help="Change U to T.")
@click.option("-r", "--reverse",is_flag=True, default=False, help="Reverse sequence.")
@click.option("-o", "--output", required=True,
              type=click.Path(),)
@click.argument("fastq", type=click.Path("r"))
def generate_fastq(fastq, u2t, reverse, output):
    """Generate FASTQ for parasail"""

    encoding = guess_type(fastq)[1]
    open_in = partial(gzip.open, mode="rt") if encoding == "gzip" else open
    open_out = partial(gzip.open, mode="wt") if output.endswith("gz") else open

    def do_reverse(record):
        record.id = f"{record.id}_rev"
        record.seq = record.seq[::-1]

        return record

    def do_u2t(record):
        record.seq = record.seq.replace("U", "T")

        return record

    tasks = []
    if u2t:
        tasks.append(do_u2t)
    if reverse:
        tasks.append(do_reverse)

    with open_in(fastq) as fin:
        with open_out(output) as fout:
            for record in SeqIO.parse(fin, 'fastq'):
                for task in tasks:
                    record = task(record)
                SeqIO.write(record, fout, "fastq")


if __name__ == '__main__':
    generate_fastq()
