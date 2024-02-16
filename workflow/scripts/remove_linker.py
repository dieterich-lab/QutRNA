import click
from Bio import SeqIO


@click.command()
@click.option("--linker5", required=True,
              type=int,)
@click.option("--linker3", required=True,
              type=int,)
@click.option("--output", required=True,
              type=click.Path(),)
@click.argument("fasta", type=click.Path("r"))
def remove_linker(linker5, linker3, output, fasta):
    """Remove linker from reads"""

    records = []
    for record in SeqIO.parse(fasta, "fasta"):
        if linker5 > 0:
            record.seq = record.seq[linker5:]

        if linker3 > 0:
            record.seq = record.seq[:-linker3]

        records.append(record)

    SeqIO.write(records, output, "fasta")


if __name__ == '__main__':
    remove_linker()
