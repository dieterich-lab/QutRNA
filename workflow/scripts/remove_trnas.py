import click

from Bio import SeqIO


@click.command()
@click.option("--output", type=click.Path(exists=False))
@click.option("--trna", "-t", type=str, multiple=True)
@click.argument("fasta", type=click.Path(exists=True))
def filter(fasta, output, trna):
    filtered_seqs = []
    records = SeqIO.parse(fasta, "fasta")

    if trna:
        filtered_seqs = [record for record in records if record.id not in trna]
    else:
        filtered_seqs = records

    SeqIO.write(filtered_seqs, output, "fasta")


if __name__ == '__main__':
    filter()
