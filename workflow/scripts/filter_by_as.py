import click
import pysam
import sys


@click.command()
@click.option("--min-alignment-score", required=True, type=int)
@click.argument("bam", type=click.Path(exists=True))
def process(min_alignment_score, bam):
  in_samfile = pysam.AlignmentFile(bam, "rb")
  out_samfile = pysam.AlignmentFile(sys.stdout, "wb", template=in_samfile)

  for read in in_samfile:
    alignment_score = int(read.get_tag("AS"))
    if alignment_score >= min_alignment_score:
      out_samfile.write(read)
  out_samfile.close()
  in_samfile.close()


if __name__ == '__main__':
    process()
