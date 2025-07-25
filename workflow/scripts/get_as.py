import click
import pysam
from collections import Counter


@click.command()
@click.argument("BAM", type=click.Path())
def process(bam):
    counts = Counter()
    in_samfile = pysam.AlignmentFile(bam, "rb")

    for read in in_samfile:

        alignment_score = read.get_tag("AS")
        counts[alignment_score] += 1

    print("\t".join(["alignment_score", "count"]))
    for alignment_score, count in counts.items():
        print("\t".join([str(alignment_score), str(count)]))


if __name__ == "__main__":
    process()