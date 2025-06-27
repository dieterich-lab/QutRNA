import argparse
import pysam

# parser related
parser = argparse.ArgumentParser(prog="get_as.py",)
parser.add_argument("bam", type=str)
args = parser.parse_args()

# input <-> output
in_samfile = pysam.AlignmentFile(args.bam, 'rb')

for read in in_samfile.fetch(until_eof=True):
    alignment_score = get_tag("AS")
    print(f"{alignment_score}\n")
