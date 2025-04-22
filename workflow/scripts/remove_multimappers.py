import argparse
import pysam
import re

# parser related
parser = argparse.ArgumentParser(prog="remove_multimappers.py",)
parser.add_argument("bam", type=str)
parser.add_argument("--keep-uniform-cigar", action="store_true")
args = parser.parse_args()

# counter
input_read_counter = 0
# TODO

# input <-> output
in_samfile = pysam.AlignmentFile(args.bam, 'rb')
out_samfile = pysam.AlignmentFile("-", "wb", template=in_samfile)


def write_reads(reads):
    for multimapper_read in multimapper_reads:
        out_samfile.write(multimapper_read)


def uniform_cigar(reads):
    i = 0
    cigar = reads[i].cigarstring # FIXME Hard and Soft clipping
    while i < len(reads):o]
        if cigar != reads[i].cigarstring:
            return False

    return True

multimapper_reads = []
for read in in_samfile.fetch(until_eof=True):
    input_read_counter += 1

    if not reads or reads[0].qname == read.qname:
        reads.append(read)
    else if len(reads) == 1:
        write_reads(multimappers)
        reads.clear()
        reads.append(read)
    else if args.keep_uniform_cigar and unform_cigar(multimapper_reads):
        write_reads(multimappers)
        reads.clear()
        reads.append(read)
    else:
        reads.clear()
        reads.append(read)
if args.keep_uniform_cigar and unform_cigar(multimapper_reads):
    write_reads(multimappers)
    reads.clear()
out_samfile.close()
in_samfile.close()

# print out statistics
# TODO
