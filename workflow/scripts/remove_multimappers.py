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
    for read in reads:
        out_samfile.write(read)


def uniform_cigar(reads):
    cigar = reads[0].cigarstring # FIXME Hard and Soft clipping
    for read in reads:
        if cigar != read.cigarstring:
            return False

    return True

multimapper_reads = []
for read in in_samfile.fetch(until_eof=True):
    input_read_counter += 1

    if not multimapper_reads or multimapper_reads[0].qname == read.qname:
        multimapper_reads.append(read)
    elif len(multimapper_reads) == 1:
        write_reads(multimapper_reads)
        multimapper_reads.clear()
        multimapper_reads.append(read)
    elif args.keep_uniform_cigar and uniform_cigar(multimapper_reads):
        write_reads(multimapper_reads)
        multimapper_reads.clear()
        multimapper_reads.append(read)
    else:
        multimapper_reads.clear()
        multimapper_reads.append(read)
if args.keep_uniform_cigar and uniform_cigar(multimapper_reads):
    write_reads(multimapper_reads)
    multimapper_reads.clear()
out_samfile.close()
in_samfile.close()

# print out statistics
# TODO
