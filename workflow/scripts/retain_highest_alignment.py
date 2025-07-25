import click
import pysam
import sys


def write_reads(reads, out_samfile):
  max_as = max(reads.keys())
  for read in reads[max_as]:
    out_samfile.write(read)


@click.command()
@click.option("--min-as", type=int)
@click.argument("bam", type=click.Path(exists=True))
def process(min_as, bam):
  # counter
  input_read_counter = 0
  in_samfile = pysam.AlignmentFile(bam, "rb")
  out_samfile = pysam.AlignmentFile(sys.stdout, "wb", template=in_samfile)

  reads = {}
  qname = ""
  for read in in_samfile:
    input_read_counter += 1

    alignment_score = int(read.get_tag("AS"))
    if min_as and alignment_score < min_as:
        continue

    if not qname or qname == read.query_name:
        reads.setdefault(alignment_score, [].append(read))
        qname = read.query_name
    else:
        write_reads(reads, out_samfile)
        reads.clear()
        reads.setdefault(alignment_score, [].append(read))
        qname = read.query_name
  if reads:
    write_reads(reads, out_samfile)
  out_samfile.close()
  in_samfile.close()


if __name__ == '__main__':
  process()
