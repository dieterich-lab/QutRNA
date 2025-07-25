from os.path import exists

import pysam
import click


@click.command()
@click.argument("BAM", type=click.Path(exists=True))
def process(bam):
  # counter
  input_read_counter = 0

  # input <-> output
  in_samfile = pysam.AlignmentFile(bam, "rb")
  out_samfile = pysam.AlignmentFile("-", "wb", template=in_samfile)

  def write_reads(reads):
      n = len(reads)
      for read in reads:
          read.tags += [('NH', n)]
          out_samfile.write(read)

  multimapper_reads = []
  for read in in_samfile:
    input_read_counter += 1

    if not multimapper_reads or multimapper_reads[0].qname == read.query_name:
      multimapper_reads.append(read)
    else:
      write_reads(multimapper_reads)
      multimapper_reads.clear()
      multimapper_reads.append(read)

  if multimapper_reads:
    write_reads(multimapper_reads)
    multimapper_reads.clear()
  out_samfile.close()
  in_samfile.close()


if __name__ == "__main__":
  process()
