import click
import pysam


def uniform_cigar(reads):
    cigar = reads[0].cigarstring
    for read in reads:
        if cigar != read.cigarstring:
            return False

    return True


@click.command()
@click.option("--keep-uniform-cigar", is_flag=True)
@click.argument("BAM", type=click.Path(exists=True))
def process(bam, keep_uniform_cigar):
  def write_reads(reads):
    for read in reads:
      out_samfile.write(read)
  # counter
  input_read_counter = 0

  # input <-> output
  in_samfile = pysam.AlignmentFile(bam, "rb")
  out_samfile = pysam.AlignmentFile("-", "wb", template=in_samfile)

  multimapper_reads = []
  for read in in_samfile:
    input_read_counter += 1

    if not multimapper_reads or multimapper_reads[0].qname == read.query_name:
      multimapper_reads.append(read)
    elif len(multimapper_reads) == 1:
      write_reads(multimapper_reads)
      multimapper_reads.clear()
      multimapper_reads.append(read)
    elif keep_uniform_cigar and uniform_cigar(multimapper_reads):
      write_reads(multimapper_reads)
      multimapper_reads.clear()
      multimapper_reads.append(read)
    else:
      multimapper_reads.clear()
      multimapper_reads.append(read)

  if keep_uniform_cigar and uniform_cigar(multimapper_reads):
    write_reads(multimapper_reads)
    multimapper_reads.clear()
  out_samfile.close()
  in_samfile.close()


if __name__ == "__main__":
  process()