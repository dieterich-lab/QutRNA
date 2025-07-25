import click
import pysam


def _trim_cigar(cigar, read):
    aln_offset = 0
    for i, (op, l) in enumerate(cigar):
        if op in (pysam.CEQUAL, pysam.CMATCH, pysam.CHARD_CLIP, pysam.CSOFT_CLIP):
            break

        match op:
            case pysam.CBACK:
                raise Exception(read)
            case pysam.CDIFF:
                cigar[i] = (pysam.CSOFT_CLIP, l)
                aln_offset += l
            case pysam.CREF_SKIP:
                raise Exception(read)
            case pysam.CDEL:
                cigar[i] = ()
                aln_offset += l
            case pysam.CINS:
                cigar[i] = (pysam.CSOFT_CLIP, l)
                #aln_offset += l
            case pysam.CPAD:
                raise Exception(read)

    return cigar, aln_offset


def trim(read, _):

    new_cigar, aln_offset = _trim_cigar(list(read.cigartuples), read)
    new_cigar, _ = _trim_cigar(list(reversed(new_cigar)), read)
    new_cigar = tuple(c for c in reversed(new_cigar) if c)

    # merge consecutive OPs, e.g.:10S10S
    merged_cigar = []
    current = [new_cigar[0][0], new_cigar[0][1]]
    for (op, l) in new_cigar[1:]:
        if current[0] == op:
            current[1] += l
        else:
            merged_cigar.append(current)
            current = [op, l]
    merged_cigar.append(current)
    new_cigar = tuple((c[0], c[1]) for c in merged_cigar)
    if new_cigar != read.cigartuples:
        read.cigartuples = new_cigar
        if aln_offset:
            read.reference_start += aln_offset

        return True

    return False


# parser related
@click.command()
@click.option("--trim-cigar", is_flag=True)
@click.option("--min-read-length", type=int)
@click.option("--max-read-length", type=int)
@click.option("--min-alignment-length", type=int)
@click.option("--max-alignment-length", type=int)
@click.option("--stats", type=str)
@click.argument("bam", type=str)
def process(trim_cigar, min_read_length, max_read_length, min_alignment_length, max_alignment_length, stats, bam):
  # input <-> output
  in_samfile = pysam.AlignmentFile(bam, "rb")
  out_samfile = pysam.AlignmentFile("-", "wb", template=in_samfile)

  # counter
  input_read_counter = 0
  trim_cigar_counter = 0
  read_length_counter = 0
  aln_length_counter = 0
  output_read_counter = 0
  for read in in_samfile:
    input_read_counter += 1
    if min_read_length or max_read_length:
      read_len = len(read.query_sequence)
      if min_read_length and min_read_length > read_len :
        read_length_counter += 1
        continue
      if max_read_length and max_read_length < read_len:
        read_length_counter += 1
        continue
    if min_alignment_length or max_alignment_length:
      aln_len = read.query_alignment_length
      if min_alignment_length and min_alignment_length > aln_len :
        aln_length_counter += 1
        continue
      if max_alignment_length and max_alignment_length < aln_len:
        aln_length_counter += 1
        continue

    if trim_cigar:
      if trim(read, out_samfile):
        trim_cigar_counter += 1
      out_samfile.write(read)
      output_read_counter += 1
    else:
      out_samfile.write(read)
      output_read_counter += 1
  out_samfile.close()
  in_samfile.close()

  # print out statistics
  if stats:
    with open(stats, "w") as f:
      f.write(f"input_reads\t{input_read_counter}\n")
      f.write(f"trimmed_cigars\t{trim_cigar_counter}\n")
      f.write(f"read_length_filtered\t{read_length_counter}\n")
      f.write(f"alignment_length_filtered\t{aln_length_counter}\n")
      f.write(f"output_reads\t{output_read_counter}\n")


if __name__ == '__main__':
    process()
