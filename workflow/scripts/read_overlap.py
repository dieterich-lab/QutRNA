import click
import gzip
import pysam
import sys

from Bio import SeqIO


def get_overlap(start1, end1, start2, end2):
    start = max(start1, start2)
    end = min(end1, end2)
    return max(end - start + 1, 0)


@click.command()
@click.option("--fasta", required=True, type=click.Path(exists=True))
@click.option("--five-adapter", type=int, default=0)
@click.option("--stats", required=True, type=click.Path())
@click.option("--three-adapter", type=int, default=0)
@click.option("--five-adapter-overlap", type=float, default=0.0)
@click.option("--trna-overlap", type=float, default=0.0)
@click.option("--three-adapter-overlap", type=float, default=0.0)
@click.argument("bam", type=click.Path(exists=True))
def process(fasta,
            five_adapter, three_adapter,
            five_adapter_overlap, trna_overlap, three_adapter_overlap,
            stats,
            bam):
    trnas = {trna.name: trna for trna in SeqIO.parse(fasta, "fasta")}

    in_samfile = pysam.AlignmentFile(bam, "rb")
    out_samfile = pysam.AlignmentFile(sys.stdout, "wb", template=in_samfile)

    line = ["read_id", "ref", "aln_start", "aln_end", "aln_score", "cigar", "ov_5", "ov_trna", "ov_3"]
    stats_file = gzip.open(stats, "wb")
    if five_adapter > 0 and five_adapter_overlap > 0:
        line.append("five_adapter_ov")
    if three_adapter > 0 and three_adapter_overlap > 0:
        line.append("three_adapter_ov")
    if trna_overlap > 0:
        line.append("trna_ov")
    stats_file.write("\t".join(line).encode())
    stats_file.write("\n".encode())

    for read in in_samfile:
        aln_start = read.query_alignment_start
        aln_end = read.query_alignment_end - 1

        line = []
        failed = False
        line.append(read.query_name)
        line.append(read.reference_name)
        line.append(str(aln_start))
        line.append(str(aln_end))
        line.append(str(read.get_tag("AS")))
        line.append(read.cigarstring)

        ov = get_overlap(0, five_adapter - 1, aln_start, aln_end)
        if five_adapter > 0 and five_adapter_overlap > 0:
            if ov / five_adapter < five_adapter_overlap:
                failed = True
        line.append(str(ov / five_adapter))

        trna_length = len(trnas[read.reference_name].seq) - five_adapter - three_adapter
        ov = get_overlap(five_adapter, five_adapter + trna_length - 1, aln_start, aln_end)
        if trna_overlap > 0:
            if ov / trna_length < five_adapter_overlap:
                failed = True
        line.append(str(ov / trna_length))

        ov = get_overlap(five_adapter + trna_length, five_adapter + trna_length + three_adapter - 1, aln_start, aln_end)
        if three_adapter > 0 and three_adapter_overlap > 0:
            if ov / three_adapter < three_adapter_overlap:
                failed = True
        line.append(str(ov / three_adapter))

        stats_file.write("\t".join(line).encode())
        stats_file.write("\n".encode())

        if not failed:
            out_samfile.write(read)
    stats_file.close()
    out_samfile.close()
    in_samfile.close()


if __name__ == '__main__':
    process()
