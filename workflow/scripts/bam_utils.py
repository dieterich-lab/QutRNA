import click
import enum
import pysam
import sys
import gzip
from contextlib import nullcontext
from collections import Counter
from Bio import SeqIO


@click.group()
def cli():
    pass


class BaseChange(enum.Enum):
    U2T = enum.auto()
    T2U = enum.auto()


def do_reverse(record):
    record.seq = record.seq[::-1]

    return record


def do_base_change(record, base1, base2):
    record.seq = record.seq.replace(base1, base2)

    return record


def pysam_stdin(fname):
    if fname == "/dev/stdin":
        return "-"

    return fname


@cli.command()
@click.option("-c", "--base-change", type = click.Choice(BaseChange, case_sensitive=False), help="Base change: U2T or T2U.")
@click.option("-r", "--reverse", is_flag=True, default=False, help="Reverse sequence.")
@click.option("-o", "--output", default="-", type=click.Path())
@click.argument("BAM", type=click.Path())
def transform(bam, base_change, reverse, output):
    """Transform read sequence"""

    # io
    in_samfile = pysam.AlignmentFile(pysam_stdin(bam), "rb")
    out_samfile = pysam.AlignmentFile(output, "wb", template=in_samfile)

    # container for tasks for each record
    tasks = []
    if base_change:
        bc = base_change.name.split("2")
        def helper(record):
            do_base_change(record, bc[0], bc[1])
        tasks.append(helper)
    if reverse:
        tasks.append(do_reverse)

    # apply all tasks to all record
    for record in in_samfile:
        for task in tasks:
            record = task(record)
        out_samfile.write(record)

    # close handlers
    out_samfile.close()
    in_samfile.close()


def check_coordinate_sorted(samfile):
    if samfile.header["HD"]["SO"] != "queryname":
        raise Exception("BAM must be sorted by coordinate!")


def pysam_stding(fname):
    if fname == "/dev/stdin":
        return "_"

    return fname


@cli.command()
@click.option("-t", "--tag", default="NH", type=str)
@click.option("-o", "--output", default="-", type=click.Path())
@click.argument("BAM", type=click.Path(exists=True))
def add_hits(bam, tag, output):
    """Add number of hits to a tag"""

    # io
    in_samfile = pysam.AlignmentFile(pysam_stdin(bam), "rb")
    check_coordinate_sorted(in_samfile)
    out_samfile = pysam.AlignmentFile(output, "wb", template=in_samfile)

    def write_records(records):
        """Add number of hits and write reads to output"""

        n = len(records)
        for record in records:
            record.tags += [(tag, n)]
            out_samfile.write(record)

    # container for multimappers
    multimappers = []
    # collect multimappers by counting consecutive records with identical read names
    for record in in_samfile:
        if not multimappers or multimappers[0].qname == record.query_name:
            multimappers.append(record)
        else:
            write_records(multimappers)
            multimappers.clear()
            multimappers.append(record)
    if multimappers:
        write_records(multimappers)
        multimappers.clear()

    # close handler
    out_samfile.close()
    in_samfile.close()


@cli.command()
@click.option("-t", "--tag", required=True, type=str)
@click.option("-c", "--column", required=True, type=str)
@click.option("-o", "--output", type=click.Path())
@click.argument("BAM", type=click.Path(exists=True))
def count_tag(bam, tag, column, output):
    """Count number of values of tag"""

    counts = Counter()
    for read in pysam.AlignmentFile(pysam_stdin(bam), "rb"):
        counts[read.get_tag(tag)] += 1

    with open(output, "w") if output else nullcontext(output) as f:
        line = "\t".join([column, "count"])
        f.write(f"{line}\n")
        for key, count in counts.items():
            line = "\t".join([str(key), str(count)])
            f.write(f"{line}\n")


@cli.command()
@click.option("--min-as", type=int)
@click.argument("bam", type=click.Path(exists=True))
def best_alignments(min_as, bam):
    """Pick the best alignments for each read"""

    # io
    in_samfile = pysam.AlignmentFile(pysam_stdin(bam), "rb")
    check_coordinate_sorted(in_samfile)
    out_samfile = pysam.AlignmentFile(sys.stdout, "wb", template=in_samfile)

    def write_records(as2records, out_samfile):
        max_as = max(as2records.keys())
        for record in as2records[max_as]:
            out_samfile.write(record)

    # container for records
    records = {}
    # current read id
    qname = ""
    # collect all alignments for a mapped read and pick best alignments
    for record in in_samfile:
        # filter by alignment score
        alignment_score = int(record.get_tag("AS"))
        if min_as and alignment_score < min_as:
            continue

        if not qname or qname == record.query_name:
            records.setdefault(alignment_score, [].append(record))
            qname = record.query_name
        else:
            write_records(records, out_samfile)
            records.clear()
            records.setdefault(alignment_score, [].append(record))
            qname = record.query_name
    if records:
        write_records(records, out_samfile)

    # close handlers
    out_samfile.close()
    in_samfile.close()


@cli.command()
@click.option("-t", "--use-tag", type=str)
@click.option("-o", "--output", default="-", type=click.Path())
@click.argument("BAM", type=click.Path(exists=True))
def filter_multimapper(bam, use_tag, output):
    """Remove multimapping reads"""

    # io
    in_samfile = pysam.AlignmentFile(pysam_stdin(bam), "rb")
    out_samfile = pysam.AlignmentFile(output, "wb", template=in_samfile)

    if use_tag:
        for record in in_samfile:
            if record.get_tag(use_tag) == "1":
                out_samfile.write(record)
    else:
        check_coordinate_sorted(in_samfile)

        def write_records(records):
            for record in records:
                out_samfile.write(record)

        multimappers = []
        for record in in_samfile:
            if not multimappers or multimappers[0].qname == record.query_name:
                multimappers.append(record)
            elif len(multimappers) == 1:
                write_records(multimappers)
                multimappers.clear()
                multimappers.append(record)
            else:
                multimappers.clear()
                multimappers.append(record)

        if len(multimappers) == 1:
            write_records(multimappers)
            multimappers.clear()

        # close handlers
        out_samfile.close()
        in_samfile.close()


def get_overlap(start1, end1, start2, end2):
    start = max(start1, start2)
    end = min(end1, end2)
    return max(end - start + 1, 0)


@cli.command()
@click.option("--fasta", type=click.Path(exists=True))
@click.option("--five-adapter", type=int, default=0)
@click.option("--stats", required=True, type=click.Path())
@click.option("--three-adapter", type=int, default=0)
@click.option("--five-adapter-overlap", type=float, default=0.0)
@click.option("--trna-overlap", type=float, default=0.0)
@click.option("--three-adapter-overlap", type=float, default=0.0)
@click.option("-o", "--output", default="-", type=click.Path())
@click.argument("BAM", type=click.Path(exists=True))
def adapter_overlap(bam,
                    fasta,
                    five_adapter, three_adapter,
                    five_adapter_overlap, trna_overlap, three_adapter_overlap,
                    stats,
                    output):
    """Filter by read-adapter overlap"""

    # trna 2 fasta
    trnas = {trna.name: trna for trna in SeqIO.parse(fasta, "fasta")}

    # io
    in_samfile = pysam.AlignmentFile(pysam_stdin(bam), "rb")
    out_samfile = pysam.AlignmentFile(output, "wb", template=in_samfile)

    # create header
    stats_file = gzip.open(stats, "wb")
    header = ["read_id", "ref", "aln_start", "aln_end", "aln_score", "cigar", "ov_5", "ov_trna", "ov_3"]
    if five_adapter > 0 and five_adapter_overlap > 0:
        header.append("five_adapter_ov")
    if three_adapter > 0 and three_adapter_overlap > 0:
        header.append("three_adapter_ov")
    if trna_overlap > 0:
        header.append("trna_ov")
    line = "\t".join(header)
    stats_file.write(f"{line}\n".encode())

    for record in in_samfile:
        aln_start = record.query_alignment_start
        aln_end = record.query_alignment_end - 1

        # container for values
        values = [record.query_name,
                  record.reference_name,
                  str(aln_start),
                  str(aln_end),
                  str(record.get_tag("AS")),
                  record.cigarstring]

        # indicator if constrains are valid
        failed = False

        # overlap with the five prime adapter
        ov = get_overlap(0, five_adapter - 1, aln_start, aln_end)
        if five_adapter > 0 and five_adapter_overlap > 0:
            if ov / five_adapter < five_adapter_overlap:
                failed = True
            values.append(str(ov / five_adapter))

        # overlap with the trna
        trna_length = len(trnas[record.reference_name].seq) - five_adapter - three_adapter
        ov = get_overlap(five_adapter, five_adapter + trna_length - 1, aln_start, aln_end)
        if trna_overlap > 0:
            if ov / trna_length < five_adapter_overlap:
                failed = True
            values.append(str(ov / trna_length))

        # overlap with the three prime adapter
        ov = get_overlap(five_adapter + trna_length, five_adapter + trna_length + three_adapter - 1, aln_start, aln_end)
        if three_adapter > 0 and three_adapter_overlap > 0:
            if ov / three_adapter < three_adapter_overlap:
                failed = True
            values.append(str(ov / three_adapter))

        # create a line from the values
        line = "\t".join(values)
        stats_file.write(f"{line}\n".encode())

        # keep record if valid
        if not failed:
            out_samfile.write(record)

    # close handlers
    stats_file.close()
    out_samfile.close()
    in_samfile.close()


def trim_cigar_helper(cigar):
    aln_offset = 0
    for i, (op, l) in enumerate(cigar):
        if op in (pysam.CEQUAL, pysam.CMATCH, pysam.CHARD_CLIP, pysam.CSOFT_CLIP):
            break

        match op:
            case pysam.CBACK:
                raise Exception()
            case pysam.CDIFF:
                cigar[i] = (pysam.CSOFT_CLIP, l)
                aln_offset += l
            case pysam.CREF_SKIP:
                raise Exception()
            case pysam.CDEL:
                cigar[i] = ()
                aln_offset += l
            case pysam.CINS:
                cigar[i] = (pysam.CSOFT_CLIP, l)
                #aln_offset += l
            case pysam.CPAD:
                raise Exception()

    return cigar, aln_offset


def do_trim_cigar(record):
    try:
        new_cigar, aln_offset = trim_cigar_helper(list(record.cigartuples))
        new_cigar, _ = trim_cigar_helper(list(reversed(new_cigar)))
    except:
        raise Exception(record)
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
    if new_cigar != record.cigartuples:
        record.cigartuples = new_cigar
        if aln_offset:
            record.reference_start += aln_offset

        return True

    return False



@cli.command()
@click.option("--trim-cigar", is_flag=True)
@click.option("--min-as", type=int)
@click.option("--min-read-length", type=int)
@click.option("--max-read-length", type=int)
@click.option("--min-alignment-length", type=int)
@click.option("--max-alignment-length", type=int)
@click.option("--stats", type=str)
@click.option("-o", "--output", type=click.Path())
@click.argument("BAM", type=str)
def filter(bam,
           trim_cigar,
           min_as,
           min_read_length, max_read_length,
           min_alignment_length, max_alignment_length,
           stats, output):
    """Filter or process reads"""

    # io
    in_samfile = pysam.AlignmentFile(pysam_stdin(bam), "rb")
    out_samfile = pysam.AlignmentFile(output, "wb", template=in_samfile)

    # counter
    counter = Counter()

    def do_filter_as(record, min_as):
        alignment_score = int(record.get_tag("AS"))
        if alignment_score >= min_as:
            counter["alignment_score"] += 1

            return True

        return False

    def do_check_read_length(record, min_length, max_length):
        read_length = len(record.query_sequence)

        if min_length and min_length < read_length:
            counter["read_length"] += 1

            return False

        if max_length and max_length > read_length:
            counter["read_length"] += 1

            return False

        return True

    def do_check_alignment_length(record, min_length, max_length):
        aln_length = record.query_alignment_length

        if min_length and min_length < aln_length:
            counter["aln_length"] += 1

            return False

        if max_length and max_length > aln_length:
            counter["aln_length"] += 1

            return False

        return True


    tasks = []
    if min_as:
        def helper(record):
            return do_filter_as(record, min_as)
        tasks.append(helper)

    if min_read_length or max_read_length:
        def helper(record):
            return do_check_read_length(record, min_read_length, max_read_length)
        tasks.append(helper)

    if min_alignment_length or max_alignment_length:
        def helper(record):
            return do_check_alignment_length(record, min_alignment_length, max_alignment_length)
        tasks.append(helper)

    if trim_cigar:
        def helper(record):
            if do_trim_cigar(record):
                counter["trim"] += 1

            return True

        tasks.append(helper)

    for record in in_samfile:
        counter["input"] += 1

        for task in tasks:
            if not task(record):
                continue

        out_samfile.write(record)
        counter["output"] += 1

    # print out statistics
    if stats:
        with open(stats, "w") as f:
            for key, value in counter.items():
                f.write(f"{key}\t{value}\n")

    # close handlers
    out_samfile.close()
    in_samfile.close()


@cli.command()
@click.option("-o", "--output", type=str)
@click.option("-t", "--tag", type=str)
@click.argument("BAM", type=click.Path(exists=True))
def count_multimapper(bam, tag, output):
    """Count multimapper"""

    in_samfile = pysam.AlignmentFile(pysam_stdin(bam), "rb")
    counter = Counter()

    if tag:
        for record in in_samfile:
            counter[record.get_tag(tag)] += 1
    elif in_samfile.header["HD"]["SO"] == "coordinate":
        multimappers = []
        for record in in_samfile:
            if not multimappers or multimappers[0].qname == record.query_name:
                multimappers.append(record)
            else:
                counter[len(multimappers)] += 1
                multimappers.clear()
                multimappers.append(record)
        if multimappers:
            counter[len(multimappers)] += 1
            multimappers.clear()
    else:
        read2hits = Counter()
        for record in in_samfile:
            read2hits[record.query_name] += 1
        for hits in read2hits.values():
            counter[hits] += 1

    with open(output) as f:
        f.write(f"hits\tcount\n")
        for key, value in counter.items():
            f.write(f"{key}\t{value}\n")


if __name__ == '__main__':
    cli()
