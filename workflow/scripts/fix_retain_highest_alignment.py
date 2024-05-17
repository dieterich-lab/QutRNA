import pysam
import click
import sys

# make sure this is sorted by readname
infile = pysam.AlignmentFile(sys.stdin, "rb")
outfile = pysam.AlignmentFile(sys.stdout, "w", template=infile)


def fix_mapq(r):
    aln_score = r.get_tag("AS")
    aln_score = min(aln_score, 255)
    r.mapping_quality = aln_score

    return r


@click.command()
@click.option("-m", "--min-mapq",
              type=int,
              default=20)
def filter(min_mapq):
    t = None
    for r in infile:
        if r.has_tag("AS"):
            r = fix_mapq(r)

        if r.mapping_quality < min_mapq:
            continue

        if t:
            if t.qname == r.qname:
                if r.mapping_quality > t.mapping_quality:
                    t = r
            else:
                outfile.write(t)
                t = r

        else:
            t = r

    if t:
        outfile.write(t)

    infile.close()
    outfile.close()


if __name__ == '__main__':
    filter()
