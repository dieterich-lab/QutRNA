import click
import enum
import pysam
import sys


class BaseChange(enum.Enum):
    U2T = enum.auto()
    T2U = enum.auto()


@click.command()
@click.option("-c", "--base-change", type = click.Choice(BaseChange, case_sensitive=False), help="Base change: U2T or T2U.")
@click.option("-r", "--reverse",is_flag=True, default=False, help="Reverse sequence.")
@click.option("-o", "--output", required=True, type=click.Path())
@click.argument("bam", type=click.Path())
def process(bam, base_change, reverse, output):

    def do_reverse(read):
        breakpoint()
        read.seq = read.seq[::-1]

        return read

    def do_base_change(base1, base2):
        def helper(read):
            read.seq = read.seq.replace(base1, base2)

            return read

        return helper

    in_samfile = pysam.AlignmentFile(bam, "rb")
    out_samfile = pysam.AlignmentFile(sys.stdout, "wb", template=in_samfile)
    tasks = []
    if base_change:
        bc = base_change.name.split("2")
        tasks.append(do_base_change(bc[0], bc[1]))
    if reverse:
        tasks.append(do_reverse)

    for read in in_samfile:
        for task in tasks:
            read = task(read)
        out_samfile.write(read)
    out_samfile.close()
    in_samfile.close()


if __name__ == '__main__':
    process()
