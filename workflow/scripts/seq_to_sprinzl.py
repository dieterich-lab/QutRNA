from io import StringIO
import click
import pandas as pd

from Bio import AlignIO


@click.command()
@click.option("--output", type=click.Path(exists=False))
@click.argument("ss_cons_mapping", type=click.Path(exists=True))
@click.argument("stk", type=click.Path(exists=True))
def mapping(ss_cons_mapping, stk, output):
    align = AlignIO.read(stk, "stockholm")
    ss_cons_mapping = pd.read_csv(ss_cons_mapping, sep="\t").set_index("cons")

    dfs = []
    for record in align:
        la_sprinzl = []
        la_seq_pos = []
        seq_pos = 0

        intron_i = 1
        last_sprinzl = ""
        for i, letter in enumerate(record.seq, start=1):
            if letter in ["A", "C", "G", "U", "a", "c", "g", "u"]:
                seq_pos += 1
                la_seq_pos.append(seq_pos)
            else:
                la_seq_pos.append("-")

            if i > max(ss_cons_mapping.index):
                sprinzl_pos = "-"
            else:
                sprinzl_pos = str(ss_cons_mapping.loc[i, "sprinzl"])

            if not sprinzl_pos in ["-", "~"]:
                last_sprinzl = sprinzl_pos

            if letter in ["A", "C", "G", "U"]:
                la_sprinzl.append(sprinzl_pos)
            elif letter in ["a", "c", "g", "u"]:
                if last_sprinzl == str(37):
                    la_sprinzl.append("i" + str(intron_i))
                    intron_i += 1
                else:
                    la_sprinzl.append(sprinzl_pos)
            else:
                la_sprinzl.append("-")
        df = pd.DataFrame(
                {"id": record.id,
                 "name": record.name,
                 "seq": list(str(record.seq)),
                 "ss": list(align.column_annotations["secondary_structure"]),
                 "seq_pos": la_seq_pos,
                 "aln_pos": range(1, len(record.seq) + 1),
                 "sprinzl": la_sprinzl})
        dfs.append(df)

    df = pd.concat(dfs)
    df.to_csv(output, sep="\t", index=False, quoting=False)


if __name__ == '__main__':
    mapping()
