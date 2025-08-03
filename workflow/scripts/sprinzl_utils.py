import click
import pandas as pd
from Bio import AlignIO


@click.group()
def cli():
    pass


# TODO add
@cli.command()
@click.option("-o", "--output", required=True, type=click.Path())
@click.option("-s", "--sprinzl", required=True, type=click.Path(exists=True))
@click.argument("STK", type=click.Path(exists=True))
def annotate_consensus(stk, sprinzl, output):
    """Label consensus secondary structure alignment"""

    # load ss consensus annotation
    sprinzl_df = pd.read_csv(sprinzl, header=None)[0].to_list()

    # load alignments
    align = AlignIO.read(stk, "stockholm")
    mapping = {}
    sprinzl_i = 0
    # pick ss consensus
    ss_cons = align.column_annotations["secondary_structure"]
    # assign sprinzl coordinates to ss consensus
    for ss_i, ss in enumerate(ss_cons, start=1):
        if not mapping:
            if ss == ":":
                mapping[ss_i] = sprinzl_df[sprinzl_i]
                sprinzl_i += 1
            elif ss == "(":
                sprinzl_i += 1
                mapping[ss_i] = sprinzl_df[sprinzl_i]
                sprinzl_i += 1
            else:
                mapping[ss_i] = "-"
        else:
            if ss == ".":
                mapping[ss_i] = "-"
            elif ss == "~":
                mapping[ss_i] = "~"
            else:
                try:
                    sprinzl = sprinzl_df[sprinzl_i]
                except IndexError:
                    sprinzl = "-"
                mapping[ss_i] = sprinzl
                sprinzl_i += 1

    df = pd.DataFrame(mapping.items(), columns=["cons", "sprinzl"])
    df["ss"] = list(ss_cons)
    df.to_csv(output, sep="\t", index=False, quoting=False)

@cli.command()
@click.option("--output", required=True, type=click.Path())
@click.option("--cons-ss-annotation", required=True, type=click.Path(exists=True))
@click.argument("STK", type=click.Path(exists=True))
def map_seq_sprinzl(stk, cons_ss_annotation, output):
    """Map sequence to sprinzl coorindates"""

    # read cmalign
    align = AlignIO.read(stk, "stockholm")
    # read conserved secondary structure with annotation (sprinzl)
    cons_ss_annotation = pd.read_csv(cons_ss_annotation, sep="\t").set_index("cons")

    dfs = []
    for record in align:
        la_sprinzl = []
        la_seq_pos = []
        seq_pos = 0

        intron_i = 1
        last_sprinzl = ""
        for i, letter in enumerate(record.seq, start=1):
            if letter in ["A", "C", "G", "U", "T", "a", "c", "g", "u", "t"]:
                seq_pos += 1
                la_seq_pos.append(seq_pos)
            else:
                la_seq_pos.append("-")

            if i > max(cons_ss_annotation.index):
                sprinzl_pos = "-"
            else:
                sprinzl_pos = str(cons_ss_annotation.loc[i, "sprinzl"])

            if not sprinzl_pos in ["-", "~"]:
                last_sprinzl = sprinzl_pos

            if letter in ["A", "C", "G", "U", "T"]:
                la_sprinzl.append(sprinzl_pos)
            elif letter in ["a", "c", "g", "u", "t"]:
                if last_sprinzl == str(37):
                    la_sprinzl.append("i" + str(intron_i))
                    intron_i += 1
                else:
                    la_sprinzl.append(sprinzl_pos)
            else:
                la_sprinzl.append("-")
        df = pd.DataFrame(
            {"id": record.id,
             "seq": list(str(record.seq)),
             "ss": list(align.column_annotations["secondary_structure"]),
             "seq_pos": la_seq_pos,
             "aln_pos": range(1, len(record.seq) + 1),
             "sprinzl": la_sprinzl})
        dfs.append(df)

    df = pd.concat(dfs)
    df.to_csv(output, sep="\t", index=False, quoting=False)


@cli.command()
@click.option("--sprinzl", required=True, help="Sequence to Sprinzl.")
@click.option("--output", required=True, help="Output FNAME")
@click.option("--linker5", default=0, help="Length of 5' linker sequence")
@click.argument("jacusa2", type=click.Path(exists=True))
def transform(sprinzl, output, linker5, jacusa2):
    """Add Sprinzl coordinates to JACUS2A output"""

    sprinzl = pd.read_csv(sprinzl, sep="\t")
    jacusa = pd.read_csv(jacusa2, sep="\t")

    i = sprinzl["id"].isin(jacusa["trna"].unique())
    sprinzl = sprinzl.loc[i, ["id", "seq_pos", "sprinzl"]]

    jacusa["n_pos"] = jacusa["seq_position"] - linker5
    jacusa = (jacusa.merge(sprinzl,
                           how="left",
                           left_on=("trna", "n_pos" ),
                           right_on=("id", "seq_pos"),
                           indicator=True)
              .drop(columns=["n_pos", "seq_pos", "id"]))
    jacusa["sprinzl"] = jacusa["sprinzl"].fillna(".")

    jacusa.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    cli()
