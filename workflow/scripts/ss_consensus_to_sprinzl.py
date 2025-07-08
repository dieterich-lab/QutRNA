from io import StringIO
import click
import pandas as pd

from Bio import AlignIO


@click.command()
@click.option("--output", type=click.Path())
@click.option("--sprinzl", type=click.Path(exists=True))
@click.argument("stk", type=click.Path(exists=True))
def annotate(sprinzl, stk, output):
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
    if not output:
        output = StringIO()
        df.to_csv(output, sep="\t", index=False, quoting=False)
        print(output.getvalue())
    else:
        df.to_csv(output, sep="\t", index=False, quoting=False)


if __name__ == '__main__':
    annotate()
