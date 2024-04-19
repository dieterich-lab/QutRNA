from io import StringIO
import click
import pandas as pd

from Bio import AlignIO

# for EUK
SPRINZL = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
           11, 12, 13, 14, 15, 16, 17, "17a", 18, 19, 20, "20a", "20b",
           21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
           31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
           41, 42, 43, 44, 45,
           "e11", "e12", "e13", "e14", "e15", "e16", "e17",
           "e1", "e2", "e3", "e4", "e5",
           "e27", "e26", "e25", "e24", "e23", "e22", "e21",
           46, 47, 48, 49, 50,
           51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
           61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
           71, 72, 73, 74, 75, 76]


@click.command()
@click.option("--output", type=click.Path(exists=False))
@click.argument("stk", type=click.Path(exists=True))
def mapping(stk, output):
    align = AlignIO.read(stk, "stockholm")
    mapping = {}
    sprinzl_i = 0
    ss_cons = align.column_annotations["secondary_structure"]
    for ss_i, ss in enumerate(ss_cons, start=1):
        if not mapping:
            if ss == ":":
                mapping[ss_i] = SPRINZL[sprinzl_i]
                sprinzl_i += 1
            elif ss == "(":
                sprinzl_i += 1
                mapping[ss_i] = SPRINZL[sprinzl_i]
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
                    sprinzl = SPRINZL[sprinzl_i]
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
    mapping()
