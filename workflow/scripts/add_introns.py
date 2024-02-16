import click
import pandas as pd
import re
import yaml
from itertools import repeat
from Bio import SeqIO
from io import StringIO

@click.command()
@click.option("--ssaln", required=True,
              type=click.File("r"),
              help="ss file from tRNAscan-SE.")
@click.option("--fasta",
              required=True,
              type=click.File("r"),
              help="FASTA file. Mature sequences having '_exononly' suffix")
@click.argument("coord_map", type=click.Path("r"))
def add_introns(ssaln, fasta, coord_map):
    """Transform secondary structure alignment for tRNA with introns"""

    fasta = {record.id: record for record in SeqIO.parse(fasta, "fasta")}
    ssaln = read_ssaln(ssaln)
    coord_map = pd.read_csv(coord_map, sep="\t")
    coord_map["trna_id"] = coord_map["trna_id"].str.replace("-.+$", "", regex=True)

    trna_names = list(fasta.keys())
    trnas_with_intron = [trna.id for trna in ssaln.values() if trna.has_intron()]
    trnas_exononly = [id for id in fasta if id.endswith("_exononly")]
    new_coords = None
    for trna_name in trnas_exononly:
        name = trna_name.replace("_exononly", "")
        coords = coord_map.loc[coord_map["ref"] == name].reset_index(drop=True)
        coords.loc[:, "ref"] = trna_name
        if new_coords is None:
            new_coords = coords
        else:
            new_coords = pd.concat([new_coords, coords], ignore_index=True)

    for trna_id in trnas_with_intron:
        coords = coord_map.loc[coord_map["trna_id"] == trna_id]
        trna = ssaln[trna_id]
        length = len(trna.intron_seq)
        trna_name = coords["ref"].unique()[0]
        pre_intron = coords[coords["seq_pos"] < trna.intron_start - 1].copy()
        post_intron = coords[coords["seq_pos"] >= trna.intron_start - 1].copy()
        post_intron["seq_pos"] += length
        intron = pd.DataFrame({
            "trna_id": length * [trna.id],
            "seq_pos": range(trna.intron_start - 1, trna.intron_end),
            "u_pos": ["i" + str(i) for i in range(1, length + 1)],
            "ref": length * [trna_name]
        })
        df = pd.concat([pre_intron, intron, post_intron]).reset_index(drop=True)
        if new_coords is None:
            new_coords = df
        else:
            new_coords = pd.concat([new_coords, df])

    for trna_name in set(trna_names) - set(new_coords["ref"].to_list()):
        coords = coord_map.loc[coord_map["ref"] == trna_name].reset_index(drop=True)
        if new_coords is None:
            new_coords = coords
        else:
            new_coords = pd.concat([new_coords, coords], ignore_index=True)

    s = StringIO()
    new_coords.to_csv(s, sep="\t", index=False)
    print(s.getvalue())


class tRNA:
    def __init__(self, id):
        self.id = id
        self.intron_start = -1
        self.intron_end = -1
        self.seq = ""
        self.str = ""

    def set_intron(self, start, end):
        """1-based [] inclusive"""
        self.intron_start = start
        self.intron_end = end

    @property
    def intron_seq(self):
        if not self.has_intron():
            return ""

        return self.seq[self.intron_start - 1:self.intron_end]

    def has_intron(self):
        return self.intron_start >= 0 and self.intron_end >= 0 and self.seq != ""

    @property
    def property_mature_seq(self):
        if not self.has_intron():
            return self.seq

        return self.seq[0: self.intron_start] + self.seq[self.intron_end - 1: ]

    @property
    def mature_str(self):
        if not self.has_intron():
            return self.str

        return self.str[0: self.intron_start] + self.str[self.intron_end - 1: ]

    def check(self):
        if not self.seq or not self.str:
            raise Exception

    @property
    def column_annotation(self):
        return {
            "secondary_structure": self.str,}

def read_ssaln(ssaln):
    trnas = {}
    trna = None
    for line in ssaln:
        line = line.strip()
        if not line:
            continue

        res = re.search(r"^(\S+)\s+(\((\d+\-\d+,?)+\))\s+Length: (\d+)\sbp", line)
        if res:
            trnascan_id = res.group(1)
            trna = tRNA(trnascan_id)
            continue

        res = re.search(r"^Possible intron: (\d+)-(\d+)", line)
        if res:
            trna.set_intron(
                int(res.group(1)),
                int(res.group(2)))

        res = re.search(r"^Seq: (\S+)", line)
        if res:
            trna.seq = res.group(1)

        res = re.search(r"^Str: (\S+)", line)
        if res:
            trna.str = res.group(1)
            trna.check()
            trnas[trna.id] = trna
            continue

    return trnas


if __name__ == '__main__':
    add_introns()
