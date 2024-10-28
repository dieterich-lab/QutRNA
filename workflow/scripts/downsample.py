from os.path import basename
import random
import string
import os

import numpy as np
import pysam
import pandas as pd
import random

dname = "/beegfs/prj/tRNA_Francesca_Tuorto/data/FT_tRNA_custom_2_biol_reps/downsampling_12-100reads/"
fnames = [dname + fname for fname in ["DNMT2_1.bam", "NSUN2_1.bam", "QTRT1_1.bam", "wt_1.bam", "DNMT2_2.bam", "NSUN2_2.bam", "QTRT1_2.bam", "wt_2.bam"]]
seeds = [random.randrange(1, 10000) for _ in range(4)]
# fnames = ""

# colect coverages per file and trna
coverage_info = []
bfnames = []
for fname in fnames:
    bfname = basename(fname)
    bfnames.append(bfname)

    infile = pysam.AlignmentFile(fname, "rb")
    for ref in infile.references:
        coverage_info.append({"ref": ref, "fname": bfname, "coverage": infile.count(ref)})
df = pd.DataFrame.from_records(coverage_info)
df = df.pivot(index="ref", columns="fname", values="coverage")
df["min_reads"] = df.min(axis="columns")

# remove trnas with zero reads
#zero_reads_i = df["min_reads"] == 0
CUT = 100
zero_reads_i = df["min_reads"] < CUT
if any(zero_reads_i):
    zero_reads_df = df[zero_reads_i]
    df = df[~zero_reads_i]
df["min_reads"] = CUT

# not needed
## calculate target fractions for sampling
#for bfname in bfnames:
#    df["target_" + bfname] = df["min_reads"] / df[bfname]

df = df.reset_index()

# create output dirs
for seed in seeds:
    try:
        os.mkdir(dname + "downsampling")
    except FileExistsError:
        pass

    try:
        os.mkdir(dname + "downsampling/" + str(seed))
    except FileExistsError:
        pass


for bfname in bfnames:
    input_file = f"{dname}{bfname}"
    input_bam = pysam.AlignmentFile(input_file, "rb")
    for seed in seeds:
        np.random.seed(seed)
        tmp_file = output_file = f"{dname}downsampling/{seed}/tmp.{bfname}"
        output_file = f"{dname}downsampling/{seed}/{bfname}"
        output_bam = pysam.AlignmentFile(tmp_file, "wb", template=input_bam)
        for idx, row in df[["ref", bfname, "min_reads"]].iterrows():
            read_count = row[bfname]
            idx = np.random.choice(read_count, row["min_reads"])
            reads = list(input_bam.fetch(row["ref"]))
            for i in idx:
               output_bam.write(reads[i])
        output_bam.close()
        pysam.sort("-o", f"{output_file}", tmp_file)
    input_bam.close()
