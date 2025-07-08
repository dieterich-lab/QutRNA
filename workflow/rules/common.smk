import re
import pandas as pd
from snakemake import shell


global FILTERS_APPLIED


REF_FASTA = "data/ref.fasta"
REF_NO_LINKER_FASTA = "results/data/no_linker.fasta"
REF_FILTERED_TRNAS_FASTA = "results/data/filtered_trnas.fasta"
CM = "data/cm.stk"
SPRINZL = "data/sprinzl.txt"
SEQ_TO_SPRINZL = "results/seq_to_sprinzl_filtered.tsv"
MODS = "data/mods.tsv"
MOD_ABBREVS = "data/mod_abbrevs.tsv"
DEFAULT_SCORE = "MDI::mismatch_score+deletion_score+insertion_score"


# FASTQ or BAM
READS = ""
if "bam" in pep.sample_table.columns and "fastq" in pep.sample_table.columns:
  raise Exception("Only column bam or fastq can be set in sample table.")
if "bam" in pep.sample_table.columns:
  READS = "bam"
if "fastq" in pep.sample_table.columns:
  READS = "fastq"

# seq or sprinzl coordinates use optional mods annotation
SCORES = "scores"
if pep.config["qutrna"]["coords"] == "sprinzl":
  SCORES += "_sprinzl"
else:
  SCORES += "_seq"
if "mods" in pep.config["qutrna"]:
  SCORES += "-mods"
SCORES += ".tsv"

# FIXME
# sprinzl CM or seq2sprinzl mapping

# FIXME nested cols
__nested_cols = []
TBL = pep.sample_table
for c in ["subsample_name", "base_calling", READS]:
  i = [isinstance(o, list) for o in TBL[c]]
  if any(i):
    __nested_cols.append(c)
if __nested_cols:
  TBL = TBL.explode(__nested_cols)
SAMPLES = TBL["sample_name"].unique().tolist()
SUBSAMPLES = TBL["subsample_name"].tolist()

# FIXME - does not work as expected
wildcard_constraints:
  SAMPLE="|".join(SAMPLES),
  SUBSAMPLE="|".join(SUBSAMPLES),
  ORIENT="|".join(["fwd", "rev"]),
  BC="|".join(["pass", "fail", "merged", "unknown"])


def fname2sample(fname):
  r = re.search("/sample~([^~/]+)", fname)
  return r.group(1)


def fname2subsample(fname):
  r = re.search("/subsample~([^~/]+)", fname)
  return r.group(1)


def create_include(name_, input_, output_, params_):
  rule:
    name: f"dyn_include_{name_}"
    input: input_
    output: output_
    params: include=params_
    run:
        if params.include == "copy":
          cmd = "cp"
        elif params.include == "link":
          cmd = "ln -s"

        shell(cmd + " {input:q} {output:q}")


create_include("ref_fasta",
               pep.config["qutrna"]["ref_fasta"],
               REF_FASTA,
               config["include"]["ref_fasta"])


if pep.config["qutrna"]["coords"] == "sprinzl":
  create_include("sprinzl",
                 pep.config["qutrna"]["sprinzl"],
                 SPRINZL,
                 config["include"].get("sprinzl", "copy"))
  if "cm" in pep.config["qutrna"]:
    create_include("cm",
                   pep.config["qutrna"]["cm"],
                   CM,
                   config["include"]["cm"])
  elif "seq_to_sprinzl" in pep.config["qutrna"]:
    SEQ_TO_SPRINZL = "data/seq_to_sprinzl.tsv"
    create_include("seq_to_sprinzl",
                   pep.config["qutrna"]["seq_to_sprinzl"],
                   SEQ_TO_SPRINZL,
                   config["include"].get("seq_to_sprinzl", "copy"))
  else:
    raise Exception("Must provide other 'cm' or 'seq_to_sprinzl'!")

if "mods" in pep.config["qutrna"]:
  create_include("mods",
                 pep.config["qutrna"]["mods"]["file"],
                 MODS,
                 config["include"]["mods"]["file"])

  if "abbrevs" in pep.config["qutrna"]:
    create_include("abbrevs",
                   pep.config["qutrna"]["mods"]["abbrevs"],
                   MODS,
                   config["include"]["mods"]["abbrevs"])


def _include_fnames(ftype):
  def helper(wildcards):
    tbl = TBL.loc[[wildcards.SAMPLE]]
    i = (tbl.subsample_name == wildcards.SUBSAMPLE) & (tbl.base_calling == wildcards.BC)
    if sum(i) != 1:
      raise Exception("")

    return tbl.loc[i, ftype].to_list()[0]

  return helper


create_include("bams",
               _include_fnames("bam"),
               "data/bams/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam",
               config["include"]["bam"])

if READS != "bams":
  create_include("fastq",
               _include_fnames("fastq"),
               "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz",
               config["include"]["fastq"])


rule remove_linker:
  input: REF_FASTA
  output: REF_NO_LINKER_FASTA
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/remove_linker.log"
  params: linker5=pep.config["qutrna"]["linker5"],
          linker3=pep.config["qutrna"]["linker3"]
  shell: """
    python {workflow.basedir}/scripts/remove_linker.py \
        --linker5 {params.linker5} \
        --linker3 {params.linker3} \
        --output {output:q} \
        {input:q} \
        2> "{log}"
  """


def _remove_trans_opts(_):
  trnas = pep.config["qutrna"].get("remove_trnas", [])

  opts = []
  for trna in trnas:
    opts.append(f"-t {trna}")

  return " ".join(opts)


rule remove_trnas:
  input: REF_NO_LINKER_FASTA
  output: REF_FILTERED_TRNAS_FASTA
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/remove_trnas.log"
  params: opts=_remove_trans_opts
  shell: """
    python {workflow.basedir}/scripts/remove_trnas.py \
        {params.opts} \
        --output {output:q} \
        {input:q} \
        2> "{log}"
  """


def _collect_input(suffix):
  def helper(wildcards):
    parsed = suffix.format(wildcards=wildcards)
    t2fnames = {}
    for sample in SAMPLES:
      df = TBL.loc[[sample]]
      # collect subsamples
      for row in df.itertuples(index=False):
        if not hasattr(row, "bam"):
          raise Exception("Not implemented yet") # TODO implement, when BAM not available
        t2fnames.setdefault("raw", []).append(f"data/bams/sample~{sample}/subsample~{row.subsample_name}/{row.base_calling}.{parsed}") # TODO what if no base calling
        for read_type in FILTERS_APPLIED:
          t2fnames.setdefault(read_type, []).append(f"results/bams/preprocessed/{read_type}/sample~{sample}/subsample~{row.subsample_name}/{row.base_calling}.{parsed}") # TODO what if no base calling
    return t2fnames

  return helper


rule collect_read_counts:
  input: unpack(_collect_input("sorted_read_count.txt"))
  output: "results/read_counts.tsv"
  run:
    dfs = []
    for read_type, fnames in input.items():
      for fname in fnames:
        df = pd.read_csv(fname, sep="\t", header=None, names=["read_count"])

        df["read_type"] = read_type
        df[["sample", "subsample", "base_calling"]] = ""
        sample = fname2sample(fname)
        df["sample"] = sample
        df["subsample"] = fname2subsample(fname)
        df["base_calling"] = re.search("/([^/]+).sorted_read_count.txt$", fname).group(1)
        df["condition"] = TBL.loc[[sample]]["condition"].unique()[0]
        df["fname"] = fname

        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(output[0], sep="\t", index=False)


rule collect_as:
  input: unpack(_collect_input("sorted_alignment_score.tsv"))
  output: "results/alignment_scores.tsv"
  run:
    dfs = []
    for read_type, fnames in input.items():
      for fname in fnames:
        df = pd.read_csv(fname, sep="\t", header=None, names=["as", "count"])

        df["read_type"] = read_type
        df[["sample", "subsample", "base_calling"]] = ""
        sample = fname2sample(fname)
        df["sample"] = sample
        df["subsample"] = fname2subsample(fname)
        df["base_calling"] = re.search("/([^/]+).sorted_alignment_scores.tsv$", fname).group(1)
        df["condition"] = TBL.loc[[sample]]["condition"].unique()[0]
        df["fname"] = fname

        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(output[0], sep="\t", index=False)


rule cut_samtools_stats:
  input: "{prefix}.stats"
  output: temp("{prefix}_stats_{type}.tsv")
  conda: "qutrna"
  resources:
    mem_mb=2000
  shell: """
    grep "^{wildcards.type}" {input:q} | cut -f 2- > {output:q}
  """


rule collect_samtools_stats:
  input: unpack(_collect_input("sorted_stats_{wildcards.type}.tsv"))
  output: "results/samtools/stats/{type}.tsv"
  run:
    type2cols = {
        "RL": ["read_length", "count"],}
    dfs = []
    for read_type, fnames in input.items():
      for fname in fnames:
        df = pd.read_csv(fname, sep="\t", header=None, names=type2cols[wildcards.type])
        sample = fname2sample(fname)
        df["read_type"] = read_type
        df["sample"] = sample
        df["subsample"] = fname2subsample(fname)
        df["base_calling"] = re.search("^(.+).sorted.+$", os.path.basename(fname)).group(1)
        df["condition"] = TBL.loc[[sample]]["condition"].unique()[0]
        df["fname"] = fname
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(output[0], sep="\t", index=False)
