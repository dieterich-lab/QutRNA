REF_FASTA = "data/ref.fasta"
REF_NO_LINKER_FASTA = "results/data/no_linker.fasta"
REF_FILTERED_TRNAS_FASTA = "results/data/filtered_trnas.fasta"
CM = "data/cm.stk"
SPRINZL = "data/sprinzl.txt"
SEQ_TO_SPRINZL = "results/seq_to_sprinzl_filtered.tsv"
MODS = "data/mods.tsv"
MOD_ABBREVS = "data/mod_abbrevs.tsv"


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


__nested_cols = []
TBL = pep.sample_table
for c in ["subsample_name", "base_calling", READS]:
  if any(isinstance(o, list) for o in TBL[c]):
    __nested_cols.append(c)
if __nested_cols:
  TBL = TBL.explode(__nested_cols)
SAMPLES = TBL["sample_name"].unique().tolist()
SUBSAMPLES = TBL["subsample_name"].tolist()

wildcard_constraints:
  SAMPLE="|".join(SAMPLES),
  SUBSAMPLE="|".join(SUBSAMPLES),
  ORIENT="|".join(["fwd", "rev"]),
  BC="|".join(["pass", "fail", "merged", "unknown"]),


def create_include(name, input, output, params):
  rule:
    name: f"dyn_include_{name}"
    input: input
    output: output
    params: include=params
    run:
        if params.include == "copy":
          cmd = "cp"
        else:
          cmd = "ln -s"

        cmd = "cp" # FIXME

        shell(cmd + " {input:q} {output:q}")


create_include("ref_fasta",
               pep.config["qutrna"]["ref_fasta"],
               REF_FASTA,
               config["include"]["ref_fasta"])


if pep.config["qutrna"]["coords"] == "sprinzl":
  create_include("sprinzl",
                 pep.config["qutrna"]["sprinzl"],
                 SPRINZL,
                 config["include"].get("sprinzl", True))
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
                   config["include"].get("seq_to_sprinzl", True))
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
  input: REF_FASTA,
  output: REF_NO_LINKER_FASTA,
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: "logs/remove_linker.log",
  params: linker5=pep.config["qutrna"]["linker5"],
          linker3=pep.config["qutrna"]["linker3"],
  shell: """
    python {workflow.basedir}/scripts/remove_linker.py \
        --linker5 {params.linker5} \
        --linker3 {params.linker3} \
        --output {output:q} \
        {input:q} \
        2> {log:q}
  """


def _remove_trans_opts(wildcards):
  trnas = pep.config["qutrna"].get("remove_trnas", [])

  opts = []
  for trna in trnas:
    opts.append(f"-t {trna}")

  return " ".join(opts)


rule remove_trnas:
  input: REF_NO_LINKER_FASTA,
  output: REF_FILTERED_TRNAS_FASTA,
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: "logs/remove_trnas.log",
  params: opts=_remove_trans_opts,
  shell: """
    python {workflow.basedir}/scripts/remove_trnas.py \
        {params.opts} \
        --output {output:q} \
        {input:q} \
        2> {log:q}
  """
