REF_FASTA = "data/ref.fasta"
REF_NO_LINKER_FASTA = "results/data/no_linker.fasta"
REF_FILTERED_TRNAS_FASTA = "results/data/filtered_trnas.fasta"
CM = "data/cm.stk"
MODS = "data/mods.tsv"
MOD_ABBREVS = "data/mod_abbrevs.tsv"

TBL = pep.sample_table.explode("subsample_name")
SAMPLES = TBL["sample_name"].unique().tolist()
SUBSAMPLES = TBL["subsample_name"].tolist()


wildcard_constraints:
  SAMPLE="|".join(SAMPLES),
  SUBSAMPLE="|".join(SUBSAMPLES),
  ORIENT="|".join(["fwd", "rev"]),
  BC="|".join(["pass", "fail", "merged", "unknown"]),

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
  SCORES += "scores_seq"
if "mods" in pep.config["qutrna"]:
  SCORES += "-mods"
SCORES += ".tsv"


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

        shell(cmd + " {input:q} {output:q}")


create_include("ref_fasta",
               pep.config["qutrna"]["ref_fasta"],
               REF_FASTA,
               config["include"]["ref_fasta"])

create_include("cm",
               pep.config["qutrna"]["cm"],
               CM,
               config["include"]["cm"])

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
    tbl = pep.sample_table.explode(["subsample_name", ftype])
    tbl = tbl.set_index("subsample_name")
    fnames = tbl.loc[wildcards.SUBSAMPLE, ftype]

    return fnames

  return helper

create_include("bams",
               _include_fnames("bam"),
               "data/bams/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam",
               config["include"]["bam"])

create_include("fastq",
             _include_fnames("fastq"),
             "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz",
             config["include"]["fastq"])


rule remove_linker:
  input: REF_FASTA,
  output: REF_NO_LINKER_FASTA,
  conda: "qutrna",
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
  log: "logs/remove_trnas.log",
  params: opts=_remove_trans_opts,
  shell: """
    python {workflow.basedir}/scripts/remove_trnas.py \
        {params.opts} \
        --output {output:q} \
        {input:q} \
        2> {log:q}
  """