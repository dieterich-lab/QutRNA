from snakemake import shell
from snakemake.io import unpack


global FILTERS_APPLIED

TRNAS_AVAILABLE = "results/data/trnas/available.txt"
TRNAS_SELECTED = "results/data/trnas/selected.txt"
TRNA_ANNOTATION = "results/data/trna_annotation.tsv"

wildcard_constraints:
  COND1="[^/~]+",
  COND2="[^/~]+",
  BAM_TYPE="[^/]+",
  SAMPLE="[^~/]+",
  SUBSAMPLE="[^~/]+",
  ORIENT="[^~/]+",
  BC="[^~/]+"


def create_include(name_, input_, output_, params_):
  rule:
    name: f"include_{name_}"
    input: input_
    output: output_
    params: include=params_
    run:
      if params.include == "copy":
        cmd = "cp"
      elif params.include == "link":
        cmd = "ln -s"

      shell(cmd + " {input:q} {output:q}")


# coordinate system: sequence or sprinzl
# in case of sprinzl handle if model or precalculate mapping
if pep.config["qutrna2"]["coords"] == "sprinzl":
  COORDS = "sprinzl"
  SPRINZL = "data/sprinzl.txt"
  create_include("sprinzl",
    pep.config["qutrna2"]["sprinzl"],
    SPRINZL,
    config["include"].get("sprinzl","copy"))
  if "cm" in pep.config["qutrna2"]:
    SPRINZL_MODE = "cm"
    CM = "data/cm.stk"
    create_include("cm",
      pep.config["qutrna2"]["cm"],
      CM,
      config["include"]["cm"])
  elif "seq_to_sprinzl" in pep.config["qutrna2"]:
    SPRINZL_MODE = "seq2sprinzl"
    create_include("seq_to_sprinzl",
      pep.config["qutrna2"]["seq_to_sprinzl"],
      "data/seq_to_sprinzl.tsv",
      config["include"].get("seq_to_sprinzl","copy"))
  else:
    raise Exception("Must provide other 'cm' or 'seq_to_sprinzl' in config!")
  SEQ_TO_SPRINZL = "results/ss/seq_to_sprinzl_filtered.tsv"
else:
  COORDS = "sequence"


# include and transform reference
REF_FASTA = "data/ref.fasta"
create_include("ref_fasta",
  pep.config["qutrna2"]["ref_fasta"],
  REF_FASTA,
  config["include"]["ref_fasta"])
REF_NO_LINKER_FASTA = "results/data/ref_no_linker.fasta"
REF_FILTERED_TRNAS_FASTA = "results/data/ref_filtered_trnas.fasta"
DEFAULT_SCORE = "MDI::mismatch_score+deletion_score+insertion_score"


# include bams or fastq from sample table
def _include_fnames(ftype):
  def helper(wildcards):
    tbl = TBL.loc[[wildcards.SAMPLE]]
    j = (tbl.subsample_name == wildcards.SUBSAMPLE) & (tbl.base_calling == wildcards.BC)
    if sum(j.tolist()) != 1:
      raise Exception("Housten, we have a problem!")

    return tbl.loc[j, ftype].to_list()[0]

  return helper

# FASTQ or BAM
if "bam" in pep.sample_table.columns and "fastq" in pep.sample_table.columns:
  raise Exception("Only column 'bam' or 'fastq' can be set in sample table.")

READS = "" # will be bam or fastq
READS_INPUT = "" # will point to a path where bam or fastq can be found
if "bam" in pep.sample_table.columns:
  READS = "bam"
  READS_INPUT = "data/bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam"

  create_include("bams",
    _include_fnames("bam"),
    READS_INPUT,
    config["include"]["bam"])

  if "bam" in config["transform"]:
    base_change = config["transform"]["bam"].get("base_change")
    reverse_seq = config["transform"]["bam"].get("reverse")
    if base_change or reverse_seq:
      _old_reads_input = READS_INPUT
      READS_INPUT = "results/bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam",

      rule transform_bam:
        input: _old_reads_input
        output: READS_INPUT
        params: base_change=lambda wildcards: "--base-change {base_change} " if base_change else "",
                reverse=lambda wildcards: "--reverse " if reverse_seq else ""
        conda: "qutrna2"
        log: "logs/transform_bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
        shell: """
          python {workflow.basedir}/scripts/bam_utils.py transform {params.reverse} {params.base_change} --output {output:q} {input:q} 2> {log:q}
        """
elif "fastq" in pep.sample_table.columns:
  READS = "fastq"
  READS_INPUT = "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz"
  REF_FASTA_REVERSED = "results/data/ref_reversed.fasta"

  create_include("fastq",
    _include_fnames("fastq"),
    READS_INPUT,
    config["include"]["fastq"])

  if config["transform"].get("fastq", {}).get("base_change"):
    _old_reads_input = "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz"
    READS_INPUT = "results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz"

    rule transform_fastq:
      input: _old_reads_input
      output: READS_INPUT
      log: "logs/transform_fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
      params: base_change=config["transform"]["fastq"]["base_change"]
      conda: "qutrna2"
      shell: """
        python {workflow.basedir}/scripts/fastq_utils.py transform --output {output:q} --base-change {params.base_change} {input:q} 2> {log:q}
      """

  rule reverse_fasta:
    input: REF_FASTA
    output: REF_FASTA_REVERSED
    log: "logs/reverse_fasta.log"
    conda: "qutrna2"
    shell: """
        python {workflow.basedir}/scripts/fasta_utils.py transform --output {output:q} --reverse {input:q} 2> {log:q}
      """
else:
  raise Exception("Missing ('bam' or ' 'fastq') column in sample table.")


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


# seq or sprinzl coordinates use optional mods annotation
SCORES = "scores"
if COORDS == "sprinzl":
  SCORES += "_sprinzl"
else:
  SCORES += "_seq"
if "mods" in pep.config["qutrna2"]:
  SCORES += "-mods"
  MODS = "data/mods.tsv"
  MOD_ABBREVS = "data/mod_abbrevs.tsv"
SCORES += ".tsv"

if "mods" in pep.config["qutrna2"]:
  create_include("mods",
                 pep.config["qutrna2"]["mods"]["file"],
                 MODS,
                 config["include"]["mods"]["file"])

  if "abbrevs" in pep.config["qutrna2"]:
    create_include("abbrevs",
                   pep.config["qutrna2"]["mods"]["abbrevs"],
                   MODS,
                   config["include"]["mods"]["abbrevs"])


rule remove_linker:
  input: REF_FASTA
  output: REF_NO_LINKER_FASTA
  conda: "qutrna2"
  log: "logs/remove_linker.log"
  params: linker5=pep.config["qutrna2"]["linker5"],
          linker3=pep.config["qutrna2"]["linker3"]
  shell: """
    python {workflow.basedir}/scripts/fasta_utils.py transform \
        --remove-linker5 {params.linker5} \
        --remove-linker3 {params.linker3} \
        --output {output:q} \
        {input:q} \
        2> {log:q}
  """


def _ignore_trans_opts(_):
  trnas = pep.config["qutrna2"].get("ignore_trnas", [])

  opts = []
  for trna in trnas:
    opts.append(f"--ignore '{trna}'")

  return " ".join(opts)


rule ignore_trnas:
  input: REF_NO_LINKER_FASTA
  output: REF_FILTERED_TRNAS_FASTA
  conda: "qutrna2"
  log: "logs/ignore_trnas.log"
  params: opts=_ignore_trans_opts
  shell: """
    python {workflow.basedir}/scripts/fasta_utils.py \
        transform  \
        {params.opts} \
        --output {output:q} \
        {input:q} \
        2> {log:q}
  """


def _aggregate_stats_helper(read_type, sample, row):
  key = [row.condition, sample, row.subsample_name, row.base_calling]
  return "|".join([read_type] + key)

# FIXME gzip
def _aggregate_stats_input(wildcards):
  if wildcards.feature == "alignment_score" and READS == "fastq":
    t2fnames = _aggregate_stats_general_input(wildcards)
    for sample in SAMPLES:
      df = TBL.loc[[sample]]
      for row in df.itertuples(index=False):
        t2fnames[_aggregate_stats_helper("mapped-rev", sample, row)] = f"results/bam/mapped/sample~{sample}/subsample~{row.subsample_name}/orient~rev/{row.base_calling}_stats/{wildcards.feature}.txt"
  elif wildcards.feature == "cutoff" and READS == "fastq":
    t2fnames = {}
    for sample in SAMPLES:
      df = TBL.loc[[sample]]
      for row in df.itertuples(index=False):
        t2fnames[_aggregate_stats_helper("mapped",sample,row)] = f"results/bam/mapped/sample~{sample}/subsample~{row.subsample_name}/{row.base_calling}_stats/{wildcards.feature}.txt"
  else:
    t2fnames = _aggregate_stats_general_input(wildcards)

  return t2fnames

def _aggregate_stats_general_input(wildcards):
  t2fnames = {}
  for sample in SAMPLES:
    df = TBL.loc[[sample]]

    # collect subsamples
    for row in df.itertuples(index=False):
      if hasattr(row, "bam"):
        t2fnames[_aggregate_stats_helper("raw", sample, row)] = f"data/bam/sample~{sample}/subsample~{row.subsample_name}/{row.base_calling}_stats/{wildcards.feature}.txt"
      elif hasattr(row, "fastq"):
        t2fnames[_aggregate_stats_helper("mapped", sample, row)] = f"results/bam/mapped/sample~{sample}/subsample~{row.subsample_name}/orient~fwd/{row.base_calling}_stats/{wildcards.feature}.txt"
      else:
        raise Exception(f"READS must be ('bam' or 'fastq')")

      for read_type in FILTERS_APPLIED:
        t2fnames[_aggregate_stats_helper(read_type, sample, row)] = f"results/bam/filtered-{read_type}/sample~{sample}/subsample~{row.subsample_name}/{row.base_calling}_stats/{wildcards.feature}.txt"

  return t2fnames

def _aggregate_stats_params(_, input):
  opts = []
  for key in input.keys():
    opts.append(f"--data {' '.join(key.split('|'))}")

  return " ".join(opts)

rule aggregate_stats:
  input: unpack(_aggregate_stats_input)
  output: "results/stats/{feature}.txt"
  conda: "qutrna2"
  log: "logs/stats/{feature}.log"
  params: opts=_aggregate_stats_params
  shell: """
    python {workflow.basedir}/scripts/aggregate_feature.py \
      --output {output:q} \
      {params.opts} \
      {input:q} \
      2> {log:q}
  """



rule trnas_available:
  input: REF_FASTA
  output: TRNAS_AVAILABLE
  conda: "qutrna2"
  log: "logs/trnas/available.log"
  shell: """
    python {workflow.basedir}/scripts/fasta_utils.py extract-seqids {intput:q} > {output:q} 2> {log:q}
  """


rule trnas_selected:
  input: REF_FILTERED_TRNAS_FASTA
  output: TRNAS_SELECTED
  conda: "qutrna2"
  log: "logs/trnas/selected.log"
  shell: """
    python {workflow.basedir}/scripts/fasta_utils.py extract-seqids {input:q} > {output:q} 2> {log:q}
  """


if "trna_annotation" in pep.config["qutrna2"]:
  create_include("trna_annotation",
                 pep.config["qutrna2"]["trna_annotation"],
                 TRNA_ANNOTATION,
                 config["include"]["trna_annotation"])
else:
  rule infer_trna_annotation:
    input: REF_FILTERED_TRNAS_FASTA
    output: TRNA_ANNOTATION
    conda: "qutrna2"
    log: "logs/infer_trna_annotation.log"
    shell: """
      python {workflow.basedir}/scripts/fasta_utils.py infer-annotation --output {output:q} {input:q} 2> {log:q}
    """