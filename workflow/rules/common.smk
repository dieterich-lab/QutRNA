from snakemake import shell
import itertools

wildcard_constraints:
  COND1="[^/~]+",
  COND2="[^/~]+",
  BAM_TYPE="[^/]+",
  SAMPLE="[^~/]+",
  SUBSAMPLE="[^~/]+",
  ALIGNMENT="[^~/]+",
  BC="(pass|fail|unknown)", #[^~/]+",
  FEATURE="[^~/\\.]+"


# Default JACUSA2 score to use: label::score1+score2+... see add_scores.R
DEFAULT_SCORE = "MDI_subsampled::norm_mismatch_score_subsampled+norm_deletion_score_subsampled+norm_insertion_score_subsampled"


########################################################################################################################

# Helper function to create include rules that will copy or link files into workspace
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



########################################################################################################################
# Include and optionally transform reference

# Path of raw reference sequence - NO CHANGES
REF_FASTA = "data/ref.fasta"
create_include("ref_fasta",
  pep.config["qutrna2"]["ref_fasta"]["file"],
  REF_FASTA,
  config["include"]["ref_fasta"])
# Reference sequence with linker removed
REF_NO_LINKER_FASTA = "results/data/ref_no_linker.fasta"
# Reference sequence without linker and trnas removed
REF_FILTERED_TRNAS_FASTA = "results/data/ref_filtered_trnas.fasta"


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


########################################################################################################################
# Include reads: bam or fastq
# transform reads on demand.

# include bams or fastq from sample table
def get_read_fnames(ftype):
  def helper(wildcards):
    tbl = TBL.loc[[wildcards.SAMPLE]]
    j = (tbl.subsample_name == wildcards.SUBSAMPLE) & (tbl.base_calling == wildcards.BC)
    if sum(j.to_list()) != 1:
      raise Exception("Houston, we have a problem!")

    return tbl.loc[j, ftype].to_list()[0]

  return helper

# FASTQ or BAM
if "bam" in pep.sample_table.columns and "fastq" in pep.sample_table.columns:
  raise Exception("Only column 'bam' or 'fastq' can be set in sample table.")

# global variable to indicate what type of reads are used: bam or fastq
READS = "" # will be bam or fastq
if "bam" in pep.sample_table.columns:
  READS = "bam"
  BAM_READS_INPUT = "data/bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam"

  create_include("bams",
    get_read_fnames("bam"),
    BAM_READS_INPUT,
    config["include"]["bam"])

  # optionally transform BAM files
  if "bam" in config["transform"]:
    base_change = config["transform"]["bam"].get("base_change")
    reverse_seq = config["transform"]["bam"].get("reverse")
    if base_change or reverse_seq:
      _old_reads_input = BAM_READS_INPUT
      BAM_READS_INPUT = "results/bam/transformed/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam",

      rule transform_bam:
        input: _old_reads_input
        # noinspection SmkNotSameWildcardsSet
        output: BAM_READS_INPUT
        params: base_change=lambda wildcards: "--base-change {base_change} " if base_change else "",
                reverse=lambda wildcards: "--reverse " if reverse_seq else ""
        conda: "qutrna2"
        log: "logs/transform_bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
        shell: """
          python {workflow.basedir}/scripts/bam_utils.py transform {params.reverse} {params.base_change} --output {output:q} {input:q} 2> {log:q}
        """
elif "fastq" in pep.sample_table.columns:
  READS = "fastq"
  FASTQ_READS_INPUT = "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz"
  BAM_READS_INPUT = "results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~real/{BC}.sorted.bam"
  REF_FASTA_RANDOM = "results/data/ref_random.fasta"

  create_include("fastq",
    get_read_fnames("fastq"),
    FASTQ_READS_INPUT,
    config["include"]["fastq"])

  # optionally transform FASTQ files
  if config["transform"].get("fastq", {}).get("base_change"):
    _old_reads_input = "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz"
    FASTQ_READS_INPUT = "results/fastq/transformed/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz"

    rule transform_fastq:
      input: _old_reads_input
      # noinspection SmkNotSameWildcardsSet
      output: FASTQ_READS_INPUT
      log: "logs/transform_fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
      params: base_change=config["transform"]["fastq"]["base_change"]
      conda: "qutrna2"
      shell: """
        python {workflow.basedir}/scripts/fastq_utils.py transform --output {output:q} --base-change {params.base_change} {input:q} 2> {log:q}
      """

  # How to calculate random alignment: rev_linker+trna, rev_trna
  if config["alignment"]["random"] == "rev_linker+trna":
    rule fasta_reverse_linker_trna:
      input: REF_FASTA
      output: REF_FASTA_RANDOM
      log: "logs/fasta_reverse.log"
      conda: "qutrna2"
      shell: """
          python {workflow.basedir}/scripts/fasta_utils.py transform --output {output:q} --reverse {input:q} 2> {log:q}
        """
  elif config["alignment"]["random"] == "rev_trna":
    rule fasta_reverse_trna_only:
      input: REF_FASTA
      output: REF_FASTA_RANDOM
      log: "logs/fasta_reverse_trna_only.log"
      conda: "qutrna2"
      params: linker5=pep.config["qutrna2"]["linker5"],
              linker3=pep.config["qutrna2"]["linker3"]
      shell: """
          python {workflow.basedir}/scripts/fasta_utils.py reverse \
            --linker5 {params.linker5} \
            --linker3 {params.linker3} \
            --output {output:q} \
            {input:q} 2> {log:q}
        """
  else:
    raise Exception("Invalid value for config['alignment']['random']!")


########################################################################################################################

__nested_cols = []
TBL = pep.sample_table
CONDITIONS = [[contrast["cond1"], contrast["cond2"]] for contrast in pep.config["qutrna2"]["contrasts"]]
CONDITIONS = set(itertools.chain(*CONDITIONS))
TBL = TBL[TBL["condition"].isin(CONDITIONS)]
for c in ["subsample_name", "base_calling", READS]:
  i = [isinstance(o, list) for o in TBL[c]]
  if any(i):
    __nested_cols.append(c)
if __nested_cols:
  TBL = TBL.explode(__nested_cols)
SAMPLES = TBL["sample_name"].unique().tolist()
SUBSAMPLES = TBL["subsample_name"].tolist()

########################################################################################################################
# Define JACUSA2 scores file name

# seq or sprinzl coordinates and use optional mods annotation
SCORES = "scores"
if pep.config["qutrna2"]["coords"] == "sprinzl":
  SCORES += "_sprinzl"
else:
  SCORES += "_seq"

if "mods" in pep.config["qutrna2"]:
  # Add mods to JACUSA2 scores
  SCORES += "-mods"
SCORES += ".tsv"

########################################################################################################################
# Optionally, include mods and mod abbreviations

if "mods" in pep.config["qutrna2"]:
  # Internal path for mods
  MODS = "data/mods.tsv"

  # Include mods file to workspace
  create_include("mods",
    pep.config["qutrna2"]["mods"]["file"],
    MODS,
    config["include"]["mods"]["file"])

  # Include mod abbreviations, if available
  if "abbrevs" in pep.config["qutrna2"]["mods"]:
    MOD_ABBREVS = "data/mod_abbrevs.tsv"
    create_include("abbrevs",
      pep.config["qutrna2"]["mods"]["abbrevs"],
      MOD_ABBREVS,
      config["include"]["mods"]["abbrevs"])
