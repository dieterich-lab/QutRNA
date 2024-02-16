REF_FASTA = "data/ref.fasta"
MODS = "data/mods.tsv"

TBL = pep.sample_table.explode("subsample_name")
SAMPLES = TBL.loc["sample_name"].to_list()
SUBSAMPLES = TBL.loc["subsample_name"].to_list()


wildcard_constraints:
  SAMPLE="|".join(SAMPLES),
  SUBSAMPLE="|".join(SUBSAMPLE),
  ORIENT="|".join(["fwd", "rev"]),
  RTYPE="|".join(["pass", "fail", "merged"]),


rule include_ref_fasta:
  input: pep.config["qutrna"]["ref_fasta"],
  output: REF_FASTA,
  params: include=config["include"]["ref_fasta"],
  run:
      if params.include == "copy":
       cmd = "cp"
      else:
       cmd = "ln -s"

      shell(cmd + " {input:q} {output:q}")


rule include_mods:
  input: pep.config["qutrna"]["mods"],
  output: MODS,
  params: include=config["include"]["mods"],
  run:
      if params.include == "copy":
       cmd = "cp"
      else:
       cmd = "ln -s"

      shell(cmd + " {input:q} {output:q}")


def _include_fnames(ftype):
  tbl = pep.sample_table.explode(["subsample_name", ftype])
  tbl = tbl.set_index("subsample_name")

  def helper(wildcards):
    return tbl.loc[wildcards.SUBSAMPLE, ftype]

  return helper


rule include_bams:
  input: _include_fnames("bam"),
  output: "data/bams/{SUBSAMPLE}.bam",
  params: include=config["include"]["bams"],
  run:
      if params.include == "copy":
       cmd = "cp"
      else:
       cmd = "ln -s"

      shell(cmd + " {input:q} {output:q}")


rule include_fastq:
  input: _include_fnames("fastq"),
  output: "data/fastq/{SUBSAMPLE}.bam",
  params: include=config["include"]["fastq"],
  run:
      if params.include == "copy":
       cmd = "cp"
      else:
       cmd = "ln -s"

      shell(cmd + " {input:q} {output:q}")


rule remove_linker:
  input: "data/ref.fasta",
  output: "results/data/no_linkers.fasta",
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
