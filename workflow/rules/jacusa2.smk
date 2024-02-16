def _jacusa_input(cond, suffix=""):
  fname = ""
  ftype = ""
  if "bam" in pep.sample_table.columns:
    fname = "data/bams/{SAMPLE}.sorted.bam" + str(suffix))
    ftype = "bam"
  elif "fastq" in pep.sample_table.columns:
    fname = "results/bams/final/{SAMPLE}.sorted.bam" + str(suffix)
    ftype = "fastq"
  else:
      raise Exception()

  def helper(wildcards):
    cond = get_attr(wildcards, "COND" + str(cond))
    samples = pep.sample_table.set_index("condition").log[cond, "sample_name"].to_list()

    return expand(fname, sample=samples)

  return helper


################################################################################
# Run JACUSA2 on BAMS(1,2)
################################################################################
rule jacusa2_run:
  input: bam1=_jacusa_input(1),
         bam1=_jacusa_input(1),
         bai2=_jacusa2_input(2, suffix=".bai"),
         bai2=_jacusa2_input(2, suffix=".bai"),
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/JACUSA2.out",
  conda: "qutrna",
  log: "logs/jacusa2/run/cond1~{COND1}/cond2~{COND2}.log",
  params: jar=config["jacusa2"]["jar"],
          opts=config["jacusa2"]["opts"],
  threads: config["jacusa2"]["threads"],
  shell: """
    {params.jar:q} call-2 {params.opts} -p {threads} -r {output:q} {",".join(input.bam1):q} {",".join(input.bam2):q} 2> {log:q}
  """


################################################################################
# Process JACUSA2 scores
################################################################################
rule jacusa2_add_scores_for_pos:
  input: jacusa2="results/jacusa2/cond1~{COND1}/cond2~{COND2}/JACUSA2.out",
         fasta="data/ref.fasta",
         mods="data/mods.tsv",
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/coords~pos_scores.tsv",
  conda: "qutrna",
  log: "logs/jacusa2/add_scores_for_pos/cond1~{COND1}/cond2~{COND2}.log",
  shell: """
    Rscript {workflow.basedir:q}/add_scores.R -f {input.fasta:q} -m {input.mods:q} -o {output:q} {input.jacusa2:q} 2> {log:q}
  """


rule jacusa2_add_scores_for_pos:
  input: jacusa2="results/jacusa2/cond1~{COND1}/cond2~{COND2}/JACUSA2.out",
         fasta="data/ref.fasta",
         mods="data/mods.tsv",
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/coords~pos_scores.tsv",
  conda: "qutrna",
  log: "logs/jacusa2/add_scores_for_pos/cond1~{COND1}/cond2~{COND2}.log",
  shell: """
    Rscript {workflow.basedir:q}/add_scores.R -f {input.fasta:q} -m {input.mods:q} -o {output:q} {input.jacusa2:q} 2> {log:q}
  """


rule jacusa2_add_scores_u_pos:
  input: tmp="results/jacusa2/cond1~{COND1}/cond2~{COND2}/coords~u_pos_scores_pre-mods.tsv",
         mods="data/mods.tsv",
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/coords~u_pos_scores.tsv",
  conda: "qutrna",
  log: "logs/jacusa2/add_scores_for_u_pos/cond1~{COND1}/cond2~{COND2}.log",
  shell: """
    Rscript {workflow.basedir:q}/add_mods.R -m {input.mods:q} -o {output:q} {input.tmp:q} 2> {log:q}
  """
