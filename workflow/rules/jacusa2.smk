def _jacusa2_input(cond_i, suffix=""):
  fname = f"results/bams/final/{{sample}}.sorted.bam{suffix}"
  cond = f"COND{cond_i}"

  def helper(wildcards):
    condition = wildcards[cond]
    tbl = pep.sample_table.set_index("condition")
    samples = tbl.loc[condition, "sample_name"]
    # FIXME
    if not isinstance(samples, str):
      samples = samples.to_list()
    else:
      samples = [samples, ]

    return expand(fname, sample=samples)

  return helper


################################################################################
# Run JACUSA2 on BAMS(1,2)
################################################################################
rule jacusa2_run:
  input: bam1=_jacusa2_input(1),
         bam2=_jacusa2_input(2),
         bai1=_jacusa2_input(1, suffix=".bai"),
         bai2=_jacusa2_input(2, suffix=".bai"),
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/JACUSA2.out",
  conda: "qutrna",
  resources:
    mem_mb=20000
  log: "logs/jacusa2/run/cond1~{COND1}/cond2~{COND2}.log",
  params: jar=config["jacusa2"]["jar"],
          opts=config["jacusa2"]["opts"],
          min_cov=config["jacusa2"]["min_cov"],
          bams1=lambda wildcards, input: ",".join(input.bam1),
          bams2=lambda wildcards, input: ",".join(input.bam2),
  threads: config["jacusa2"]["threads"],
  shell: """
    {params.jar} call-2 \
        {params.opts} \
        -c {params.min_cov} \
        -p {threads} \
        -r {output:q} \
        {params.bams1} \
        {params.bams2} \
        2> {log:q}
  """


################################################################################
# Process JACUSA2 scores
################################################################################
rule jacusa2_add_scores:
  input: jacusa2="results/jacusa2/cond1~{COND1}/cond2~{COND2}/JACUSA2.out",
         fasta=REF_FASTA,
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/scores.tsv",
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: "logs/jacusa2/add_scores/cond1~{COND1}/cond2~{COND2}.log",
  shell: """
    Rscript {workflow.basedir:q}/scripts/add_scores.R -f {input.fasta:q} -o {output:q} {input.jacusa2:q} 2> {log:q}
  """


rule jacusa2_add_seq_mods:
  input: scores="results/jacusa2/cond1~{COND1}/cond2~{COND2}/scores.tsv",
         mods=MODS,
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/scores_seq-mods.tsv",
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: "logs/jacusa2/add_seq_mods/cond1~{COND1}/cond2~{COND2}.log",
  shell: """
    Rscript {workflow.basedir:q}/scripts/add_mods.R \
        -m {input.mods:q} \
        -o {output:q} \
        {input.scores:q} \
        2> {log:q}
  """


rule jacusa2_add_sprinzl_mods:
  input: scores="results/jacusa2/cond1~{COND1}/cond2~{COND2}/scores_sprinzl.tsv",
         mods=MODS,
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/scores_sprinzl-mods.tsv",
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: "logs/jacusa2/add_sprinzl_mods/cond1~{COND1}/cond2~{COND2}.log",
  shell: """
    Rscript {workflow.basedir:q}/scripts/add_mods.R \
        --sprinzl \
        -m {input.mods:q} \
        -o {output:q} \
        {input.scores:q} \
        2> {log:q}
  """
