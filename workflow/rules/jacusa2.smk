import pandas as pd
import re


global DEFAULT_SCORE
global FILTERS_APPLIED
global MODS
global REF_FASTA
global SCORES


def _jacusa2_input(cond_i, suffix=""):
  fname = f"results/bams/{{bam_type}}/{{sample}}.sorted.bam{suffix}"
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

    return expand(fname, sample=samples, allow_missing=True)

  return helper


################################################################################
# Run JACUSA2 on BAMS(1 vs. 2)
################################################################################
rule jacusa2_run:
  input: bam1=_jacusa2_input(1),
         bam2=_jacusa2_input(2),
         bai1=_jacusa2_input(1, suffix=".bai"),
         bai2=_jacusa2_input(2, suffix=".bai"),
         ref=REF_FASTA
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/bams~{bam_type}/JACUSA2.out"
  conda: "qutrna"
  resources:
    mem_mb=20000
  log: "logs/jacusa2/run/cond1~{COND1}/cond2~{COND2}/bams~{bam_type}.log"
  params: jar=config["jacusa2"]["jar"],
          opts=config["jacusa2"]["opts"],
          min_cov=config["jacusa2"]["min_cov"],
          bams1=lambda wildcards, input: ",".join(input.bam1),
          bams2=lambda wildcards, input: ",".join(input.bam2)
  threads: config["jacusa2"]["threads"]
  shell: """
    {params.jar} call-2 \
        -R {input.ref} \
        {params.opts} \
        -c {params.min_cov} \
        -p {threads} \
        -r {output:q} \
        {params.bams1} \
        {params.bams2} \
        2> {log}
  """


################################################################################
# Process JACUSA2 scores
################################################################################
def _jacusa_add_scores(_):
  scores = {plot.get("score", DEFAULT_SCORE) for plot in config["plots"]}

  return "-s " + ",".join(scores)

rule jacusa2_add_scores:
  input: jacusa2="results/jacusa2/cond1~{COND1}/cond2~{COND2}/bams~{bam_type}/JACUSA2.out"
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/bams~{bam_type}/scores_seq.tsv"
  conda: "qutrna"
  resources:
    mem_mb=10000
  log: "logs/jacusa2/add_scores/cond1~{COND1}/cond2~{COND2}/bams~{bam_type}.log"
  params: stats=_jacusa_add_scores,
          remove_exon_only=config.get("_remove_exon_only", False)
  shell: """
    Rscript {workflow.basedir:q}/scripts/add_scores.R --remove {params.remove_exon_only} {params.stats} -o {output:q} {input.jacusa2:q} 2> "{log}"
  """


rule jacusa2_add_mods:
  input: scores="results/jacusa2/cond1~{COND1}/cond2~{COND2}/bams~{bam_type}/scores_{coord_type}.tsv",
         mods=MODS
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/bams~{bam_type}/scores_{coord_type}-mods.tsv"
  conda: "qutrna"
  resources:
    mem_mb=10000
  params:
    coords=lambda wildcards: "--sprinzl" if wildcards.coord_type == "sprinzl" else ""
  log: "logs/jacusa2/add_seq_mods/cond1~{COND1}/cond2~{COND2}/bams~{bam_type}/{coord_type}.log"
  shell: """
    Rscript {workflow.basedir:q}/scripts/add_mods.R \
        {params.coords} \
        -m {input.mods:q} \
        -o {output:q} \
        {input.scores:q} \
        2> "{log}"
  """


def _input_jacusa2_max_scores(_):
  targets = []

  bam_types = ["final",]
  if config["plot_preprocessed"] and FILTERS_APPLIED:
    bam_types += FILTERS_APPLIED
  for conds in pep.config["qutrna"]["contrasts"]:
    fname =  "results/jacusa2/cond1~{cond1}/cond2~{cond2}/bams~{bam_type}/" + SCORES
    targets.extend(expand(fname,
               cond1=conds["cond1"], cond2=conds["cond2"],
               bam_type=bam_types))

  return targets


def _params_config_scores():
  scores = [plot["score"] for plot in config["plots"]]

  return list(set([score.split("::")[0] for score in scores]))


rule jacusa2_max_scores:
  input: _input_jacusa2_max_scores
  output: "results/jacusa2/max_scores.tsv"
  params: scores=_params_config_scores()
  run:
    dfs = []
    for f in input:
      df = pd.read_csv(f, sep="\t")[params.scores].apply(max).to_frame().T
      cond1 = re.search("/cond1~([^/]+)/", f).group(1)
      cond2 = re.search("/cond2~([^/]+)/", f).group(1)
      bam_type = re.search("/bams~([^/]+)/", f).group(1)
      df["cond1"] = cond1
      df["cond2"] = cond2
      df["bam_type"] = bam_type
      dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(output[0], sep="\t", index=False)
