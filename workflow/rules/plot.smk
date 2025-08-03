from os.path import split

from snakemake.io import directory, unpack

global DEFAULT_SCORE
global REF_FASTA
global SCORES
global SPRINZL
global TRNA_ANNOTATION


################################################################################
# Plot JACUSA2 score: Mis+Del+Ins
################################################################################


PLOTS = {plot["id"]: plot for plot in config["heatmap_plots"]}


def _plot_heatmap_opts(wildcards, input):
  plot = PLOTS[wildcards.plot_id]

  opts = {}

  col_mapping = {
      "seq": "seq_position",
      "sprinzl": "sprinzl"}
  opts["--position_column"] = col_mapping[pep.config["qutrna2"]["coords"]]

  opts[f"--group"] =  plot["group"]
  if plot.get("title"):
    opts["--title"] = f"'{plot['title']}'"

  abbrevs = pep.config["qutrna2"].get("mods", {}).get("abbrevs")
  if abbrevs:
    opts["--modmap"] = abbrevs

  linker5 = str(pep.config["qutrna2"]["linker5"])
  opts["--five_adapter"] = linker5
  if "ref_fasta_prefix" in pep.config["qutrna2"]:
    opts["--remove_prefix"] = pep.config["qutrna2"]["ref_fasta_prefix"]

  # fixme
  for contrast in pep.config["qutrna2"]["contrasts"]:
    if wildcards.COND1 == contrast["cond1"] and wildcards.COND2 == contrast["cond2"] and "flag" in contrast:
      flag_positions = ",".join(contrast["flag"])
      opts["--flag_positions"] =  f"'{flag_positions}'"

  score = DEFAULT_SCORE
  if "score" in plot:
    score = plot["score"]
  opts["--score_column"] = score.split("::")[0]
  if "sprinzl" in pep.config["qutrna2"]:
    opts["--sprinzl"] = input.sprinzl # FIXME

  if plot["harmonize_scaling"] != "NONE":
    if plot["harmonize_scaling"] == "CONTRAST":
      harmonize = "-1"
    elif plot["harmonize_scaling"] == "FILTER":
      harmonize = "-2"
      opts["--max_scores"] = input.max_scores
    elif plot["harmonize_scaling"] == "CONTRAST_FILTER":
      harmonize = "-3"
      opts["--max_scores"] = input.max_scores
    else:
      harmonize = str(plot["harmonize_scaling"])
    opts["--harmonize_scaling"] = harmonize

  s = " ".join([f"{key}={value}" for key, value in opts.items()])
  if "opts" in plot and plot["opts"]:
    s = f"{s} {plot['opts']}"

  return s


def _plot_heatmap_input(wildcards):
  plot = PLOTS[wildcards.plot_id]

  d = {
    "scores": f"results/jacusa2/cond1~{{COND1}}/cond2~{{COND2}}/bam~{{bam_type}}/{SCORES}",
    "fasta": REF_FASTA,
    "trna_annotation": TRNA_ANNOTATION
  }

  if "sprinzl" in pep.config["qutrna2"]:
    d["sprinzl"] = SPRINZL

  if plot["harmonize_scaling"] == "FILTER" or plot["harmonize_scaling"] == "CONTRAST_FILTER":
    d["max_scores"] = "results/jacusa2/max_scores.tsv"

  return d


rule plot_heatmap:
  input: unpack(_plot_heatmap_input)
  output: directory("results/plots/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}")
  conda: "qutrna2"
  log: "logs/plot/heatmap/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}.log"
  params: opts=_plot_heatmap_opts,
          basedir=workflow.basedir
  shell: """
    ( export QUTRNA2="{params.basedir}"
      mkdir -p {output} && \
      Rscript {workflow.basedir}/scripts/plot_heatmap.R \
          --condition1 '{wildcards.COND1}' \
          --condition2 '{wildcards.COND2}' \
          --bam_type '{wildcards.bam_type}' \
          --ref_fasta {input.fasta:q} \
          --trna_annotation {input.trna_annotation:q} \
          --output_dir {output:q} \
          {params.opts} \
          {input.scores:q} \
          ) 2> {log:q}
  """


rule plot_record_count:
  input: "results/stats/record_count.txt"
  output: "results/plots/record_count/{type}.pdf"
  conda: "qutrna2"
  log: "logs/plot/record_count/{type}.log"
  shell: """
    Rscript {workflow.basedir}/scripts/plot_record_count.R \
         --type {wildcards.type} \
         --output {output:q} {input:q} \
         2> {log:q}
  """
rule plot_record_count_custom:
  input: "results/stats/record_count.txt"
  output: "results/plots/record_count/custom/{plot_id}.pdf"
  conda: "qutrna2"
  log: "logs/plot/read_record_custom/{plot_id}.log"
  params: opts=lambda wildcards: [plot.opts for plot in config["record_count_plots"] if plot["id"] == wildcards.plot_id][0]
  shell: """
    Rscript {workflow.basedir}/scripts/plot_record_count.R \
         {params.opts} \
         --output {output:q} {input:q} \
         2> {log:q}
  """


rule plot_read_length:
  input: "results/stats/samtools_RL.txt"
  output: "results/plots/read_length/{type}.pdf"
  conda: "qutrna2"
  log: "logs/plot/read_length/{type}.log"
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_length.R \
         --type {wildcards.type} \
         --output {output:q} {input:q} \
         2> {log:q}
  """
rule plot_read_length_custom:
  input: "results/stats/samtools_RL.txt"
  output: "results/plots/read_length/custom/{plot_id}.pdf"
  conda: "qutrna2"
  log: "logs/plot/read_length_custom/{plot_id}.log"
  params: opts=lambda wildcards: [plot.opts for plot in config["read_length_plots"] if plot["id"] == wildcards.plot_id][0]
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_length.R \
         {params.opts} \
         --output {output:q} {input:q} \
         2> {log:q}
  """


rule plot_multimapper:
  input: "results/stats/multimapper.txt"
  output: "results/plots/multimapper.pdf"
  conda: "qutrna2"
  log: "logs/plot/multimapper.log"
  shell: """
    Rscript {workflow.basedir}/scripts/plot_multimapper.R \
         --output {output:q} {input:q} \
         2> {log:q}
  """


rule plot_threshold_summary:
  input: score="results/stats/alignment_score.txt",
         cutoff="results/stats/cutoff.txt"
  output: "results/plots/alignment/threshold_summary.pdf"
  conda: "qutrna2"
  log: "logs/plot/threshold_summary.log"
  params: bam_types=",".join(["mapped", "mapped-rev"])
  shell: """
    Rscript {workflow.basedir}/scripts/plot_threshold_summary.R \
         --type {params.bam_types:q} \
         --cutoff {input.cutoff:q} \
         --output {output:q} \
         {input.score:q} \
         2> {log:q}
  """
