from snakemake.io import directory, unpack

global DEFAULT_SCORE
global REF_FASTA
global SCORES
global SPRINZL


################################################################################
# Plot JACUSA2 score: Mis+Del+Ins
################################################################################

def _plot_heatmap_opts(wildcards, input):
  opts = []

  plots = {plot["id"]: plot for plot in config["plots"]}
  plot = plots[wildcards.plot_id]

  col_mapping = {
      "seq": "seq_position",
      "sprinzl": "sprinzl"}

  opts.append("--position_column " +
              col_mapping[pep.config["qutrna2"]["coords"]])

  opts.append("--split " + plot["trnas"])
  if plot.get("title"):
    opts.append("--title '" + plot["title"] + "'")
  if "opts" in plot:
    opts.append(plot["opts"])

  abbrevs = pep.config["qutrna2"].get("mods", {}).get("abbrevs")
  if abbrevs:
    opts.append("--modmap " + abbrevs)

  opts.append("--five_adapter " + str(pep.config["qutrna2"]["linker5"]))
  if "ref_fasta_prefix" in pep.config["qutrna2"]:
    opts.append("--remove_prefix=" + pep.config["qutrna2"]["ref_fasta_prefix"])

  for contrast in pep.config["qutrna2"]["contrasts"]:
    if wildcards.COND1 == contrast["cond1"] and wildcards.COND2 == contrast["cond2"] and "flag" in contrast:
      opts.append("--flag_positions=" + ",".join(contrast["flag"]))

  if "score" in plot:
    opts.append("--score_column=" + plot["score"].split("::")[0])
  else:
    opts.append("--score_column=" + DEFAULT_SCORE.split("::")[0])

  if "sprinzl" in pep.config["qutrna2"]:
    opts.append("--sprinzl=" + input.sprinzl)

  return " ".join(opts)


def _plot_heatmap_input(_):
  d = {"scores": "results/jacusa2/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/" + SCORES,
       "fasta": REF_FASTA,}

  if "sprinzl" in pep.config["qutrna2"]:
    d["sprinzl"] = SPRINZL

  d["max_scores"] = "results/jacusa2/max_scores.tsv"

  return d

rule plot_heatmap:
  input: unpack(_plot_heatmap_input)
  output: directory("results/plots/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}")
  conda: "qutrna2"
  log: "logs/plot/heatmap/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}.log"
  params: opts=_plot_heatmap_opts
  shell: """
    ( mkdir -p {output} && \
      Rscript {workflow.basedir}/scripts/plot_score.R \
          --condition1 {wildcards.COND1:q} \
          --condition2 {wildcards.COND2:q} \
          --bam_type {wildcards.bam_type:q} \
          --max_scores {input.max_scores:q} \
          --ref_fasta {input.fasta:q} \
          --output_dir {output:q} {input.scores:q} \
          {params.opts} ) 2> {log:q}
  """


# TODO precision-recal cutoff alignment score

rule plot_read_count:
  input: "results/stats/read_count.txt"
  output: "results/plots/read_count/{type}.pdf"
  conda: "qutrna2"
  log: "logs/plot/read_count/{type}.log"
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_count.R \
         --type {wildcards.type} \
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
         cutoff="results/stats/cutoff.txt",
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
