from snakemake.io import directory, unpack

global DEFAULT_SCORE
global REF_FASTA
global SCORES
global SPRINZL_LABELS
global TRNA_ANNOTATION
global MAX_SCORES


################################################################################
# Plot JACUSA2 score: Mis+Del+Ins
################################################################################


def get_plot(plot_type, plot_id):
  plots = config["plots"][plot_type]
  for plot in plots:
    if plot["id"] == plot_id:
      return plot

  return {}


def _plot_heatmap_opts(wildcards, input):
  plot = get_plot("heatmap", wildcards.plot_id)

  opts = {}

  col_mapping = {
      "seq": "seq_position",
      "sprinzl": "sprinzl"}
  opts["--position_column"] = col_mapping[pep.config["qutrna2"]["coords"]]

  if "group" in plot:
    opts[f"--group"] =  plot["group"]
  if plot.get("title"):
    opts["--title"] = f"'{plot['title']}'"

  linker5 = str(pep.config["qutrna2"]["linker5"])
  opts["--five_adapter"] = linker5

  for contrast in pep.config["qutrna2"]["contrasts"]:
    if wildcards.COND1 == contrast["cond1"] and wildcards.COND2 == contrast["cond2"] and "flag" in contrast:
      flag_positions = ",".join(contrast["flag"])
      opts["--flag_positions"] =  f"'{flag_positions}'"

  score = DEFAULT_SCORE
  if "score" in plot:
    score = plot.get("score", DEFAULT_SCORE)
  opts["--score_column"] = score.split("::")[0]
  if "sprinzl" in pep.config["qutrna2"]:
    opts["--sprinzl"] = input.sprinzl

  if "abbrevs" in pep.config["qutrna2"].get("mods", {}):
    opts["--abbrevs"] = input.abbrevs

  if plot.get("harmonize_scaling"):
    if plot["harmonize_scaling"] in ["CONTRASTS", "OVERALL"]:
      opts["--max_scores"] = input.max_scores
    opts["--harmonize_scaling"] = plot["harmonize_scaling"]

  s = " ".join([f"{key}={value}" for key, value in opts.items()])
  if "opts" in plot and plot["opts"]:
    s = f"{s} {plot['opts']}"

  return s


def _plot_heatmap_input(wildcards):
  plot = get_plot("heatmap", wildcards.plot_id)
  input = {
    "scores": f"results/jacusa2/cond1~{{COND1}}/cond2~{{COND2}}/bam~{{bam_type}}/{SCORES}",
    "trna_annotation": TRNA_ANNOTATION
  }

  if "sprinzl" in pep.config["qutrna2"]:
    input["sprinzl"] = SPRINZL_LABELS

  if "abbrevs" in pep.config["qutrna2"].get("mods",{}):
    global MOD_ABBREVS

    input["abbrevs"] = MOD_ABBREVS

  if plot.get("harmonize_scaling", "") in ["CONTRASTS", "OVERALL"]:
    input["max_scores"] = MAX_SCORES

  return input


rule plot_heatmap:
  input: unpack(_plot_heatmap_input)
  output: dir=directory("results/plots/scores/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}")
  conda: "qutrna2"
  log: "logs/plot/heatmap/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}.log"
  params: opts=_plot_heatmap_opts,
          basedir=workflow.basedir
  shell: """
    ( export QUTRNA2="{params.basedir}"
      mkdir -p {output.dir} && \
      Rscript {workflow.basedir}/scripts/plot_heatmap.R \
          --condition1 '{wildcards.COND1}' \
          --condition2 '{wildcards.COND2}' \
          --bam_type '{wildcards.bam_type}' \
          --trna_annotation {input.trna_annotation:q} \
          --output_dir {output.dir:q} \
          {params.opts} \
          {input.scores:q} \
          ) 2> {log:q}
  """


rule plot_record_count:
  input: "results/stats/record_count.txt"
  output: "results/plots/record_count/{FEATURE}.pdf"
  conda: "qutrna2"
  log: "logs/plot/record_count/{FEATURE}.log"
  params: basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_record_count.R \
         --type {wildcards.FEATURE} \
         --output {output:q} {input:q} \
         2> {log:q}
  """
rule plot_record_count_custom:
  input: "results/stats/record_count.txt"
  output: "results/plots/record_count/custom/{plot_id}.pdf"
  conda: "qutrna2"
  log: "logs/plot/read_record_custom/{plot_id}.log"
  params: opts= lambda wildcards: get_plot("record_count", wildcards.plot_id).get("opts",""),
          basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_record_count.R \
         {params.opts} \
         --output {output:q} {input:q} \
         2> {log:q}
  """


rule plot_read_length:
  input: "results/stats/read_length.txt"
  output: "results/plots/read_length/{FEATURE}.pdf"
  conda: "qutrna2"
  log: "logs/plot/read_length/{FEATURE}.log"
  params: basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_read_length.R \
         --type {wildcards.FEATURE} \
         --output {output:q} {input:q} \
         2> {log:q}
  """
rule plot_read_length_custom:
  input: "results/stats/read_length.txt"
  output: "results/plots/read_length/custom/{plot_id}.pdf"
  conda: "qutrna2"
  log: "logs/plot/read_length_custom/{plot_id}.log"
  params: opts=lambda wildcards: get_plot("read_length", wildcards.plot_id).get("opts", ""),
          basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_read_length.R \
         {params.opts} \
         --output {output:q} {input:q} \
         2> {log:q}
  """


rule plot_threshold_summary:
  input: score="results/stats/alignment_score.txt",
         cutoff="results/stats/cutoff.txt"
  output: "results/plots/alignment/threshold_summary.pdf"
  conda: "qutrna2"
  log: "logs/plot/threshold_summary.log"
  params: bam_types=",".join(["mapped", "mapped-random"]),
          basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_threshold_summary.R \
         --type {params.bam_types:q} \
         --cutoff {input.cutoff:q} \
         --output {output:q} \
         {input.score:q} \
         2> {log:q}
  """
rule plot_threshold_summary_custom:
  input: score="results/stats/alignment_score.txt",
    cutoff="results/stats/cutoff.txt"
  output: "results/plots/alignment/threshold_summary/custom/{plot_id}.pdf"
  conda: "qutrna2"
  log: "logs/plot/threshold_summary_custom/{plot_id}.log"
  params: opts=lambda wildcards: get_plot("threshold_summary", wildcards.plot_id).get("opts", ""),
          bam_types=",".join(["mapped", "mapped-random"]),
          basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_threshold_summary.R \
         {params.opts} \
         --type {params.bam_types:q} \
         --cutoff {input.cutoff:q} \
         --output {output:q} \
         {input.score:q} \
         2> {log:q}
  """
