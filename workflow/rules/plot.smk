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
              col_mapping[pep.config["qutrna"]["coords"]])

  opts.append("--split " + plot["trnas"])
  if plot.get("title"):
    opts.append("--title '" + plot["title"] + "'")
  if "opts" in plot:
    opts.append(plot["opts"])

  abbrevs = pep.config["qutrna"].get("mods", {}).get("abbrevs")
  if abbrevs:
    opts.append("--modmap " + abbrevs)

  opts.append("--five_adapter " + str(pep.config["qutrna"]["linker5"]))
  if "ref_fasta_prefix" in pep.config["qutrna"]:
    opts.append("--remove_prefix=" + pep.config["qutrna"]["ref_fasta_prefix"])

  for contrast in pep.config["qutrna"]["contrasts"]:
    if wildcards.COND1 == contrast["cond1"] and wildcards.COND2 == contrast["cond2"] and "flag" in contrast:
      opts.append("--flag_positions=" + ",".join(contrast["flag"]))

  # TODO
  if "coverages" in pep.config["qutrna"]:
    opts.append("--coverage_info=" + input.coverages)

  if "score" in plot:
    opts.append("--score_column=" + plot["score"].split("::")[0])
  else:
    opts.append("--score_column=" + DEFAULT_SCORE.split("::")[0])

  if "sprinzl" in pep.config["qutrna"]:
    opts.append("--sprinzl=" + input.sprinzl)

  return " ".join(opts)


def _plot_heatmap_input(wildcards):
  d = {"scores": "results/jacusa2/cond1~{COND1}/cond2~{COND2}/" + SCORES,
       "fasta": REF_FASTA,}
  # FIXME
  if "coverages" in pep.config["qutrna"]:
    d["coverages"] = pep.config["qutrna"]["coverages"]

  if "sprinzl" in pep.config["qutrna"]:
    d["sprinzl"] = SPRINZL

  return d


rule plot_heatmap:
  input: unpack(_plot_heatmap_input)
  output: directory("results/plots/cond1~{COND1}/cond2~{COND2}/{plot_id}"),
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: "logs/plot/heatmap/cond1~{COND1}/cond2~{COND2}/{plot_id}.log",
  params: opts=_plot_heatmap_opts,
  shell: """
    ( mkdir -p {output} && \
      Rscript {workflow.basedir}/scripts/plot_score.R \
          --condition1 {wildcards.COND1:q} \
          --condition2 {wildcards.COND2:q} \
          --ref_fasta {input.fasta:q} \
          --output_dir {output:q} {input.scores:q} \
          {params.opts} ) 2> {log:q}
  """


rule plot_read_counts:
  input: "results/read_counts.tsv",
  output: "results/plots/read_counts.pdf",
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: "logs/plot/read_counts.log",
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_counts.R \
         --output {output:q} {input:q} \
         2> {log:q}
  """


rule plot_read_length:
  input: "results/samtools/stats/RL.tsv",
  output: "results/plots/read_length.pdf",
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: "logs/plot/read_length.log",
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_length.R \
         --output {output:q} {input:q} \
         2> {log:q}
  """


rule plot_multimapper:
  input: "results/stats/multimapper.tsv",
  output: "results/plots/multimapper.pdf",
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: "logs/plot/multimapper.log",
  shell: """
    Rscript {workflow.basedir}/scripts/plot_multimapper.R \
         --output {output:q} {input:q} \
         2> {log:q}
  """
