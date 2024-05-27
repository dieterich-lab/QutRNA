################################################################################
# Plot JACUSA2 score: Mis+Del+Ins
################################################################################

def _plot_heatmap_opts(wildcards):
  opts = []

  plots = {plot["id"]: plot for plot in config["plots"]}
  plot = plots[wildcards.plot_id]

  col_mapping = {
      "seq": "Pos3",
      "sprinzl": "sprinzl"}

  opts.append("--column " +
              col_mapping[pep.config["qutrna"]["coords"]])

  opts.append("--split " + plot["trnas"])
  if plot.get("title"):
    opts.append("--title '" + plot["title"] + "'")
  if "opts" in plot:
    opts.append(plot["opts"])

  abbrevs = pep.config["qutrna"].get("mods", {}).get("abbrevs")
  if abbrevs:
    opts.append("--modmap " + abbrevs)

  opts.append("--left " + str(pep.config["qutrna"]["linker5"]))
  if "ref_fasta_prefix" in pep.config["qutrna"]:
    opts.append("--remove_prefix=" + pep.config["qutrna"]["ref_fasta_prefix"])

  for contrast in pep.config["qutrna"]["contrasts"]:
    if wildcards.COND1 == contrast["cond1"] and wildcards.COND2 == contrast["cond2"] and "flag" in contrast:
      opts.append("--flag=" + ",".join(contrast["flag"]))

  return " ".join(opts)


rule plot_heatmap:
  input: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/" + SCORES,
  output: directory("results/plots/cond1~{COND1}/cond2~{COND2}/{plot_id}"),
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: "logs/plot/heatmap/cond1~{COND1}/cond2~{COND2}/{plot_id}.log",
  params: opts=_plot_heatmap_opts,
  shell: """
    ( mkdir -p {output} && \
      Rscript {workflow.basedir}/scripts/plot_score.R \
          --cond1 {wildcards.COND1:q} \
          --cond2 {wildcards.COND2:q} \
          --output {output:q} {input:q} \
          {params.opts} ) 2> {log:q}
  """
