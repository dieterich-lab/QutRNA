################################################################################
# Plot JACUSA2 score: Mis+Del+Ins
################################################################################
# TODO add tRNA family
rule plot_heatmap:
  input: "results/jacusa2/",
  output: directory("results/plots/heatmap/cond1~{COND1}/cond2~{COND2}/"),
  conda: "qutrna",
  log: "logs/plot/heatmap/cond1~{COND1}/cond2~{COND2}.log",
  shell: """
    Rscript {workflow.basedir}/plot_score.R --output {output:q} {input:q} 2> {log:q}
  """
