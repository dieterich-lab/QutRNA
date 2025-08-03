global REF_FASTA


##############################################################################
# Use GPU-assisted or only parasail to map reads
##############################################################################
if config["alignment"]["method"] == "gpu":
  include: "parasail_gpu.smk"
else:
  if config["parasail"]["lines"] > 0:
    include: "parasail_split.smk"
  else:
    include: "parasail_map.smk"


##############################################################################
# Filter alignments by random score distribution
##############################################################################

rule parasail_infer_cutoff:
  input: fwd="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd/{BC}_stats/alignment_score.txt",
         rev="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~rev/{BC}_stats/alignment_score.txt"
  output: score_plot="results/plots/alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}/alignment_score.pdf",
          cutoff="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_stats/cutoff.txt"
  params: precision=config["alignment"]["precision"]
  log: "logs/parasail/infer_cutoff/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
  conda: "qutrna2"
  shell: """
    Rscript --vanilla {workflow.basedir}/scripts/alignment_score_cutoff.R \
      -S {output.score_plot:q} \
      -C {output.cutoff:q} \
      -p {params.precision:q} \
      --forward {input.fwd:q} --reverse {input.rev:q} 2> {log:q}
  """
