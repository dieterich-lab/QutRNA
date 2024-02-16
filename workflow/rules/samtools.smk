rule samtools_index:
  input: "{prefix}.sorted.bam",
  output: "{prefix}.sorted.bam.bai",
  conda: "qutrna",
  log: "logs/samtools/index/{prefix}.log",
  shell: """
    samtools index {input:q} 2> {log:q}
  """


rule samtools_get_score:
  input: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}sorted.bam"
  output: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}_score.txt"
  conda: "qutrna",
  log: "logs/samtools/get_score/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}.log",
  shell: """
    (samtools view {input:q} | cut -f 5 > {output:q}) 2> {log:q}
  """


rule samtools_filter_by_cutoff:
  input: bam="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd_{RTYPE}.sorted.bam",
         cutoff="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{RTYPE}_cutoff.txt",
  output: "results/bams/final/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{RTYPE}.sorted.bam",
  conda: "qutrna",
  log: "logs/samtools/filter_by_cutoff/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{RTYPE}.log",
  shell: """
    samtools view -b -q `cat {input.cutoff:q}` \
      -o {output:q} {input.bam:q} 2> {log:q}
  """


# FIXME pass, fail, merged
def _samtools_merge_input(rtype):
  tbl =  pep.sample_table.explode("subsample_name")
  fname = "results/bams/final/sample~{sample}/subsample~{subsample}/{rtype}.sorted.bam"

  def helper(wildcards):
    subsamples = tbl.loc[wildcards.SAMPLE, "subsample_name"].to_list()

    return expand(fname, sample=wildcards.SAMPLE, subsample=subsamples, rtype=rtype)

  return helper


def _samtools_merge_input(wildcards):
  fnames = []
  for row in tbl.loc[wildcards.SAMPLE].iterrows():
    fname = "results/bams/final/sample~{sample}/subsample~{subsample}/{rtype}.sorted.bam"

  return fnames


rule samtools_merge:
  input: _samtools_merge_input,
  output: "results/bams/final/{SAMPLE}.sorted.bam",
  conda: "qutrna",
  log: "logs/samtools/merge/{SAMPLE}.log",
  shell: """
    samtools merge {output:q} {input.pass:q} {input.fail:q} 2> {log}
  """
