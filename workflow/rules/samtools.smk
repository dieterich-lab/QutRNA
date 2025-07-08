global FILTERS_APPLIED
global REF_FASTA
global READS
global TBL


rule samtools_index:
  input: "{prefix}.sorted.bam"
  output: "{prefix}.sorted.bam.bai"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools/index/{prefix}.log"
  shell: """
    samtools index {input:q} 2> "{log}"
  """


rule samtools_fa_index:
  input: "{prefix}.fa"
  output: "{prefix}.fa.fai"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools_fa_index/{prefix}.log"
  shell: """
    samtools faindex {input:q} 2> "{log}"
  """


rule samtools_stats:
  input: "{prefix}.bam"
  output: "{prefix}.stats"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools/stats/{prefix}.log"
  shell: """
    samtools stats {input:q} > {output:q} 2> "{log}"
  """


rule samtools_coverage:
  input: "{prefix}.bam"
  output: "{prefix}_coverage.tsv"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools/coverage/{prefix}.log"
  shell: """
    samtools coverage {input:q} > {output:q} 2> "{log}"
  """


rule samtools_read_count:
  input: "{prefix}.bam"
  output: "{prefix}_read_count.txt"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools/read_count/{prefix}.log"
  shell: """
    samtools view -c {input:q} > {output:q} 2> "{log}"
  """


rule samtools_multimappers:
  input: "{prefix}.bam"
  output: "{prefix}_multimappers.tsv"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools/multimappers/{prefix}.log"
  shell: """
    ( samtools view {input:q} | \
      cut -f1,3 | \
      sort -k1,1 | \
      python  {workflow.basedir}/count_multimapper.py ) > {output:q} 2> "{log}"
  """


rule samtools_get_as:
  input: "{prefix}.bam"
  output: "{prefix}_as.tsv"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools/get_as/{prefix}.log"
  shell: """
    ( python {workflow.basedir}/get_as.py ) {input:q} | sort | uniq -c ) > {output:q} 2> "{log}"
  """


# TODO
rule samtools_get_score:
  input: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sorted.bam"
  output: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_score.txt"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools/get_score/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log"
  shell: """
    (samtools view {input:q} | cut -f 5 > {output:q}) 2> "{log}"
  """


# FIXME
rule samtools_filter_by_cutoff:
  input: bam="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd/{BC}.sorted.bam",
         cutoff="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_cutoff.txt"
  output: "results/bams/final/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam"
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/samtools/filter_by_cutoff/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
  shell: """
    samtools view -b -q `cat {input.cutoff:q}` \
      -o {output:q} {input.bam:q} 2> "{log}"
  """


def _samtools_merge_input(wildcards):
  fnames = []
  for row in TBL.loc[[wildcards.SAMPLE]].itertuples():
    if READS == "fastq":
        fname = f"results/bams/final/sample~{row.sample_name}/subsample~{row.subsample_name}/{row.base_calling}.sorted.bam"
    elif READS == "bam":
      if wildcards.bam_type == "final":
        filter_ = FILTERS_APPLIED[-1]
        fname = f"results/bams/preprocessed/{filter_}/sample~{row.sample_name}/subsample~{row.subsample_name}/{row.base_calling}.sorted.bam"
      else:
        fname = f"results/bams/preprocessed/{{bam_type}}/sample~{row.sample_name}/subsample~{row.subsample_name}/{row.base_calling}.sorted.bam"
    else:
      raise Exception()
    fnames.append(fname)

  return fnames


rule samtools_merge:
  input: _samtools_merge_input
  output: "results/bams/{bam_type}/{SAMPLE}.sorted.bam"
  conda: "qutrna"
  resources:
    mem_mb=10000
  log: "logs/samtools/merge/{SAMPLE}/{bam_type}.log"
  shell: """
    samtools merge {output:q} {input:q} 2> "{log}"
  """


