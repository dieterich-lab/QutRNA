rule samtools_index:
  input: "{prefix}.sorted.bam",
  output: "{prefix}.sorted.bam.bai",
  conda: "qutrna",
  log: "logs/samtools/index/{prefix}.log",
  shell: """
    samtools index {input:q} 2> {log:q}
  """


rule samtools_get_score:
  input: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.sorted.bam"
  output: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}_score.txt"
  conda: "qutrna",
  log: "logs/samtools/get_score/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.log",
  shell: """
    (samtools view {input:q} | cut -f 5 > {output:q}) 2> {log:q}
  """


rule samtools_filter_by_cutoff:
  input: bam="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd_{BC}.sorted.bam",
         cutoff="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_cutoff.txt",
  output: "results/bams/final/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam",
  conda: "qutrna",
  log: "logs/samtools/filter_by_cutoff/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log",
  shell: """
    samtools view -b -q `cat {input.cutoff:q}` \
      -o {output:q} {input.bam:q} 2> {log:q}
  """


def _samtools_merge_input(wildcards):
  fnames = []
  for row in pep.sample_table.loc[[wildcards.SAMPLE]].itertuples():
    if READS == "fastq":
        fname = f"results/bams/final/sample~{row.sample_name}/subsample~{row.subsample_name}/{row.base_calling}.sorted.bam"
    elif READS == "bam":
        fname = f"data/bams/sample~{row.sample_name}/subsample~{row.subsample_name}/{row.base_calling}.sorted.bam"
    else:
      raise Exception()
    fnames.append(fname)

  return fnames


rule samtools_merge:
  input: _samtools_merge_input,
  output: "results/bams/final/{SAMPLE}.sorted.bam",
  conda: "qutrna",
  log: "logs/samtools/merge/{SAMPLE}.log",
  shell: """
    samtools merge {output:q} {input:q} {input:q} 2> {log:q}
  """


def _samtools_merge_reads_input(wildcards):
  output_dir = checkpoints.parasail_split_reads.get(**wildcards).output[0]
  fnames, = glob_wildcards(os.path.join(ck_output, "part_{fname}.fastq.gz"))

  return expand(os.path.join(output_dir, "part_{fname}.fastq.gz"), fname=fnames)


rule samtools_merge_split_reads:
  input: bams=_samtools_merge_reads_input,
         fasta=REF_FASTA,
  output: bam="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.sorted.bam",
          bai="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.sorted.bam.bai",
  conda: "qutrna",
  log: "logs/samtools/merge_split_reads/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.sorted.bam",
  threads: 4
  shell: """
    ( samtools merge - {input.bams} \
        samtools calmd /dev/stdin {input.fasta:q} | \
        samtools sort -@ {threads}-o {output.bam:q} /dev/stdin && samtools index {output.bam:q} ) 2> {log:q}
  """
