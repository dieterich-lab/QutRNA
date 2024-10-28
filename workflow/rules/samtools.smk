# FIXME
#__SAMTOOLS_PROCESS_BAM_INPUT = {
#    "bam": "data/bams/raw/{filename}.bam",
#    "bai": "data/bams/raw/{filename}.bam.bai",
#  }
#if config["process"]["calmd"]:
#  d["ref"] = REF_FASTA
#
#rule samtools_process_bam:
#  input: __SAMTOOLS_PROCESS_BAM_INPUT,
#  output: "results/data/bams/processed/{filename}.bam",
#  params:
#    filter=config["process"]["filter"],
#    calmd=config["process"]["calmd"],
#  log: join_path("logs/samtools/processed/{filename}.log")
#  run:
#    if not params.filter and not params.calmd:
#      cmd = "ln -s {input.bam} {output}"
#      shell("( " + cmd + " > {output} ) 2> {log}")
#    else:
#      filter_cmd = "samtools view {params.filter} -b {input.bam}"
#
#      cmds = [filter_cmd, ]
#      if params.calmd:
#        calmd_cmd = "samtools calmd -b /dev/stdin {input.ref}"
#        cmds.append(calmd_cmd)
#      cmd = " | ".join(cmds)
#      shell("( " + cmd + " > {output} ) 2> {log}")


rule samtools_index:
  input: "{prefix}.sorted.bam",
  output: "{prefix}.sorted.bam.bai",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: "logs/samtools/index/{prefix}.log",
  shell: """
    samtools index {input:q} 2> {log:q}
  """


rule samtools_coverage:
  input: "{prefix}.bam",
  output: "{prefix}_coverage.tsv",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: "logs/samtools/coverage/{prefix}.log",
  shell: """
    samtools coverage {input:q} > {output:q} 2> {log:q}
  """


# TODO
#rule merge_coverage:
#  input: ""
#  output: "coverage.tsv"
#  conda: "qutrna",
#  resources:
#    mem_mb=4000
#  log: "logs/samtools/merge_coverage.log",
#  shell: """
#    python {workflow.basedir}/merge_coverage.py > {output:q} 2> {log:q}
#  """


rule samtools_get_score:
  input: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sorted.bam"
  output: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_score.txt"
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: "logs/samtools/get_score/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log",
  shell: """
    (samtools view {input:q} | cut -f 5 > {output:q}) 2> {log:q}
  """


rule samtools_filter_by_cutoff:
  input: bam="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd/{BC}.sorted.bam",
         cutoff="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_cutoff.txt",
  output: "results/bams/final/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: "logs/samtools/filter_by_cutoff/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log",
  shell: """
    samtools view -b -q `cat {input.cutoff:q}` \
      -o {output:q} {input.bam:q} 2> {log:q}
  """


def _samtools_merge_input(wildcards):
  fnames = []
  for row in TBL.loc[[wildcards.SAMPLE]].itertuples():
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
  resources:
    mem_mb=10000
  log: "logs/samtools/merge/{SAMPLE}.log",
  shell: """
    samtools merge {output:q} {input:q} 2> {log:q}
  """


if config["parasail"]["lines"] > 0:
  def _samtools_merge_reads_input(wildcards):
    split_reads = checkpoints.parasail_split_reads.get(**wildcards).output[0]
    fnames, = glob_wildcards(
        os.path.join(
          split_reads,
          "part_{fname}.fastq.gz"))

    output_dir = "results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split"
    return expand(os.path.join(output_dir, "part_{fname}_raw.bam"),
                  fname=fnames,
                  allow_missing=True)


  rule samtools_merge_split_reads:
    input: bams=_samtools_merge_reads_input,
           fasta=REF_FASTA,
    output: bam="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.bam",
            bai="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.bam.bai",
    conda: "qutrna",
    resources:
      mem_mb=10000
    log: "logs/samtools/merge_split_reads/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sorted.bam",
    shell: """
      ( samtools merge - {input.bams} | \
          samtools sort -o {output.bam:q} /dev/stdin && samtools index {output.bam:q} ) 2> {log:q}
    """
