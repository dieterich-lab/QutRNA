from snakemake.io import directory, temp, glob_wildcards, expand


global REF_FASTA
global REF_FASTA_REVERSED
global READS_INPUT


##############################################################################
# Use parasail to map and align reads but split reads set in to smaller
# parts
##############################################################################

checkpoint parasail_split_reads:
  input: READS_INPUT
  output: temp(directory("results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_split"))
  log: "logs/parasail/split_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
  conda: "qutrna2"
  params: lines=config["parasail"]["lines"]
  shell: """
    mkdir -p {output}
    ( gunzip -c {input:q} | \
        split -l {params.lines} --filter 'gzip -c > $FILE.fastq.gz' /dev/stdin {output:q}/part_ ) 2> {log:q}
  """


rule parasail_map_split_reads:
  input: fastq="results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_split/part_{part}.fastq.gz",
         ref_fasta= lambda wildcards: REF_FASTA if wildcards.ORIENT == "fwd" else REF_FASTA_REVERSED
  output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.sam")
  log: "logs/parasail/map_split_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.log"
  conda: "qutrna2"
  threads: 1
  params: parasail_opts=config["parasail"]["opts"]
  shell: """
    parasail_aligner {params.parasail_opts} \
                      -t {threads} \
                      -O SAMH \
                      -f {input.ref_fasta:q} \
                      -g {output:q} \
                      -q {input.fastq:q} 2> {log:q}
  """


rule parasail_map_split_postprocess:
  input: sam="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.sam",
         ref=lambda wildcards: REF_FASTA if wildcards.ORIENT == "fwd" else REF_FASTA_REVERSED
  output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.bam")
  log: "logs/parasail/map_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.bam"
  benchmark: "benchmarks/parasail/map_split_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.tsv"
  conda: "qutrna2"
  params:
    min_aln_score=config["alignment"]["parasail"]["min_aln_score"] # FIXME
  shell: """
    (
      samtools view -b -F 4 {input.sam:q} | \
      samtools calmd /dev/stdin {input.ref:q} | \
      samtools sort -n -O bam /dev/stdin | \
      python {workflow.basedir}/scripts/bam_utils.py best-alignment --min-as {params.min_aln_score} /dev/stdin | \
      python {workflow.basedir}/scripts/bam_utils.py add-nh /dev/stdin | \
      samtools sort -O bam /dev/stdin > {output:q}
    ) 2> {log:q}
  """


def _samtools_merge_reads_input(wildcards):
  split_reads = checkpoints.parasail_split_reads.get(**wildcards).output[0]
  fnames, = glob_wildcards(
    os.path.join(
      split_reads,
      "part_{fname}.fastq.gz"))

  output_dir = "results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split"
  return expand(os.path.join(output_dir, "part_{fname}_raw.bam"),
    fname=fnames,
    allow_missing=True)


rule samtools_merge_split_reads:
  input: bams=_samtools_merge_reads_input,
    fasta=REF_FASTA
  output: bam="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.bam",
          bai="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.bam.bai"
  conda: "qutrna2"
  log: "logs/samtools/merge_split_reads/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sorted.bam"
  shell: """
    (
      samtools merge - {input.bams:q} | \
      samtools sort -o {output.bam:q} /dev/stdin && samtools index {output.bam:q}
    ) 2> {log:q}
"""
