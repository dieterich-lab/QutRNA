from snakemake.io import directory, temp, glob_wildcards


global REF_FASTA


##############################################################################
# Use parasail to map and align reads but split reads set in to smaller
# parts
##############################################################################

checkpoint parasail_split_reads:
  input: "results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.fastq.gz"
  output: temp(directory("results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split"))
  log: "logs/parasail/split_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log"
  conda: "qutrna"
  resources:
    mem_mb=2000
  params: lines=config["parasail"]["lines"]
  shell: """
    mkdir -p {output}
    ( gunzip -c {input:q} | \
        split -l {params.lines} --filter 'gzip -c > $FILE.fastq.gz' /dev/stdin {output:q}/part_ ) 2> "{log}"
  """


rule parasail_map_split_reads:
  input: fastq="results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.fastq.gz",
         ref_fasta=REF_FASTA
  output: temp("results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}_raw.bam")
  log: "logs/parasail/map_split_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.log"
  conda: "qutrna"
  resources:
    mem_mb=16000
  threads: config["parasail"]["threads"]
  params: parasail_opts=config["parasail"]["opts"]
  shell: """
  (
    parasail_aligner {params.parasail_opts} \
                      -t {threads} \
                      -O SAMH \
                      -f {input.ref_fasta:q} \
                      -g /dev/stdout \
                      -q {input.fastq:q}
      samtools view -bS /dev/stdin | \
      samtools calmd --output-fmt BAM /dev/stdin {input.ref_fasta:q} > {output:q} \
    ) 2> "{log}"
  """


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
    fasta=REF_FASTA
  output: bam="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.bam",
    bai="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.bam.bai"
  conda: "qutrna"
  resources:
    mem_mb=10000
  log: "logs/samtools/merge_split_reads/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sorted.bam"
  shell: """
  ( samtools merge - {input.bams:q} | \
      samtools sort -o {output.bam:q} /dev/stdin && samtools index {output.bam:q} ) 2> "{log}"
"""
