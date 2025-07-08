global REF_FASTA


##############################################################################
# Use GPU-assisted alignment
##############################################################################
rule parasail_gpu_assisted:
  input: fastq="results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.fastq.gz",
         ref_fasta=REF_FASTA
  output: "results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sorted.bam"
  log: "logs/parasail/gpu_assisted/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.log"
  conda: "qutrna"
  resources:
    mem_mb=10000
  threads: config["gpu"]["threads"]
  benchmark: "benchmarks/parasail/gpu_assisted/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.tsv"
  params:
    gpu_opts=config["gpu"]["opts"],
    min_aln_score=config["min_aln_score"],
    pre=config["gpu"]["pre"]
  shell: """
    (
      {params.pre}

      gpu-tRNA-mapper \
        --readFileName {input.fastq} \
        --referenceFileName {input.ref_fasta} \
        --outputFileName /dev/stdout \
        {params.gpu_opts} \
        --minAignmentScore {params.min_aln_score} | \
      samtools view -b -o /dev/stdout /dev/stdin | \
      samtools sort -o {output} /dev/stdin
    ) 2> {log}
  """
