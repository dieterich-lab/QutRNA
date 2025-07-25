from snakemake.io import temp


global REF_FASTA
global REF_FASTA_REVERSED
global READS_INPUT


##############################################################################
# Use GPU-assisted alignment
##############################################################################
rule parasail_gpu_assisted:
  input: fastq=READS_INPUT,
         ref_fasta=lambda wildcards: REF_FASTA if wildcards.ORIENT == "fwd" else REF_FASTA_REVERSED
  output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sam")
  log: "logs/parasail/gpu_assisted/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log"
  benchmark: "benchmarks/parasail/gpu_assisted/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.txt"
  threads: 1
  params:
    gpu_opts=config["gpu"]["opts"],
    min_aln_score=config["params"]["min_aln_score"],
    pre=config["gpu"]["pre"]
  shell: """
    (
      {params.pre}
      export OMP_NUM_THREADS={threads}
      source  /etc/profile.d/modules.sh
      module load gpu-tRNA-mapper
      gpu-tRNA-mapper \
        --readFileName {input.fastq:q} \
        --referenceFileName {input.ref_fasta:q} \
        --outputFileName {output:q} \
        {params.gpu_opts} \
        --minAlignmentScore {params.min_aln_score:q}
    ) 2> {log:q}
  """


rule parasail_gpu_assisted_postprocess:
  input: sam="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sam",
         ref=lambda wildcards: REF_FASTA if wildcards.ORIENT == "fwd" else REF_FASTA_REVERSED
  output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sorted.bam")
  conda: "qutrna2"
  benchmark: "benchmarks/parasail/gpu_assisted_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.txt"
  log: "logs/parasail/gpu_assisted_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log"
  shell: """
    (
      samtools view -F 4 -b {input.sam:q} | \
      samtools calmd /dev/stdin {input.ref:q} | \
      python {workflow.basedir}/scripts/add_NH.py /dev/stdin | \
      samtools sort -O bam -o {output:q} /dev/stdin
    ) 2> {log:q}
  """
