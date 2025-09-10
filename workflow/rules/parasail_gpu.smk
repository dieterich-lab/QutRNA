from snakemake.io import temp


global REF_FASTA
global REF_FASTA_RANDOM
global FASTQ_READS_INPUT

##############################################################################
# Use GPU-assisted alignment

rule parasail_gpu_assisted:
  input: fastq=FASTQ_READS_INPUT,
         ref_fasta=lambda wildcards: REF_FASTA if wildcards.ALIGNMENT == "real" else REF_FASTA_RANDOM
  output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.sam")
  log: "logs/parasail/gpu_assisted/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.log"
  benchmark: repeat("benchmarks/parasail/gpu_assisted/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.txt", config.get("_benchmark_repeats", 1))
  threads: 1
  params:
    gpu_opts=config["gpu"]["opts"],
    min_aln_score=config["alignment"]["min_aln_score"],
    pre=config["gpu"]["pre"],
    bin=config["gpu"].get("bin", "gpu-tRNA-mapper")
  shell: """
    (
      {params.pre}
      export OMP_NUM_THREADS={threads}
      {params.bin} \
        --readFileName {input.fastq:q} \
        --referenceFileName {input.ref_fasta:q} \
        --outputFileName {output:q} \
        {params.gpu_opts} \
        --minAlignmentScore {params.min_aln_score:q} 
    ) 2> {log:q}
  """


rule samtools_sam_to_bam:
  input: "results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.sam",
  output: "results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.bam",
  log: "logs/samtools_sam_to_bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.log"
  benchmark: "benchmarks/sam_to_bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.txt"
  conda: "qutrna2"
  shell: """
    samtools view -b -o {output:q} -F 4 {input:q} 2> {log:q}
  """


rule parasail_gpu_assisted_postprocess:
  input: bam="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.bam",
         ref=lambda wildcards: REF_FASTA if wildcards.ALIGNMENT == "real" else REF_FASTA_RANDOM
  output: "results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.sorted.bam"
  conda: "qutrna2"
  benchmark: repeat("benchmarks/parasail/gpu_assisted_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.txt", config.get("_benchmark_repeats", 1))
  log: "logs/parasail/gpu_assisted_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.log"
  shell: """
    (
      samtools calmd {input.bam:q} {input.ref:q} | \
      samtools sort -n -O bam /dev/stdin | \
      python {workflow.basedir}/scripts/bam_utils.py add-hits /dev/stdin | \
      samtools sort -O bam -o {output:q} /dev/stdin
    ) 2> {log:q}
  """
