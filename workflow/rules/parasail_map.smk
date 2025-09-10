from snakemake.io import temp


global REF_FASTA
global REF_FASTA_RANDOM
global FASTQ_READS_INPUT


# FIXME -t 8 \
rule parasail_map:
  input: fastq=FASTQ_READS_INPUT,
         ref_fasta=lambda wildcards: REF_FASTA if wildcards.ALIGNMENT == "real" else REF_FASTA_RANDOM
  output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.sam")
  log: "logs/parasail/map/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.log"
  conda: "qutrna2"
  threads: 1
  benchmark: repeat("benchmarks/parasail/map/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.txt", config.get("_benchmark_repeats", 1))
  params: parasail_opts=config["parasail"]["opts"],
          pre=config["parasail"].get("pre", "")
  shell: """
    (
      {params.pre}
      parasail_aligner {params.parasail_opts} \
                        -O SAMH \
                        -f {input.ref_fasta:q} \
                        -g {output:q} \
                        < {input.fastq:q}
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


rule parasail_map_postprocess:
    input: bam="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.bam",
           ref= lambda wildcards: REF_FASTA if wildcards.ALIGNMENT == "real" else REF_FASTA_RANDOM
    output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.sorted.bam")
    log: "logs/parasail/map_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.log"
    benchmark: repeat("benchmarks/parasail/map_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~{ALIGNMENT}/{BC}.txt", config.get("_benchmark_repeats", 1))
    conda: "qutrna2"
    params:
        min_aln_score=config["alignment"]["min_aln_score"]
    shell: """
    (
      python {workflow.basedir}/scripts/bam_utils.py best-alignment --min-as {params.min_aln_score} {input.bam:q} | \
      samtools calmd - {input.ref:q} | \
      samtools sort -n -O bam /dev/stdin | \
      python {workflow.basedir}/scripts/bam_utils.py add-hits /dev/stdin | \
      samtools sort -O bam -o {output:q} /dev/stdin
    ) 2> {log:q}
  """
