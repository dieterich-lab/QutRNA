from snakemake.io import temp


global REF_FASTA
global REF_FASTA_REVERSED
global READS_INPUT


rule parasail_map:
    input: fastq=READS_INPUT,
           ref_fasta=lambda wildcards: REF_FASTA if wildcards.ORIENT == "fwd" else REF_FASTA_REVERSED
    output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sam")
    log: "logs/parasail/map/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log"
    conda: "qutrna2"
    threads: 1
    benchmark: "benchmarks/parasail/map/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.bam"
    params: parasail_opts=config["parasail"]["opts"]
    shell: """
      parasail_aligner {params.parasail_opts} \
                        -t {threads} \
                        -O SAMH \
                        -f {input.ref_fasta:q} \
                        -g {output:q} \
                        < {input.fastq:q} 2> {log:q}
    """


rule parasail_map_postprocess:
    input: sam="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.sam",
           ref= lambda wildcards: REF_FASTA if wildcards.ORIENT == "fwd" else REF_FASTA_REVERSED
    output: temp("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.bam")
    log: "logs/parasail/map_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.bam"
    benchmark: "benchmarks/parasail/map_postprocess/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.tsv"
    conda: "qutrna2"
    params:
        min_aln_score=config["params"]["min_aln_score"]
    shell: """
    (
      samtools view -b -F 4 {input.sam:q} | \
      samtools calmd /dev/stdin {input.ref:q} | \
      samtools sort -n /dev/stdin | \
      python {workflow.basedir}/scripts/retain_highest_alignment.py --min-aln-score {params.min_aln_score} /dev/stdin \
      python {workflow.basedir}/scripts/add_NH.py /dev/stdin \
      samtools sort /dev/stdin > {output:q}
    ) 2> {log:q}
  """
