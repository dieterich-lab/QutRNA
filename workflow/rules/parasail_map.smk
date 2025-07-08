from snakemake.io import temp


global REF_FASTA


rule parasail_map:
    input: fastq="results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.fastq.gz",
           ref_fasta=REF_FASTA
    output: temp("results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.bam")
    log: "logs/parasail/map/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log"
    conda: "qutrna"
    resources:
      mem_mb=16000
    threads: config["parasail"]["threads"]
    benchmark: "benchmarks/parasail/map/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.bam"
    params: parasail_opts=config["parasail"]["opts"]
    shell: """
      (
      parasail_aligner {params.parasail_opts} \
                        -t {threads} \
                        -O SAMH \
                        -f {input.ref_fasta:q} \
                        -g /dev/stdout \
                        < {input.fastq:q} 2> {log:q}
        samtools view -bS /dev/stdout | \
        samtools calmd --output-fmt BAM /dev/stdin {input.ref_fasta:q} > {output:q}
      ) 2> "{log}"
    """
