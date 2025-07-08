from snakemake.io import temp


global REF_FASTA


_PRE_PROCESS_FASTQ_PARAMS = {
    "fwd": "",
    "rev": "-r",
}


##############################################################################
# Merge nanopore FASTQ files, transform Us to Ts, and
#  reverse nucleotide sequence to obtain random score distribution
##############################################################################
rule parasail_pre_process_fastq:
  input: "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.fastq.gz"
  output: temp("results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.fastq.gz")
  conda: "qutrna"
  resources:
    mem_mb=2000
  log: "logs/parasail/pre_process_fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log"
  params:
    opts=lambda wildcards: _PRE_PROCESS_FASTQ_PARAMS[wildcards.ORIENT]
  benchmark:
    "benchmarks/parasail_pre_process_fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.fastq.gz"
  shell: """
    python {workflow.basedir}/scripts/generate_fastq.py \
        -t {params.opts} \
        -o {output:q} \
        {input:q} 2> "{log}"
  """


##############################################################################
# Use GPU-assisted or only parasail to map reads
##############################################################################
if config["params"]["alignment"] == "gpu":
  include: "parasail_gpu.smk"
else:
  ############################################################################
  # Retain highest scoring alignments with a minimum alignment
  ############################################################################
  rule parasail_retain_highest_scoring_alignment:
    input: bam="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{prefix}_raw.bam",
           ref_fasta=REF_FASTA
    output: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{prefix}.sorted.bam"
    log: "logs/parasail/retain_highest_scoring_alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{prefix}.log"
    conda: "qutrna"
    resources:
      mem_mb=6000
    params: min_aln_score=config["params"]["min_aln_score"]
    shell: """
      (
        samtools sort -n -m 4G {input.bam:q} | \
          python {workflow.basedir}/scripts/fix_retain_highest_alignment.py --min-mapq {params.min_aln_score} | \
          samtools sort -o {output:q} -m 4G /dev/stdin && samtools index {output:q} 
      ) 2> "{log}"
    """

  if config["parasail"]["lines"] > 0:
    include: "parasail_split.smk"
  else:
    include: "parasail_map.smk"


##############################################################################
# Filter alignments by random score distribution
##############################################################################
# TODO nice plot
rule parasail_filter_by_random_score:
  input: fwd="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd/{BC}_score.txt",
         rev="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~rev/{BC}_score.txt"
  output: prc_plot="results/plots/alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_prc_reads.pdf",
          score_plot="results/plots/alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_score.pdf",
          cutoff="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_cutoff.txt"
  params: precision=config["params"]["precision"],
          title=lambda wildcards: f"tRNA Alignment Score {wildcards.BC} distributions"
  log: "logs/parasail/filter_by_random_score/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
  conda: "qutrna"
  resources:
    mem_mb=2000
  shell: """
    Rscript --vanilla {workflow.basedir}/scripts/alignment_score_cutoff.R \
      -P {output.prc_plot:q} \
      -S {output.score_plot:q} \
      -C {output.cutoff:q} \
      -p {params.precision} \
      -t {params.title:q} \
      --forward {input.fwd:q} --reverse {input.rev:q} 2> "{log}"
  """
