_PRE_PROCESS_FASTQ_PARAMS = {
    "fwd": "",
    "rev": "-r",
}


##############################################################################
# Merge nanopore FASTQ files, transform Us to Ts, and
#  reverse nucleotide sequence to obtain random score distribution
##############################################################################
rule parasail_pre_process_fastq:
  input: "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}_{RYPTE}.fastq.gz"
  output: "results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}.fastq.gz"
  conda: "qutrna",
  log: "logs/parasail/pre_process_fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}.log",
  params:
    opts=lambda wildcards: _PRE_PROCESS_FASTQ_PARAMS[wildcards.DIR]
  shell: """
    ( {workflow.basedir}/scripts/generate_fastq.sh -t {params.opts} {input:q} | \
        gzip -c > {output:q} ) 2> {log:q}
  """


##############################################################################
# Use parasail to map and align reads
##############################################################################
rule parasail_map_reads:
  input: fastq="results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}.fastq.gz",
         ref_fasta=REF_FASTA,
  output: bam="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_RTYPE.sorted.bam",
          bai="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_RTYPE.sorted.bam.bai",
  log: "logs/parasail/map_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_RTYPE.log",
  threads: config["parasail"]["threads"]
  shell: """
    ( {workflow.basedir}/parasail_align.sh -b {output:q} \
        -t {THREADS} -f {input.ref_fasta:q} {input.fastq:q} ) 2> {log:q}
  """


##############################################################################
# Retain highest scoring alignments with a minimum alignment
##############################################################################
rule parasail_retain_highest_scoring_alignment:
  input: "results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}.sorted.bam",
  output: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}.sorted.bam",
  log: "logs/parasail/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{RTYPE}.log",
  threads: config["samtools"]["threads"]
  params: min_core=config["min_score"],
  shell: """
    ( {workflow.basedir}/retain_highest_alignments.sh -b {output:q} \
        -t {threads} -s {params.min_score} {input:q} ) 2> {log:q}
  """


##############################################################################
# Filter alignment by random score distribution
##############################################################################
rule parasail_filter_by_random_score:
  input: fwd="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd_{RTYPE}.sorted.bam",
         rev="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~rev_{RTYPE}.sorted.bam",
  output: plot="results/plots/alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}.pdf",
          cutoff="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{RTYPE}_cutoff.txt",
  params: precision=config["precision"],
          title=lambda wildcards: f"tRNA Alignment Score {wildcards.RTYPE} distributions",
  log: "logs/parasail/filter_by_random_score/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{RTYPE}.log",
  shell: """
    Rscript --vanilla {workflow.basedir}/scripts/alignment_score_cutoff.R \
      -o {output.plot:q} \
      -p {params.precision:q} \
      -t {params.title:q} \
      --forward {input.fwd:q} --reverse {input.rev:q} 2> {log:q}
  """
