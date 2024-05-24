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
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: "logs/parasail/pre_process_fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log",
  params:
    opts=lambda wildcards: _PRE_PROCESS_FASTQ_PARAMS[wildcards.ORIENT]
  shell: """
    python {workflow.basedir}/scripts/generate_fastq.py \
        -t {params.opts} \
        -o {output:q} \
        {input:q} 2> {log:q}
  """


##############################################################################
# Use parasail to map and align reads
##############################################################################

def _batch_size():
  if config["parasail"]["batch_size"]:
    return "-b " + str(config["parasail"]["batch_size"])

  return ""

if config["parasail"]["lines"] > 0:
  checkpoint parasail_split_reads:
    input: "results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.fastq.gz",
    output: temp(directory("results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split")),
    log: "logs/parasail/split_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log",
    conda: "qutrna"
    resources:
      mem_mb=2000
    params: lines=config["parasail"]["lines"],
    shell: """
      mkdir -p {output}
      ( gunzip -c {input:q} | \
          split -l {params.lines} --filter 'gzip -c > $FILE.fastq.gz' /dev/stdin {output:q}/part_ ) 2> {log:q}
    """


  rule parasail_map_split_reads:
    input: fastq="results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.fastq.gz",
           ref_fasta=REF_FASTA,
    output: temp("results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}_raw.bam"),
    log: "logs/parasail/map_split_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_split/part_{part}.log",
    conda: "qutrna"
    resources:
      mem_mb=16000
    threads: config["parasail"]["threads"]
    params: parasail_opts=config["parasail"]["opts"],
            parasail_batch_size=_batch_size(),
    shell: """
    (
      parasail_aligner {params.parasail_opts} {params.parasail_batch_size} \
                        -t {threads} \
                        -O SAMH \
                        -f {input.ref_fasta:q} \
                        -g {output:q}.tmp \
                        -q {input.fastq:q}
      samtools view -bS {output}.tmp | \
      samtools calmd --output-fmt BAM /dev/stdin {input.ref_fasta:q} > {output:q} && \
      rm {output}.tmp ) 2> {log:q}
    """
else:
  rule parasail_map:
    input: fastq="results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.fastq.gz",
           ref_fasta=REF_FASTA,
    output: temp("results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}_raw.bam"),
    log: "logs/parasail/map/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{BC}.log",
    conda: "qutrna"
    resources:
      mem_mb=16000
    threads: config["parasail"]["threads"]
    params: parasail_opts=config["parasail"]["opts"],
            parasail_batch_size=_batch_size(),
    shell: """
      parasail_aligner {params.parasail_opts} {params.parasail_batch_size} \
                        -t {threads} \
                        -O SAMH \
                        -f {input.ref_fasta:q} \
                        -g {output:q}.tmp \
                        < {input.fastq:q} 2> {log:q}
      samtools view -bS {output}.tmp | \
        samtools calmd --output-fmt BAM /dev/stdin {input.ref_fasta:q} > {output:q}
      rm {output}.tmp
    """

##############################################################################
# Retain highest scoring alignments with a minimum alignment
##############################################################################

rule parasail_retain_highest_scoring_alignment:
  input: bam="results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{prefix}_raw.bam",
         ref_fasta=REF_FASTA,
  output: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{prefix}.sorted.bam",
  log: "logs/parasail/retain_highest_scoring_alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}/{prefix}.log",
  conda: "qutrna"
  resources:
    mem_mb=6000
  params: min_aln_score=config["params"]["min_aln_score"],
  shell: """
    (
      samtools sort -n -m 4G {input.bam:q} | \
      python {workflow.basedir}/scripts/fix_retain_highest_alignment.py --min-mapq {params.min_aln_score} | \
      samtools sort -o {output:q} -m 4G /dev/stdin && samtools index {output:q} ) 2> {log:q}
  """


##############################################################################
# Filter alignment by random score distribution
##############################################################################
rule parasail_filter_by_random_score:
  input: fwd="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd/{BC}_score.txt",
         rev="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~rev/{BC}_score.txt",
  output: prc_plot="results/plots/alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_prc_reads.pdf",
          score_plot="results/plots/alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_score.pdf",
          cutoff="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_cutoff.txt",
  params: precision=config["params"]["precision"],
          title=lambda wildcards: f"tRNA Alignment Score {wildcards.BC} distributions",
  log: "logs/parasail/filter_by_random_score/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log",
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
      --forward {input.fwd:q} --reverse {input.rev:q} 2> {log:q}
  """
