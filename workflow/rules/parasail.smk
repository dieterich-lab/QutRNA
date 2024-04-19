_PRE_PROCESS_FASTQ_PARAMS = {
    "fwd": "",
    "rev": "-r",
}


##############################################################################
# Merge nanopore FASTQ files, transform Us to Ts, and
#  reverse nucleotide sequence to obtain random score distribution
##############################################################################
rule parasail_pre_process_fastq:
  input: "data/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}_{BC}.fastq.gz"
  output: temp("results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.fastq.gz")
  conda: "qutrna",
  log: "logs/parasail/pre_process_fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.log",
  params:
    opts=lambda wildcards: _PRE_PROCESS_FASTQ_PARAMS[wildcards.ORIENT]
  shell: """
    ( {workflow.basedir}/scripts/generate_fastq.sh -t {params.opts} {input:q} | \
        gzip -c > {output:q} ) 2> {log:q}

  """


##############################################################################
# Use parasail to map and align reads
##############################################################################

checkpoint parasail_split_reads:
  input: "results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.fastq.gz",
  output: temp(directory("results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}_split")),
  log: "logs/parasail/split_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.log",
  params: lines=config["parasail"]["lines"],
  shell: """
    ( gunzip -c {input:q} | \
        split -l {params.lines} --filter 'gzip -c > $FILE.fastq.gz' /dev/stdin {output:q}/part_ ) 2 > {log:q}
  """


rule parasail_map_split_reads:
  input: fastq="results/fastq/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}_split/part_{part}.fastq.gz",
         ref_fasta=REF_FASTA,
  output: bam=temp("results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}_split/part_{part}.sorted.bam"),
          bai=temp("results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}_split/part_{part}.sorted.bam.bai"),
  log: "logs/parasail/map_split_reads/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}_split/part_{part}.log",
  conda: "qutrna" # TODO load samtools, parasail
  threads: config["parasail"]["threads"]
  params: parasail_opts=config["parasail"]["opts"],
  shell: """
    # TODO -g ${PART%.fastq.gz}.sam \ # TODO stdout
    ( parasail_aligner {params.parasail_opts} \
                       -t {threads} \
                       -O SAMH \
                       -f {input.ref_fasta:q} \
                       -q {input.fastq:q} | \
      gawk -v OFS="\t" \
        ' BEGIN {{ HEADER=0 }}
          $$0 ~ /^@/ {{ if (HEADER==0) {{ print }} ; next }}
          $$0 !~ /^@/ {{ HEADER=1 }}
          $$0 ~ /AS:i:([0-9]+)/ {{
                                 AS=int(gensub(/.+AS:i:([0-9]+).+/, "\\1", "g")) ;
                                 (AS <= 255) ? $$5=AS : $$5=255 ;
                                 print ;
                                 next }} ;
                               {{ print }} ' | \
      samtools view -bS /dev/stdin | \
      samtools calmd /dev/stdin {input.ref_fasta} | \
      samtools sort -@ {threads} -m 2G -o {output.bam} /dev/stdin && samtools index {output.bam} ) 2> {log:q}
  """


##############################################################################
# Retain highest scoring alignments with a minimum alignment
##############################################################################

rule parasail_retain_highest_scoring_alignment:
  input: "results/bams/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.sorted.bam",
  output: "results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.sorted.bam",
  log: "logs/parasail/retain_highest_scoring_alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~{ORIENT}_{BC}.log",
  params: min_aln_core=config["params"]["min_aln_score"],
  shell: """
    {workflow.basedir}/retain_highest_alignments.sh \
        -b {output:q} \
        -s {params.min_aln_score} {input:q} 2> {log:q}
  """


##############################################################################
# Filter alignment by random score distribution
##############################################################################
rule parasail_filter_by_random_score:
  input: fwd="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd_{BC}.sorted.bam",
         rev="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~rev_{BC}.sorted.bam",
  output: plot="results/plots/alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.pdf",
          cutoff="results/bams/filtered/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_cutoff.txt",
  params: precision=config["params"]["precision"],
          title=lambda wildcards: f"tRNA Alignment Score {wildcards.BC} distributions",
  log: "logs/parasail/filter_by_random_score/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log",
  shell: """
    Rscript --vanilla {workflow.basedir}/scripts/alignment_score_cutoff.R \
      -o {output.plot:q} \
      -p {params.precision} \
      -t {params.title:q} \
      --forward {input.fwd:q} --reverse {input.rev:q} 2> {log:q}
  """
