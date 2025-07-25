global FILTERS_APPLIED
global REF_FASTA
global READS
global TBL


rule samtools_bam_index:
  input: "{prefix}.sorted.bam"
  output: "{prefix}.sorted.bam.bai"
  conda: "qutrna2"
  log: "logs/samtools/index_bam/{prefix}.log"
  shell: """
    samtools index {input:q} 2> {log:q}
  """


rule samtools_fasta_index:
  input: "{prefix}.fasta"
  output: "{prefix}.fasta.fai"
  conda: "qutrna2"
  log: "logs/samtools_fasta_index/{prefix}.log"
  shell: """
    samtools faidx {input:q} 2> {log:q}
  """


rule samtools_stats:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/samtools.txt"
  conda: "qutrna2"
  log: "logs/samtools/stats/{prefix}.log"
  shell: """
    samtools stats {input:q} > {output:q} 2> {log:q}
  """


SAMTOOLS_STATS_TYPE2COLS = {
  "RL": ["read_length", "count"],
}

rule samtools_cut_stats:
  input: "{prefix}_stats/samtools.txt"
  output: "{prefix}_stats/samtools_{type}.txt"
  conda: "qutrna2"
  params: cols=lambda wildcards: "\t".join(SAMTOOLS_STATS_TYPE2COLS[wildcards.type])
  log: "logs/samtools_cut_stats/{prefix}/{type}.log"
  shell: """
    (
      echo -e "{params.cols}" > {output:q}
      grep "^{wildcards.type}" {input:q} | cut -f 2- >> {output:q}
    ) 2> {log:q}
  """


rule samtools_coverage:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/samtools_coverage.txt"
  conda: "qutrna2"
  log: "logs/samtools/coverage/{prefix}.log"
  shell: """
    samtools coverage {input:q} > {output:q} 2> {log:q}
  """


rule samtools_read_count:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/read_count.txt"
  conda: "qutrna2"
  log: "logs/samtools/read_count/{prefix}.log"
  shell: """
    samtools view -c {input:q} | awk 'BEGIN {{ print "reads" }} ; {{ print }} '> {output:q} 2> {log:q}
  """


rule samtools_multimapper:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/multimapper.txt"
  conda: "qutrna2"
  log: "logs/samtools/multimapper/{prefix}.log"
  shell: """
    ( samtools view {input:q} | \
      cut -f1,3 | \
      sort -k1,1 | \
      python  {workflow.basedir}/count_multimapper.py ) > {output:q} 2> {log:q}
  """


rule samtools_alignment_score:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/alignment_score.txt"
  conda: "qutrna2"
  log: "logs/samtools/alignment_score/{prefix}.log"
  shell: """
    python {workflow.basedir}/scripts/get_as.py {input:q} > {output:q} 2> {log:q}
  """


def _samtools_merge_input(wildcards):
  fnames = []
  for row in TBL.loc[[wildcards.SAMPLE]].itertuples():
    fname = f"results/bam/{{BAM_TYPE}}/sample~{row.sample_name}/subsample~{row.subsample_name}/{row.base_calling}.sorted.bam"
    fnames.append(fname)

  return fnames


rule samtools_merge:
  input: _samtools_merge_input
  output: "results/bam/{BAM_TYPE}/{SAMPLE}.sorted.bam"
  conda: "qutrna2"
  log: "logs/samtools/merge/{BAM_TYPE}/{SAMPLE}.log"
  threads: 1
  shell: """
    samtools merge -@ {threads} {output:q} {input:q} 2> {log:q}
  """