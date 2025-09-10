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


SAMTOOLS_FEATURES_TYPE2COLS = {
  "RL": ["read_length", "count"],
}

rule samtools_cut_stats:
  input: "{prefix}_stats/samtools.txt"
  output: "{prefix}_stats/samtools_{type}.txt"
  conda: "qutrna2"
  params: cols=lambda wildcards: "\t".join(SAMTOOLS_FEATURES_TYPE2COLS[wildcards.type])
  log: "logs/samtools_cut_stats/{prefix}/{type}.log"
  shell: """
    (
      echo -e "{params.cols}" > {output:q}
      grep "^{wildcards.type}" {input:q} | cut -f 2- >> {output:q}
    ) 2> {log:q}
  """


rule samtools_read_length:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/read_length.txt"
  conda: "qutrna2"
  log: "logs/samtools/read_length/{prefix}.log"
  shell: """
    python {workflow.basedir}/scripts/bam_utils.py count-record-length --output {output:q} {input:q} 2> {log:q}
  """


rule samtools_record_count:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/record_count.txt"
  conda: "qutrna2"
  log: "logs/samtools/record_count/{prefix}.log"
  shell: """
    python {workflow.basedir}/scripts/bam_utils.py count-records --simple-count --output {output:q} {input:q} 2> {log:q}
  """


rule samtools_multimapper:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/multimapper.txt"
  conda: "qutrna2"
  log: "logs/samtools/multimapper/{prefix}.log"
  shell: """
    python  {workflow.basedir}/scripts/bam_utils.py count-multimapper --output {output:q} {input:q} 2> {log:q}
  """


rule samtools_alignment_score:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_stats/alignment_score.txt"
  conda: "qutrna2"
  log: "logs/samtools/alignment_score/{prefix}.log"
  shell: """
    python {workflow.basedir}/scripts/bam_utils.py count-tag -t AS -c alignment_score --output {output:q} {input:q} 2> {log:q}
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
