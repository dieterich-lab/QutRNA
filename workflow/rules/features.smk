from snakemake.io import unpack


global BAM_READS_INPUT
global FASTQ_READS_INPUT
global FILTERS_APPLIED
global READS
global SAMPLES
global TBL


def _aggregate_feature_helper(read_type, sample, row):
  key = [row.condition, sample, row.subsample_name, row.base_calling]
  return "|".join([read_type] + key)

def _aggregate_feature_input(wildcards):
  if wildcards.FEATURE == "alignment_score" and READS == "fastq":
    t2fnames = _aggregate_feature_general_input(wildcards)
    for sample in SAMPLES:
      df = TBL.loc[[sample]]
      for row in df.itertuples(index=False):
        t2fnames[_aggregate_feature_helper("mapped-random", sample, row)] = f"results/bam/mapped/sample~{sample}/subsample~{row.subsample_name}/alignment~random/{row.base_calling}_stats/{wildcards.FEATURE}.{wildcards.suffix}"
  elif wildcards.FEATURE == "cutoff" and READS == "fastq":
    t2fnames = {}
    for sample in SAMPLES:
      df = TBL.loc[[sample]]
      for row in df.itertuples(index=False):
        t2fnames[_aggregate_feature_helper("mapped", sample, row)] = f"results/bam/mapped/sample~{sample}/subsample~{row.subsample_name}/{row.base_calling}_stats/{wildcards.FEATURE}.{wildcards.suffix}"
  elif wildcards.FEATURE in ["record_count", "read_length"] and READS == "fastq":
    fname = FASTQ_READS_INPUT.replace(".fastq.gz", "")
    fname = f"{fname}_stats/{wildcards.FEATURE}.{wildcards.suffix}"
    t2fnames = {}
    for sample in SAMPLES:
      df = TBL.loc[[sample]]
      for row in df.itertuples(index=False):
        key = _aggregate_feature_helper("fastq", sample, row)
        t2fnames[key] = fname.format(
          SAMPLE=sample,
          SUBSAMPLE=row.subsample_name,
          BC=row.base_calling
        )
    t2fnames.update(_aggregate_feature_general_input(wildcards))
  else:
    t2fnames = _aggregate_feature_general_input(wildcards)

  return t2fnames

def _aggregate_feature_general_input(wildcards):
  if READS == "bam":
    read_type = "raw"
  elif READS == "fastq":
    read_type = "mapped"
  else:
    raise Exception()

  fname = BAM_READS_INPUT.replace(".sorted.bam","")
  fname = f"{fname}_stats/{wildcards.FEATURE}.{wildcards.suffix}"

  t2fnames = {}
  for sample in SAMPLES:
    df = TBL.loc[[sample]]

    # collect subsamples
    for row in df.itertuples(index=False):
      key = _aggregate_feature_helper(read_type, sample, row)

      t2fnames[key] = fname.format(
        SAMPLE=sample,
        SUBSAMPLE=row.subsample_name,
        BC=row.base_calling
      )
      for filter in FILTERS_APPLIED:
        t2fnames[_aggregate_feature_helper(filter, sample, row)] = f"results/bam/filtered-{filter}/sample~{sample}/subsample~{row.subsample_name}/{row.base_calling}_stats/{wildcards.FEATURE}.{wildcards.suffix}"

  return t2fnames

#
def _aggregate_feature_params(_, input):
  opts = []
  for key in input.keys():
    opts.append(f"--data {' '.join(key.split('|'))}")

  return " ".join(opts)

rule aggregate_feature:
  input: unpack(_aggregate_feature_input)
  output: "results/stats/{FEATURE}.{suffix}"
  conda: "qutrna2"
  log: "logs/stats/{FEATURE}_{suffix}.log"
  params: opts=_aggregate_feature_params
  shell: """
    python {workflow.basedir}/scripts/aggregate_feature.py \
      --output {output:q} \
      {params.opts} \
      {input:q} \
      2> {log:q}
  """

########################################################################################################################
# FASTQ features

rule fastq_read_count:
  input: "{prefix}.fastq.gz"
  output: "{prefix}_stats/record_count.txt"
  conda: "qutrna2"
  log: "logs/fastq/read_count/{prefix}.log"
  shell: """
    (
      gzip -cd {input:q} | \
      awk -v OFS="\t" ' END {{ print "unique_reads","multimapper_reads","reads","total_records" ; print "NA","NA",NR/4,NR/4 }} ' \
    ) > {output:q} 2> {log:q}
  """


rule fastq_read_length:
  input: "{prefix}.fastq.gz"
  output: "{prefix}_stats/read_length.txt"
  conda: "qutrna2"
  log: "logs/fastq/read_length/{prefix}.log"
  shell: """
    (
      gzip -cd {input:q} | \
        awk ' NR % 4 == 2 {{ print(length($0)) }} ' | \
        sort | \
        uniq -c | \
        sort -k2,2 -V | \
        awk -v OFS="\t" ' BEGIN {{ print "read_length","unique_reads","multimapper_reads","reads","records" }} ; {{ print $2,"NA","NA",$1,"NA" }} ' \
    ) > {output:q} 2> {log:q}
  """
