import pandas as pd

global create_include
global REF_FILTERED_TRNAS_FASTA


########################################################################################################################
# Coordinate system: seq(uence) or sprinzl

# final sequence to sprinzl mapping
SEQ_TO_SPRINZL_FINAL = "results/ss/seq_to_sprinzl_filtered.tsv"

# consensus annotation
SPRINZL_LABELS = "data/sprinzl_labels.txt"
create_include("sprinzl_labels",
  pep.config["qutrna2"]["sprinzl"]["labels"],
  SPRINZL_LABELS,
  config["include"].get("sprinzl_labels","copy"))

# in case of sprinzl, handle if model or precalculated mapping is to be used
if pep.config["qutrna2"]["coords"] == "sprinzl" and "cm" in pep.config["qutrna2"]["sprinzl"]:

  # How mapping to sprinzl coordinates are calculate, by alignment or from existing mapping
  # cm -> covarince model -> calculate secondary structure alignment and create mapping
  # seq_to_sprinzl -> use exisiting mapping
  if "cm" in pep.config["qutrna2"]["sprinzl"]:
    SPRINZL_MODE = "cm"
    # covariance model destination
    CM = "data/cm.stk"
    create_include("cm",
      pep.config["qutrna2"]["sprinzl"]["cm"],
      CM,
      config["include"]["cm"])
    SEQ_TO_SPRINZL_INIT = "results/ss/seq_to_sprinzl.tsv"
  elif "seq_to_sprinzl" in pep.config["qutrna2"]["sprinzl"]:
    SPRINZL_MODE = "seq2sprinzl"
    SEQ_TO_SPRINZL_INIT = "data/seq_to_sprinzl.tsv"
    create_include("seq_to_sprinzl",
      pep.config["qutrna2"]["sprinzl"]["seq_to_sprinzl"],
      SEQ_TO_SPRINZL_INIT,
      config["include"].get("seq_to_sprinzl", "copy"))
  else:
    raise Exception("Must provide either 'cm' or 'seq_to_sprinzl' in pep!")

  rule cmalign_run:
    input: cm=CM,
           fasta=REF_FILTERED_TRNAS_FASTA
    output: "results/cmalign/align.stk"
    log: "logs/cmalign/run.log"
    threads: 1
    params: opts=config["cmalign"]["opts"]
    conda: "qutrna2"
    shell: """
      cmalign {params.opts} \
        --cpu {threads} \
        -o {output:q} \
        {input.cm:q} \
        {input.fasta:q} \
        2> {log}
   """


  rule ss_annotate_consensus:
    input: stk="results/cmalign/align.stk",
           sprinzl=SPRINZL_LABELS
    output: "results/ss/consensus_annotated.tsv"
    conda: "qutrna2"
    log: "logs/ss/annotate_consensus.log"
    shell: """
      python {workflow.basedir}/scripts/sprinzl_utils.py annotate-consensus \
        --output {output:q} \
        --sprinzl {input.sprinzl:q} \
        {input.stk:q} \
        2> {log:q}
    """


  rule ss_map_seq_to_sprinzl:
    input: align="results/cmalign/align.stk",
           consensus_annotated="results/ss/consensus_annotated.tsv"
    output: SEQ_TO_SPRINZL_INIT
    conda: "qutrna2"
    log: "logs/ss/map_seq_to_sprinzl.log"
    params: map_introns=lambda wildcards: "--map-introns" if pep.config["qutrna2"]["sprinzl"].get("map_introns", False) else ""
    shell: """
      python {workflow.basedir}/scripts/sprinzl_utils.py map-seq-sprinzl \
        {params.map_introns} \
        --output {output:q} \
        --cons-ss-annotation {input.consensus_annotated:q} \
        {input.align:q} \
        2> {log:q}
    """
else:
  SEQ_TO_SPRINZL_INIT = "data/seq_to_sprinzl.tsv"

  create_include("seq_to_sprinzl",
    pep.config["qutrna2"]["sprinzl"]["seq_to_sprinzl"],
    SEQ_TO_SPRINZL_INIT,
    config["include"].get("seq_to_sprinzl", "copy"))

########################################################################################################################
# Filter sprinzl mapping and apply

rule ss_seq_to_sprinzl_final:
  input: SEQ_TO_SPRINZL_INIT
  output: SEQ_TO_SPRINZL_FINAL
  run:
    df = pd.read_csv(input[0], sep="\t")
    df = df[df["sprinzl"] != "-"]
    df.to_csv(str(output[0]), sep="\t", index=False, quoting=False)


rule ss_transform:
  input: jacusa2="results/jacusa2/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/scores_seq.tsv",
         seq_sprinzl=SEQ_TO_SPRINZL_FINAL
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/scores_sprinzl.tsv"
  conda: "qutrna2"
  log: "logs/ss/transform/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}.log"
  params: linker5=pep.config["qutrna2"]["linker5"]
  shell: """
    python {workflow.basedir:q}/scripts/sprinzl_utils.py transform \
        --sprinzl {input.seq_sprinzl:q} \
        --output {output:q} \
        --linker5 {params.linker5} \
        {input.jacusa2:q} 2> {log:q}
  """
