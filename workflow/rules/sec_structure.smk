import pandas as pd

global SPRINZL_MODE
global SEQ_TO_SPRINZL_INIT
global SEQ_TO_SPRINZL_FINAL
global REF_FILTERED_TRNAS_FASTA
global SPRINZL


if SPRINZL_MODE == "cm":
  global CM

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
           sprinzl=SPRINZL
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
    shell: """
      python {workflow.basedir}/scripts/sprinzl_utils.py map-seq-sprinzl \
        --output {output:q} \
        --cons-ss-annotation {input.consensus_annotated:q} \
        {input.align:q} \
        2> {log:q}
    """

else:
  _SS_SEQ_TO_SPRINZL = "data/seq_to_sprinzl.tsv"

########################################################################################################################

rule ss_seq_to_sprinzl_final:
  input: SEQ_TO_SPRINZL_INIT
  output: SEQ_TO_SPRINZL_FINAL
  run:
    df = pd.read_csv(input[0],sep="\t")
    df = df[df["sprinzl"] != "-"]
    df.to_csv(output[0],sep="\t",index=False,quoting=False)

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
