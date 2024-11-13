rule ss_transform:
  input: jacusa2="results/jacusa2/cond1~{COND1}/cond2~{COND2}/scores_seq.tsv",
         sprinzl=SEQ_TO_SPRINZL,
  output: temp("results/jacusa2/cond1~{COND1}/cond2~{COND2}/scores_sprinzl.tsv"),
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: "logs/ss/transform/cond1~{COND1}/cond2~{COND2}.log",
  params: linker5=pep.config["qutrna"]["linker5"],
  shell: """
    python {workflow.basedir:q}/scripts/transform.py \
        --sprinzl {input.sprinzl:q} \
        --output {output:q} \
        --linker5 {params.linker5} \
        {input.jacusa2:q} 2> {log:q}
  """


if "cm" in pep.config["qutrna"]:
  rule cmalign_run:
    input: cm=CM,
           fasta=REF_FILTERED_TRNAS_FASTA,
    output: "results/cmalign/align.stk",
    log: "logs/cmalign/run.log",
    threads: config["cmalign"]["threads"],
    params: opts=config["cmalign"]["opts"],
    conda: "qutrna",
    resources:
      mem_mb=20000
    shell: """
      cmalign {params.opts} --cpu {threads} -o {output:q} \
          {input.cm:q} {input.fasta:q} \
          2> {log}
    """


  rule ss_consensus_add_sprinzl:
    input: stk="results/cmalign/align.stk",
           sprinzl=SPRINZL,
    output: "results/ss_consensus_to_sprinzl.tsv",
    conda: "qutrna",
    resources:
      mem_mb=2000
    log: "logs/ss/ss_consensus_add_sprinzl.log",
    shell: """
      python {workflow.basedir}/scripts/ss_consensus_add_sprinzl.py \
          --output {output:q} \
          --sprinzl {input.sprinzl:q} \
          {input.stk:q} \
          2> {log:q}
    """


  rule ss_seq_to_sprinzl:
    input: align="results/cmalign/align.stk",
           ss_to_sprinzl="results/ss_consensus_to_sprinzl.tsv",
    output: "results/seq_to_sprinzl.tsv",
    conda: "qutrna",
    resources:
      mem_mb=2000
    log: "logs/ss/seq_to_sprinzl.log",
    shell: """
      python {workflow.basedir}/scripts/seq_to_sprinzl.py \
          --output {output:q} {input.ss_to_sprinzl:q} {input.align:q} \
          2> {log:q}
    """


  rule ss_seq_to_filtered_sprinzl:
    input: "results/seq_to_sprinzl.tsv",
    output: "results/seq_to_sprinzl_filtered.tsv",
    log: "logs/ss/seq_to_filtered_sprinzl.log",
    run:
      df = pd.read_csv(input[0], sep="\t")
      df = df[df["sprinzl"] != "-"]
      df.to_csv(output[0], sep="\t", index=False, quoting=False)
