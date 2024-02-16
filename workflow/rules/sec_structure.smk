rule transform:
  input: jacusa2="results/jacusa2/cond1~{COND1}/cond2~{COND2}/coords~pos_scores.tsv",
         numbering="results/coord_trans.tsv",
  output: temp("results/jacusa2/cond1~{COND1}/cond2~{COND2}/coords~u_pos_scores_pre-mods.tsv"),
  conda: "qutrna",
  log: "logs/transform/cond1~{COND1}/cond2~{COND2}.log",
  params: linker5=pep.config["qutrna"]["linker5"],
  shell: """
    python {workflow.basedir:q}/scripts/transform.py \
        --numbering {input.numbering:q} \
        --output {output:q} \
        --linker5 {params.linker5} \
        --aln_prefix "" \
        {input.jacusa2:q} 2> {log:q}
  """


rule trnascan_se_run:
  input: fasta="results/data/no_linkers.fasta",
  output: result="results/trnascan_se/result.txt",
          structure="results/trnascan_se/structure.ss",
          isospecific="results/trnascan_se/isospecific.txt",
          stats="results/trnascan_se/stats.txt",
          bed="results/trnascan_se/result.bed",
          gff="results/trnascan_se/result.gff",
          fasta="results/trnascan_se/result.fasta",
          log="results/trnascan_se/log.txt",
  conda: "qutrna",
  log: "logs/trnascan_se.log",
  params: opts=config["trnascan_se"]["opts"],
  threads: config["trnascan_se"]["threads"],
  shell: """
    tRNAscan-SE {params.opts} \
        -o {output.result:q} \
        -f {output.structure:q} \
        -s {output.isospecific:q} \
        -m {output.stats:q} \
        -b {output.bed:q} \
        -j {output.gff:q} \
        -a {output.fasta:q} \
        -l {output.log:q} \
        {input.fasta:q} 2> {log:q}
  """
