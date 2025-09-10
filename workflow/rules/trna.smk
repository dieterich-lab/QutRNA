global create_include
global REF_FASTA
global REF_NO_LINKER_FASTA
global REF_FILTERED_TRNAS_FASTA


## Global paths to files
# Set of trans in reference fasta
TRNAS_AVAILABLE = "results/data/trnas/available.txt"
# Final set of trnas of filtering with "ignore_trna".
# This set will be used for variant calling and secondary structure alignment
TRNAS_SELECTED = "results/data/trnas/selected.txt"
# Provides: trna type, raw sequence (without adpters)
TRNA_ANNOTATION = "results/data/trna_annotation.tsv"
# Helper to define tRNA to omit from secondary structure alignment


def _ignore_trans_opts(_):
  opts = "--base-change T2U"
  pattern = pep.config["qutrna2"].get("sprinzl", {}).get("ignore_trnas", "")
  if pattern:
    opts += f" --ignore '{pattern}'"

  return opts

rule ignore_trnas:
  input: REF_NO_LINKER_FASTA
  output: REF_FILTERED_TRNAS_FASTA
  conda: "qutrna2"
  log: "logs/ignore_trnas.log"
  params: opts=_ignore_trans_opts
  shell: """
    python {workflow.basedir}/scripts/fasta_utils.py \
        transform  \
        {params.opts} \
        --output {output:q} \
        {input:q} \
        2> {log:q}
  """


rule trnas_available:
  input: REF_FASTA
  output: TRNAS_AVAILABLE
  conda: "qutrna2"
  log: "logs/trnas/available.log"
  shell: """
    python {workflow.basedir}/scripts/fasta_utils.py extract-seqids {input:q} > {output:q} 2> {log:q}
  """


rule trnas_selected:
  input: REF_FILTERED_TRNAS_FASTA
  output: TRNAS_SELECTED
  conda: "qutrna2"
  log: "logs/trnas/selected.log"
  shell: """
    python {workflow.basedir}/scripts/fasta_utils.py extract-seqids --output {output:q} {input:q} 2> {log:q}
  """


# Infer tRNA annotation from FASTA or use provided
if "trna_annotation" in pep.config["qutrna2"]:
  create_include("trna_annotation",
    pep.config["qutrna2"]["trna_annotation"],
    TRNA_ANNOTATION,
    config["include"]["trna_annotation"])
else:
  rule infer_trna_annotation:
    input: REF_FILTERED_TRNAS_FASTA
    output: TRNA_ANNOTATION
    conda: "qutrna2"
    log: "logs/infer_trna_annotation.log"
    shell: """
      python {workflow.basedir}/scripts/fasta_utils.py infer-annotation --output {output:q} {input:q} 2> {log:q}
    """
