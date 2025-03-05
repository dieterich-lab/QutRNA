# FIXME return value vs. rule

def _filter_input_helper(i):
   bam = FILTER2OUTPUT[i - 1]
   d = {"bam": bam,
        "bai": bam + ".bai",}

   return d


def _filter_samtools(i):
  config_filter = config["preprocess"]["samtools"]["filter"]
  config_calmd = config["preprocess"]["samtools"]["calmd"]
  if not config_filter and not config_calmd:
    return

  input_ = _filter_input_helper(i)
  if config_calmd:
     input_["ref"] = REF_FASTA

  output_ = "results/bams/preprocessed/samtools/{filename}.bam"
  FILTER2OUTPUT.append(output_)
  FILTERS_APPLIED.append("samtools")

  rule:
    name: "dny_filter_samtools",
    input: **input_,
    output: temp(output_),
    log: f"logs/preprocessed/{i}_samtools/{{filename}}.log"
    params:
      filter=config_filter,
      calmd=config_calmd,
    run:
      cmds = ["samtools view {params.filter} -b {input.bam}", ]
      if params.calmd:
        cmds.append("samtools calmd -b /dev/stdin {input.ref}")
      cmd = " | ".join(cmds)
      cmd = "( " + cmd + " > {output[0]} ) 2> {log}"

      shell(cmd)


def _filter_trim_cigar(i):
  config_trim_cigar = config["preprocess"]["trim_cigar"]
  if not config_trim_cigar:
    return

  input_ = _filter_input_helper(i)
  bam = "results/bams/preprocessed/trim_cigar/{filename}.bam"
  output_ = {"bam": temp(bam),
             "stats": "results/bams/preprocessed/trim_cigar/{filename}_stats.tsv", }
  FILTER2OUTPUT.append(bam)
  FILTERS_APPLIED.append("trim_cigar")

  rule:
    name: "dyn_filter_trim_cigar"
    input: **input_,
    output: **output_,
    log: f"logs/preprocessed/{i}_trim_cigar/{{filename}}.log"
    params:
      trim_cigar=config_trim_cigar,
    run:
      cmds = [
        "python " + os.path.join(workflow.basedir, "scripts", "process_read.py") + " --stats {output.stats} --trim-cigar {input.bam} ",
        "samtools sort"]

      cmd = " | ".join(cmds)
      cmd = "( " + cmd + " > {output.bam} ) 2> {log}"
      shell(cmd)


def _filter_read_length(i):
  config_min = config["preprocess"]["read_length"]["min"]
  config_max = config["preprocess"]["read_length"]["max"]
  if not config_min and not config_max:
    return

  input_ =_filter_input_helper(i)
  bam = "results/bams/preprocessed/read_length/{filename}.bam"
  output_ = {"bam": temp(bam),
             "stats": "results/bams/preprocessed/read_length/{filename}_stats.tsv", }
  FILTER2OUTPUT.append(bam)
  FILTERS_APPLIED.append("read_length")

  rule:
    name: "dyn_filter_read_length",
    input: **input_,
    output: **output_,
    log: f"logs/preprocessed/{i}_read_length/{{filename}}.log"
    params:
      min=config["preprocess"]["read_length"]["min"],
      max=config["preprocess"]["read_length"]["max"],
    run:
      process_read_opts = []
      if params.min:
        process_read_opts.append(f"--min-read-length {params.min}")
      if params.max:
        process_read_opts.append(f"--max-read-length {params.max}")
      cmd = "python " + os.path.join(workflow.basedir, "scripts", "process_read.py") + " " + " ".join(process_read_opts) + " --stats {output.stats} {input.bam} > {output.bam} 2> {log}"
      shell(cmd)


def _filter_alignment_length(i):
  config_min = config["preprocess"]["alignment_length"]["min"]
  config_max = config["preprocess"]["alignment_length"]["max"]
  if not config_min and not config_max:
    return

  input_ =_filter_input_helper(i)
  bam = "results/bams/preprocessed/alignment_length/{filename}.bam"
  output_ = {"bam": temp(bam),
             "stats": "results/bams/preprocessed/alignment_length/{filename}_stats.tsv", }
  FILTER2OUTPUT.append(bam)
  FILTERS_APPLIED.append("alignment_length")

  rule:
    name: "dyn_filter_alignment_length",
    input: **input_,
    output: **output_,
    log: f"logs/preprocessed/{i}_alignment_length/{{filename}}.log"
    params:
      min=config["preprocess"]["alignment_length"]["min"],
      max=config["preprocess"]["alignment_length"]["max"],
    run:
      process_read_opts = []
      if params.min:
        process_read_opts.append(f"--min-alignment-length {params.min}")
      if params.max:
        process_read_opts.append(f"--max-alignment-length {params.max}")
      cmd = "python " + os.path.join(workflow.basedir, "scripts", "process_read.py") + " " + " ".join(process_read_opts) + " --stats {output.stats} {input.bam} > {output.bam} 2> {log}"
      shell(cmd)


# preprocess -> filter function
_filters = {
    "samtools": _filter_samtools,
    "trim_cigar": _filter_trim_cigar,
    "read_length": _filter_read_length,
    "alignment_length": _filter_alignment_length,
}


FILTERS_APPLIED = []
# apply filters sequentially
FILTER2OUTPUT = ["data/bams/{filename}.bam", ]
if "preprocess" in config:
  for i, filter_name in enumerate(config["preprocess"], start=1):
    filter_rule = _filters[filter_name]
    filter_rule(i)
