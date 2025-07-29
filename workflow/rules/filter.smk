from snakemake.io import temp, unpack

global READS
global REF_FASTA

FILTERS_APPLIED = []

class Filter:

  OUTPUT = []

  def __init__(self, filter_name):
    self.filter_name = filter_name
    self.input = {
      "bam": self.OUTPUT[-1],
      "bai": self.OUTPUT[-1] + ".bai"
    }
    self.output = {
      "bam": temp(
        f"results/bam/filtered-{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}.sorted.bam")
    }
    self.cmds = []
    self.params = {}

  def process(self) -> bool:
    if not self.filter_name in config["filter"]:
      return False

    return self._process()

  @property
  def config(self):
    return config["filter"][self.filter_name]

  def _process(self) -> bool:
    raise NotImplementedError

  def shell_formatted(self):
    cmd = "| ".join(self.cmds)
    cmd = f"( {cmd} ) 2> {{log:q}}"

    return cmd

  def create_rule(self):
    if self.process():
      self.OUTPUT.append(self.output["bam"])
      FILTERS_APPLIED.append(self.filter_name)
      rule:
        name: f"filter_{self.filter_name}"
        input: unpack(lambda _: self.input)
        output: **self.output
        benchmark: f"benchmarks/filter/{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}.txt"
        conda: "qutrna2"
        log: f"logs/filter/{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}.log"
        params: **self.params
        shell: self.shell_formatted()


class FilterRandomAlignment(Filter):
  def __init__(self):
    super().__init__("random_alignment")

  def _process(self):
    self.input["cutoff"] = "results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_stats/cutoff.txt"
    self.cmds.append(
      f"python {workflow.basedir}/scripts/filter_by_as.py --min-alignment-score `sed -n '2p' {{input.cutoff:q}}` {{input.bam:q}} > {{output.bam:q}}")

    return True

class FilterSamtools(Filter):
  def __init__(self):
    super().__init__("samtools")

  def _process(self):
    input_ = "{input.bam:q}"

    if self.config["filter"]:
      self.params["filter"] = self.config["filter"]
      self.cmds.append(f"samtools view {{params.filter}} -b {input_}")
      input_ = "/dev/stdin"

    if "calmd" in config and config["calmd"]:
      self.input["ref"] = REF_FASTA
      self.cmds.append(f"samtools calmd -b {input_}")

    if self.cmds:
      self.cmds[-1] = f"{self.cmds[-1]} > {{output.bam:q}}"

    return self.cmds

class FilterTrimCigar(Filter):
  def __init__(self):
    super().__init__("trim_cigar")

  def _process(self):
    if not self.config:
      return False

    self.output["stats"] = f"results/bam/filtered-{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}_stats/{self.filter_name}.txt"
    self.cmds.append(
      f"python {workflow.basedir}/scripts/process_read.py --stats {{output.stats:q}} --trim-cigar {{input.bam:q}} > {{output.bam:q}}")

    return True

class FilterReadLength(Filter):
  def __init__(self):
    super().__init__("read_length")

  def _process(self):
    opts = []
    for k in ("min", "max"):
      if k in self.config and self.config[k]:
        opts.append(f"--{k}-read-length {self.config[k]}")
    if not opts:
      return False

    self.params["opts"] = opts
    self.output["stats"] = f"results/bam/filtered-{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}_stats/{self.filter_name}.txt"
    self.cmds.append(
      f"python {workflow.basedir}/scripts/process_read.py {{params.opts:q}} {{input.bam:q}} > {{output.bam}}")

    return True


class FilterMultimapper(Filter):
  def __init__(self):
    super().__init__("multimapper")

  def _process(self):
    if not self.config:
      return False

    self.cmds.append(
      f"python {workflow.basedir}/scripts/remove_multimapper.py {{input.bam:q}} > {{output.bam:q}}")

    return True


class FilterAdapterOverlap(Filter):
  def __init__(self):
    super().__init__("adapter_overlap")

  def _process(self):
    self.input["ref_fasta"] = REF_FASTA
    opts = []
    conf = self.config
    pep_conf = pep.config["qutrna2"]
    pep_pairs = {
      "linker5": "--five-adapter",
      "linker3": "--three-adapter"
    }
    for key, opt in pep_pairs.items():
      if pep_conf.get(key):
        opts.append(f"{opt} {pep_conf[key]}")
    conf_pairs = {
      "five_linker_overlap": "--five-adapter-overlap",
      "trna_overlap": "--trna-overlap",
      "three_linker_overlap": "--three-adapter-overlap",
    }
    for key, opt in conf_pairs.items():
      if conf.get(key):
        opts.append(f"{opt} {conf[key]}")
    if not opts:
      return False

    self.params["opts"] = opts
    self.output["stats"] = f"results/bam/filtered-{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}_stats/{self.filter_name}.txt.gzip"
    self.cmds.append(
      f"python {workflow.basedir}/scripts/read_overlap.py --fasta {{input.ref_fasta}} {{params.opts}} --stats {{output.stats:q}} {{input.bam:q}} > {{output.bam:q}}")

    return True

if READS == "bam":
  Filter.OUTPUT.append("data/bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam")
elif READS == "fastq":
  Filter.OUTPUT.append("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd/{BC}.sorted.bam")

if "filter" in config:
  _FILTERS = {i.filter_name: i for i in [c() for c in Filter.__subclasses__()]}
  for _filter_name in config["filter"]:
    _FILTERS[_filter_name].create_rule()

def _filter_final_input(_):
  return Filter.OUTPUT[-1]

rule filter_final:
  input: _filter_final_input
  output: "results/bam/final/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam"
  shell: """
    cp {input} {output}
  """
