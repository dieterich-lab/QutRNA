from abc import ABC
from snakemake.io import temp

global READS
global REF_FASTA

# Container for applied filters (filter name).
FILTERS_APPLIED = []

class Filter(ABC):
  """
  Represents a filter class that is successively applied on BAM files.
  """

  # container for chain of BAM output files.
  OUTPUT = []

  def __init__(self, filter_name):
    self.filter_name = filter_name
    # add more dependencies if necessary
    self.input = {}
    # add more output if necessary
    self.output = {
      "bam": f"results/bam/filtered-{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}.sorted.bam"
    }
    # implementing class must populate this with a command that uses snakemake: input.bam, output.bam and log
    self.cmds = []
    # parameters to be populated to shell command
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
      self.input["bam"] = self.OUTPUT[-1]
      self.input["bai"] = self.input["bam"] + ".bai"
      self.OUTPUT.append(self.output["bam"])
      self.output["bam"] = temp(self.output["bam"])
      FILTERS_APPLIED.append(self.filter_name)
      rule:
        name: f"filter_{self.filter_name}"
        input: **self.input
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
      f"python {workflow.basedir}/scripts/bam_utils.py filter --min-as `sed -n '2p' {{input.cutoff:q}}` --output {{output.bam:q}} {{input.bam:q}}")

    return True

class FilterSamtools(Filter):
  def __init__(self):
    super().__init__("samtools")

  def _process(self):
    input_ = "{input.bam:q}"

    if self.config["filter"]:
      self.params["filter"] = self.config["filter"]
      self.cmds.append(f"samtools view {{params.filter}} -b {input_}")
      # any subsequent command will use the pipe and stdin
      input_ = "/dev/stdin"

    if "calmd" in config and config["calmd"]:
      self.input["ref"] = REF_FASTA
      self.cmds.append(f"samtools calmd -b {input_}")

    if self.cmds:
      # if any commands defined, add the output redirect
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
      f"python {workflow.basedir}/scripts/bam_utils.py filter --stats {{output.stats:q}} --trim-cigar --output {{output.bam:q}} {{input.bam:q}}")

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
      f"python {workflow.basedir}/scripts/bam_utils.py filter --stats {{output.stats:q}} {{params.opts}} --output {{output.bam:q}} {{input.bam:q}}")

    return True


class FilterAlignmentLength(Filter):
  def __init__(self):
    super().__init__("alignment_length")

  def _process(self):
    opts = []
    for k in ("min", "max"):
      if k in self.config and self.config[k]:
        opts.append(f"--{k}-alignment-length {self.config[k]}")
    if not opts:
      return False

    self.params["opts"] = opts
    self.output["stats"] = f"results/bam/filtered-{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}_stats/{self.filter_name}.txt"
    self.cmds.append(
      f"python {workflow.basedir}/scripts/bam_utils.py filter --stats {{output.stats:q}} {{params.opts}} --output {{output.bam:q}} {{input.bam:q}}")

    return True


class FilterMultimapper(Filter):
  def __init__(self):
    super().__init__("multimapper")

  def _process(self):
    if not self.config:
      return False

    self.cmds.append("samtools sort -n {input.bam:q}")
    self.cmds.append(f"python {workflow.basedir}/scripts/bam_utils.py filter-multimapper /dev/stdin")
    self.cmds.append("samtools sort -o {output.bam:q} /dev/stdin")

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
    # FIXME gzip
    self.output["stats"] = f"results/bam/filtered-{self.filter_name}/sample~{{SAMPLE}}/subsample~{{SUBSAMPLE}}/{{BC}}_stats/{self.filter_name}.txt.gzip"
    self.cmds.append(
      f"python {workflow.basedir}/scripts/bam_utils.py adapter-overlap --fasta {{input.ref_fasta}} {{params.opts}} --stats {{output.stats:q}} --output {{output.bam:q}} {{input.bam:q}}")

    return True


# Starting BAMs for filtering pipeline.
# Depends on read input BAMs or FASTQ.
if READS == "bam":
  Filter.OUTPUT.append("data/bam/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam")
elif READS == "fastq":
  Filter.OUTPUT.append("results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/orient~fwd/{BC}.sorted.bam")


# apply filters as defined in the config
if "filter" in config:
  _FILTERS = {i.filter_name: i for i in [c() for c in Filter.__subclasses__()]}
  for _filter_name in config["filter"]:
    if config["filter"][_filter_name]["apply"]:
      _FILTERS[_filter_name].create_rule()


# Final BAMs coorespond to last filter output.
# Intermediate BAMs are temporary.
rule filter_final:
  input: lambda wildcards: Filter.OUTPUT[-1]
  output: "results/bam/final/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.sorted.bam"
  shell: """
    cp {input:q} {output:q}
  """