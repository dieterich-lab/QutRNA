# QutRNA

A pipeline for the alignment and detection of RNA modifications in tRNAs from Oxford Nanopore direct RNA sequencing reads.

The pipeline consists of three steps:

1. Local alignment of tRNAs reads with [parasail](https://github.com/jeffdaily/parasail/)
2. Detection of modifications with [JACUSA2](https://github.com/dieterich-lab/JACUSA2/)
3. Visualization of RNA modifications

The most convenient way to use [QutRNA](https://github.com/dieterich-lab/QutRNA) 
is to clone the entire repository and install dependencies with [conda](https://docs.conda.io/en/latest/).

Go to your desired <LOCAL-DIRECTORY> and clone the repository:
```bash
cd <LOCAL-DIRECTORY>
git clone https://github.com/dieterich-lab/QutRNA
```

## Alignment

Our alignment workflow employs an optimal local alignment strategy using the implementation as provided by the [parasail](https://github.com/jeffdaily/parasail/) software.
In summary, optimal alignments may produce drastically different tRNA read mappings and are superior to heuristic alignments.

Our strategy to assess the statistical significance is rooted in a simulation-based approach, which produces random alignments.
Briefly, we reverse input sequences.
Then, we compute alignments in **forward** orientation with all original reads and **reverse** orientation. We classify alignments in the **forward** orientation as true and the ones with **reverse** read orientation as false.

Alignment precision is then defined by TP/(TP+FP) and alignment recall by TP/(TP+FN) according to some score threshold *t*, where TP: true positive, FP: false positive, and FN: false negative.

We calculate an optimal threshold for a given precision and filter mapped reads accordingly.

We employ the alignment strategy separately on reads that are designated by Guppy as *fastq_pass* or *fastq_fail* and merge the results subsequently.

## Installation & Requirements:

Currently, the pipeline requires an X86\_64 architecture. (We are working on supporting OSX.)

* [R 4.3](https://www.r-project.org/)
* Java >= 11
* [JACUSA2 2.0.4](https://github.com/dieterich-lab/JACUSA2)
* [JACUSA2helper 1.9.9600](https://github.com/dieterich-lab/JACUSA2helper)
* [samtools](https://www.htslib.org/)
* [parasail v2.6.2](https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz)
* (check `conda.yaml` in the repository for a complete)

We provide a YAML file to create a [conda](https://docs.conda.io/en/latest/) environment with all necessary software with the exception of [parasail](https://github.com/jeffdaily/parasail/) (no package in conda).

Create a conda environment with:
```console
conda env create -n qutrna -f <LOCAL-REPOSITORY>/conda.yaml
conda activate qutrna
```

If removing white space from final plots is desired, a workking TEX environment is required. Unfortunatelly, the TEX environment available via conda is not applicable.

### parasail

Follow [parasail compiling](https://github.com/jeffdaily/parasail?tab=readme-ov-file#autotools-build) instructions. 

If you install and compile from sources and you use conda, it is imperative to compile parasail within the conda environment. 

1. Create and activate conda environment (see above)
2. Download parasail [v2.6.2.tar.gz](https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz)
3. Compile and install. Replace `<DESTINATION>` with your desired path for parasail:

```console
tar -zvpf v2.6.2.tar.gz
cd v2.6.2
autoreconf -fi
configure --prefix=<DESTINATION>
make
make install
```

Make sure, [parasail](https://github.com/jeffdaily/parasail) is installed in `$PATH`. 


### Snakemake workflow

The workflow is implemented with [snakemake](https://github.com/snakemake/snakemake) and encompasses:

* tRNA alignment with [parasail](https://github.com/jeffdaily/parasail), 
* RNA modification detection with [JACUSA2](https://github.com/dieterich-lab/JACUSA2),
* and visualization in Sprinzl coordinates.

The workflow can be configered with YAML files:
* analysis.yaml : analysis specific config, e.g.: parameters of tools.
* data.yaml : data specific config, e.g.: reference sequence, sample description.


#### Config: analysis

Custom parameters and plots can be defined in a custom `analysis.yaml`.
Otherwise, default values are used for the tools.

A minimal custom `analysis.yaml` is required to set the `precision` and `minimal alignment score`:
```yaml
params:
  precision: 0.95
  min_aln_score: 10
```

##### parasail

[parasail](https://github.com/jeffdaily/parasail) is used to perform fast local alignment of reads against reference sequences of tRNAs. 
The following defaults values for parasail are set in QutRNA and can be overritten in a custom `analysis.yaml`: 

```yaml
[...]
parasail:
  opts: -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d
  batch_size: 1000
  threads: 1
  lines: 0
[...]
```
`opts` defines parasail specific command line options (check [parasail](https://github.com/jeffdaily/parasail) for details)
`batch_size` defines the batch size of reads, the parameter influences main memory requirements (check [parasail](https://github.com/jeffdaily/parasail) for details)
`threads` sets the number of parallel threads to use. Adjust to your local computing machine
If `lines` is > 0, each FASTQ will be split in files with the number of lines. Make sure that the number is divisible by 4! If splitting input is desired, choose depending on the number of raw reads and pick a reasonably high number of lines (> 10000).


##### JACUSA2

[JACUSA2](https://github.com/dieterich-lab/JACUSA2) is used to detect RNA modifications by means of scoring mismatches and INDELs.

The following default values are defined for JACUSA2 and can be overwritten in a custom `analysis.yaml`:

```yaml
[...]
jacusa2:
  opts:-m 1 -q 1 -p 1 -D -i -a D,Y -P1 FR-SECONDSTRAND -P2 FR-SECONDSTRAND
  min_cov: 10
  threads: 2
[...]
```

`opts` defines JACUSA2 specific command line options (check [JACUASA2](https://github.com/dieterich-lab/JACUSA2) for details)
`min_cov` defines the minimum number of reads at a given position in EACH BAM file to consider for RNA modification detection
`threads` sets the number of parallel threads to use. Adjust to your local computing machine

##### cmalign

[cmalign] is used to perform secondary structure alignment of tRNAs which is required to for Sprinzl coordinates.

```yaml
[...]
cmalign:
  opts: --notrunc --nonbanded -g
  threads: 2
```

`opts` defines cmalign specific command line options (check [cmalign]() for details )
`threads` sets the number of parallel threads to use. Adjust to your local computing machine

The underlying covariance model is data specific and therefore defined in `data.yaml`.


#### Config: data

Every analysis requries a `data.yaml` where details about the underlying data are defined.

In the following an extensive example of `data.yaml` with descriptions is presented:

```yaml
pep_version: 2.0.0             # [Required] by Snakemake

sample_table: sample_table.tsv # [Required] Filename of sample description

qutrna:
  output_dir: <OUTPUT-DIR> # [Required] Where QutRNA output will be written to
  cm: <PATH-TO-CM>         # [Required] Path to custom covariance model 
  ref_fasta: <PATH-TO-REF-FASTA>    # [Required] Path to reference sequence
  ref_fasta_prefix: "Homo_sapiens_" # (Optional) Prefix to remove from sequence ID in visualization
  coords: sprinzl                   # [Required] Possible values are 'seq' or 'sprinzl'.
                                    #            WARNING! Coordinates of RNA modifications must be compatible
  mods:
    file: <PATH-TO-MODS>     # (Optional) Path to RNA modifications
    abbrev: <PATH-TO-ABBREV> # (Optional) Path to RNA modification abbreviations
  linker5: <INTEGER> # [Required] length of 5' linker in nt
  linker3: <INTEGER> # [Required] lentgh of 3' linker in nt
  remove_trnas:
    [seq1, ]         # (Optional) Sequence IDs to ignore for secondary structure alignment and visualization
                     #            However, those Sequence IDs will be used for the alignment.
  contrasts:         # [Required] Define combinations of conditions to analyse.
                     #            Multiple comparisons are possible.
    - cond1: <CONDITION1> # must match to condition in sample_table.tsv
      cond2: <CONDITION2> # must match to condition in sample_table.tsv
```

#### Sample table

Sample description `sample_table.tsv` must be TAB-separated and contain the following columns:


| condition | sample_name | subsample_name | base_calling | fastq\|bam |
| --------- | ----------- | -------------- | ------------ | ---------- |
| ...       | ...         | ...            | ...          | ...        |

`condition`: Name of the respective condition. Will be used in `data.yaml` to define contrasts.
`sample_name`: Name of the sample. A sample can consist of muliple FASTQ or BAM files. Samples with the same `sample_name` will be merged before RNA modification detection. 
`subsample_name`: Name of the subsample (see above). A sample can consist of multipe subsamples, e.g.: tech. replicates.
`base_calling`: Metainformation describing base calling for the respective row.
                Possible values are: 'pass', 'fail', 'merged' or 'unknown'.
`fastq` or `bam`: Only one column is permitted. In either case, the absolute path of the sequencing data is expected. If `fastq` is used, then the path to GZIPPED FASTQ sequencing reads is expected. If column `bam` is provided, then mapped reads in the BAM file format are expected.

#### RNA modifications

The file with RNA modification information is expected to be TAB-separated and contain the following columns:

| trna                      | pos         | mod            |
| ------------------------- | ----------- | -------------- |
| (should match ref. fasta) | ...         | ...            |

The file with abbreviations with RNA modifications is expected to be TAB-separated and contain the following columns:

| short_name                  | abbrev |
| --------------------------- | ------ |
| (should match column 'mod') | ...    |

#### Options for visualization

JACUAS2 scores are visualized with the script `workflow/scripts/plot_score.R`.
Plots can be added to `analysis.yaml` with custom plot options:

```yaml
[...]
plots:
  - id: <PLOT-ID>                 # [Required] name of the directory for the plot, in results/plots/<PLOT-ID>
    trnas: isoacceptor|isodecoder # [Required] How to group tRNA in the output
    opts: "--sort"                # (Optional) Command line options for plot_score.R
                                  #            Here: sort tRNAs by median read coverage.
[...]
```

The script `workflow/scripts/plot_score.R` supports the following options that can be added to `opts`:
```
--title=TITLE         Title for each plot. The following patterns can be used: {anti_codon} or {amnio_acid}.

--hide_varm           Hide variable arm
--hide_mods           Hide RNA modification annotation, if available.
                      
--show_introns        Show introns
--show_coverage       Show a barplot with median read coverage for each tRNA

--positions=POSITIONS Restrict output to POSITIONS separated by ","

--crop                Crop final pdf to remove white space (requires TEX environment)
```

Multiple plots, each with unique `<PLOT-ID>` and different options are supported.

##### Crop output

The final heatmap plot will be surrounded by white space that can not be removed by adjusting parameters in ggplot.
Although, it is possible to use `pdfcrop.pl` to remove the white space, a working TEX environment is required.
The script is included in the conda environment of QutRNA but unfortunatelly the requirements cannot be satisfied within conda.
It is beyond this introduction, to go in to details how to setup up TEX. Please check your OS or distribution how to setup TEX.

If you manage to setup TEX, add `--crop` to plot options to remove white space from final heatmap plots.

### Executing the workflow

Setup `analysis.yaml`, `data.yaml`, and `sample_table.tsv`.
Replace `<QUTRNA>` with the directory where you cloned the repository, add paths to the YAML files and start the workflow with:

```console
snakemake -c 1 -f <QUTRNA>/workflow/Snakefile --pep data.yaml --configfile=analysis.yaml
```

In case, you are only interested in the parasail alignment, add the intermediate target `alignment`:

```console
snakemake -c 1 -f <QUTRNA>/workflow/Snakefile --pep data.yaml --configfile=analysis.yaml alignment
```

#### Output

The output of the pipeline can  be found in the output directorty that you provided in `data.yaml`.

TODO results/plots/

### Examples

We provide a data set to explore the pipeline. Necessary data for *S.pombe* was deposited in the repository in the `data`directory:

* *S.pombe* tRNAAsp IVT and
* *S.pombe* tRNAAsp IVT Q.

Run the pipeline with: `example/run_example.sh`.

The output will be in `example/output`.

# How to cite

...

# License

See LICENSE for details

# References

* [RNA modification mapping with JACUSA2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02676-0)

## Software
* [parasail](https://github.com/jeffdaily/parasail/)
* [JACUSA2](https://github.com/dieterich-lab/JACUSA2)
* [JACUSA2helper](https://github.com/dieterich-lab/JACUSA2helper)
