# QutRNA

A pipeline for the alignment and detection of RNA modifications in tRNAs from Oxford Nanopore direct RNA sequencing reads.

The pipeline consists of two main steps:

1. Local alignment of tRNAs reads with parasail 
2. Detection of modifications with JACUSA2

The most convenient way to use [QutRNA](https://github.com/dieterich-lab/QutRNA) 
is to clone the entire repository.
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
Then, we compute alignments in fwd orientation with all original reads and reverse orientation. We classify alignments in the fwd orientation as true and the ones with reverse read orientation as false.

Alignment precision is then defined by TP/(TP+FP) and alignment recall by TP/(TP+FN) according to some score threshold *t*. TP: true positive, FP: false positive and FN: false negative.

We calculate an optimal threshold for a given precision and filter mapped reads accordingly.

We employ the alignment strategy separately on reads that are designated by Guppy as *fastq_pass* or *fastq_fail* and merge the results subsequently.

### Installation & Requirements:

* [R 4.3](https://www.r-project.org/)
* pheatmap
* yardstick
* bash, awk
* [samtools 1.16](https://www.htslib.org/)
* ([parasail v2.5](https://github.com/jeffdaily/parasail/archive/refs/tags/v2.5.tar.gz))

We provide a YAML file to create a [conda](https://docs.conda.io/en/latest/) environment with all necessary software except [parasail](https://github.com/jeffdaily/parasail/).

Create a conda environment with :
```bash
conda env create -n qutrna-align -f alignment_conda.yaml
conda activate qutrna-align
```

#### parasail

Check installation instructions at (https://github.com/jeffdaily/parasail/). 
If you install and compile from sources,
we recommend to add `-DBUILD_SHARED_LIBS=OFF` to your `cmake` command.
Otherwise, make sure that `parasail` is in the `$PATH` and the dynamic library can be found.


### 1. Generate FASTQ

Process ONT raw FASTQ files and reverse sequence and/or transform **U**s to **T**s. Use script on raw *fastq_pass* or *fastq_fail* ONT directory.

```bash
Usage: ./generate_fastq.sh [ -r ] [ -t ] FASTQ-DIR-NAME
-r    Reverse sequence of called bases
-t    Tranform Us to Ts
```

### 2. Align reads

Map processed FASTQ reads to a FASTA reference sequence.
To save memory, the script will split FASTQ files into chunks of 2500 reads and merge all chunks at the end.

```bash
Usage: ./parasail_align.sh  -b <OUTPUT-BAM> -f <FASTA> [ -t <THREADS> ] FASTQ
-b <OUTPUT-BAM>   Filename for aligned reads
-f <FASTA>        Filename of fasta sequence reference
-t <THREADS>      Number of threads, default=1
```

### 3. Retain highest scoring alignments

Given a BAM file, this script will filter mapped reads with a minimal alignment score and 
retain only the highest scoring alignments for each read - per read there might be more than one high scoring alignment.

```
Usage: ./retain_highest_alignments.sh  -b <FILTERED-BAM> [ -s <MIN_SCORE> ] [ -t <THREADS> ] BAM
-b <OUTPUT-BAM>   Filename for filtered reads
-s <MIN-SCORE>    Minimal alignment score, default=10
-t <THREADS>      Number of threads, default=1
```

### 4. Plot: Optimial alignment score

Evaluate and determine optimal alignment score by random alignment of reverse sequences at a given precision.

```bash
Usage: alignment_score_cutoff.R [options] JACUSA2_SCORES
-f FORWARD, --forward=FORWARD
        Scores of reads mapping to forward

-r REVERSE, --reverse=REVERSE
        Scores of reads mapping to reverse

-o OUTPUT, --output=OUTPUT
        Output directory

-p PRECISION, --precision=PRECISION
        Precision cutoff, default=0.95

-t TITLE, --title=TITLE
        Alignment score plot title
```

### Run alignment pipeline

This script will execute all the aforementioned steps and produce the final BAM file that can be used in the [JACUSA2 analysis](#markdown-header-jacusa2-analysis) part.

```bash
Usage: ./optimize_alignment.sh  -f <FASTA> -o <OUT-DIR> [ -s <MIN-SCORE> ] [ -p <PRECISION ] [ -t <THREADS> ] IN_DIR
-f <FASTA>        Filename of fasta sequence reference
-b <OUT-DIR>      Directory for output
-s <MIN-SCORE>    Minimal alignment score, default=10
-s <PRECISION>    Precision cutoff, default=0.95
-t <THREADS>      Number of threads
```

### Examples

We provide a data set to explore the pipeline. Necessary data for *S.pombe* was deposited in the repository in the `data`directory:

* *S.pombe* tRNAAsp IVT and
* *S.pombe* tRNAAsp IVT Q.

Run the exemplary alignment analysis with: `alignment/run_examples.sh`.

The output will be in `output/alignment`.

## JACUSA2 analysis

[JACUSA2](https://github.com/dieterich-lab/JACUSA2) compares paired samples and detects positions, where 2 conditions differ in the base composition, or the number of INDELs.
[JACUSA2](https://github.com/dieterich-lab/JACUSA2) employs DirichletMultinomial and BetaBinomial distributions to model base compositions and INDELs, respectively. 

We calculate mismatch, insertion, and deletion scores with [JACUSA2](https://github.com/dieterich-lab/JACUSA2) and combine them to create the Mis+Del+Ins score.

Finally,we provide diagnostic plots where the score is plotted as a heatmap and known modifications are displayed to verify the predictions.

### Installation & Requirements

* Java >= 11
* [JACUSA2 2.0.4](https://github.com/dieterich-lab/JACUSA2)
* [JACUSA2helper 1.9.9600](https://github.com/dieterich-lab/JACUSA2helper)

Create a conda environment with all requirements with:
```bash
conda env create -n qutrna-jacusa2 -f JACUSA2_conda.yaml
conda activate qutrna-jacusa2
```

### Run JACUSA2 analysis and plot

```bash
Usage: ./analysis.sh  -f <FASTA> -o <OUT-DIR> -m <MODS> BAMS1(,) BAMS2(,)
-o <OUT-DIR>      Directory for output
-f <FASTA>        Filename of fasta sequence reference
-m <mods>         CSV file with known modifications (position: 0-index)
```
The results of JACUSA2 are in `<OUT-DIR>/JACUSA2.out`.
The processed scores can be explored in `<OUT-DIR>/scores.csv`.
Finally, PDFs`<OUT-DIR>/main.pdf` and `<OUT-DIR>/small.pdf` provide a consise visualization of scores and known modifications.

### Examples

It is required to [Run alignment pipeline](#markdown-header-run-alignment-pipeline) first, before continuing with JACUAS2 analysis: `JACUSA2/run_examples.sh`.

The results are in `output/JACUSA/...`.

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
