#!/bin/bash

# TODO
set -e

IN_DIR=""
FASTA=""
OUT_DIR=""
PRECISION=0.95
MIN_SCORE=10
THREADS=1

usage() {
  echo "Usage: $0  -f <FASTA> -o <OUT-DIR> [ -s <MIN-SCORE> ] [ -p <PRECISION ] [ -t <THREADS> ] IN_DIR" 1>&2
}

exit_with_error() {
  MSG="$1"
  if [ "$MSG" != "" ]
  then
    echo "ERROR: $MSG" 1>&2
  fi
  usage
  exit 1
}

while getopts "f:o:p:s:t:" options
do
  case "${options}" in
    f)
      FASTA="${OPTARG}"
      ;;
    o)
      OUT_DIR="${OPTARG}"
      ;;
    s)
      MIN_SCORE="${OPTARG}"
      ;;
    p)
      PRECISION="${OPTARG}"
      ;;
    t)
      THREADS="${OPTARG}"
      ;;
    *)
      exit_with_error
      ;;
  esac
done

IN_DIR=${@:$OPTIND:1}
if [ "$FASTA" = "" ]
then
  exit_with_error "FASTA is required!"
fi
if [ ! -e "$FASTA" ]
then
  exit_with_error "FASTA '$FASTA' cannot be accessed!"
fi
if [ "$OUT_DIR" = "" ]
then
  exit_with_error "OUT-DIR is required!"
fi
mkdir -p $OUT_DIR

for RTYPE in pass fail
do
  ##############################################################################
  # Merge nanopore FASTQ files, transform Us to Ts, and
  #  reverse nucleotide sequence to obtain random score distribution
  ##############################################################################
  echo "Merge fwd $RTYPE FASTQ" 1>&2
  ./generate_fastq.sh -t $IN_DIR/fastq_$RTYPE | gzip -c > $OUT_DIR/fwd_$RTYPE.fastq.gz
  echo "Merge rev $RTYPE FASTQ" 1>&2
  ./generate_fastq.sh -r -t $IN_DIR/fastq_$RTYPE | gzip -c > $OUT_DIR/rev_$RTYPE.fastq.gz

  ##############################################################################
  # Use parasail to map and align reads
  ##############################################################################
  echo "Map fwd $RTYPE FASTQ" 1>&2
  ./parasail_align.sh -b $OUT_DIR/mapped_fwd_$RTYPE.bam \
    -t $THREADS -f $FASTA $OUT_DIR/fwd_$RTYPE.fastq.gz
  echo "Map rev $RTYPE FASTQ" 1>&2
  ./parasail_align.sh -b $OUT_DIR/mapped_rev_$RTYPE.bam \
    -t $THREADS -f $FASTA $OUT_DIR/rev_$RTYPE.fastq.gz

  ##############################################################################
  # Retain highest scoring alignments with a minimum alignment score of $MIN_SCORE
  ##############################################################################
  echo "Filter fwd $RTYPE BAM" 1>&2
  ./retain_highest_alignments.sh -b $OUT_DIR/mapped_fwd_$RTYPE.filtered.bam \
    -t $THREADS -s $MIN_SCORE $OUT_DIR/mapped_fwd_$RTYPE.bam
  echo "Filter rev $RTYPE BAM" 1>&2
  ./retain_highest_alignments.sh -b $OUT_DIR/mapped_rev_$RTYPE.filtered.bam \
    -t $THREADS -s $MIN_SCORE $OUT_DIR/mapped_rev_$RTYPE.bam

  ##############################################################################
  # Filter alignment by random score distribution
  ##############################################################################
  echo "Filter by random score $RTYPE $BAM" 1>&2
  samtools view $OUT_DIR/mapped_fwd_$RTYPE.filtered.bam | \
    cut -f 5 > $OUT_DIR/score_fwd_$RTYPE.txt
  samtools view $OUT_DIR/mapped_rev_$RTYPE.filtered.bam | \
    cut -f 5 > $OUT_DIR/score_rev_$RTYPE.txt

  Rscript --vanilla alignment_score_cutoff.R \
    -o $OUT_DIR/score_$RTYPE \
    -p $PRECISION \
    -t "tRNA Alignment Score $RTYPE distributions" \
    --forward $OUT_DIR/score_fwd_$RTYPE.txt --reverse $OUT_DIR/score_rev_$RTYPE.txt


  samtools view --write-index -b -q `cat $OUT_DIR/score_$RTYPE/CutOff.txt` \
    -o $OUT_DIR/final_$RTYPE.bam $OUT_DIR/mapped_fwd_$RTYPE.filtered.bam
done

echo "Merge reads (pass, fail) BAM" 1>&2
samtools merge -r --write-index $OUT_DIR/final_merged.bam \
  $OUT_DIR/final_pass.bam $OUT_DIR/final_fail.bam
