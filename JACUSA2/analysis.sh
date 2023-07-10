#!/bin/bash
set -e

# Run JACUSA2 analysis pipeline.

SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

# JAR or wrapper to use for analysis
JACUSA2="jacusa2"
# Default JACUSA2 options
JACUSA2_OPTS=" -m 1 -q 1 -c 4 -p 1 -D -i -a D,Y -P1 FR-SECONDSTRAND -P2 FR-SECONDSTRAND "

FASTA=""   # FASTA sequence reference
MODS=""    # CSV file with known modifications (position: 0-index)
OUT_DIR="" # output directory
THREADS="1" # number of threads to use for computation
BAMS1=""   # ","-separated set of BAM files for condition1
BAMS2=""   # ","-separated set of BAM files for condition2

usage() {
  echo "Usage: $0  -f <FASTA> -o <OUT-DIR> -m <MODS> BAMS1(,) BAMS2(,)" 1>&2
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

while getopts "f:o:m:t:" options
do
  case "${options}" in
    f)
      FASTA="${OPTARG}"
      ;;
    o)
      OUT_DIR="${OPTARG}"
      ;;
    m)
      MODS="${OPTARG}"
      ;;
    m)
      THREADS="2"
      ;;
    *)
      exit_with_error
      ;;
  esac
done

BAMS1=${@:$OPTIND:1}
if [ "$BAMS1" = "" ]
then
  exit_with_error "BAMS1 is required!"
fi
BAMS2=${@:(($OPTIND+1)):1}
if [ "$BAMS2" = "" ]
then
  exit_with_error "BAMS2 is required!"
fi
if [ "$FASTA" = "" ]
then
  exit_with_error "FASTA is required!"
fi
if [ ! -e "$FASTA" ]
then
  exit_with_error "FASTA '$FASTA' cannot be accessed!"
fi
if [ "$MODS" = "" ]
then
  exit_with_error "MODS is required!"
fi
if [ ! -e "$MODS" ]
then
  exit_with_error "FASTA '$MODS' cannot be accessed!"
fi
if [ "$OUT_DIR" = "" ]
then
  exit_with_error "OUT-DIR is required!"
fi
mkdir -p $OUT_DIR

################################################################################
# Run JACUSA2 on BAMS(1,2)
################################################################################
JACUSA2_OUT="$OUT_DIR/JACUSA2.out"
echo "Run JACUSA2" 1>&2
$JACUSA2 call-2 $JACUSA2_OPTS -p $THREADS -r $JACUSA2_OUT $BAMS1 $BAMS2

################################################################################
# Process JACUSA2 scores
################################################################################
JACUSA2_SCORE_OUT="$OUT_DIR/scores.csv"
echo "Add scores" 1>&2
Rscript $DIR/add_scores.R -f $FASTA -m $MODS -o $JACUSA2_SCORE_OUT $JACUSA2_OUT

################################################################################
# Plot JACUSA2 score: Mis+Del+Ins
################################################################################
PLOT_OUT="$OUT_DIR/plots"
echo "Create plots" 1>&2
Rscript $DIR/plot_score.R --output $PLOT_OUT $JACUSA2_SCORE_OUT
