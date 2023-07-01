#!/bin/bash
set -e

IN_BAM=""
OUT_BAM=""
MIN_SCORE=10
THREADS=1

usage() {
  echo "Usage: $0  -b <FILTERED-BAM> [ -s <MIN_SCORE> ] [ -t <THREADS> ] BAM" 1>&2
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

while getopts "b:s:t:" options
do
  case "${options}" in
    b)
      OUT_BAM="${OPTARG}"
      ;;
    s)
      MIN_SCORE="${OPTARG}"
      ;;
    t)
      THREADS="${OPTARG}"
      ;;
    *)
      exit_with_error
      ;;
  esac
done


IN_BAM=${@:$OPTIND:1}
if [ "$IN_BAM" = "" ]
then
  exit_with_error "BAM is required!"
fi
if [ ! -e "$IN_BAM" ]
then
  exit_with_error "BAM '$IN_BAM' cannot be accessed!"
fi
if [ "$OUT_BAM" = "" ]
then
  exit_with_error "FILTERED-BAM is required!"
fi
if [ -e "$OUT_BAM" ]
then
  exit_with_error "Output already exists '$OUT_BAM'!"
fi
if [ "$IN_BAM" = "$OUT_BAM" ]
then
  exit_with_error "Input and output cannot be the same!"
fi

{
  samtools view -H $IN_BAM && \
  samtools view $IN_BAM | \
    sort -k1,1 -k5,5r | \
    awk -v MIN_SCORE="$MIN_SCORE" \
      ' BEGIN { LAST_READ="" ; SCORE="" } ;
        $5>=MIN_SCORE {
                        if (LAST_READ!=$1) {
                          LAST_READ=$1 ;
                          SCORE=$5 ;
                          print ;
                        } else if ($5==SCORE) {
                          print
                        }
                      } '
} | \
samtools view -bS /dev/stdin | \
samtools sort -@ $THREADS -m 2G -o $OUT_BAM /dev/stdin
samtools index $OUT_BAM
