#!/bin/bash

# Align reads with parasail

set -e

FASTQ=""  # FASTQ reads
FASTA=""  # FASTA sequence reference (indexed)
BAM=""    # Output BAM file
THREADS=1 # Number of threads to be used

usage() {
  echo "Usage: $0  -b <OUTPUT-BAM> -f <FASTA> [ -t <THREADS> ] FASTQ" 1>&2
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

while getopts "b:f:t:" options
do
  case "${options}" in
    b)
      BAM="${OPTARG}"
      ;;
    f)
      FASTA="${OPTARG}"
      ;;
    t)
      THREADS="${OPTARG}"
      ;;
    *)
      exit_with_error
      ;;
  esac
done

FASTQ=${@:$OPTIND:1}
if [ "$FASTQ" = "" ]
then
  exit_with_error "FASTQ required!"
fi
if [ ! -e "$FASTQ" ]
then
  exit_with_error "Missing FASTQ '$FASTQ'!"
fi
if [ "$FASTA" = "" ]
then
  exit_with_error "FASTA required!"
fi
if [ ! -e "$FASTA" ]
then
  exit_with_error "Missing FASTA '$FASTA'!"
fi
if [ "$BAM" = "" ]
then
  exit_with_error "OUTPUT-BAM required!"
fi
if [ -e "$BAM" ]
then
  exit_with_error "BAM already exists '$BAM'!"
fi

TMP_OUT="$(mktemp $BAM.XXXXXXXXXX)" || { exit_with_error "Failed to create temp file!" ; }

gunzip -c $FASTQ | split -l 10000 --filter 'gzip -c > $FILE.fastq.gz' /dev/stdin $TMP_OUT.part_

for PART in $TMP_OUT.part_*.fastq.gz
do
  parasail_aligner -a sw_trace_striped_sse41_128_16 \
                   -M 1 \
                   -X 1 \
                   -e 1 \
                   -o 1 \
                   -c 20 \
                   -x \
                   -d \
                   -t $THREADS \
                   -O SAMH \
                   -f ${FASTA} \
                   -g ${PART%.fastq.gz}.sam \
                   -q $PART
  rm $PART
done

cat $TMP_OUT.part_*.sam | \
  gawk -v OFS="\t" \
    ' BEGIN { HEADER=0 }
      $0 ~ /^@/ { if (HEADER==0) { print } ; next }
      $0 !~ /^@/ { HEADER=1 }
      $0 ~ /AS:i:([0-9]+)/ {
                             AS=int(gensub(/.+AS:i:([0-9]+).+/, "\\1", "g")) ;
                             (AS <= 255) ? $5=AS : $5=255 ;
                             print ;
                             next } ;
                           { print } ' | \
  samtools view -bS /dev/stdin | \
  samtools calmd /dev/stdin $FASTA | \
  samtools sort -@ $THREADS -m 2G -o $BAM /dev/stdin && samtools index $BAM

rm $TMP_OUT.part_*.sam
rm $TMP_OUT
