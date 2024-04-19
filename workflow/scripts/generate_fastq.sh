#!/bin/bash

# Reads FASTQ reads, transform Us to Ts and optionally reverse sequence.

set -e

FNAME=""  # FASTQ file
U2T=0     # (-t) transform Us to Ts
REVERSE=0 # (-r) Reverse FASTQ bases

usage() {
  echo "Usage: $0 [ -r ] [ -t ] FASTQ" 1>&2
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

while getopts "rt" options
do
  case "${options}" in
    r)
      REVERSE=1
      ;;
    t)
      U2T=1
      ;;
    *)
      exit_abnormal
      ;;
  esac
done

FNAME=${@:$OPTIND:1}
if [ "$FNAME" = "" ]
then
  exit_with_error "You must provide a FASTQ.gz !"
fi
if [ ! -e "$FNAME" ]
then
  exit_with_error "File '$FNAME' does not exist!"
fi
if [ "$REVERSE" -eq 1 ]
then
  cat $FNAME | \
    gunzip -c | \
    while read LINE
    do
      echo `echo $LINE | cut -d" " -f1`"_rev" && \
      read LINE && echo $LINE | rev && \
      read LINE && echo $LINE && \
      read LINE && echo $LINE
    done | \
    gawk -v U2T="$U2T" \
      ' NR%4==2 { if (U2T==1) { gsub("U", "T", $0) } ; print ; next } ;
                { print } '
else
  cat $FNAME | \
    gunzip -c | \
    gawk -v U2T="$U2T" \
    ' NR%4==1 { $0=gensub(/^([^ ]+).+/, "\\1", "g") ; print ; next } ;
      NR%4==2 { if (U2T==1) { gsub("U", "T", $0) } ; print ; next } ;
              { print } '
fi
