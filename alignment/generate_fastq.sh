#!/bin/bash

# Concatenate FASTQ reads, transform Us to Ts and optionally reverse sequence.

set -e

DNAME=""  # FASTQ directory
U2T=0     # transform Us to Ts
REVERSE=0 # Reverse FASTQ bases

usage() {
  echo "Usage: $0 [ -r ] [ -t ] FASTQ-DIR-NAME" 1>&2
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

DNAME=${@:$OPTIND:1}
if [ "$DNAME" = "" ]
then
  exit_with_error "You must provide a directory with FASTQ.gz files!"
fi
if [ ! -d "$DNAME" ]
then
  exit_with_error "Directory '$DNAME' does not exist!"
fi

if [ "$REVERSE" -eq 1 ]
then
  cat $DNAME/*.gz | \
    gunzip -c | \
    while read LINE
    do
      echo `echo $LINE | cut -d" " -f1`"_rev" && \
      read LINE && echo $LINE | rev && \
      read LINE && echo $LINE && \
      read LINE && echo $LINE
    done | \
    awk -v U2T="$U2T" \
      ' NR%4==2 { if (U2T==1) { gsub("U", "T", $0) } ; print ; next } ;
                { print } '
else
  cat $DNAME/*.gz | \
    gunzip -c | \
    awk -v U2T="$U2T" \
    ' NR%4==1 { $0=gensub(/^([^ ]+).+/, "\\1", "g") ; print ; next } ;
      NR%4==2 { if (U2T==1) { gsub("U", "T", $0) } ; print ; next } ;
              { print } '
fi
