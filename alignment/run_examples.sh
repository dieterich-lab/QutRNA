#!/bin/bash
set -e

# run S.pombe TODO
#./optimize_alignment.sh \
#  -f "../data/S.pombe.fasta" \
#  -o "S.pombe_tRNAAsp_IVT" \
#  "../data/S.pombe_tRNAAsp_IVT"

./optimize_alignment.sh \
  -f "../data/S.pombe.fasta" \
  -o "S.pombe_tRNAAsp_IVT_Q" \
  "../data/S.pombe_tRNAAsp_IVT_Q"
