#!/bin/bash
set -e

# Run example alignment for S.Pombe

# optimize alignment for S.pombe tRNAAsp IVT
./optimize_alignment.sh \
  -f "../data/S.pombe.fasta" \
  -o "S.pombe_tRNAAsp_IVT" \
  "../data/S.pombe_tRNAAsp_IVT"

# optimize alignment for S.pombe tRNAAsp IVT Q
./optimize_alignment.sh \
  -f "../data/S.pombe.fasta" \
  -o "S.pombe_tRNAAsp_IVT_Q" \
  "../data/S.pombe_tRNAAsp_IVT_Q"
