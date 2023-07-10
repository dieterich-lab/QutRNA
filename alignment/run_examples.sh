#!/bin/bash
set -e

# Run example alignment for S.Pombe
mkdir -p ../output/S.pombe_tRNAASP_IVT
mkdir -p ../output/S.pombe_tRNAASP_IVT_Q

# optimize alignment for S.pombe tRNAAsp IVT
./optimize_alignment.sh \
  -f "../data/S.pombe.fasta" \
  -o "../output/S.pombe_tRNAAsp_IVT" \
  "../data/S.pombe_tRNAAsp_IVT"

# optimize alignment for S.pombe tRNAAsp IVT Q
./optimize_alignment.sh \
  -f "../data/S.pombe.fasta" \
  -o "../output/S.pombe_tRNAAsp_IVT_Q" \
  "../data/S.pombe_tRNAAsp_IVT_Q"
