#!/bin/bash

# Run analysis for S.pombe tRNAAsp IVT vs IVT Q.
# MAKE SURE YOU RUN "alignment/run_examples.sh" FIRST!

./analysis.sh -f ../data/S.pombe.fasta \
              -o ../JACUSA2/S.pombe_tRNAAsp_IVT-Q \
              -m ../data/S.pombe_mods.csv \
              ../output/alignment/S.pombe_tRNAAsp_IVT/final_merged.bam \
              ../output/alignment/S.pombe_tRNAAsp_IVT_Q/final_merged.bam
