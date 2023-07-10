#!/bin/bash

# Run analysis for S.pombe tRNAAsp IVT vs IVT Q.
# MAKE SURE YOU RUN "alignment/run_examples.sh" FIRST!

./analysis.sh -f ../data/S.pombe.fasta \
              -o S.pombe_tRNAAsp_IVT-Q \
              -m ../data/S.pombe_mods.csv \
              ../output/S.pombe_tRNAAsp_IVT/final_merged.bam \
              ../output/S.pombe_tRNAAsp_IVT_Q/final_merged.bam
