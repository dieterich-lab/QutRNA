#!/bin/bash

snakemake -c 1 -f ../workflow/Snakefile --pep data.yaml --configfile=analysis
