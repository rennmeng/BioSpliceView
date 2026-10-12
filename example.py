#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run BioSpliceView on the bundled mini dataset.

Equivalent command line (run from the project root):

    python main.py \
        -r examples/mini_dataset/rmats \
        -g examples/mini_dataset/mm10_mini.fa \
        -t examples/mini_dataset/refGene_mini.gtf \
        -o results_example --job-id mini \
        --fdr 0.05 --psi 0.1 --minjc 10

Defaults below point to the bundled example files; edit them to analyse your
own data.
"""

from scripts.cli import main

# --- inputs (defaults: bundled example dataset) ---------------------------
RMATS_DIR = "examples/mini_dataset/rmats"
GENOME = "examples/mini_dataset/mm10_mini.fa"
GTF = "examples/mini_dataset/refGene_mini.gtf"

# --- output & parameters ---------------------------------------------------
OUTPUT_DIR = "results_example"
JOB_ID = "mini"
FDR = "0.05"
PSI = "0.1"
MINJC = "10"

ARGS = [
    "-r", RMATS_DIR,
    "-g", GENOME,
    "-t", GTF,
    "-o", OUTPUT_DIR,
    "--job-id", JOB_ID,
    "--fdr", FDR,
    "--psi", PSI,
    "--minjc", MINJC,
]

if __name__ == "__main__":
    main(ARGS)
