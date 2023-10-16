#!/usr/bin/env bash
set -euo pipefail

python3 simulate_eSTR_data.py \
    -r ../data/dummy_str_regions_filt.csv \
    -g ../data/dummy_somatic_mutation_calls_filt.csv \
    -s 42 \
    -o ../data/gene_expression/dummy_tumor_gene_expression.csv
