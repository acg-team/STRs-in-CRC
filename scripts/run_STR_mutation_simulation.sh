#!/usr/bin/env bash
set -euo pipefail

UNFILTERED_REGIONS="https://drive.google.com/uc?id=1b9zk34yIZypr_9sJRReAf0PMBaj9DWJP"
REGIONS="../data/dummy_str_regions_filt.csv"
MSS="../data/dummy_data_MSS.csv"
MSI="../data/dummy_data_MSI.csv"
OUTPUT="../data/dummy_STR_mutation_data.csv"

python3 filter_STR_regions.py \
    -r ${UNFILTERED_REGIONS} \
    -o ${REGIONS} \
    -n 1000

python3 simulate_STR_mutation_data.py \
    -r ${REGIONS} \
    -o ${MSS} \
    -n 25 \
    -p .0142 \
    -b .5 \
    -z 0.555 \
    -s 42

python3 simulate_STR_mutation_data.py \
    -r ${REGIONS} \
    -o ${MSI} \
    -n 25 \
    -p .0518 \
    -b .91 \
    -z 0.354 \
    -s 42

python3 merge_dummy_datasets.py \
    --mss ${MSS} \
    --msi ${MSI} \
    -o ${OUTPUT}

rm ${MSS}
rm ${MSI}