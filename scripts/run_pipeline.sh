#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/processed data/raw

python src/build_chembl_egfr_dataset.py \
  --target CHEMBL203 \
  --out data/processed/chembl_egfr_activities.tsv

python src/build_bindingdb_egfr_dataset.py \
  --uniprot P00533 \
  --cutoff 10000 \
  --out data/processed/bindingdb_egfr_activities.tsv \
  --raw-json-out data/raw/bindingdb_egfr_10k_raw.json

python src/make_conplex_files.py \
  --chembl data/processed/chembl_egfr_activities.tsv \
  --bindingdb data/processed/bindingdb_egfr_activities.tsv \
  --mutants configs/egfr_mutants.json \
  --out-affinities data/processed/egfr_affinities_long.tsv \
  --out-pairs data/processed/egfr_conplex_pairs.tsv
