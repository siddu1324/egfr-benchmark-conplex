#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/processed

python src/make_conplex_files.py \
  --chembl data/processed/chembl_egfr_activities.tsv \
  --bindingdb data/processed/bindingdb_egfr_activities.tsv \
  --mutants configs/egfr_mutants_v3.json \
  --out-affinities data/processed/egfr_affinities_long.tsv \
  --out-pairs data/processed/egfr_conplex_pairs.tsv \
  --review-out data/processed/egfr_unknown_or_skipped_rows.tsv

python scripts/inspect_mutant_dataset.py
