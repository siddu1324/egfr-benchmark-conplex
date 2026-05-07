# EGFR ConPLex Dataset Automation

This project builds EGFR / EGFR-mutant ligand datasets from public bioactivity sources and exports ConPLex-ready `pairs.tsv`.

## Outputs

- `data/processed/chembl_egfr_activities.tsv`
- `data/processed/bindingdb_egfr_activities.tsv`
- `data/processed/egfr_affinities_long.tsv`
- `data/processed/egfr_conplex_pairs.tsv`

## Dataset idea

Each bioactivity row keeps:

- molecule_id
- canonical_smiles
- target_id
- mutation_label
- protein_sequence
- affinity_type
- affinity_value_nm
- pAffinity
- source_db
- assay_description
- document/source reference

ConPLex prediction file format:

```tsv
protein_id    molecule_id    protein_sequence    molecule_smiles
```

No header.

## Setup

```bash
conda create -n egfr-conplex python=3.10 -y
conda activate egfr-conplex

pip install -r requirements.txt
```

Optional ConPLex install:

```bash
conda create -n conplex-dti python=3.9 -y
conda activate conplex-dti
pip install conplex-dti
conplex-dti download --to . --models ConPLex_v1_BindingDB
```

## Run

```bash
python src/build_chembl_egfr_dataset.py \
  --target CHEMBL203 \
  --out data/processed/chembl_egfr_activities.tsv

python src/build_bindingdb_egfr_dataset.py \
  --uniprot P00533 \
  --cutoff 100000000 \
  --out data/processed/bindingdb_egfr_activities.tsv

python src/make_conplex_files.py \
  --chembl data/processed/chembl_egfr_activities.tsv \
  --bindingdb data/processed/bindingdb_egfr_activities.tsv \
  --mutants configs/egfr_mutants.json \
  --out-affinities data/processed/egfr_affinities_long.tsv \
  --out-pairs data/processed/egfr_conplex_pairs.tsv
```

Run ConPLex:

```bash
conda activate conplex-dti

conplex-dti predict \
  --data-file data/processed/egfr_conplex_pairs.tsv \
  --model-path ./models/ConPLex_v1_BindingDB.pt \
  --outfile data/processed/egfr_conplex_predictions.tsv
```

## Important validation notes

1. Do not treat IC50, Ki, and Kd as the same experimentally.
2. Keep `affinity_type` as a separate column.
3. Use pAffinity only as a standardized numerical scale.
4. Verify mutant sequence numbering before generating final sequences.
5. LLMs should extract/label data from assay descriptions and papers, not invent values.
