#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

aff = Path("data/processed/egfr_affinities_long.tsv")
rev = Path("data/processed/egfr_unknown_or_skipped_rows.tsv")

if aff.exists():
    df = pd.read_csv(aff, sep="\t")
    print("\n=== AFFINITY MUTATION COUNTS ===")
    if len(df):
        print(df["mutation_label"].value_counts(dropna=False).to_string())
        print("\n=== SAMPLE MUTANT ROWS ===")
        cols = ["source_db", "molecule_id", "mutation_label", "affinity_type", "affinity_value_nm", "assay_description"]
        print(df[cols].head(20).to_string(index=False))
    else:
        print("Affinity file exists but has 0 rows.")
else:
    print("No egfr_affinities_long.tsv found.")

if rev.exists():
    df = pd.read_csv(rev, sep="\t")
    print("\n=== UNKNOWN/SKIPPED SAMPLE ASSAYS ===")
    if len(df):
        cols = ["source_db", "target_pref_name", "affinity_type", "assay_description"]
        print(df[cols].head(40).to_string(index=False))
    else:
        print("Review file exists but has 0 rows.")
else:
    print("No review file found.")
