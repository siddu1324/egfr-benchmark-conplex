#!/usr/bin/env python3
"""
Benchmark ConPLex predictions against experimental EGFR affinity data.

Patch v2:
- Handles ConPLex outputs where columns are reversed:
  prediction file may be molecule_id, target_id, score instead of target_id, molecule_id, score.
- If normal merge gives 0 rows, it automatically swaps target_id/molecule_id and retries.

Run:
python scripts/benchmark_conplex_results.py \
  --affinities data/processed/egfr_affinities_long.tsv \
  --predictions data/processed/conplex_results/egfr_conplex_predictions_full.tsv \
  --outdir data/processed/conplex_results
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def find_prediction_score_col(df: pd.DataFrame) -> str:
    preferred = [
        "prediction", "predictions", "score", "scores", "pred_score",
        "binding_score", "probability", "y_pred", "pred", "conplex_score"
    ]
    lower_map = {c.lower(): c for c in df.columns}
    for name in preferred:
        if name in lower_map:
            return lower_map[name]

    numeric_cols = []
    for c in df.columns:
        if c.lower() in {"target_id", "protein_id", "molecule_id", "drug_id"}:
            continue
        converted = pd.to_numeric(df[c], errors="coerce")
        if converted.notna().sum() > 0:
            numeric_cols.append(c)

    if not numeric_cols:
        raise ValueError(
            "Could not auto-detect prediction score column. "
            f"Prediction columns are: {list(df.columns)}"
        )
    return numeric_cols[-1]


def read_predictions(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)

    # First try headered TSV.
    df = pd.read_csv(p, sep="\t")

    # If no useful header, reread headerless.
    expected_ids = {"target_id", "protein_id", "molecule_id", "drug_id"}
    has_known_cols = any(c.lower() in expected_ids for c in df.columns)

    if not has_known_cols:
        raw = pd.read_csv(p, sep="\t", header=None)

        if raw.shape[1] == 3:
            raw.columns = ["col0", "col1", "prediction"]
        elif raw.shape[1] == 5:
            raw.columns = ["col0", "col1", "protein_sequence", "canonical_smiles", "prediction"]
        elif raw.shape[1] >= 4:
            cols = [f"col_{i}" for i in range(raw.shape[1])]
            raw.columns = cols
            raw = raw.rename(columns={"col_0": "col0", "col_1": "col1", cols[-1]: "prediction"})
        else:
            raise ValueError(f"Unexpected prediction file shape: {raw.shape}")

        # Guess direction from values.
        # EGFR_* is target_id; CHEMBL* is molecule_id.
        col0 = raw["col0"].astype(str)
        col1 = raw["col1"].astype(str)

        if col0.str.startswith("EGFR_").mean() > col1.str.startswith("EGFR_").mean():
            raw = raw.rename(columns={"col0": "target_id", "col1": "molecule_id"})
        else:
            raw = raw.rename(columns={"col0": "molecule_id", "col1": "target_id"})

        df = raw

    # Normalize possible named columns.
    rename = {}
    for c in df.columns:
        cl = c.lower()
        if cl in {"protein_id", "target", "target_name"}:
            rename[c] = "target_id"
        if cl in {"drug_id", "compound_id", "ligand_id", "molecule"}:
            rename[c] = "molecule_id"
    df = df.rename(columns=rename)

    if "target_id" not in df.columns or "molecule_id" not in df.columns:
        # Fallback: first two columns. Direction guessed by EGFR_ prefix.
        c0, c1 = df.columns[0], df.columns[1]
        col0 = df[c0].astype(str)
        col1 = df[c1].astype(str)
        if col0.str.startswith("EGFR_").mean() > col1.str.startswith("EGFR_").mean():
            df = df.rename(columns={c0: "target_id", c1: "molecule_id"})
        else:
            df = df.rename(columns={c0: "molecule_id", c1: "target_id"})

    score_col = find_prediction_score_col(df)
    df["conplex_score"] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=["target_id", "molecule_id", "conplex_score"])

    df["target_id"] = df["target_id"].astype(str)
    df["molecule_id"] = df["molecule_id"].astype(str)

    return df[["target_id", "molecule_id", "conplex_score"]].drop_duplicates()


def safe_corr(x, y, method: str):
    if len(x) < 3:
        return np.nan
    return pd.Series(x).corr(pd.Series(y), method=method)


def rmse(y_true, y_pred):
    if len(y_true) == 0:
        return np.nan
    return float(np.sqrt(np.mean((np.asarray(y_true) - np.asarray(y_pred)) ** 2)))


def mae(y_true, y_pred):
    if len(y_true) == 0:
        return np.nan
    return float(np.mean(np.abs(np.asarray(y_true) - np.asarray(y_pred))))


def summarize_metrics(df: pd.DataFrame, group_name="ALL") -> dict:
    valid = df.dropna(subset=["conplex_score", "pAffinity_exp"]).copy()

    pearson = safe_corr(valid["conplex_score"], valid["pAffinity_exp"], "pearson")
    spearman = safe_corr(valid["conplex_score"], valid["pAffinity_exp"], "spearman")

    # Also compute flipped correlations in case lower ConPLex score means stronger binding.
    pearson_flipped = safe_corr(-valid["conplex_score"], valid["pAffinity_exp"], "pearson")
    spearman_flipped = safe_corr(-valid["conplex_score"], valid["pAffinity_exp"], "spearman")

    if len(valid) >= 3 and valid["conplex_score"].nunique() > 1:
        coef = np.polyfit(valid["conplex_score"], valid["pAffinity_exp"], deg=1)
        pred_scaled = np.polyval(coef, valid["conplex_score"])
        scaled_rmse = rmse(valid["pAffinity_exp"], pred_scaled)
        scaled_mae = mae(valid["pAffinity_exp"], pred_scaled)
    else:
        scaled_rmse = np.nan
        scaled_mae = np.nan

    if len(valid) >= 10:
        q90 = valid["conplex_score"].quantile(0.90)
        q10 = valid["conplex_score"].quantile(0.10)
        top_mean = valid.loc[valid["conplex_score"] >= q90, "pAffinity_exp"].mean()
        bottom_mean = valid.loc[valid["conplex_score"] <= q10, "pAffinity_exp"].mean()
        top_minus_bottom = top_mean - bottom_mean

        # Flipped top/bottom: if lower score is better.
        low_score_mean = valid.loc[valid["conplex_score"] <= q10, "pAffinity_exp"].mean()
        high_score_mean = valid.loc[valid["conplex_score"] >= q90, "pAffinity_exp"].mean()
        flipped_top_minus_bottom = low_score_mean - high_score_mean
    else:
        top_mean = bottom_mean = top_minus_bottom = flipped_top_minus_bottom = np.nan

    return {
        "group": group_name,
        "n_pairs": len(valid),
        "n_targets": valid["target_id"].nunique(),
        "n_molecules": valid["molecule_id"].nunique(),
        "pearson_conplex_vs_pAffinity": pearson,
        "spearman_conplex_vs_pAffinity": spearman,
        "pearson_flipped_score_vs_pAffinity": pearson_flipped,
        "spearman_flipped_score_vs_pAffinity": spearman_flipped,
        "linear_scaled_RMSE_pAffinity": scaled_rmse,
        "linear_scaled_MAE_pAffinity": scaled_mae,
        "top10pct_pred_mean_pAffinity": top_mean,
        "bottom10pct_pred_mean_pAffinity": bottom_mean,
        "top_minus_bottom_pAffinity": top_minus_bottom,
        "flipped_top_minus_bottom_pAffinity": flipped_top_minus_bottom,
        "mean_pAffinity": valid["pAffinity_exp"].mean(),
        "median_pAffinity": valid["pAffinity_exp"].median(),
    }


def build_grouped_affinities(aff: pd.DataFrame) -> pd.DataFrame:
    return (
        aff.dropna(subset=["pAffinity"])
        .groupby(["target_id", "molecule_id"], as_index=False)
        .agg(
            mutation_label=("mutation_label", lambda x: x.mode().iat[0] if len(x.mode()) else x.iloc[0]),
            pAffinity_exp=("pAffinity", "median"),
            affinity_value_nm_median=("affinity_value_nm", "median"),
            n_exp_measurements=("pAffinity", "size"),
            affinity_types=("affinity_type", lambda x: ",".join(sorted(set(map(str, x))))),
            assay_examples=("assay_description", lambda x: " | ".join(list(pd.Series(x).dropna().astype(str).drop_duplicates().head(3))))
        )
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--affinities", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--affinity-type", default=None, choices=["IC50", "Ki", "Kd"])
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    aff = pd.read_csv(args.affinities, sep="\t")
    pred = read_predictions(args.predictions)

    required = {"target_id", "molecule_id", "mutation_label", "pAffinity", "affinity_value_nm", "affinity_type"}
    missing = required - set(aff.columns)
    if missing:
        raise ValueError(f"Affinity file missing columns: {missing}")

    aff["target_id"] = aff["target_id"].astype(str)
    aff["molecule_id"] = aff["molecule_id"].astype(str)
    aff["pAffinity"] = pd.to_numeric(aff["pAffinity"], errors="coerce")
    aff["affinity_value_nm"] = pd.to_numeric(aff["affinity_value_nm"], errors="coerce")

    if args.affinity_type:
        aff = aff[aff["affinity_type"].astype(str).str.upper() == args.affinity_type.upper()].copy()

    grouped = build_grouped_affinities(aff)

    merged = pred.merge(grouped, on=["target_id", "molecule_id"], how="inner")
    merge_mode = "normal"

    # Automatic repair: if no rows merged, swap the prediction ID columns and retry.
    if merged.empty:
        pred_swapped = pred.rename(columns={"target_id": "molecule_id", "molecule_id": "target_id"})
        merged_swapped = pred_swapped.merge(grouped, on=["target_id", "molecule_id"], how="inner")
        if not merged_swapped.empty:
            pred = pred_swapped
            merged = merged_swapped
            merge_mode = "swapped_prediction_columns"

    merged = merged.dropna(subset=["conplex_score", "pAffinity_exp"])

    if merged.empty:
        print("No rows merged even after trying swapped prediction columns.")
        print("\nPrediction sample:")
        print(pred.head().to_string(index=False))
        print("\nAffinity sample:")
        print(grouped.head().to_string(index=False))
        raise SystemExit(1)

    merged_path = outdir / "benchmark_merged_predictions.tsv"
    merged.to_csv(merged_path, sep="\t", index=False)

    overall = pd.DataFrame([summarize_metrics(merged, "ALL")])
    overall.insert(1, "merge_mode", merge_mode)
    overall_path = outdir / "benchmark_overall_metrics.tsv"
    overall.to_csv(overall_path, sep="\t", index=False)

    by_mut_rows = []
    for mut, sub in merged.groupby("mutation_label"):
        by_mut_rows.append(summarize_metrics(sub, mut))
    by_mut = pd.DataFrame(by_mut_rows).sort_values("n_pairs", ascending=False)
    by_mut.insert(1, "merge_mode", merge_mode)
    by_mut_path = outdir / "benchmark_by_mutation.tsv"
    by_mut.to_csv(by_mut_path, sep="\t", index=False)

    top = merged.sort_values("conplex_score", ascending=False).head(100)
    top_path = outdir / "benchmark_top_predictions.tsv"
    top.to_csv(top_path, sep="\t", index=False)

    print("\n=== BENCHMARK COMPLETE ===")
    print(f"Merge mode: {merge_mode}")
    print(f"Merged rows: {len(merged):,}")
    print(f"Wrote: {merged_path}")
    print(f"Wrote: {overall_path}")
    print(f"Wrote: {by_mut_path}")
    print(f"Wrote: {top_path}")

    print("\n=== OVERALL METRICS ===")
    print(overall.to_string(index=False))

    print("\n=== BY MUTATION ===")
    show_cols = [
        "group", "n_pairs",
        "pearson_conplex_vs_pAffinity", "spearman_conplex_vs_pAffinity",
        "spearman_flipped_score_vs_pAffinity",
        "top_minus_bottom_pAffinity", "flipped_top_minus_bottom_pAffinity"
    ]
    print(by_mut[show_cols].to_string(index=False))


if __name__ == "__main__":
    main()
