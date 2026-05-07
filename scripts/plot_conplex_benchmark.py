#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.merged, sep="\t")
    plt.figure(figsize=(8,6))
    plt.scatter(df["conplex_score"], df["pAffinity_exp"], alpha=0.35)
    plt.xlabel("ConPLex prediction score")
    plt.ylabel("Experimental pAffinity, median per pair")
    plt.title("ConPLex score vs experimental EGFR affinity")
    scatter_path = outdir / "plot_conplex_vs_pAffinity.png"
    plt.tight_layout(); plt.savefig(scatter_path, dpi=200); plt.close()
    counts = df["mutation_label"].value_counts()
    keep = counts[counts >= 20].index
    sub = df[df["mutation_label"].isin(keep)].copy()
    order = sub.groupby("mutation_label")["conplex_score"].median().sort_values(ascending=False).index.tolist()
    plt.figure(figsize=(10,6))
    data = [sub[sub["mutation_label"] == m]["conplex_score"] for m in order]
    plt.boxplot(data, labels=order, showfliers=False)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("ConPLex prediction score")
    plt.title("ConPLex scores by EGFR mutant")
    box_path = outdir / "plot_scores_by_mutation.png"
    plt.tight_layout(); plt.savefig(box_path, dpi=200); plt.close()
    print(f"Wrote {scatter_path}")
    print(f"Wrote {box_path}")

if __name__ == "__main__":
    main()
