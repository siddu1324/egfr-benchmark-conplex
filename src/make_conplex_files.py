#!/usr/bin/env python3
"""
Merge EGFR ChEMBL + BindingDB data, infer mutant labels conservatively,
generate EGFR mutant sequences, and export ConPLex-ready pairs.tsv.

V3 fixes:
1. Does NOT treat plain "EGFR" as WT.
2. Does NOT parse random mutation-like strings by default.
3. Uses UniProt P00533 canonical full EGFR sequence by default.
4. Supports only configured EGFR mutants unless --allow-any-mutations is passed.
5. Writes unknown/unmapped rows to a review file.
"""

import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd
import requests

SIMPLE_MUT_RE = re.compile(r"\b([A-Z])(\d{2,4})([A-Z])\b")


def fetch_uniprot_fasta(uniprot_id: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    r = requests.get(url, timeout=90)
    r.raise_for_status()
    lines = [x.strip() for x in r.text.splitlines() if x.strip() and not x.startswith(">")]
    seq = "".join(lines)
    if not seq:
        raise RuntimeError(f"Could not fetch UniProt sequence for {uniprot_id}")
    return seq


def apply_mutation_sequence(reference_sequence: str, sequence_start: int, label: str) -> str:
    if label == "WT":
        return reference_sequence

    seq = list(reference_sequence)
    for mut in label.split("_"):
        m = SIMPLE_MUT_RE.fullmatch(mut.strip())
        if not m:
            raise ValueError(f"Unsupported mutation label: {label}. Only point mutations are supported here.")
        ref_aa, pos_str, alt_aa = m.groups()
        pos = int(pos_str)
        idx = pos - sequence_start
        if idx < 0 or idx >= len(seq):
            raise ValueError(f"{mut} maps outside sequence. Check sequence_start={sequence_start}.")
        observed = seq[idx]
        if observed != ref_aa:
            raise ValueError(
                f"Reference mismatch for {mut}: expected {ref_aa} at UniProt {pos}, "
                f"but sequence has {observed}. Check whether you are using full EGFR P00533 numbering."
            )
        seq[idx] = alt_aa
    return "".join(seq)


def compile_alias_patterns(aliases: Dict[str, List[str]]) -> List[tuple]:
    compiled = []
    labels = sorted(aliases.keys(), key=lambda x: (x.count("_"), len(x)), reverse=True)
    for label in labels:
        for alias in aliases[label]:
            alias = alias.strip()
            if not alias:
                continue
            parts = re.split(r"[/_\-\s]+", alias)
            if len(parts) > 1:
                pattern = r"\b" + r"[\s/_\-]+".join(map(re.escape, parts)) + r"\b"
            else:
                pattern = r"\b" + re.escape(alias) + r"\b"
            compiled.append((label, re.compile(pattern, flags=re.IGNORECASE)))
    return compiled


def infer_mutation_label(text: str, alias_patterns: List[tuple], allow_any_mutations: bool = False) -> str:
    if not isinstance(text, str):
        return "UNKNOWN"
    for label, pattern in alias_patterns:
        if pattern.search(text):
            return label
    if allow_any_mutations:
        found = SIMPLE_MUT_RE.findall(text.upper())
        if found:
            return "_".join([f"{a}{p}{b}" for a, p, b in found])
    return "UNKNOWN"


def molecule_safe_id(x):
    s = str(x)
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s[:80]


def load_optional(path):
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists():
        print(f"WARNING: input file not found: {p}")
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def p_affinity(value_nm):
    try:
        value_nm = float(value_nm)
        if value_nm <= 0:
            return None
        return -math.log10(value_nm * 1e-9)
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chembl")
    parser.add_argument("--bindingdb")
    parser.add_argument("--mutants", required=True)
    parser.add_argument("--out-affinities", required=True)
    parser.add_argument("--out-pairs", required=True)
    parser.add_argument("--review-out", default="data/processed/egfr_unknown_or_skipped_rows.tsv")
    parser.add_argument("--include-unknown", action="store_true")
    parser.add_argument("--allow-any-mutations", action="store_true")
    args = parser.parse_args()

    with open(args.mutants) as f:
        cfg = json.load(f)

    use_uniprot = cfg.get("use_uniprot_sequence", True)
    uniprot_id = cfg.get("reference", {}).get("uniprot", "P00533")
    sequence_start = int(cfg.get("sequence_start", 1))

    if use_uniprot:
        print(f"Fetching canonical EGFR sequence from UniProt {uniprot_id}...")
        reference_sequence = fetch_uniprot_fasta(uniprot_id)
        sequence_start = 1
    else:
        reference_sequence = cfg["sequence"]

    print(f"Reference sequence length: {len(reference_sequence)} aa")
    print(f"Sequence start numbering: {sequence_start}")

    alias_patterns = compile_alias_patterns(cfg["mutation_aliases"])

    dfs = []
    for one in [load_optional(args.chembl), load_optional(args.bindingdb)]:
        if not one.empty:
            dfs.append(one)
    if not dfs:
        raise SystemExit("No input data loaded.")

    df = pd.concat(dfs, ignore_index=True)
    if "assay_description" not in df.columns:
        df["assay_description"] = ""
    df["assay_description"] = df["assay_description"].fillna("").astype(str)

    df["mutation_label"] = df["assay_description"].apply(
        lambda x: infer_mutation_label(x, alias_patterns, args.allow_any_mutations)
    )

    review_cols = [
        "source_db", "molecule_id", "target_source_id", "target_pref_name",
        "affinity_type", "affinity_value_nm", "assay_description", "mutation_label"
    ]
    for c in review_cols:
        if c not in df.columns:
            df[c] = None

    skipped = df[df["mutation_label"] == "UNKNOWN"].copy()
    Path(args.review_out).parent.mkdir(parents=True, exist_ok=True)
    skipped[review_cols].drop_duplicates().to_csv(args.review_out, sep="\t", index=False)

    if not args.include_unknown:
        df = df[df["mutation_label"] != "UNKNOWN"].copy()

    sequences = {}
    bad_labels = {}
    for label in sorted(df["mutation_label"].dropna().unique()):
        if label == "UNKNOWN":
            continue
        try:
            sequences[label] = apply_mutation_sequence(reference_sequence, sequence_start, label)
        except Exception as e:
            bad_labels[label] = str(e)
            sequences[label] = None

    if bad_labels:
        print("\nWARNING: some labels could not be mapped:")
        for k, v in bad_labels.items():
            print(f"  {k}: {v}")

    df["target_id"] = "EGFR_" + df["mutation_label"].astype(str)
    df["protein_sequence"] = df["mutation_label"].map(sequences)

    before = len(df)
    df = df.dropna(subset=["canonical_smiles", "protein_sequence"])
    print(f"Dropped {before - len(df):,} rows without SMILES or valid protein sequence.")

    if "pAffinity" not in df.columns:
        df["pAffinity"] = df["affinity_value_nm"].apply(p_affinity)

    affinity_cols = [
        "molecule_id", "canonical_smiles", "target_id", "mutation_label", "protein_sequence",
        "affinity_type", "affinity_relation", "affinity_value_nm", "unit", "pAffinity",
        "source_db", "target_source_id", "target_pref_name", "assay_description",
        "assay_chembl_id", "document_chembl_id", "document_id", "activity_id", "source_url"
    ]
    for c in affinity_cols:
        if c not in df.columns:
            df[c] = None

    affinities = df[affinity_cols].drop_duplicates()
    Path(args.out_affinities).parent.mkdir(parents=True, exist_ok=True)
    affinities.to_csv(args.out_affinities, sep="\t", index=False)

    pairs = affinities[["target_id", "molecule_id", "protein_sequence", "canonical_smiles"]].copy()
    pairs["molecule_id"] = pairs["molecule_id"].apply(molecule_safe_id)
    pairs = pairs.drop_duplicates()
    pairs.to_csv(args.out_pairs, sep="\t", index=False, header=False)

    print(f"\nSaved {len(affinities):,} affinity rows to {args.out_affinities}")
    print(f"Saved {len(pairs):,} ConPLex pairs to {args.out_pairs}")
    print(f"Saved {len(skipped):,} UNKNOWN/skipped rows to {args.review_out}")
    print("\nMutation-label counts:")
    if len(affinities):
        print(affinities["mutation_label"].value_counts(dropna=False).to_string())
    else:
        print("No mutation-labeled rows found. Inspect review file and add aliases to config.")


if __name__ == "__main__":
    main()
