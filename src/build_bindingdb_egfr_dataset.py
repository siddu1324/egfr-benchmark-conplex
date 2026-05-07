#!/usr/bin/env python3
"""
Build EGFR bioactivity data from BindingDB.

Important:
- A huge cutoff like 100000000 can trigger 504 Gateway Time-out.
- Start with smaller cutoffs: 10000, 100000, or 1000000.
- This script uses retries and supports both BindingDB's singular and plural UniProt endpoints.

Recommended run:
python src/build_bindingdb_egfr_dataset.py \
  --uniprot P00533 \
  --cutoff 10000 \
  --out data/processed/bindingdb_egfr_activities.tsv
"""

import argparse
import json
import math
import re
import time
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
import requests


def p_affinity(value_nm):
    try:
        value_nm = float(value_nm)
        if value_nm <= 0:
            return None
        return -math.log10(value_nm * 1e-9)
    except Exception:
        return None


def parse_number(x):
    if x is None:
        return None
    s = str(x).strip()
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    return float(m.group(0))


def get_any(d: Dict, keys: List[str]):
    for k in keys:
        if k in d and d[k] not in [None, ""]:
            return d[k]
    lower = {str(k).lower(): v for k, v in d.items()}
    for k in keys:
        v = lower.get(k.lower())
        if v not in [None, ""]:
            return v
    return None


def flatten_json(obj: Any) -> List[Dict]:
    rows = []

    def walk(x):
        if isinstance(x, dict):
            keys = {str(k).lower() for k in x.keys()}
            if (
                any("smiles" in k for k in keys)
                or any(k in keys for k in ["ic50", "ki", "kd", "monomerid"])
            ):
                rows.append(x)
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)

    walk(obj)
    return rows


def normalize_bindingdb_json(data: Any, uniprot: str) -> List[Dict]:
    candidates = flatten_json(data)
    rows = []

    for item in candidates:
        if not isinstance(item, dict):
            continue

        smiles = get_any(item, ["smiles", "SMILES", "canonical_smiles", "Ligand SMILES", "Ligand_SMILES"])
        mol_id = get_any(item, ["monomerid", "monomerId", "MonomerID", "Monomer ID", "BDBM", "bdbm", "ligand"])
        target_name = get_any(item, ["targetName", "Target Name", "target_name", "Target", "Protein", "Name"])
        assay = get_any(item, ["assay", "Assay Description", "assay_description", "Assay", "Target Chain Sequence"])
        doc = get_any(item, ["pmid", "PubMed", "DOI", "doi", "Article DOI", "Institution"])

        for typ in ["IC50", "Ki", "Kd"]:
            val = get_any(item, [typ, typ.lower(), f"{typ} (nM)", f"{typ}_nM"])
            num = parse_number(val)
            if num is not None:
                rows.append({
                    "source_db": "BindingDB",
                    "molecule_id": f"BDBM{mol_id}" if mol_id and not str(mol_id).startswith("BDB") else mol_id,
                    "canonical_smiles": smiles,
                    "target_source_id": uniprot,
                    "target_pref_name": target_name,
                    "affinity_type": typ,
                    "affinity_relation": None,
                    "affinity_value_nm": num,
                    "unit": "nM",
                    "assay_description": assay,
                    "document_id": doc,
                    "raw_keys": ",".join(map(str, item.keys())),
                })

    return rows


def request_with_retries(url: str, params: dict, retries: int = 3, timeout: int = 180) -> requests.Response:
    last = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 504:
                print(f"Attempt {attempt}: BindingDB timed out with 504. Retrying after wait...")
                time.sleep(5 * attempt)
                last = r
                continue
            r.raise_for_status()
            return r
        except Exception as e:
            print(f"Attempt {attempt} failed: {e}")
            last = e
            time.sleep(5 * attempt)
    if isinstance(last, requests.Response):
        last.raise_for_status()
    raise RuntimeError(f"BindingDB request failed after {retries} attempts: {last}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--uniprot", default="P00533")
    parser.add_argument("--cutoff", type=int, default=10000, help="nM cutoff. Avoid massive values first.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--raw-json-out", default=None, help="Optional debug path to save raw BindingDB JSON")
    parser.add_argument("--endpoint", choices=["plural", "singular"], default="plural")
    args = parser.parse_args()

    if args.cutoff > 1000000:
        print("WARNING: cutoff > 1,000,000 nM can be very large and may cause BindingDB 504 timeouts.")

    if args.endpoint == "plural":
        url = "https://bindingdb.org/rest/getLigandsByUniprots"
        params = {"uniprot": args.uniprot, "cutoff": str(args.cutoff), "response": "application/json"}
    else:
        url = "https://bindingdb.org/rest/getLigandsByUniprot"
        params = {"uniprot": f"{args.uniprot};{args.cutoff}", "response": "application/json"}

    r = request_with_retries(url, params=params)

    try:
        data = r.json()
    except Exception as e:
        raise SystemExit(f"BindingDB did not return JSON. Response starts:\n{r.text[:2000]}") from e

    if args.raw_json_out:
        Path(args.raw_json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.raw_json_out).write_text(json.dumps(data, indent=2)[:10_000_000])

    rows = normalize_bindingdb_json(data, args.uniprot)
    df = pd.DataFrame(rows)

    if df.empty:
        raise SystemExit(
            "No rows parsed from BindingDB response. Save raw JSON and inspect structure:\n"
            f"python src/build_bindingdb_egfr_dataset.py --uniprot {args.uniprot} --cutoff {args.cutoff} "
            "--raw-json-out data/raw/bindingdb_debug.json --out data/processed/tmp.tsv"
        )

    df["affinity_value_nm"] = pd.to_numeric(df["affinity_value_nm"], errors="coerce")
    df["pAffinity"] = df["affinity_value_nm"].apply(p_affinity)
    df = df.dropna(subset=["canonical_smiles", "affinity_value_nm"])

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"Saved {len(df):,} BindingDB rows to {args.out}")
    print("Endpoint counts:")
    print(df["affinity_type"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
