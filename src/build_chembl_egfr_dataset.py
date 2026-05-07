#!/usr/bin/env python3
"""
Build EGFR bioactivity data from ChEMBL.

Fixes:
1. ChEMBL activity JSON key is "activities", not "activitys".
2. Adds safer pagination.
3. Fetches molecule SMILES in batches using molecule_chembl_id__in.
4. Filters IC50/Ki/Kd + nM again in Python for safety.

Run:
python src/build_chembl_egfr_dataset.py \
  --target CHEMBL203 \
  --out data/processed/chembl_egfr_activities.tsv
"""

import argparse
import math
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
import requests
from tqdm import tqdm


BASE = "https://www.ebi.ac.uk/chembl/api/data"
RESOURCE_KEYS = {
    "activity": "activities",
    "molecule": "molecules",
    "assay": "assays",
    "target": "targets",
}


def get_json(url: str, params: Optional[dict] = None, sleep: float = 0.05) -> dict:
    r = requests.get(url, params=params, timeout=90)
    r.raise_for_status()
    time.sleep(sleep)
    return r.json()


def paginate(endpoint: str, params: dict, limit: int = 1000) -> Iterable[dict]:
    offset = 0
    total = None
    key = RESOURCE_KEYS.get(endpoint, endpoint + "s")

    while True:
        p = dict(params)
        p.update({"limit": limit, "offset": offset, "format": "json"})
        data = get_json(f"{BASE}/{endpoint}.json", params=p)

        records = data.get(key, [])
        if total is None:
            total = data.get("page_meta", {}).get("total_count", 0)
            print(f"{endpoint}: {total:,} records reported by ChEMBL")

        if not records:
            break

        for row in records:
            yield row

        offset += limit
        if total is not None and offset >= total:
            break


def chunks(items: List[str], size: int = 100):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def fetch_molecule_smiles_batched(molecule_ids: List[str]) -> Dict[str, str]:
    smiles = {}
    unique_ids = sorted(set([m for m in molecule_ids if isinstance(m, str) and m]))

    for batch in tqdm(list(chunks(unique_ids, 100)), desc="Fetching ChEMBL SMILES batches"):
        params = {
            "molecule_chembl_id__in": ",".join(batch),
            "limit": 1000,
            "format": "json",
        }
        try:
            data = get_json(f"{BASE}/molecule.json", params=params)
            for mol in data.get("molecules", []):
                mid = mol.get("molecule_chembl_id")
                structs = mol.get("molecule_structures") or {}
                smi = structs.get("canonical_smiles")
                if mid and smi:
                    smiles[mid] = smi
        except Exception as e:
            print(f"WARNING: SMILES batch failed for {batch[:3]}... {e}")

    return smiles


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
    parser.add_argument("--target", default="CHEMBL203", help="ChEMBL target ID for EGFR")
    parser.add_argument("--out", required=True)
    parser.add_argument("--limit", type=int, default=1000)
    parser.add_argument("--no-smiles", action="store_true", help="Skip molecule SMILES fetching for quick debugging")
    args = parser.parse_args()

    params = {
        "target_chembl_id": args.target,
        "standard_units__iexact": "nM",
        "standard_type__in": "IC50,Ki,Kd",
    }

    rows = []
    allowed_types = {"IC50", "KI", "KD"}

    for a in paginate("activity", params=params, limit=args.limit):
        standard_type = str(a.get("standard_type") or "").upper()
        standard_units = str(a.get("standard_units") or "").lower()
        standard_value = a.get("standard_value")
        molecule_chembl_id = a.get("molecule_chembl_id")

        if standard_type not in allowed_types:
            continue
        if standard_units != "nm":
            continue
        if standard_value in [None, "", "None"]:
            continue
        if not molecule_chembl_id:
            continue

        rows.append({
            "source_db": "ChEMBL",
            "molecule_id": molecule_chembl_id,
            "canonical_smiles": None,
            "target_source_id": a.get("target_chembl_id"),
            "target_pref_name": a.get("target_pref_name"),
            "affinity_type": a.get("standard_type"),
            "affinity_relation": a.get("standard_relation"),
            "affinity_value_nm": standard_value,
            "unit": a.get("standard_units"),
            "assay_chembl_id": a.get("assay_chembl_id"),
            "assay_description": a.get("assay_description"),
            "document_chembl_id": a.get("document_chembl_id"),
            "activity_id": a.get("activity_id"),
            "bao_endpoint": a.get("bao_endpoint"),
            "pchembl_value": a.get("pchembl_value"),
            "data_validity_comment": a.get("data_validity_comment"),
        })

    df = pd.DataFrame(rows)

    if df.empty:
        raise SystemExit(
            "No usable ChEMBL rows after filtering. Try checking the first API page manually:\n"
            "https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=CHEMBL203&limit=5"
        )

    if not args.no_smiles:
        smiles_map = fetch_molecule_smiles_batched(df["molecule_id"].dropna().unique().tolist())
        df["canonical_smiles"] = df["molecule_id"].map(smiles_map)
        before = len(df)
        df = df.dropna(subset=["canonical_smiles"])
        print(f"Dropped {before - len(df):,} rows without SMILES.")
    else:
        print("WARNING: --no-smiles used; canonical_smiles will be empty.")

    df["affinity_value_nm"] = pd.to_numeric(df["affinity_value_nm"], errors="coerce")
    df["pAffinity"] = df["affinity_value_nm"].apply(p_affinity)
    df = df.dropna(subset=["affinity_value_nm"])

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"Saved {len(df):,} ChEMBL rows to {args.out}")
    print("Endpoint counts:")
    print(df["affinity_type"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
