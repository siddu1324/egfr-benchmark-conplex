#!/usr/bin/env python3
"""
Quick API sanity checks.

Run:
python scripts/debug_api_shapes.py
"""

import json
import requests

print("ChEMBL activity first page...")
r = requests.get(
    "https://www.ebi.ac.uk/chembl/api/data/activity.json",
    params={"target_chembl_id": "CHEMBL203", "limit": 2, "format": "json"},
    timeout=60,
)
r.raise_for_status()
data = r.json()
print("Top-level keys:", data.keys())
print("Record key exists?", "activities" in data)
print("First activity keys:", list(data["activities"][0].keys())[:30])
print(json.dumps(data["activities"][0], indent=2)[:1500])

print("\nBindingDB small request...")
r = requests.get(
    "https://bindingdb.org/rest/getLigandsByUniprots",
    params={"uniprot": "P00533", "cutoff": "100", "response": "application/json"},
    timeout=90,
)
print("Status:", r.status_code)
print("Text starts:", r.text[:1000])
