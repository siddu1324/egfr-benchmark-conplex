# LLM Paper/Table Extraction Prompt for EGFR Binding Data

You are extracting EGFR ligand-binding or inhibition data from a scientific paper, patent, or assay table.

Return ONLY valid JSON with this schema:

[
  {
    "compound_name": "",
    "compound_id_in_paper": "",
    "smiles": null,
    "target": "EGFR",
    "mutation_label": "WT | L858R | T790M | C797S | L858R_T790M | L858R_T790M_C797S | G719C | G719S | L861Q | S768I | UNKNOWN",
    "affinity_type": "IC50 | Ki | Kd | EC50 | UNKNOWN",
    "affinity_value": null,
    "affinity_unit": "nM | uM | M | UNKNOWN",
    "assay_description": "",
    "cell_line_or_system": "",
    "source_table_or_figure": "",
    "source_quote": "",
    "confidence": "high | medium | low",
    "notes": ""
  }
]

Rules:
- Do not invent values.
- If a SMILES string is not present, return null.
- Preserve exact values and units.
- If the value is written as >1000 nM or <1 nM, preserve the relation in notes.
- Mutation label must come from the assay description, table heading, target name, or surrounding text.
- If target says only EGFR without wild-type or mutant wording, set mutation_label to UNKNOWN, not WT.
