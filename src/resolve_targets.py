"""Resolve ChEMBL target IDs from UniProt accessions for the pain-target panel.

Targets (the paracetamol / AM404 central analgesia axis):
  TRPV1 - transient receptor potential vanilloid 1
  CB1   - cannabinoid receptor 1
  FAAH  - fatty acid amide hydrolase

Human is primary; rat/mouse orthologs are collected separately because a large
fraction of published TRPV1 and FAAH pharmacology is rodent. They are kept with
an explicit organism tag rather than silently pooled.
"""

import json
import time
from pathlib import Path

from chembl_client import get_json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = {
    "TRPV1": {
        "human": "Q8NER1",
        "rat": "O35433",
        "mouse": "Q704Y3",
    },
    "CB1": {
        "human": "P21554",
        "rat": "P20272",
        "mouse": "P47746",
    },
    "FAAH": {
        "human": "O00519",
        "rat": "P97612",
        "mouse": "O08914",
    },
}


def resolve(accession: str):
    """Return SINGLE PROTEIN ChEMBL targets matching a UniProt accession."""
    payload = get_json(
        "target.json",
        {"target_components__accession": accession, "limit": 50},
    )
    hits = payload.get("targets", [])
    return [
        {
            "target_chembl_id": h["target_chembl_id"],
            "pref_name": h.get("pref_name"),
            "organism": h.get("organism"),
            "target_type": h.get("target_type"),
            "accession": accession,
        }
        for h in hits
        if h.get("target_type") == "SINGLE PROTEIN"
    ]


def main():
    resolved = {}
    for target, orgs in PANEL.items():
        resolved[target] = {}
        for org, acc in orgs.items():
            hits = resolve(acc)
            resolved[target][org] = hits
            label = hits[0]["target_chembl_id"] if hits else "NOT FOUND"
            name = hits[0]["pref_name"] if hits else ""
            print(f"{target:6s} {org:6s} {acc:8s} -> {label:14s} {name}")
            time.sleep(0.3)

    path = OUT / "targets.json"
    path.write_text(json.dumps(resolved, indent=2), encoding="utf-8")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
