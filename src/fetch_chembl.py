"""Harvest raw bioactivity records from ChEMBL for the pain-target panel.

Design notes
------------
* No `assay_type` filter is applied server-side. TRPV1 is an ion channel and a
  large share of its potency data is *functional* (assay_type 'F': calcium
  imaging, patch clamp) rather than binding ('B'). Filtering to 'B' alone, as
  is conventional for GPCR/enzyme panels, would silently discard most TRPV1
  pharmacology. We pull everything and filter in `build_dataset.py`, where the
  decision is visible and auditable.
* `action_type` is retained because it distinguishes AGONIST from ANTAGONIST.
  A naive "active/inactive" label collapses capsaicin (agonist) and capsazepine
  (antagonist) into the same class; keeping the field lets us model mechanism.
* `data_validity_comment` and `potential_duplicate` are retained as quality
  flags rather than being dropped at fetch time.

Output: data/raw/chembl_activities_raw.csv (one row per activity record)
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from chembl_client import fetch_all

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
CKPT = RAW / "chembl_by_target"
RAW.mkdir(parents=True, exist_ok=True)
CKPT.mkdir(parents=True, exist_ok=True)

KEEP_FIELDS = [
    "activity_id",
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "molecule_pref_name",
    "canonical_smiles",
    "target_chembl_id",
    "target_organism",
    "target_pref_name",
    "assay_chembl_id",
    "assay_type",
    "assay_description",
    "bao_label",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "action_type",
    "activity_comment",
    "data_validity_comment",
    "potential_duplicate",
    "document_chembl_id",
    "document_journal",
    "document_year",
    "src_id",
]


def load_panel() -> list[tuple[str, str, str]]:
    """Return (target_name, organism, chembl_id) for every resolved target."""
    targets = json.loads((RAW / "targets.json").read_text(encoding="utf-8"))
    panel = []
    for name, orgs in targets.items():
        for org, hits in orgs.items():
            for hit in hits:
                panel.append((name, org, hit["target_chembl_id"]))
    return panel


def harvest_target(chembl_id: str) -> list[dict]:
    recs = fetch_all(
        "activity.json",
        {"target_chembl_id": chembl_id},
        collection="activities",
        page_size=1000,
        workers=4,
    )
    return [{k: r.get(k) for k in KEEP_FIELDS} for r in recs]


def main():
    panel = load_panel()
    frames = []
    for name, org, cid in panel:
        ckpt = CKPT / f"{name}_{org}_{cid}.csv"
        if ckpt.exists():
            df = pd.read_csv(ckpt, low_memory=False)
            print(f"\n== {name} / {org} / {cid}  [cached {len(df)} records]")
            frames.append(df)
            continue

        print(f"\n== {name} / {org} / {cid}", flush=True)
        rows = harvest_target(cid)
        print(f"   -> {len(rows)} activity records", flush=True)
        if not rows:
            continue
        df = pd.DataFrame(rows)
        df["target_name"] = name
        df["target_organism_group"] = org
        # checkpoint immediately: this harvest takes ~an hour against a flaky
        # endpoint and must survive an interruption
        df.to_csv(ckpt, index=False)
        frames.append(df)

    full = pd.concat(frames, ignore_index=True)
    out = RAW / "chembl_activities_raw.csv"
    full.to_csv(out, index=False)
    print(f"\nwrote {out}  shape={full.shape}")
    print("\nrecords per target/organism:")
    print(full.groupby(["target_name", "target_organism_group"]).size())


if __name__ == "__main__":
    main()
