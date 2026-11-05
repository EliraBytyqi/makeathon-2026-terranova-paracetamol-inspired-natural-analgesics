"""Build the plant-derived compound library from the LOTUS frozen metadata dump.

LOTUS gives referenced structure-organism pairs: every (compound, organism) edge
carries a literature DOI. That provenance is what makes the knowledge graph
defensible rather than decorative.

Two tiers are produced:

  tier = "pain_seed"   compounds from the hand-curated ethnopharmacology seed
                       list (data/reference/pain_plants_seed.csv), matched at
                       species level where possible, else genus level.
  tier = "plant_background"  compounds from any Archaeplastida (land plant)
                       source. Needed so the screening model is not evaluated
                       only on compounds we already believe in - without a
                       background set, "top ranked" is meaningless.

Fungal/bacterial/marine sources are excluded: the project premise is
plant-derived analgesics.

Output: data/interim/np_library_raw.csv  (one row per compound-organism pair)
"""

from __future__ import annotations

import gzip
import csv
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LOTUS = ROOT / "data" / "external" / "lotus" / "260413_frozen_metadata.csv.gz"
SEED = ROOT / "data" / "reference" / "pain_plants_seed.csv"
OUT = ROOT / "data" / "interim"
OUT.mkdir(parents=True, exist_ok=True)

# Columns we carry forward out of LOTUS's 39.
COLS = [
    "structure_inchikey",
    "structure_smiles",
    "structure_smiles_2D",
    "structure_molecular_formula",
    "structure_exact_mass",
    "structure_xlogp",
    "structure_cid",
    "structure_nameTraditional",
    "structure_taxonomy_npclassifier_01pathway",
    "structure_taxonomy_npclassifier_02superclass",
    "structure_taxonomy_npclassifier_03class",
    "structure_taxonomy_classyfire_03class",
    "organism_name",
    "organism_taxonomy_02kingdom",
    "organism_taxonomy_06family",
    "organism_taxonomy_08genus",
    "organism_taxonomy_09species",
    "reference_doi",
    "structure_wikidata",
    "organism_wikidata",
]

PLANT_KINGDOM = "Archaeplastida"


def load_seed():
    seed = pd.read_csv(SEED)
    species = {s.strip().lower() for s in seed["species"].dropna()}
    genera = {g.strip().lower() for g in seed["genus"].dropna()}
    meta = {}
    for _, r in seed.iterrows():
        meta[r["species"].strip().lower()] = {
            "traditional_pain_use": r["traditional_pain_use"],
            "relevant_target_hypothesis": r["relevant_target_hypothesis"],
            "common_name": r["common_name"],
        }
    return species, genera, meta


def main():
    species_set, genus_set, seed_meta = load_seed()
    print(f"seed list: {len(species_set)} species across {len(genus_set)} genera")

    rows = []
    n_total = 0
    n_plant = 0
    kingdom_counts = defaultdict(int)

    csv.field_size_limit(10_000_000)
    with gzip.open(LOTUS, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh)
        for rec in reader:
            n_total += 1
            if n_total % 250_000 == 0:
                print(f"  scanned {n_total:,} pairs, kept {len(rows):,}", flush=True)

            kingdom = (rec.get("organism_taxonomy_02kingdom") or "").strip()
            kingdom_counts[kingdom] += 1
            if kingdom != PLANT_KINGDOM:
                continue
            n_plant += 1

            smiles = (rec.get("structure_smiles") or "").strip()
            ikey = (rec.get("structure_inchikey") or "").strip()
            if not smiles or not ikey:
                continue

            sp = (rec.get("organism_taxonomy_09species") or "").strip().lower()
            gn = (rec.get("organism_taxonomy_08genus") or "").strip().lower()

            if sp in species_set:
                tier, match = "pain_seed", "species"
            elif gn in genus_set:
                tier, match = "pain_seed", "genus"
            else:
                tier, match = "plant_background", "none"

            out = {c: rec.get(c) for c in COLS}
            out["tier"] = tier
            out["seed_match_level"] = match
            info = seed_meta.get(sp, {})
            out["traditional_pain_use"] = info.get("traditional_pain_use", "")
            out["relevant_target_hypothesis"] = info.get("relevant_target_hypothesis", "")
            rows.append(out)

    df = pd.DataFrame(rows)
    out_path = OUT / "np_library_raw.csv"
    df.to_csv(out_path, index=False)

    print(f"\nscanned {n_total:,} LOTUS pairs")
    print(f"plant (Archaeplastida) pairs: {n_plant:,}")
    print(f"kept rows: {len(df):,}")
    print("\ntop source kingdoms in LOTUS:")
    for k, v in sorted(kingdom_counts.items(), key=lambda x: -x[1])[:8]:
        print(f"   {k or '(blank)':20s} {v:,}")
    print("\nby tier (pairs):")
    print(df["tier"].value_counts())
    print("\nunique compounds by tier:")
    print(df.groupby("tier")["structure_inchikey"].nunique())
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
