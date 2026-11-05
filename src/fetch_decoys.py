"""Fetch a presumed-inactive decoy pool from ChEMBL.

Why this is necessary
---------------------
The attrition audit exposed a serious problem: after curation, the active rate
is 84.4% for TRPV1, 75.1% for FAAH and 59.7% for CB1. That is not biology, it
is publication bias. ChEMBL is assembled from medicinal-chemistry papers, and
people publish compounds that worked. Genuine inactives are under-reported.

A classifier trained on 84% actives learns that nearly everything is active.
Applied to 150,000 plant natural products it will return high probabilities
almost everywhere, and the resulting "ranking" is driven by noise rather than
by signal. Class weighting does not fix this: the problem is not the loss
function, it is that the model has never been shown what a non-ligand looks
like.

The standard remedy in virtual screening is decoys: molecules presumed
inactive because nothing suggests otherwise. Real screening hit rates for a
given target are well under 1%, so a random drug-like molecule is
overwhelmingly likely to be inactive against TRPV1.

Property matching matters
-------------------------
Decoys are matched to each target's actives on molecular weight and cLogP
(DUD-E style). Without matching, the actives and decoys differ in gross
physical properties and the model learns a trivial shortcut - "heavy and
greasy means active" - which inflates measured performance and transfers
nothing. Matching forces the model to use structure.

Output: data/raw/decoy_pool.csv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import random
from concurrent.futures import ThreadPoolExecutor

from chembl_client import get_json, total_count

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

TARGET_POOL_SIZE = 60_000


def main():
    out = RAW / "decoy_pool.csv"
    if out.exists():
        print(f"{out} already exists; delete it to re-fetch")
        return

    # Drug-like molecules from ChEMBL, drawn without reference to our targets.
    # max_phase is deliberately NOT filtered: restricting to approved drugs
    # would make the decoys a biased, tiny, structurally distinctive set.
    params = {
        "molecule_properties__mw_freebase__gte": 150,
        "molecule_properties__mw_freebase__lte": 700,
        "molecule_type": "Small molecule",
    }
    total = total_count("molecule.json", params)
    print(f"{total:,} candidate molecules in ChEMBL")

    # Sample pages spread across the whole ID range rather than taking the
    # first N. ChEMBL IDs are assigned roughly chronologically, so the first
    # 60k molecules would over-represent older, pre-2000 medicinal chemistry.
    page_size = 1000
    n_pages = TARGET_POOL_SIZE // page_size
    max_offset = max(total - page_size, 0)
    rng = random.Random(42)
    offsets = sorted(rng.sample(range(0, max_offset, page_size),
                                min(n_pages, max_offset // page_size)))
    print(f"fetching {len(offsets)} pages of {page_size} spread across the range")

    recs = []
    done = [0]

    def work(off):
        p = dict(params, limit=page_size, offset=off)
        got = get_json("molecule.json", p).get("molecules", [])
        done[0] += 1
        if done[0] % 10 == 0:
            print(f"   {done[0]}/{len(offsets)} pages", flush=True)
        return got

    with ThreadPoolExecutor(max_workers=4) as pool:
        for got in pool.map(work, offsets):
            recs.extend(got)
    print(f"fetched {len(recs):,} molecule records")

    rows = []
    for r in recs:
        struct = r.get("molecule_structures") or {}
        props = r.get("molecule_properties") or {}
        smi = struct.get("canonical_smiles")
        if not smi:
            continue
        rows.append({
            "molecule_chembl_id": r.get("molecule_chembl_id"),
            "canonical_smiles": smi,
            "mw": props.get("mw_freebase"),
            "alogp": props.get("alogp"),
            "max_phase": r.get("max_phase"),
        })

    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"wrote {out}  shape={df.shape}")


if __name__ == "__main__":
    main()
