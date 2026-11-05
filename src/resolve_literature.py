"""Resolve the curated literature compound set to structures via PubChem.

The literature set (data/reference/literature_actives.csv) is deliberately NOT
training data. It is an external validation set: compounds whose activity on
TRPV1 / CB1 / FAAH is reported in the primary literature, many of which are
natural products, plus synthetic reference ligands used as positive controls.

Keeping it separate matters. If these compounds were folded into training, the
project's headline claim - "the model recovers known natural analgesic
chemistry" - would be circular.
"""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd
import requests

from standardize import standardize

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data" / "reference"
OUT = ROOT / "data" / "interim"
OUT.mkdir(parents=True, exist_ok=True)

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
HEADERS = {"User-Agent": "terranova-np-screening/0.1 (academic cheminformatics)"}

# PubChem name lookup is unreliable for a few of these; pin them explicitly.
MANUAL_CID = {
    "AM404": 4666781,
    "N-benzyl-linoleamide": 6442400,
    "Delta-9-THC": 16078,
    "PF-3845": 25021165,
    "7-Hydroxyflavone": 5281894,
}


def cid_for_name(name: str) -> int | None:
    if name in MANUAL_CID:
        return MANUAL_CID[name]
    url = f"{PUBCHEM}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    for _ in range(5):
        try:
            r = requests.get(url, headers=HEADERS, timeout=45)
            if r.status_code == 200:
                cids = r.json().get("IdentifierList", {}).get("CID", [])
                return cids[0] if cids else None
            if r.status_code == 404:
                return None
        except requests.exceptions.RequestException:
            pass
        time.sleep(1.5)
    return None


def smiles_for_cid(cid: int) -> str | None:
    """Fetch an isomeric SMILES for a CID.

    PubChem renamed the SMILES properties (IsomericSMILES/CanonicalSMILES ->
    SMILES/ConnectivitySMILES) in 2025, so try the modern names first and fall
    back to the legacy ones rather than assuming either.
    """
    for prop in ("SMILES", "IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES"):
        url = f"{PUBCHEM}/compound/cid/{cid}/property/{prop}/JSON"
        for _ in range(3):
            try:
                r = requests.get(url, headers=HEADERS, timeout=45)
                if r.status_code == 200:
                    props = r.json()["PropertyTable"]["Properties"][0]
                    val = props.get(prop)
                    if val:
                        return val
                    break
                if r.status_code in (400, 404):
                    break
            except (requests.exceptions.RequestException, KeyError, IndexError):
                pass
            time.sleep(1.5)
    return None


def main():
    lit = pd.read_csv(REF / "literature_actives.csv")
    names = sorted(lit["compound_name"].unique())
    print(f"resolving {len(names)} unique compounds via PubChem\n")

    rows = []
    for name in names:
        cid = cid_for_name(name)
        smi = smiles_for_cid(cid) if cid else None
        std = standardize(smi) if smi else None
        ok = bool(std and std.ok)
        rows.append(
            {
                "compound_name": name,
                "pubchem_cid": cid,
                "raw_smiles": smi,
                "std_smiles": std.smiles if ok else None,
                "inchikey": std.inchikey if ok else None,
                "resolve_status": "ok" if ok else (std.reason if std else "no_structure"),
            }
        )
        print(f"  {name:24s} CID={str(cid):>10s}  {rows[-1]['resolve_status']}")
        time.sleep(0.3)

    struct = pd.DataFrame(rows)
    merged = lit.merge(struct, on="compound_name", how="left")
    out = OUT / "literature_actives_resolved.csv"
    merged.to_csv(out, index=False)

    n_ok = struct["resolve_status"].eq("ok").sum()
    print(f"\nresolved {n_ok}/{len(struct)} compounds")
    failed = struct[struct["resolve_status"] != "ok"]["compound_name"].tolist()
    if failed:
        print(f"UNRESOLVED (need manual CIDs): {failed}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
