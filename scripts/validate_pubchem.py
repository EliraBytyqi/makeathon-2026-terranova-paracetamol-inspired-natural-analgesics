"""Retrospective validation of the screen against PubChem BioAssay.

This is the closest thing to a real experimental test the project can run
without a laboratory, and it is the only analysis here with a proper
denominator.

Every other validation in this project suffers the same flaw: the comparison
sets contain only compounds somebody bothered to publish as active. PubChem
BioAssay records the Activity Outcome explicitly - Active AND Inactive - so a
compound that was screened and did nothing is in the data. That means we can
ask the question that actually matters:

    Of the library compounds that have since been measured against our
    targets, are the ones our model ranked highly more likely to be
    genuinely active than the ones it ranked low?

That is a real hit-rate comparison on real experimental outcomes, with real
negatives. Compounds present in any production training set are excluded, so
nothing the model was taught counts as a prediction.

    python scripts/validate_pubchem.py
"""
from __future__ import annotations

import glob
import json
import os
import pickle
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RAW = ROOT / "data" / "external" / "pubchem"
RAW.mkdir(parents=True, exist_ok=True)

UA = "terranova-academic-research/1.0"
PANEL = {
    "TRPV1": ["Q8NER1", "O35433", "Q704Y3"],
    "CB1": ["P21554", "P20272", "P47746"],
    "FAAH": ["O00519", "P97612", "O08914"],
}


def fetch(url: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 200:
        return True
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=180) as r:
            dest.write_bytes(r.read())
        time.sleep(0.4)          # be polite to PubChem
        return True
    except Exception as exc:
        print(f"    {dest.name}: {exc}")
        return False


def load_bioactivity() -> pd.DataFrame:
    frames = []
    for target, accs in PANEL.items():
        for acc in accs:
            dest = RAW / f"{acc}.json"
            url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/protein/"
                   f"accession/{acc}/concise/JSON")
            if not fetch(url, dest):
                continue
            try:
                blob = json.loads(dest.read_text(encoding="utf-8"))
                cols = blob["Table"]["Columns"]["Column"]
                rows = [r["Cell"] for r in blob["Table"]["Row"]]
            except Exception as exc:
                print(f"    {acc}: unreadable ({exc})")
                continue
            df = pd.DataFrame(rows, columns=cols)
            df["target"] = target
            df["accession"] = acc
            frames.append(df)
            n_act = (df["Activity Outcome"] == "Active").sum()
            n_in = (df["Activity Outcome"] == "Inactive").sum()
            print(f"    {target:6s} {acc}: {len(df):7,} records "
                  f"({n_act:,} active / {n_in:,} inactive)")
    if not frames:
        raise SystemExit("no PubChem data retrieved")
    return pd.concat(frames, ignore_index=True)


def cids_to_inchikey(cids: list[int]) -> dict[int, str]:
    """Batch CID -> InChIKey through PUG REST."""
    out: dict[int, str] = {}
    cids = [int(c) for c in cids]
    CHUNK = 400
    for i in range(0, len(cids), CHUNK):
        block = cids[i:i + CHUNK]
        dest = RAW / f"keys_{i:06d}.json"
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
               + ",".join(map(str, block)) + "/property/InChIKey/JSON")
        if not fetch(url, dest):
            continue
        try:
            blob = json.loads(dest.read_text(encoding="utf-8"))
            for p in blob["PropertyTable"]["Properties"]:
                if "InChIKey" in p:
                    out[int(p["CID"])] = p["InChIKey"]
        except Exception:
            pass
        if (i // CHUNK) % 5 == 0:
            print(f"      resolved {len(out):,} / {len(cids):,} CIDs", flush=True)
    return out


def main() -> None:
    print("=== 1. pulling PubChem BioAssay records for the panel ===")
    bio = load_bioactivity()
    bio = bio[bio["Activity Outcome"].isin(["Active", "Inactive"])].copy()
    bio["CID"] = pd.to_numeric(bio["CID"], errors="coerce")
    bio = bio[bio["CID"].notna()]
    print(f"  total usable outcome records: {len(bio):,} "
          f"({(bio['Activity Outcome']=='Active').sum():,} active / "
          f"{(bio['Activity Outcome']=='Inactive').sum():,} inactive)")
    print(f"  unique compounds: {bio['CID'].nunique():,}")

    print("\n=== 2. resolving CIDs to InChIKeys ===")
    key_map = cids_to_inchikey(sorted(bio["CID"].unique()))
    bio["inchikey"] = bio["CID"].map(key_map)
    bio = bio[bio["inchikey"].notna()]
    bio["skeleton"] = bio["inchikey"].str[:14]
    print(f"  resolved {bio['inchikey'].nunique():,} structures")

    print("\n=== 3. joining to the screen ===")
    screen = pd.read_csv(RESULTS / "screening_all.csv",
                         usecols=["std_inchikey", "compound_name", "best_target",
                                  "target_score", "composite_score", "in_domain",
                                  "p_TRPV1", "p_CB1", "p_FAAH", "tier"],
                         low_memory=False)
    screen["skeleton"] = screen["std_inchikey"].str[:14]

    train_keys = set()
    for p in glob.glob(str(ROOT / "models" / "*.pkl")):
        if os.path.basename(p) == "dili.pkl":
            continue
        with open(p, "rb") as fh:
            train_keys.update(pickle.load(fh).get("train_inchikeys", []))
    train_skel = {k[:14] for k in train_keys if isinstance(k, str)}
    print(f"  training compounds to exclude: {len(train_skel):,} skeletons")

    rows = []
    for target in PANEL:
        b = bio[bio["target"] == target]
        # one verdict per compound: Active if ANY assay called it active
        verdict = (b.groupby("skeleton")["Activity Outcome"]
                     .apply(lambda s: "Active" if (s == "Active").any() else "Inactive"))
        s = screen[~screen["skeleton"].isin(train_skel)].copy()
        s["outcome"] = s["skeleton"].map(verdict)
        s = s[s["outcome"].notna()]
        s["target"] = target
        s["model_score"] = s[f"p_{target}"]
        rows.append(s)

    joined = pd.concat(rows, ignore_index=True)
    joined.to_csv(RESULTS / "pubchem_validation.csv", index=False)
    print(f"  library compounds with an independent measured outcome, "
          f"never trained on: {len(joined):,}")

    if len(joined) < 10:
        print("\n  too few overlaps for a hit-rate analysis")
        return

    print("\n=== 4. DOES OUR RANKING ENRICH FOR REAL ACTIVES? ===")
    summary = {"n_total": int(len(joined))}
    for target in PANEL:
        s = joined[joined["target"] == target]
        if len(s) < 10:
            print(f"  {target}: only {len(s)} compounds, skipped")
            continue
        s = s.sort_values("model_score", ascending=False)
        y = (s["outcome"] == "Active").astype(int).to_numpy()
        base = y.mean()
        n = len(y)
        res = {"n": int(n), "base_rate": round(float(base), 4)}
        print(f"\n  {target}: {n:,} measured compounds, "
              f"{int(y.sum()):,} active ({base:.1%} base rate)")
        for frac in (0.01, 0.05, 0.10, 0.25):
            k = max(1, int(round(n * frac)))
            hit = y[:k].mean()
            ef = hit / base if base else float("nan")
            res[f"ef_{int(frac*100)}pct"] = round(float(ef), 2)
            res[f"hitrate_{int(frac*100)}pct"] = round(float(hit), 4)
            print(f"     top {frac:4.0%}: {int(y[:k].sum()):4d}/{k:4d} active "
                  f"= {hit:6.1%}   enrichment {ef:5.2f}x")
        try:
            auc = float(sps.mannwhitneyu(
                s[s["outcome"] == "Active"]["model_score"],
                s[s["outcome"] == "Inactive"]["model_score"],
                alternative="greater").pvalue)
            from sklearn.metrics import roc_auc_score, average_precision_score
            res["roc_auc"] = round(float(roc_auc_score(y, s["model_score"])), 4)
            res["pr_auc"] = round(float(average_precision_score(y, s["model_score"])), 4)
            res["mannwhitney_p"] = round(auc, 6)
            print(f"     ROC-AUC {res['roc_auc']:.3f}   PR-AUC {res['pr_auc']:.3f}   "
                  f"Mann-Whitney p = {auc:.2e}")
        except Exception as exc:
            print(f"     stats unavailable: {exc}")
        summary[target] = res

    (RESULTS / "pubchem_validation.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nwrote {RESULTS/'pubchem_validation.csv'}")


if __name__ == "__main__":
    main()
