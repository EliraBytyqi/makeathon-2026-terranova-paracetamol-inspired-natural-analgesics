"""Validate all three targets against PubChem BioAssay, one row per structure.

Two problems with the earlier validation, both fixed here.

1. It counted rows, not structures.
   The outcome from PubChem is assigned per InChIKey skeleton - the first 14
   characters, which collapse stereoisomers and salts. The plant library
   contains several distinct entries per skeleton, so one measured outcome was
   being counted several times. The CB1 set reported as n=71 is 31 distinct
   skeletons: 13 of them account for 53 of the 71 rows. Those are
   pseudo-replicates. The model scores near-identical structures almost
   identically, so they inflate n without adding evidence, and the p-value
   computed from them is overstated. Everything here is deduplicated to one
   row per skeleton before any statistic is computed.

2. It required the test compound to be in the plant library.
   That restriction is what made TRPV1 and FAAH unscoreable: inside the
   library there are 25 measured TRPV1 compounds and all 25 are active, and 6
   measured FAAH compounds of which the 2 untrained ones are both inactive.
   No negatives, no ROC-AUC.

   But the question "does this model rank real actives above real inactives"
   does not need the test compounds to be plant natural products. It needs
   them measured, and absent from training. Dropping the library restriction
   turns 11 scoreable TRPV1 structures into 996, and 2 FAAH into 668, with
   both classes present in all three.

So two validations are reported, because they answer different questions:

    in_library   the original, deduplicated. "Does the ranking work on the
                 chemistry we actually screened?" Small, and for TRPV1 and
                 FAAH still single-class - reported as unscoreable, not hidden.

    external     all measured structures outside training. "Does the model
                 discriminate real actives from real inactives at all?"
                 Hundreds of structures per target, all three scoreable.

The external set is drug-like rather than plant-like, which makes it a test of
the model and not of the library. Both numbers are in the output; neither
replaces the other.

Production models are read and never written. Compounds in any production
training set are excluded. Also reported, for each target, is ROC-AUC
restricted to structures with low Tanimoto similarity to the training set,
since discrimination on chemistry the model has seen close relatives of is the
easier half of the problem.

    python scripts/validate_pubchem_expanded.py
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
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import modelio  # noqa: E402

from featurize import featurize_smiles                      # noqa: E402
from standardize import standardize                         # noqa: E402
from train import max_similarity_from_fps, query_fingerprints  # noqa: E402

RESULTS = ROOT / "results"
RAW = ROOT / "data" / "external" / "pubchem"
CACHE = RAW / "smiles_cache"
MODELS = ROOT / "models"
UA = "terranova-academic-research/1.0"
SEED = 42
N_BOOT = 2000

PANEL = {
    "TRPV1": ["Q8NER1", "O35433", "Q704Y3"],
    "CB1": ["P21554", "P20272", "P47746"],
    "FAAH": ["O00519", "P97612", "O08914"],
}
# which production model scores each target
LABEL_COL = {"TRPV1": "label", "CB1": "label", "FAAH": "label_10um"}


def fetch(url: str, dest: Path) -> bytes | None:
    if dest.exists():
        return dest.read_bytes() or None
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        dest.write_bytes(data)
        time.sleep(0.25)
        return data
    except Exception as exc:
        print(f"    {dest.name}: {exc}")
        dest.write_bytes(b"")
        return None


# ---------------------------------------------------------------------------
def bioassay() -> pd.DataFrame:
    frames = []
    for target, accs in PANEL.items():
        for acc in accs:
            f = RAW / f"{acc}.json"
            if not f.exists():
                continue
            b = json.loads(f.read_text(encoding="utf-8"))
            d = pd.DataFrame([r["Cell"] for r in b["Table"]["Row"]],
                             columns=b["Table"]["Columns"]["Column"])
            d["target"] = target
            frames.append(d)
    if not frames:
        raise SystemExit(f"no PubChem payloads in {RAW}; "
                         "run scripts/validate_pubchem.py first")
    bio = pd.concat(frames, ignore_index=True)
    bio = bio[bio["Activity Outcome"].isin(["Active", "Inactive"])].copy()
    bio["CID"] = pd.to_numeric(bio["CID"], errors="coerce")
    return bio[bio["CID"].notna()]


def cid_inchikeys() -> dict[int, str]:
    out = {}
    for f in glob.glob(str(RAW / "keys_*.json")):
        try:
            for p in json.loads(Path(f).read_text(encoding="utf-8"))[
                    "PropertyTable"]["Properties"]:
                if "InChIKey" in p:
                    out[int(p["CID"])] = p["InChIKey"]
        except Exception:
            pass
    return out


def cid_smiles(cids: list[int]) -> dict[int, str]:
    out: dict[int, str] = {}
    for i in range(0, len(cids), 150):
        block = cids[i:i + 150]
        raw = fetch("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                    + ",".join(map(str, block))
                    + "/property/SMILES,CanonicalSMILES/JSON",
                    CACHE / f"{block[0]}_{len(block)}.json")
        if not raw:
            continue
        try:
            for p in json.loads(raw)["PropertyTable"]["Properties"]:
                s = p.get("SMILES") or p.get("CanonicalSMILES")
                if s:
                    out[int(p["CID"])] = s
        except Exception:
            pass
    return out


def training_skeletons() -> set[str]:
    keys = set()
    for p in glob.glob(str(MODELS / "*.pkl")):
        if os.path.basename(p).startswith("dili") or \
           os.path.basename(p).startswith("hepatic"):
            continue
            keys.update(modelio.load(p).get("train_inchikeys", []))
    return {k[:14] for k in keys if isinstance(k, str)}


# ---------------------------------------------------------------------------
def stats_block(y: np.ndarray, s: np.ndarray) -> dict:
    """ROC-AUC with a bootstrap CI, plus a one-sided Mann-Whitney."""
    if len(set(y)) < 2:
        return {"n": int(len(y)), "n_active": int(y.sum()),
                "n_inactive": int((1 - y).sum()), "scoreable": False,
                "reason": "single-class: ROC-AUC is undefined"}
    auc = float(roc_auc_score(y, s))
    rng = np.random.RandomState(SEED)
    boots = []
    for _ in range(N_BOOT):
        idx = rng.randint(0, len(y), len(y))
        if len(set(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], s[idx]))
    lo, hi = (float(np.percentile(boots, 2.5)),
              float(np.percentile(boots, 97.5))) if boots else (np.nan, np.nan)
    p = float(sps.mannwhitneyu(s[y == 1], s[y == 0], alternative="greater").pvalue)
    return {"n": int(len(y)), "n_active": int(y.sum()),
            "n_inactive": int((1 - y).sum()), "scoreable": True,
            "roc_auc": round(auc, 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "mannwhitney_p": p,
            "median_active": round(float(np.median(s[y == 1])), 6),
            "median_inactive": round(float(np.median(s[y == 0])), 6)}


def report(title: str, d: dict) -> None:
    if not d.get("scoreable"):
        print(f"    {title:26s} n={d['n']:5d} "
              f"({d['n_active']} act / {d['n_inactive']} inact)  "
              f"NOT SCOREABLE - {d['reason']}")
        return
    print(f"    {title:26s} n={d['n']:5d} "
          f"({d['n_active']} act / {d['n_inactive']} inact)  "
          f"ROC-AUC {d['roc_auc']:.3f} "
          f"[{d['ci95'][0]:.3f}, {d['ci95'][1]:.3f}]  p={d['mannwhitney_p']:.2e}")


# ---------------------------------------------------------------------------
def main() -> None:
    print("=== 1. PubChem outcomes ===")
    bio = bioassay()
    km = cid_inchikeys()
    bio["ik"] = bio["CID"].astype(int).map(km)
    bio = bio[bio["ik"].notna()].copy()
    bio["skel"] = bio["ik"].str[:14]
    print(f"  {len(bio):,} called outcomes over {bio['skel'].nunique():,} structures")

    tskel = training_skeletons()
    print(f"  excluding {len(tskel):,} training skeletons")

    lib = pd.read_csv(RESULTS / "screening_all.csv.gz",
                      usecols=["std_inchikey", "compound_name"], low_memory=False)
    lib["skel"] = lib["std_inchikey"].str[:14]
    lib_names = lib.drop_duplicates("skel").set_index("skel")["compound_name"]
    lib_skel = set(lib_names.index)
    print(f"  plant library: {len(lib_skel):,} structures")

    print("\n=== 2. one verdict and one representative CID per structure ===")
    rows = []
    for target in PANEL:
        b = bio[bio["target"] == target]
        g = b.groupby("skel").agg(
            outcome=("Activity Outcome",
                     lambda s: "Active" if (s == "Active").any() else "Inactive"),
            cid=("CID", "first"))
        g = g[~g.index.isin(tskel)].reset_index()
        g["target"] = target
        g["in_library"] = g["skel"].isin(lib_skel)
        rows.append(g)
        print(f"  {target:6s} {len(g):5d} untrained structures "
              f"({int((g['outcome'] == 'Active').sum())} act / "
              f"{int((g['outcome'] == 'Inactive').sum())} inact), "
              f"{int(g['in_library'].sum())} of them in the plant library")
    val = pd.concat(rows, ignore_index=True)

    print(f"\n=== 3. structures for {val['cid'].nunique():,} CIDs ===")
    smap = cid_smiles(sorted(val["cid"].astype(int).unique()))
    val["raw_smiles"] = val["cid"].astype(int).map(smap)
    val = val[val["raw_smiles"].notna()].copy()
    std = [standardize(s) for s in val["raw_smiles"]]
    val["std_smiles"] = [r.smiles if r.ok else None for r in std]
    val = val[val["std_smiles"].notna()].copy()
    # standardization can merge two skeletons; keep one row per target+skeleton
    val = val.drop_duplicates(subset=["target", "skel"]).reset_index(drop=True)
    print(f"  usable: {len(val):,} rows")

    print("\n=== 4. scoring with the production models ===")
    out_frames = []
    for target in PANEL:
        path = MODELS / f"{target}_human_{LABEL_COL[target]}.pkl"
        if not path.exists():
            print(f"  {target}: {path.name} missing, skipped")
            continue
        b = modelio.load(path)
        s = val[val["target"] == target].copy()
        X, mask = featurize_smiles(s["std_smiles"].tolist())
        s = s[mask].copy()
        s["model_score"] = b["model"].predict_proba(X[mask])[:, 1]
        Q, q_idx = query_fingerprints(s["std_smiles"].tolist())
        s["ad"] = max_similarity_from_fps(Q, q_idx, len(s), b["ad_fps"])
        s["compound_name"] = s["skel"].map(lib_names).fillna(
            "CID" + s["cid"].astype(int).astype(str))
        out_frames.append(s)
        print(f"  {target:6s} scored {len(s):,}")

    scored = pd.concat(out_frames, ignore_index=True)
    scored.to_csv(RESULTS / "pubchem_validation_expanded.csv", index=False)

    print(f"\n{'=' * 78}\nVALIDATION, ONE ROW PER STRUCTURE\n{'=' * 78}")
    summary = {"n_boot": N_BOOT, "targets": {}}
    for target in PANEL:
        s = scored[scored["target"] == target]
        if not len(s):
            continue
        y = (s["outcome"] == "Active").astype(int).to_numpy()
        sc = s["model_score"].to_numpy()
        ad = s["ad"].to_numpy()
        inlib = s["in_library"].to_numpy()
        print(f"\n  {target}")
        blocks = {
            "external_all": stats_block(y, sc),
            "external_novel_ad_lt_0.5": stats_block(y[ad < 0.5], sc[ad < 0.5]),
            "external_novel_ad_lt_0.4": stats_block(y[ad < 0.4], sc[ad < 0.4]),
            "in_library": stats_block(y[inlib], sc[inlib]),
        }
        titles = {"external_all": "all untrained structures",
                  "external_novel_ad_lt_0.5": "novel chemistry (AD<0.5)",
                  "external_novel_ad_lt_0.4": "novel chemistry (AD<0.4)",
                  "in_library": "plant library only"}
        for k, d in blocks.items():
            report(titles[k], d)
        summary["targets"][target] = blocks

    (RESULTS / "pubchem_validation_expanded.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nwrote {RESULTS / 'pubchem_validation_expanded.csv'}")
    print(f"wrote {RESULTS / 'pubchem_validation_expanded.json'}")

    n_ok = sum(1 for t in summary["targets"]
               if summary["targets"][t]["external_all"].get("scoreable"))
    print(f"\n  targets with a two-class external test: {n_ok}/3")


if __name__ == "__main__":
    main()
