"""Where does the hepatic panel actually stop working? Measure, don't assume.

The gate on the hepatic term needs a threshold. Two ways to pick one:

  linear shrinkage   p_eff = 0.5 + min(AD/0.4, 1) * (p - 0.5).
                     Assumes confidence rises linearly with similarity, which
                     nothing here has shown.
  hard gate          use the model where it has demonstrated skill, substitute
                     a neutral 0.5 everywhere else.

The second is the more honest shape, but only if the cut-off comes from data.
So this bins the plant-natural-product validation by applicability domain and
reports discrimination in each bin, per endpoint, against that endpoint's own
measured label. If performance decays smoothly a shrinkage function is
defensible; if it falls off a cliff, a hard gate at the cliff is right.

Each head is scored on its own endpoint, not on a pooled "positive in any
endpoint" label - that pooled label was a mistake in the earlier run, because
most structures carry only one measured endpoint, so the maximum pulled in
predictions for assays the compound was never put through.

Caveat that cannot be fixed here: the 1,367 plant compounds with hepatic
labels sit at median AD 0.562, while the screening library sits at 0.226.
Even the lowest bin below is less alien than most of what gets screened, so
these numbers are an upper bound on out-of-domain performance.

    python scripts/hepatic_ad_bins.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.metrics import brier_score_loss, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from featurize import featurize_smiles                       # noqa: E402
from train import build_ad_fps, make_models, max_similarity_from_fps, \
    query_fingerprints                                       # noqa: E402

DATA = ROOT / "data" / "processed" / "hepatic_endpoints.csv"
RESULTS = ROOT / "results"
SEED = 42
N_BOOT = 2000
HEADS = ("mito", "hepg2", "clinical")
BINS = [(0.0, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 1.01)]


def block(y, s) -> dict:
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    if len(y) < 15 or len(set(y)) < 2:
        return {"n": int(len(y)), "n_pos": int(y.sum()), "scoreable": False}
    auc = float(roc_auc_score(y, s))
    rng = np.random.RandomState(SEED)
    b = [roc_auc_score(y[i], s[i]) for i in
         (rng.randint(0, len(y), len(y)) for _ in range(N_BOOT))
         if len(set(y[i])) > 1]
    lo, hi = (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))) if b else (np.nan, np.nan)
    return {"n": int(len(y)), "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()),
            "scoreable": True, "roc_auc": round(auc, 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "ci_excludes_chance": bool(lo > 0.5),
            "brier": round(float(brier_score_loss(y, s)), 4),
            "p": float(sps.mannwhitneyu(s[y == 1], s[y == 0],
                                        alternative="greater").pvalue)}


def main() -> None:
    df = pd.read_csv(DATA)
    df["skel"] = df["inchikey"].str[:14]
    lib = pd.read_csv(RESULTS / "screening_all.csv.gz",
                      usecols=["std_inchikey"], low_memory=False)
    lib_skel = set(lib["std_inchikey"].str[:14])
    is_np = df["skel"].isin(lib_skel)
    test, train = df[is_np].copy(), df[~is_np].copy()

    out: dict = {"bins": {}, "note": (
        "plant natural products only; each head scored on its own measured "
        "endpoint; heads retrained with every library structure removed")}

    for name in HEADS:
        tr = train[train["endpoint"] == name]
        Xtr, mtr = featurize_smiles(tr["std_smiles"].tolist())
        tr = tr[mtr]
        ytr = tr["y"].to_numpy().astype(int)
        if len(set(ytr)) < 2:
            continue
        model = make_models()["random_forest"].fit(Xtr[mtr], ytr)
        ad_fps = build_ad_fps(tr["std_smiles"].tolist())

        te = test[test["endpoint"] == name].drop_duplicates("skel")
        Xte, mte = featurize_smiles(te["std_smiles"].tolist())
        te = te[mte].reset_index(drop=True)
        y = te["y"].to_numpy().astype(int)
        s = model.predict_proba(Xte[mte])[:, 1]
        Q, qi = query_fingerprints(te["std_smiles"].tolist())
        ad = max_similarity_from_fps(Q, qi, len(te), ad_fps)

        print(f"\n{'=' * 78}\n{name}  (trained on {len(ytr):,}, tested on "
              f"{len(y):,} plant structures)\n{'=' * 78}")
        overall = block(y, s)
        print(f"  {'ALL':14s} n={overall['n']:5d} "
              f"({overall['n_pos']:4d} pos) ROC {overall['roc_auc']:.3f} "
              f"[{overall['ci95'][0]:.3f}, {overall['ci95'][1]:.3f}]")
        print(f"  {'AD bin':14s} {'n':>5s} {'pos':>5s} {'ROC-AUC':>8s} "
              f"{'95% CI':>16s} {'Brier':>7s}  above chance?")
        rows = {"overall": overall}
        for lo, hi in BINS:
            sel = (ad >= lo) & (ad < hi)
            d = block(y[sel], s[sel])
            key = f"{lo:.1f}-{hi:.1f}"
            rows[key] = d
            if not d.get("scoreable"):
                print(f"  {key:14s} {d['n']:5d} {d['n_pos']:5d}  "
                      f"{'too few / single-class':>33s}")
                continue
            flag = "YES" if d["ci_excludes_chance"] else "no"
            print(f"  {key:14s} {d['n']:5d} {d['n_pos']:5d} {d['roc_auc']:8.3f} "
                  f"[{d['ci95'][0]:.3f}, {d['ci95'][1]:.3f}] {d['brier']:7.3f}  {flag}")
        out["bins"][name] = rows

    # where is the library, relative to these bins?
    print(f"\n{'=' * 78}\nWHERE THE SCREENING LIBRARY ACTUALLY SITS\n{'=' * 78}")
    s = pd.read_csv(RESULTS / "screening_all.csv.gz",
                    usecols=["ad_dili"], low_memory=False)["ad_dili"].dropna()
    lib_bins = {}
    for lo, hi in BINS:
        n = int(((s >= lo) & (s < hi)).sum())
        lib_bins[f"{lo:.1f}-{hi:.1f}"] = {"n": n, "frac": round(n / len(s), 4)}
        print(f"  AD {lo:.1f}-{hi:.1f}  {n:8,d}  {n / len(s):6.1%}")
    out["library_ad_distribution"] = lib_bins

    (RESULTS / "hepatic_ad_bins.json").write_text(json.dumps(out, indent=2),
                                                  encoding="utf-8")
    print(f"\nwrote {RESULTS / 'hepatic_ad_bins.json'}")


if __name__ == "__main__":
    main()
