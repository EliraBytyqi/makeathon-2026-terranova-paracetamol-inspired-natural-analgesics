"""Does the hepatic panel work on plant chemistry? The test that was missing.

Every hepatic number so far comes from marketed drugs, while the screen is
applied to 145,164 plant natural products. Marketed-drug models failing on
natural-product chemistry is a well-known failure mode, so a panel validated
only on drugs is a panel validated on the wrong thing.

It turns out the test is available. 1,367 structures in the screening library
also carry a measured hepatic outcome - 116 in DILIrank, 844 in the HepG2
cytotoxicity assay, 1,328 in the mitochondrial assay - and all three are
two-class. They are currently inside the training data, which is why they were
never a test. So: remove every plant-library structure from training, retrain
all three heads on what is left, and score them.

The combined label mirrors how the screen combines the heads. The screen takes
the maximum over heads, so a structure counts as positive if any endpoint
measured it positive, and the prediction it is scored against is the maximum
of the three predictions.

Reported alongside is the same split by applicability domain used everywhere
else, because 65.7% of the library sits below AD 0.4 against hepatic training
data and the number for that region is the one that governs the screen.

    python scripts/validate_hepatic_naturalproducts.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from featurize import featurize_smiles                       # noqa: E402
from train import build_ad_fps, make_models, max_similarity_from_fps, \
    query_fingerprints                                       # noqa: E402

DATA = ROOT / "data" / "processed" / "hepatic_endpoints.csv"
RESULTS = ROOT / "results"
SEED = 42
N_BOOT = 2000
HEADS = ("clinical", "mito", "hepg2")


def stats_block(y, s) -> dict:
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    if len(y) < 10 or len(set(y)) < 2:
        return {"n": int(len(y)), "scoreable": False, "reason": "too few or single-class"}
    auc = float(roc_auc_score(y, s))
    rng = np.random.RandomState(SEED)
    b = []
    for _ in range(N_BOOT):
        i = rng.randint(0, len(y), len(y))
        if len(set(y[i])) > 1:
            b.append(roc_auc_score(y[i], s[i]))
    return {"n": int(len(y)), "n_positive": int(y.sum()),
            "n_negative": int((1 - y).sum()), "scoreable": True,
            "roc_auc": round(auc, 4),
            "ci95": [round(float(np.percentile(b, 2.5)), 4),
                     round(float(np.percentile(b, 97.5)), 4)] if b else None,
            "mannwhitney_p": float(sps.mannwhitneyu(
                s[y == 1], s[y == 0], alternative="greater").pvalue)}


def show(label, d):
    if not d.get("scoreable"):
        print(f"  {label:34s} n={d['n']:5d}  not scoreable ({d['reason']})")
        return
    ci = d["ci95"]
    print(f"  {label:34s} n={d['n']:5d} ({d['n_positive']:4d} pos/{d['n_negative']:4d} neg)"
          f"  ROC {d['roc_auc']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]  p={d['mannwhitney_p']:.2e}")


def main() -> None:
    df = pd.read_csv(DATA)
    df["skel"] = df["inchikey"].str[:14]

    lib = pd.read_csv(RESULTS / "screening_all.csv.gz",
                      usecols=["std_inchikey"], low_memory=False)
    lib_skel = set(lib["std_inchikey"].str[:14])

    is_np = df["skel"].isin(lib_skel)
    print(f"hepatic endpoint rows: {len(df):,}")
    print(f"  of which are plant-library structures: {int(is_np.sum()):,} "
          f"({df[is_np]['skel'].nunique():,} unique)\n")

    test = df[is_np].copy()
    train = df[~is_np].copy()

    preds, train_smiles = {}, []
    print("retraining each head with ALL plant-library structures removed")
    for name in HEADS:
        tr = train[train["endpoint"] == name]
        X, m = featurize_smiles(tr["std_smiles"].tolist())
        tr = tr[m]
        y = tr["y"].to_numpy().astype(int)
        if len(set(y)) < 2:
            print(f"  {name:9s} single-class after exclusion, skipped")
            continue
        model = make_models()["random_forest"].fit(X[m], y)
        preds[name] = model
        train_smiles.extend(tr["std_smiles"].tolist())
        print(f"  {name:9s} trained on {len(y):5d} ({int(y.sum())} positive)")

    ad_fps = build_ad_fps(train_smiles)
    out = {"design": ("all three heads retrained with every screening-library "
                      "structure removed, then scored on those structures"),
           "per_endpoint": {}, "combined": {}}

    print("\n=== each head, on plant natural products it has never seen ===")
    for name in HEADS:
        if name not in preds:
            continue
        sub = test[test["endpoint"] == name].drop_duplicates("skel")
        X, m = featurize_smiles(sub["std_smiles"].tolist())
        sub = sub[m]
        y = sub["y"].to_numpy().astype(int)
        s = preds[name].predict_proba(X[m])[:, 1]
        d = stats_block(y, s)
        out["per_endpoint"][name] = d
        show(name, d)

    # ---- the combined score, exactly as the screen forms it ---------------
    uni = test.drop_duplicates("skel")[["skel", "std_smiles"]].reset_index(drop=True)
    pos = test.groupby("skel")["y"].max()
    uni["y"] = uni["skel"].map(pos).astype(int)
    X, m = featurize_smiles(uni["std_smiles"].tolist())
    uni = uni[m].reset_index(drop=True)
    X = X[m]
    stack = np.column_stack([preds[n].predict_proba(X)[:, 1] for n in preds])
    combined = stack.max(axis=1)
    y = uni["y"].to_numpy()

    Q, qi = query_fingerprints(uni["std_smiles"].tolist())
    ad = max_similarity_from_fps(Q, qi, len(uni), ad_fps)

    print("\n=== the quantity the screen uses: max over heads ===")
    d = stats_block(y, combined)
    out["combined"]["all"] = d
    show("COMBINED, all plant structures", d)
    print(f"  median applicability domain: {np.median(ad):.3f}")

    print("\n=== combined, by applicability domain ===")
    for lo, hi, label in [(0.0, 0.4, "AD < 0.4  (novel)"),
                          (0.4, 1.01, "AD >= 0.4 (familiar)")]:
        sel = (ad >= lo) & (ad < hi)
        d = stats_block(y[sel], combined[sel])
        out["combined"][label.split()[0] + label.split()[1]] = d
        show(label, d)

    out["ad_median"] = round(float(np.median(ad)), 4)
    (RESULTS / "hepatic_validation_naturalproducts.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {RESULTS / 'hepatic_validation_naturalproducts.json'}")


if __name__ == "__main__":
    main()
