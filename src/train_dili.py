"""Train the drug-induced liver injury (DILI) model - safety layer B.

Data: TDC's DILI benchmark (475 marketed drugs, 236 DILI-positive / 239
negative), aggregated by the FDA National Center for Toxicological Research
and derived from Xu et al. 2015.

What this model is
------------------
A weak, literature-derived prior over hepatic liability, trained on 475
compounds. That is a very small dataset by any standard, and the chemical
space is marketed small-molecule drugs - not plant natural products. Scores
on flavonoid glycosides and terpenoids are extrapolation.

We therefore:
  * evaluate under a scaffold split, where the number is noticeably worse
    than under a random split, and report both;
  * carry an applicability-domain similarity alongside every prediction;
  * treat the output as one input to a ranking, never as a safety verdict.

A model like this is useful for deprioritising, not for clearing. Saying so
explicitly is more scientifically defensible than quoting an AUROC and moving
on, and DILI is notoriously hard to predict from structure alone because much
of the clinical risk is idiosyncratic, immune-mediated and dose-dependent.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupShuffleSplit, StratifiedKFold

from featurize import featurize_smiles, scaffold_groups
from standardize import standardize
from train import build_ad_fps, evaluate

ROOT = Path(__file__).resolve().parents[1]
DILI_CSV = ROOT / "data" / "external" / "dili" / "tdc_dili.csv"
MODELS = ROOT / "models"
RESULTS = ROOT / "results"
MODELS.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

RANDOM_STATE = 42


def main():
    df = pd.read_csv(DILI_CSV)
    print(f"TDC DILI: {len(df)} compounds, {int(df['Y'].sum())} positive")

    # identical standardization to every other structure in the project
    std = [standardize(s) for s in df["SMILES"]]
    df["std_smiles"] = [r.smiles for r in std]
    ok = df["std_smiles"].notna()
    print(f"standardized OK: {int(ok.sum())}/{len(df)}")
    df = df[ok].reset_index(drop=True)

    X, mask = featurize_smiles(df["std_smiles"].tolist())
    X = X[mask]
    y = df["Y"].values[mask]
    smiles = [s for s, m in zip(df["std_smiles"], mask) if m]
    groups = scaffold_groups(smiles)
    print(f"{len(y)} usable compounds, {len(set(groups))} scaffolds, "
          f"positive rate {y.mean():.1%}")

    def model():
        return RandomForestClassifier(
            n_estimators=700, min_samples_leaf=2, max_features="sqrt",
            class_weight="balanced_subsample", n_jobs=-1,
            random_state=RANDOM_STATE,
        )

    report = {"n": int(len(y)), "positive_rate": round(float(y.mean()), 4),
              "n_scaffolds": int(len(set(groups))), "splits": {}}

    # scaffold split: repeated, because 475 compounds makes a single split noisy
    scaff_metrics = []
    for seed in range(5):
        gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed)
        tr, te = next(gss.split(X, y, groups=groups))
        if len(set(y[te])) < 2:
            continue
        m = model().fit(X[tr], y[tr])
        scaff_metrics.append(evaluate(y[te], m.predict_proba(X[te])[:, 1]))

    def mean_metrics(ms):
        keys = ["pr_auc", "roc_auc", "mcc", "balanced_acc"]
        return {k: round(float(np.mean([m[k] for m in ms])), 4) for k in keys} | {
            f"{k}_sd": round(float(np.std([m[k] for m in ms])), 4) for k in keys
        }

    report["splits"]["scaffold_5x"] = mean_metrics(scaff_metrics)
    print("\nscaffold split (5 repeats, mean +/- sd):")
    s = report["splits"]["scaffold_5x"]
    for k in ["pr_auc", "roc_auc", "mcc", "balanced_acc"]:
        print(f"   {k:14s} {s[k]:.3f} +/- {s[f'{k}_sd']:.3f}")

    # random split for contrast
    rnd_metrics = []
    for seed in range(5):
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        tr, te = next(skf.split(X, y))
        m = model().fit(X[tr], y[tr])
        rnd_metrics.append(evaluate(y[te], m.predict_proba(X[te])[:, 1]))
    report["splits"]["random_5x"] = mean_metrics(rnd_metrics)
    r = report["splits"]["random_5x"]
    print("\nrandom split (5 repeats, mean):")
    for k in ["pr_auc", "roc_auc", "mcc", "balanced_acc"]:
        print(f"   {k:14s} {r[k]:.3f} +/- {r[f'{k}_sd']:.3f}")
    print(f"\n  generalisation gap (ROC-AUC): "
          f"{r['roc_auc'] - s['roc_auc']:+.3f}")

    final = model().fit(X, y)
    bundle = {
        "model": final, "model_name": "random_forest", "task": "DILI",
        "train_smiles": smiles, "ad_fps": build_ad_fps(smiles),
        "report": report,
        "caveat": (
            "Trained on 475 marketed drugs. Natural products are largely "
            "outside this chemical space. Use for deprioritisation only; "
            "this is not a safety assessment."
        ),
    }
    path = MODELS / "dili.pkl"
    with open(path, "wb") as fh:
        pickle.dump(bundle, fh)
    (RESULTS / "dili_report.json").write_text(json.dumps(report, indent=2),
                                              encoding="utf-8")
    print(f"\nsaved {path}")

    # Sanity probe. Verified with check_dili_labels.py: paracetamol and
    # isoniazid - two of the most famous hepatotoxins in medicine - are NOT in
    # this benchmark at all, while diclofenac, troglitazone and valproate are
    # positives and aspirin is a negative. So any score this model gives
    # paracetamol is extrapolation, and it gets it wrong (~0.37, i.e. "low
    # risk"). That is the clearest possible demonstration of why the structural
    # alert layer is not redundant: the alert set flags paracetamol at maximum
    # severity on the exact mechanism that makes it hepatotoxic, and this
    # statistical model cannot, because the benchmark encodes idiosyncratic
    # DILI rather than intrinsic dose-dependent toxicity.
    probes = {
        "paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "isoniazid": "NNC(=O)c1ccncc1",
        "curcumin": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
        "kaempferol": "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
        "capsaicin": "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
        "caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    }
    Xp, mp = featurize_smiles(list(probes.values()))
    pr = final.predict_proba(Xp[mp])[:, 1]
    print("\nDILI probability probes:")
    for name, p in zip(probes, pr):
        print(f"   {name:14s} {p:.3f}")
    print("\n  NOTE: paracetamol is absent from the TDC DILI benchmark, so its")
    print("  score here is extrapolation and reads as low risk. The structural")
    print("  alert layer flags it at severity 3. Use both layers, not either.")


if __name__ == "__main__":
    main()
