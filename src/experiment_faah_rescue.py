"""Can a handful of curated natural actives rescue the FAAH blind spot?

A chemical-space diagnostic established that ChEMBL's FAAH
records contain essentially no flavonoids (2 of 4,979 rows, none active), and
that the resulting model scores kaempferol at 6e-6 and genistein at 6.5e-5.
That is extrapolation into empty space, not a judgement.

This experiment asks the obvious follow-up: **is the blind spot fixable with a
small amount of curated literature data?**

Protocol - leave-one-out over the known natural FAAH inhibitors:

    for each natural inhibitor i:
        train on  ChEMBL FAAH data + decoys + (the other natural inhibitors)
        test on   compound i, which the model has never seen

Compound i is never in its own training set, so its score is an honest
out-of-sample prediction. Comparing it against the baseline model (trained on
ChEMBL alone) measures exactly how much a handful of literature actives buys.

This is a real experiment with a real answer either way. If LOO scores stay
near zero, the flavonoid chemotype is simply unlearnable from this
representation and the right conclusion is to use similarity search instead. If
they jump, the remedy is cheap curation.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier

from decoys import match_decoys, prepare_pool
from featurize import featurize_smiles, scaffold
from standardize import standardize

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

RANDOM_STATE = 42
DECOYS_PER_ACTIVE = 3

# The natural-product FAAH inhibitors that motivate the project, with the
# potency reported in the primary literature (micromolar).
NATURAL_FAAH = {
    "kaempferol": ("O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12", 5.0),
    "genistein": ("O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12", 1.3),
    "biochanin A": ("COc1ccc(-c2coc3cc(O)cc(O)c3c2=O)cc1", 2.4),
    "daidzein": ("O=c1c(-c2ccc(O)cc2)coc2cc(O)ccc12", 14.0),
    "7-hydroxyflavone": ("O=c1cc(-c2ccccc2)oc2cc(O)ccc12", 2.04),
    "silybin": ("COc1cc(C2Oc3ccc(C4Oc5cc(O)cc(O)c5C(=O)C4O)cc3OC2CO)ccc1O", 5.08),
}


def model():
    return LGBMClassifier(
        n_estimators=600, learning_rate=0.05, num_leaves=63,
        min_child_samples=10, subsample=0.8, subsample_freq=1,
        colsample_bytree=0.5, class_weight="balanced",
        random_state=RANDOM_STATE, n_jobs=-1, verbose=-1,
    )


def main():
    agg = pd.read_csv(PROC / "bioactivity_aggregated.csv")
    faah = agg[(agg["target_name"] == "FAAH")
               & (agg["organism_group"] == "human")]
    # use the 10 uM label: these natural inhibitors sit at 1-15 uM, so the
    # 1 uM cut is the wrong question for them (see build_dataset.py)
    faah = faah.dropna(subset=["label_10um"]).copy()
    faah = faah.drop_duplicates(subset="inchikey")
    print(f"ChEMBL FAAH/human (10 uM label): {len(faah):,} compounds, "
          f"{int(faah['label_10um'].sum()):,} active")

    # standardize the natural inhibitors
    nat = []
    for name, (smi, um) in NATURAL_FAAH.items():
        r = standardize(smi)
        if not r.ok:
            print(f"  WARNING: {name} failed standardization ({r.reason})")
            continue
        nat.append({"name": name, "inchikey": r.inchikey,
                    "std_smiles": r.smiles, "reported_uM": um})
    nat = pd.DataFrame(nat)
    print(f"natural FAAH inhibitors: {len(nat)}")

    # none of them should already be in the ChEMBL set
    overlap = set(nat["inchikey"]) & set(faah["inchikey"])
    print(f"overlap with ChEMBL FAAH data: {len(overlap)} "
          f"(expected 0 - confirms the blind spot)")

    # decoy pool, shared across folds
    raw_pool = pd.read_csv(RAW / "decoy_pool.csv")
    pool = prepare_pool(raw_pool, set(agg["inchikey"]) | set(nat["inchikey"]))
    print(f"decoy pool: {len(pool):,}")

    actives = faah[faah["label_10um"] == 1]
    active_scaffolds = {scaffold(s) for s in actives["std_smiles"]}
    active_scaffolds.discard("")
    dec = match_decoys(actives, pool, n_per_active=DECOYS_PER_ACTIVE,
                       seed=RANDOM_STATE, exclude_scaffolds=active_scaffolds)
    dec = dec[~dec["inchikey"].isin(set(faah["inchikey"]))]
    print(f"decoys: {len(dec):,}")

    base = pd.concat([
        pd.DataFrame({"inchikey": faah["inchikey"],
                      "std_smiles": faah["std_smiles"],
                      "y": faah["label_10um"].astype(int)}),
        pd.DataFrame({"inchikey": dec["inchikey"],
                      "std_smiles": dec["std_smiles"], "y": 0}),
    ], ignore_index=True)
    print(f"base training set: {len(base):,} ({base['y'].mean():.1%} active)\n")

    # ---- baseline: ChEMBL only ----------------------------------------
    Xb, mb = featurize_smiles(base["std_smiles"].tolist())
    yb = base["y"].values[mb]
    Xb = Xb[mb]
    clf = model().fit(Xb, yb)
    Xn, mn = featurize_smiles(nat["std_smiles"].tolist())
    baseline = clf.predict_proba(Xn)[:, 1]

    # ---- leave-one-out with the other naturals added ------------------
    # The base feature matrix is identical in every fold, so featurize it once
    # and only stack the handful of natural-product rows that change. Rebuilding
    # 7,378 fingerprints per fold made this experiment ~10x slower than it
    # needed to be.
    loo = []
    for i in range(len(nat)):
        row = nat.iloc[i]
        other_idx = [j for j in range(len(nat)) if j != i]
        Xa = np.vstack([Xb, Xn[other_idx]])
        ya = np.concatenate([yb, np.ones(len(other_idx), dtype=int)])
        m = model().fit(Xa, ya)
        p = float(m.predict_proba(Xn[[i]])[:, 1][0])
        loo.append(p)
        print(f"  {row['name']:20s} baseline {baseline[i]:.6f} -> "
              f"LOO {p:.4f}   ({row['reported_uM']} uM reported)", flush=True)

    nat["p_baseline"] = baseline
    nat["p_leave_one_out"] = loo
    nat["fold_change"] = nat["p_leave_one_out"] / np.maximum(nat["p_baseline"], 1e-9)

    out = RESULTS / "faah_rescue_experiment.csv"
    nat.to_csv(out, index=False)

    print(f"\n{'':20s} {'baseline':>10s} {'LOO':>10s}")
    print(f"{'median':20s} {np.median(baseline):10.6f} {np.median(loo):10.4f}")
    print(f"{'n above 0.5':20s} {int((baseline>=0.5).sum()):10d} "
          f"{int((np.array(loo)>=0.5).sum()):10d}")

    summary = {
        "n_natural_inhibitors": int(len(nat)),
        "overlap_with_chembl": int(len(overlap)),
        "baseline_median_p": round(float(np.median(baseline)), 8),
        "loo_median_p": round(float(np.median(loo)), 6),
        "baseline_n_above_0.5": int((baseline >= 0.5).sum()),
        "loo_n_above_0.5": int((np.array(loo) >= 0.5).sum()),
        "interpretation": (
            "Each compound is scored by a model that never saw it. The "
            "baseline column is a model trained on ChEMBL FAAH data alone; the "
            "LOO column adds only the other five natural inhibitors. Any large "
            "gain is attributable to those few curated literature actives, and "
            "shows the blind spot is a data-coverage problem rather than a "
            "limitation of the representation."
        ),
    }
    (RESULTS / "faah_rescue_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
