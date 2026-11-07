"""Does the activity threshold explain the external-validation result?

FAAH carries a 10 uM label variant because the plant inhibitors that motivate
this project sit at 1-15 uM, and a 1 uM cut would train the model to reject
them. TRPV1 and CB1 never got that treatment - but the plant TRPV1 modulators
in our literature set sit at 2-40 uM, so they are "inactive" by those models'
own definition. Scoring them low is correct behaviour, not a miss.

This trains the missing 10 uM variants and re-scores the held-out literature
compounds under both cuts, against the SAME library background sample, so the
two are directly comparable.

It writes new model files (TRPV1_*_label_10um.pkl, CB1_*_label_10um.pkl) and
results/threshold_experiment.{csv,json}. It overwrites nothing that exists.

    python scripts/experiment_threshold.py
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from featurize import featurize_smiles          # noqa: E402
from train import train_one, max_similarity     # noqa: E402

PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
MODELS = ROOT / "models"
BACKGROUND_N = 30_000
SEED = 42


def load_decoys(agg):
    cache = ROOT / "data" / "interim" / "decoy_pool_standardized.csv"
    if not cache.exists():
        # Without decoys the comparison is worthless: the 1 uM models were
        # trained against presumed-inactives at a ~20% active rate, and a
        # 10 uM model trained on raw measured data sits at 86%. That is a
        # decoy comparison wearing a threshold comparison's clothes.
        raw = ROOT / "data" / "raw" / "decoy_pool.csv"
        if not raw.exists():
            raise SystemExit("no decoy pool on disk - refusing to run an "
                             "unmatched comparison")
        from decoys import prepare_pool
        print("standardizing the decoy pool once (slow, then cached)...",
              flush=True)
        pool = prepare_pool(pd.read_csv(raw), set(agg["inchikey"]))
        pool.to_csv(cache, index=False)
        print(f"  cached {len(pool):,} usable decoys")
        return pool
    pool = pd.read_csv(cache)
    pool["cell"] = [tuple(eval(c)) if isinstance(c, str) else c for c in pool["cell"]]
    print(f"decoy pool: {len(pool):,} standardized (cached)")
    return pool


def main():
    agg = pd.read_csv(PROC / "bioactivity_aggregated.csv")
    decoys = load_decoys(agg)

    # ---- train the two missing 10 uM variants -------------------------
    trained = {}
    for target in ("TRPV1", "CB1"):
        df = agg[agg["target_name"] == target]
        print(f"\n{'='*66}\ntraining {target} at the 10 uM cut\n{'='*66}", flush=True)
        rep = train_one(target, "human", df, "label_10um", decoys)
        if rep is None:
            print(f"  {target}: not trainable, skipping")
            continue
        path = MODELS / f"{target}_human_label_10um.pkl"
        with open(path, "rb") as fh:
            trained[target] = pickle.load(fh)

    if not trained:
        raise SystemExit("nothing trained")

    # ---- the old 1 uM models, for comparison --------------------------
    old = {}
    for target in trained:
        with open(MODELS / f"{target}_human_label.pkl", "rb") as fh:
            old[target] = pickle.load(fh)

    # ---- a common library background ----------------------------------
    print(f"\nsampling {BACKGROUND_N:,} library compounds as background...", flush=True)
    bg = pd.read_csv(RESULTS / "screening_all.csv",
                     usecols=["std_smiles"], low_memory=False)
    bg = bg.sample(n=min(BACKGROUND_N, len(bg)), random_state=SEED)
    Xbg, mbg = featurize_smiles(bg["std_smiles"].tolist())
    Xbg = Xbg[mbg]
    print(f"  featurized {Xbg.shape[0]:,}")

    bg_scores = {}
    for target in trained:
        bg_scores[(target, "10um")] = trained[target]["model"].predict_proba(Xbg)[:, 1]
        bg_scores[(target, "1um")] = old[target]["model"].predict_proba(Xbg)[:, 1]
        print(f"  scored background for {target}")

    # ---- the held-out literature compounds ----------------------------
    lit = pd.read_csv(RESULTS / "literature_validation.csv")
    held = lit[(lit["in_chembl_training"] == False)
               & (lit["target"].isin(trained))].copy()
    Xl, ml = featurize_smiles(held["std_smiles"].tolist())

    rows = []
    for i, (_, r) in enumerate(held.iterrows()):
        if not ml[i]:
            continue
        t = r["target"]
        x = Xl[i:i + 1]
        out = {
            "compound_name": r["compound_name"],
            "target": t,
            "value_uM": r["value_uM"],
            "interaction": r["interaction"],
        }
        for tag in ("1um", "10um"):
            model = (old if tag == "1um" else trained)[t]["model"]
            p = float(model.predict_proba(x)[0, 1])
            pct = float((bg_scores[(t, tag)] < p).mean() * 100)
            out[f"p_{tag}"] = p
            out[f"pct_{tag}"] = round(pct, 1)
        out["pct_gain"] = round(out["pct_10um"] - out["pct_1um"], 1)
        out["ad"] = round(float(max_similarity([r["std_smiles"]],
                                               trained[t]["ad_fps"])[0]), 3)
        rows.append(out)

    res = pd.DataFrame(rows).sort_values("pct_gain", ascending=False)
    res.to_csv(RESULTS / "threshold_experiment.csv", index=False)

    print(f"\n{'='*78}")
    print("HELD-OUT LITERATURE COMPOUNDS: 1 uM model vs 10 uM model")
    print(f"{'='*78}")
    show = res[["compound_name", "target", "value_uM", "p_1um", "pct_1um",
                "p_10um", "pct_10um", "pct_gain", "ad"]].copy()
    show["p_1um"] = show["p_1um"].map(lambda v: f"{v:.2e}")
    show["p_10um"] = show["p_10um"].map(lambda v: f"{v:.2e}")
    print(show.to_string(index=False))

    summary = {
        "background_n": int(Xbg.shape[0]),
        "n_compounds": int(len(res)),
        "median_pct_1um": round(float(res["pct_1um"].median()), 1),
        "median_pct_10um": round(float(res["pct_10um"].median()), 1),
        "median_gain": round(float(res["pct_gain"].median()), 1),
        "n_improved": int((res["pct_gain"] > 0).sum()),
        "top10pct_1um": int((res["pct_1um"] >= 90).sum()),
        "top10pct_10um": int((res["pct_10um"] >= 90).sum()),
        "above_half_1um": int((res["p_1um"] >= 0.5).sum()),
        "above_half_10um": int((res["p_10um"] >= 0.5).sum()),
    }
    (RESULTS / "threshold_experiment.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print(f"\n{'='*78}\nSUMMARY\n{'='*78}")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"\nwrote {RESULTS/'threshold_experiment.csv'}")


if __name__ == "__main__":
    main()
