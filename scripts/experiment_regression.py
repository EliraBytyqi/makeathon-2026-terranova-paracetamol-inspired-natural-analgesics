"""Can regression on potency beat thresholded classification for ranking?

The deployed models are binary classifiers built on a potency threshold, and
that threshold costs us twice:

  1. It DISCARDS the ambiguous band. 4,315 compounds sit between the active and
     inactive cuts and are dropped - 926 for TRPV1, 2,289 for CB1, 1,100 for
     FAAH. That band is 1-10 uM, which is exactly where the plant compounds we
     are hunting sit. We are throwing away the most relevant training data we
     have.

  2. It forces a binary decision on a continuous quantity. Screening does not
     need a decision, it needs a RANKING - and a regressor on pActivity ranks
     directly, with no threshold to pick wrong.

This trains LightGBM regressors on pActivity and compares them with the
deployed classifiers on the one thing that matters: where the held-out
literature actives land against the library background.

Two variants, because regression needs something to anchor the low end:
  measured_only - trained on measured potencies alone
  with_decoys   - presumed-inactive decoys pinned at pActivity 4.0 (100 uM)

    python scripts/experiment_regression.py
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn.model_selection import GroupShuffleSplit
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from featurize import featurize_smiles, scaffold_groups   # noqa: E402

PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
MODELS = ROOT / "models" / "experiments"
MODELS.mkdir(parents=True, exist_ok=True)

SEED = 42
BACKGROUND_N = 30_000
DECOY_PACTIVITY = 4.0        # 100 uM: presumed weak, not merely "not active"
MODEL_CUT_UM = {"TRPV1": 1.0, "CB1": 1.0, "FAAH": 10.0}


def regressor() -> LGBMRegressor:
    return LGBMRegressor(
        n_estimators=700, learning_rate=0.05, num_leaves=63,
        min_child_samples=10, subsample=0.8, subsample_freq=1,
        colsample_bytree=0.5, random_state=SEED, n_jobs=-1, verbose=-1,
    )


def load_decoy_pool() -> pd.DataFrame | None:
    cache = ROOT / "data" / "interim" / "decoy_pool_standardized.csv"
    if not cache.exists():
        return None
    pool = pd.read_csv(cache)
    pool["cell"] = [tuple(eval(c)) if isinstance(c, str) else c for c in pool["cell"]]
    return pool


def main() -> None:
    agg = pd.read_csv(PROC / "bioactivity_aggregated.csv")
    pool = load_decoy_pool()
    if pool is None:
        raise SystemExit("standardized decoy pool missing")

    # ---- the common background, featurized once ------------------------
    print(f"sampling {BACKGROUND_N:,} library compounds as background...", flush=True)
    bg = pd.read_csv(RESULTS / "screening_all.csv", usecols=["std_smiles"],
                     low_memory=False).sample(n=BACKGROUND_N, random_state=SEED)
    Xbg, mbg = featurize_smiles(bg["std_smiles"].tolist())
    Xbg = Xbg[mbg]
    print(f"  featurized {Xbg.shape[0]:,}")

    # ---- the held-out literature compounds -----------------------------
    lit = pd.read_csv(RESULTS / "literature_validation.csv")
    held = lit[lit["in_chembl_training"] == False].copy()
    held["fair"] = held["value_uM"] <= held["target"].map(MODEL_CUT_UM)
    Xl, ml = featurize_smiles(held["std_smiles"].tolist())

    out_rows, summaries = [], {}

    for scope in ("human", "all"):
        for variant in ("measured_only", "with_decoys"):
            tag = f"{scope}/{variant}"
            print(f"\n{'='*70}\n{tag}\n{'='*70}", flush=True)
            per_target = {}

            for target in ("TRPV1", "CB1", "FAAH"):
                sub = agg[agg["target_name"] == target]
                if scope == "human":
                    sub = sub[sub["organism_group"] == "human"]
                sub = (sub.sort_values("n_measurements", ascending=False)
                          .drop_duplicates(subset="inchikey", keep="first"))
                sub = sub.dropna(subset=["pactivity"])

                smiles = sub["std_smiles"].tolist()
                y = sub["pactivity"].astype(float).to_numpy()

                if variant == "with_decoys":
                    # anchor the weak end; sample 2 decoys per measured compound
                    n_dec = min(len(sub) * 2, len(pool))
                    dec = pool.sample(n=n_dec, random_state=SEED)
                    dec = dec[~dec["inchikey"].isin(set(sub["inchikey"]))]
                    smiles = smiles + dec["std_smiles"].tolist()
                    y = np.concatenate([y, np.full(len(dec), DECOY_PACTIVITY)])

                X, mask = featurize_smiles(smiles)
                X, y = X[mask], y[mask]
                kept = [s for s, m in zip(smiles, mask) if m]
                groups = scaffold_groups(kept)

                # scaffold-split check: can it RANK held-out scaffolds?
                gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
                tr, te = next(gss.split(X, y, groups=groups))
                m = regressor().fit(X[tr], y[tr])
                pred = m.predict(X[te])
                rho = float(sps.spearmanr(y[te], pred).statistic)

                final = regressor().fit(X, y)
                with open(MODELS / f"reg_{target}_{scope}_{variant}.pkl", "wb") as fh:
                    pickle.dump({"model": final, "target": target, "scope": scope,
                                 "variant": variant, "n": int(len(y)),
                                 "spearman_scaffold": round(rho, 4)}, fh)

                per_target[target] = final
                print(f"  {target:6s} n={len(y):6,}  scaffold Spearman rho={rho:.3f}")

            # ---- score background and the literature set ----------------
            bg_pred = {t: per_target[t].predict(Xbg) for t in per_target}
            for i, (_, r) in enumerate(held.iterrows()):
                if not ml[i]:
                    continue
                t = r["target"]
                p = float(per_target[t].predict(Xl[i:i + 1])[0])
                pct = float((bg_pred[t] < p).mean() * 100)
                out_rows.append({
                    "scope": scope, "variant": variant,
                    "compound_name": r["compound_name"], "target": t,
                    "value_uM": r["value_uM"], "fair": bool(r["fair"]),
                    "pred_pactivity": round(p, 3),
                    "percentile": round(pct, 1),
                    "classifier_percentile": r["library_percentile"],
                })

            fair_rows = [x for x in out_rows
                         if x["scope"] == scope and x["variant"] == variant and x["fair"]]
            pcts = np.array([x["percentile"] for x in fair_rows])
            base = np.array([x["classifier_percentile"] for x in fair_rows])
            w = sps.wilcoxon(pcts - 50, alternative="greater")
            summaries[tag] = {
                "n_fair": len(pcts),
                "median_percentile": round(float(np.median(pcts)), 1),
                "classifier_median": round(float(np.median(base)), 1),
                "gain": round(float(np.median(pcts) - np.median(base)), 1),
                "wilcoxon_p": round(float(w.pvalue), 4),
                "in_top10pct": int((pcts >= 90).sum()),
                "in_top25pct": int((pcts >= 75).sum()),
            }
            s = summaries[tag]
            print(f"  -> fair compounds n={s['n_fair']}  median percentile "
                  f"{s['classifier_median']} (classifier) -> {s['median_percentile']} "
                  f"(regression)  p={s['wilcoxon_p']}  top10%={s['in_top10pct']}")

    df = pd.DataFrame(out_rows)
    df.to_csv(RESULTS / "regression_experiment.csv", index=False)
    (RESULTS / "regression_experiment.json").write_text(
        json.dumps(summaries, indent=2), encoding="utf-8")

    print(f"\n{'='*70}\nSUMMARY - median percentile of held-out literature actives\n{'='*70}")
    print(f"{'variant':26s} {'n':>3s} {'classifier':>11s} {'regression':>11s} "
          f"{'gain':>6s} {'p':>7s} {'top10%':>7s}")
    for tag, s in summaries.items():
        print(f"{tag:26s} {s['n_fair']:3d} {s['classifier_median']:11.1f} "
              f"{s['median_percentile']:11.1f} {s['gain']:+6.1f} {s['wilcoxon_p']:7.3f} "
              f"{s['in_top10pct']:7d}")
    print(f"\nwrote {RESULTS/'regression_experiment.csv'}")


if __name__ == "__main__":
    main()
