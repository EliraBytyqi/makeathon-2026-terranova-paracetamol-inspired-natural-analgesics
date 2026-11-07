"""The two fixes together: potency regression on ChEMBL + NPASS.

Two independent changes each improved ranking of the held-out literature
actives. They address different problems, so they should compose:

  regression on pActivity   recovers the 4,315 compounds the classification
                            threshold discarded (the 1-10 uM band, which is
                            exactly where plant compounds sit) and removes the
                            threshold decision entirely. Best variant measured
                            55.1 -> 77.6 median percentile.

  NPASS augmentation        adds 312 real natural-product measurements, mostly
                            weak ones - hard negatives in the chemical space
                            the screen actually operates in. Measured
                            55.1 -> 67.0, p = 0.027.

Note the regression result only holds WITHOUT decoys: pinning tens of
thousands of molecules at a fabricated pActivity of 4.0 creates a point mass
the regressor learns to sit on, compressing the range it needs. Decoys help
classification (class balance) and hurt regression. So the configuration below
uses no assumed labels at all - every training row is a real measurement.

    python scripts/experiment_combined.py
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy import stats as sps
from sklearn.model_selection import GroupShuffleSplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from featurize import featurize_smiles, scaffold_groups   # noqa: E402

PROC = ROOT / "data" / "processed"
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"
MODELS = ROOT / "models" / "experiments"
MODELS.mkdir(parents=True, exist_ok=True)

SEED = 42
BACKGROUND_N = 30_000
MODEL_CUT_UM = {"TRPV1": 1.0, "CB1": 1.0, "FAAH": 10.0}


def regressor() -> LGBMRegressor:
    return LGBMRegressor(
        n_estimators=700, learning_rate=0.05, num_leaves=63,
        min_child_samples=10, subsample=0.8, subsample_freq=1,
        colsample_bytree=0.5, random_state=SEED, n_jobs=-1, verbose=-1,
    )


def main() -> None:
    merged = pd.read_csv(INTERIM / "bioactivity_with_npass.csv")
    base = pd.read_csv(PROC / "bioactivity_aggregated.csv")
    n_npass = merged["inchikey"].nunique() - base["inchikey"].nunique()
    print(f"training table: {len(merged):,} rows "
          f"({merged['inchikey'].nunique():,} compounds, ~{n_npass} from NPASS)")
    print("every row is a real measurement - no decoys, no assumed labels\n")

    lit = pd.read_csv(RESULTS / "literature_validation.csv")
    held = lit[lit["in_chembl_training"] == False].copy()
    held["fair"] = held["value_uM"] <= held["target"].map(MODEL_CUT_UM)
    # belt and braces: the merge already excluded these, verify it held
    leaked = set(merged["inchikey"]) & set(lit["inchikey"].dropna())
    leaked_npass = set(merged[merged.get("source", "") == "NPASS"]["inchikey"]) & set(
        held["inchikey"].dropna())
    print(f"validation compounds present via ChEMBL (expected, pre-existing): {len(leaked)}")
    print(f"validation compounds leaked in via NPASS (must be 0): {len(leaked_npass)}\n")
    if leaked_npass:
        raise SystemExit("LEAK: NPASS introduced validation compounds into training")

    print(f"sampling {BACKGROUND_N:,} background compounds...", flush=True)
    bg = pd.read_csv(RESULTS / "screening_all.csv", usecols=["std_smiles"],
                     low_memory=False).sample(n=BACKGROUND_N, random_state=SEED)
    Xbg, mbg = featurize_smiles(bg["std_smiles"].tolist())
    Xbg = Xbg[mbg]
    Xl, ml = featurize_smiles(held["std_smiles"].tolist())

    rows, per_target_rho = [], {}
    for target in ("TRPV1", "CB1", "FAAH"):
        sub = merged[merged["target_name"] == target].dropna(subset=["pactivity"])
        sub = (sub.sort_values("n_measurements", ascending=False)
                  .drop_duplicates(subset="inchikey", keep="first"))
        X, mask = featurize_smiles(sub["std_smiles"].tolist())
        y = sub["pactivity"].astype(float).to_numpy()[mask]
        X = X[mask]
        kept = [s for s, m in zip(sub["std_smiles"].tolist(), mask) if m]
        groups = scaffold_groups(kept)

        gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
        tr, te = next(gss.split(X, y, groups=groups))
        m = regressor().fit(X[tr], y[tr])
        rho = float(sps.spearmanr(y[te], m.predict(X[te])).statistic)
        per_target_rho[target] = round(rho, 3)

        final = regressor().fit(X, y)
        with open(MODELS / f"combined_{target}.pkl", "wb") as fh:
            pickle.dump({"model": final, "target": target, "n": int(len(y)),
                         "spearman_scaffold": round(rho, 4),
                         "recipe": "pActivity regression, ChEMBL+NPASS, no decoys"}, fh)
        print(f"  {target:6s} n={len(y):6,}  scaffold Spearman rho={rho:.3f}")

        bg_pred = final.predict(Xbg)
        for i, (_, r) in enumerate(held.iterrows()):
            if not ml[i] or r["target"] != target:
                continue
            p = float(final.predict(Xl[i:i + 1])[0])
            rows.append({
                "compound_name": r["compound_name"], "target": target,
                "value_uM": r["value_uM"], "fair": bool(r["fair"]),
                "pred_pactivity": round(p, 3),
                "pct_combined": round(float((bg_pred < p).mean() * 100), 1),
                "pct_baseline": r["library_percentile"],
            })

    res = pd.DataFrame(rows)
    res.to_csv(RESULTS / "combined_experiment.csv", index=False)

    fair = res[res["fair"]]
    a = fair["pct_combined"].to_numpy(dtype=float)
    b = fair["pct_baseline"].to_numpy(dtype=float)
    w = sps.wilcoxon(a - 50, alternative="greater")
    paired = sps.wilcoxon(a - b, alternative="greater")
    summary = {
        "recipe": "pActivity regression on ChEMBL + NPASS, all species, no decoys",
        "scaffold_spearman": per_target_rho,
        "n_fair": int(len(a)),
        "median_baseline": round(float(np.median(b)), 1),
        "median_combined": round(float(np.median(a)), 1),
        "gain": round(float(np.median(a) - np.median(b)), 1),
        "wilcoxon_vs_chance_p": round(float(w.pvalue), 4),
        "paired_vs_baseline_p": round(float(paired.pvalue), 4),
        "top10pct_baseline": int((b >= 90).sum()),
        "top10pct_combined": int((a >= 90).sum()),
        "top25pct_combined": int((a >= 75).sum()),
        "improved": int((a > b).sum()),
    }
    (RESULTS / "combined_experiment.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print(f"\n{'='*78}\nHELD-OUT LITERATURE ACTIVES\n{'='*78}")
    print(res.sort_values("pct_combined", ascending=False)[
        ["compound_name", "target", "value_uM", "fair",
         "pct_baseline", "pct_combined"]].to_string(index=False))
    print(f"\n{'='*78}\nSUMMARY\n{'='*78}")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
