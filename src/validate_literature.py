"""External validation against the curated literature set.

This is the project's real test. The literature set contains natural-product
modulators of TRPV1, CB1 and FAAH whose activities are reported in the primary
literature, plus synthetic reference ligands as positive controls. Compounds
that also appear in ChEMBL training data are identified and reported
separately - only the genuinely held-out ones constitute external validation.

The question is not "what is the AUC". With ~24 compounds that number would be
meaningless. The question is: does the model rank compounds that are known to
engage a target above the bulk of the plant library? That is exactly the
screening decision the pipeline exists to make, and it is answerable here by
reporting each compound's percentile against the full library background.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from featurize import featurize_smiles
from safety import physchem, profile
from screen import load_models, pick_model_per_target
from train import max_similarity
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)


def main():
    lit = pd.read_csv(INTERIM / "literature_actives_resolved.csv")
    lit = lit[lit["inchikey"].notna()].copy()
    # One row per (compound, target): several compounds have multiple reported
    # values that disagree across studies (kaempferol's FAAH Ki of 5 uM vs a
    # review's IC50 of 1.06 uM, for instance). Keep the highest-confidence
    # record, ordering explicitly rather than relying on alphabetical sort.
    rank = {"high": 0, "medium": 1, "low": 2}
    lit["_conf_rank"] = lit["confidence"].map(rank).fillna(9)
    lit = (lit.sort_values("_conf_rank")
              .drop_duplicates(subset=["compound_name", "target"], keep="first")
              .drop(columns=["_conf_rank"]))
    print(f"literature validation set: {len(lit)} compound-target claims, "
          f"{lit['compound_name'].nunique()} unique compounds")

    bundles = load_models()
    chosen = pick_model_per_target(bundles)
    print(f"models: {list(chosen)}")

    train_keys = set()
    for b in bundles.values():
        train_keys.update(b["train_inchikeys"])

    smiles = lit["std_smiles"].tolist()
    X, mask = featurize_smiles(smiles)

    for t, b in chosen.items():
        probs = np.full(len(lit), np.nan)
        probs[mask] = b["model"].predict_proba(X[mask])[:, 1]
        lit[f"p_{t}"] = probs
        lit[f"ad_{t}"] = max_similarity(
            [s if m else "" for s, m in zip(smiles, mask)], b["ad_fps"])

    lit["in_chembl_training"] = lit["inchikey"].isin(train_keys)

    # score on the claimed target
    lit["p_claimed_target"] = [
        r.get(f"p_{r['target']}", np.nan) for _, r in lit.iterrows()
    ]
    lit["ad_claimed_target"] = [
        r.get(f"ad_{r['target']}", np.nan) for _, r in lit.iterrows()
    ]

    # percentile against the full plant library background
    np_path = RESULTS / "screening_all.csv"
    if np_path.exists():
        bg = pd.read_csv(np_path, low_memory=False)
        pct = []
        for _, r in lit.iterrows():
            col = f"p_{r['target']}"
            if col in bg.columns and pd.notna(r["p_claimed_target"]):
                pct.append(round(float(
                    (bg[col] < r["p_claimed_target"]).mean() * 100), 1))
            else:
                pct.append(np.nan)
        lit["library_percentile"] = pct
    else:
        lit["library_percentile"] = np.nan

    # safety profile of each literature compound
    srows = []
    for smi in lit["std_smiles"]:
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        sp = profile(mol)
        pc = physchem(mol)
        srows.append({
            "n_reactive_alerts": sp.n_alerts,
            "max_alert_severity": sp.max_severity,
            "reactive_alerts": "|".join(sp.alerts),
            "n_pains": len(sp.pains),
            "cns_like": pc.get("cns_like", 0),
        })
    lit = pd.concat([lit.reset_index(drop=True),
                     pd.DataFrame(srows)], axis=1)

    # DILI
    dili_path = ROOT / "models" / "dili.pkl"
    if dili_path.exists():
        with open(dili_path, "rb") as fh:
            db = pickle.load(fh)
        dp = np.full(len(lit), np.nan)
        dp[mask] = db["model"].predict_proba(X[mask])[:, 1]
        lit["p_dili"] = dp

    out = RESULTS / "literature_validation.csv"
    lit.to_csv(out, index=False)

    cols = ["compound_name", "target", "interaction", "value_uM",
            "p_claimed_target", "library_percentile", "ad_claimed_target",
            "in_chembl_training", "n_reactive_alerts", "p_dili"]
    cols = [c for c in cols if c in lit.columns]

    held_out = lit[~lit["in_chembl_training"]]
    seen = lit[lit["in_chembl_training"]]

    print(f"\n=== HELD OUT (not in ChEMBL training) - {len(held_out)} claims ===")
    print(held_out[cols].sort_values(
        "p_claimed_target", ascending=False).to_string(index=False))

    print(f"\n=== SEEN IN TRAINING - {len(seen)} claims (sanity check only) ===")
    print(seen[cols].sort_values(
        "p_claimed_target", ascending=False).to_string(index=False))

    summary = {
        "n_claims": int(len(lit)),
        "n_held_out": int(len(held_out)),
        "n_seen_in_training": int(len(seen)),
    }
    if len(held_out):
        hp = held_out["p_claimed_target"].dropna()
        summary["held_out_median_probability"] = round(float(hp.median()), 4)
        summary["held_out_frac_above_0.5"] = round(float((hp >= 0.5).mean()), 4)
        lp = held_out["library_percentile"].dropna()
        if len(lp):
            summary["held_out_median_library_percentile"] = round(float(lp.median()), 1)
            summary["held_out_frac_top10pct_of_library"] = round(
                float((lp >= 90).mean()), 4)
    if len(seen):
        sp_ = seen["p_claimed_target"].dropna()
        summary["seen_median_probability"] = round(float(sp_.median()), 4)

    (RESULTS / "literature_validation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print("\n=== SUMMARY ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
