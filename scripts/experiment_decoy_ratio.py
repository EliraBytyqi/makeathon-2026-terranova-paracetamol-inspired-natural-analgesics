"""How many presumed-inactive decoys improve external ranking?

The original production setting used three property-matched presumed-inactive
decoys per active, justified by an assumed class-balance target. This
experiment compared ratios 0, 1, 2, 3 and 5 against external PubChem assay
outcomes. The current production setting is one decoy per active, with sample
weight 0.4; three is a historical comparison, not the current setting.

The candidate pool contains 60,000 ChEMBL molecules sampled across the ID
range. They are presumed inactive because no target activity was found, not
because every molecule was experimentally confirmed inactive.

Ratios tested: 0 (measured data only), 1, 2, 3 (historical setting), 5.

Each is judged on two numbers that disagree on purpose:

  internal   scaffold-split ROC-AUC, where most negatives ARE decoys. Adding
             decoys inflates this almost by construction, so on its own it
             cannot answer the question.
  external   the deduplicated PubChem set - real measured actives against real
             measured inactives, none of them trained on. Decoys cannot
             flatter this number, so this is the one that decides.

Production models are never touched; everything here trains in memory.

    python scripts/experiment_decoy_ratio.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from decoys import match_decoys                      # noqa: E402
from featurize import featurize_smiles, scaffold_groups, scaffold  # noqa: E402
from train import make_models                        # noqa: E402

PROC = ROOT / "data" / "processed"
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"
SEED = 42
RATIOS = [0, 1, 2, 3, 5]
TARGETS = {"TRPV1": "label", "CB1": "label", "FAAH": "label_10um"}


def build(target: str, label_col: str, agg: pd.DataFrame,
          pool: pd.DataFrame, ratio: int) -> pd.DataFrame:
    """Reproduce train_one's training frame at a given decoy ratio."""
    sub = agg[(agg["target_name"] == target)
              & (agg["organism_group"] == "human")].dropna(subset=[label_col]).copy()
    sub = sub.drop_duplicates(subset="inchikey", keep="first")
    sub = sub[["inchikey", "std_smiles", label_col]].copy()
    sub["evidence"] = "measured"
    if ratio == 0:
        return sub.reset_index(drop=True)

    actives = sub[sub[label_col] == 1]
    act_scaf = {scaffold(s) for s in actives["std_smiles"]}
    act_scaf.discard("")
    dec = match_decoys(actives, pool, n_per_active=ratio, seed=SEED,
                       exclude_scaffolds=act_scaf)
    dec = dec[~dec["inchikey"].isin(set(sub["inchikey"]))]
    rows = pd.DataFrame({"inchikey": dec["inchikey"],
                         "std_smiles": dec["std_smiles"],
                         label_col: 0, "evidence": "presumed_inactive"})
    return pd.concat([sub, rows], ignore_index=True)


def main() -> None:
    agg = pd.read_csv(PROC / "bioactivity_aggregated.csv")
    pool = pd.read_csv(INTERIM / "decoy_pool_standardized.csv")
    pool["cell"] = [tuple(eval(c)) if isinstance(c, str) else c
                    for c in pool["cell"]]

    ext = pd.read_csv(RESULTS / "pubchem_validation_expanded.csv")
    print(f"external test set: {len(ext)} structures "
          f"({ext['target'].value_counts().to_dict()})\n")

    out: dict = {}
    for target, label_col in TARGETS.items():
        e = ext[ext["target"] == target]
        Xe, me = featurize_smiles(e["std_smiles"].tolist())
        ye = (e["outcome"].to_numpy()[me] == "Active").astype(int)
        Xe = Xe[me]

        print(f"{'=' * 72}\n{target}  (external test: {len(ye)} structures, "
              f"{int(ye.sum())} active / {int((1 - ye).sum())} inactive)\n{'=' * 72}")
        print(f"  {'ratio':>5s} {'train n':>8s} {'active%':>8s} "
              f"{'internal ROC':>13s} {'EXTERNAL ROC':>13s} {'ext PR-AUC':>11s}")
        out[target] = {}

        for ratio in RATIOS:
            df = build(target, label_col, agg, pool, ratio)
            X, mask = featurize_smiles(df["std_smiles"].tolist())
            X = X[mask]
            y = df[label_col].to_numpy()[mask].astype(int)
            smiles = [s for s, m in zip(df["std_smiles"], mask) if m]
            if len(set(y)) < 2:
                continue

            # internal: scaffold split, same protocol as production
            groups = scaffold_groups(smiles)
            sp = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
            tr, te = next(sp.split(X, y, groups=groups))
            m = make_models()["lightgbm"]
            internal = (roc_auc_score(y[te], m.fit(X[tr], y[tr])
                                      .predict_proba(X[te])[:, 1])
                        if len(set(y[te])) > 1 else float("nan"))

            # external: train on everything, score the held-out measured set
            full = make_models()["lightgbm"].fit(X, y)
            s = full.predict_proba(Xe)[:, 1]
            ext_roc = float(roc_auc_score(ye, s))
            ext_pr = float(average_precision_score(ye, s))

            out[target][str(ratio)] = {
                "n_train": int(len(y)), "active_rate": round(float(y.mean()), 4),
                "internal_scaffold_roc_auc": round(float(internal), 4),
                "external_roc_auc": round(ext_roc, 4),
                "external_pr_auc": round(ext_pr, 4),
            }
            tag = "  <- production" if ratio == 3 else ""
            print(f"  {ratio:5d} {len(y):8,d} {y.mean():7.1%} "
                  f"{internal:13.3f} {ext_roc:13.3f} {ext_pr:11.3f}{tag}")

        best = max(out[target], key=lambda r: out[target][r]["external_roc_auc"])
        prod = out[target].get("3", {}).get("external_roc_auc", float("nan"))
        print(f"\n  best external: ratio {best} "
              f"({out[target][best]['external_roc_auc']:.3f}) "
              f"vs production ratio 3 ({prod:.3f})\n")

    (RESULTS / "decoy_ratio_experiment.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {RESULTS / 'decoy_ratio_experiment.json'}")


if __name__ == "__main__":
    main()
