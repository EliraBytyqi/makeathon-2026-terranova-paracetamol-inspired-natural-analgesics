"""Screen the plant compound library against the trained target models.

Pipeline
--------
  plant compounds
    -> remove anything present in the training data (leakage control)
    -> score against each target model
    -> applicability-domain similarity to each model's training set
    -> reactive-metabolite alerts + PAINS + CNS-likeness
    -> composite rank

The composite score is deliberately transparent and its weights are printed,
because an opaque ranking formula is unfalsifiable:

    composite = target_score
                * exp(-ALERT_PENALTY * alert_severity_sum)
                * ad_confidence
                * cns_bonus

Leakage control is not a formality. A natural-product hit list that silently
contains the model's own training compounds will rank capsaicin first on
TRPV1 and look spectacular while demonstrating nothing. We remove exact
InChIKey matches and report how many were removed, then re-add the known
compounds separately as positive controls so we can check the model recovers
them.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

import modelio
from featurize import featurize_smiles
from train import max_similarity_from_fps, query_fingerprints

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
MODELS = ROOT / "models"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ALERT_PENALTY = 0.35     # per unit of summed alert severity
AD_FLOOR = 0.25          # below this Tanimoto, predictions are extrapolation
CNS_BONUS = 1.15         # central analgesia is the hypothesis, so reward CNS-likeness
DILI_WEIGHT = 0.35       # raised from 0.30, but less than first set - see note

# This weight used to be 0.30, and the reason was specific: the old safety
# model was trained on TDC's 475-drug DILI benchmark, which contains neither
# paracetamol nor isoniazid, and it scored paracetamol at 0.37 - "low risk"
# for the index compound of this entire project. A model that misses that has
# not earned a large say in a ranking.
#
# That model is gone. The hepatic panel is trained on 12,444 compounds across
# FDA DILIrank 2.0 and two Tox21 HepG2 assays, and paracetamol is in it as
# most-DILI-concern. Internal scaffold ROC-AUC is 0.79 (clinical), 0.86
# (mitochondrial) and 0.77 (cytotoxicity), which first suggested a weight
# around 0.45.
#
# scripts/validate_hepatic.py then measured the panel properly, and the
# external numbers are much lower than the internal ones:
#
#   clinical head, trained on DILIrank v1 and tested on the drugs the FDA
#   added in v2:                    ROC-AUC 0.68, and 0.65 on new scaffolds
#                                   with a CI that reaches 0.50
#   mito head predicting clinical DILI (it has never seen a clinical label):
#                                   ROC-AUC 0.66
#   hepg2 head, same test:          ROC-AUC 0.57, not significant
#
# Resubstitution on the clinical head is 0.998, which is what memorisation
# looks like and why the paracetamol score is not evidence of anything.
#
# So: more than 0.30, because the panel is 27x larger and does contain the
# index compound. Well short of 0.45, because its true external
# discrimination is around 0.65. The structural alerts keep the hard veto -
# they encode mechanism rather than correlation, and they are the only layer
# that can see metabolism-dependent toxicity at all.


def load_models(scope: str = "human"):
    """Load the per-target activity models.

    dili.pkl lives in the same directory but is a safety model with a different
    bundle shape, so it is loaded separately where it is needed, not here.
    """
    bundles = {}
    for path in sorted(MODELS.glob("*.pkl")):
        if path.name == "dili.pkl":
            continue
        b = modelio.load(path)
        if "target" not in b:
            continue
        key = (b["target"], b["scope"], b["label_col"])
        bundles[key] = b
    return bundles


def pick_model_per_target(bundles):
    """One model per target: human scope preferred, 1 uM label preferred,
    except FAAH where the 10 uM cut matches the natural-product potency range."""
    chosen = {}
    targets = sorted({k[0] for k in bundles})
    for t in targets:
        pref = [
            (t, "human", "label_10um" if t == "FAAH" else "label"),
            (t, "all", "label_10um" if t == "FAAH" else "label"),
            (t, "human", "label"),
            (t, "all", "label"),
        ]
        for key in pref:
            if key in bundles:
                chosen[t] = bundles[key]
                break
    return chosen


def add_composite_columns(df: pd.DataFrame, target_cols: list[str],
                          ad_cols: list[str]) -> pd.DataFrame:
    """Attach the best-target pick, applicability domain and composite score.

    Module level rather than nested inside main() so it can be exercised by
    a small subset before committing to a full 150k-compound run.
    """
    df = df.copy()
    df["best_target"] = (df[target_cols].idxmax(axis=1)
                         .str.replace("p_", "", regex=False))
    df["target_score"] = df[target_cols].max(axis=1)
    # AD of the model that produced the best score. Vectorised: iterrows over a
    # 150k-row frame is needlessly slow here.
    ad_matrix = df[ad_cols].to_numpy()
    best_idx = np.argmax(df[target_cols].to_numpy(), axis=1)
    df["ad_best"] = ad_matrix[np.arange(len(df)), best_idx]
    df["ad_confidence"] = np.clip(df["ad_best"] / 0.5, 0.0, 1.0)
    df["in_domain"] = (df["ad_best"] >= AD_FLOOR).astype(int)
    df["safety_multiplier"] = np.exp(
        -ALERT_PENALTY * df["alert_severity_sum"].fillna(0))
    df["cns_multiplier"] = np.where(df["cns_like"] == 1, CNS_BONUS, 1.0)
    # p_hepatic is the max over the three hepatic heads; p_dili is the
    # clinical head alone, kept as a fallback for results produced before the
    # panel existed.
    hep = ("p_hepatic" if "p_hepatic" in df.columns
           else "p_dili" if "p_dili" in df.columns else None)
    if hep:
        df["dili_multiplier"] = 1.0 - DILI_WEIGHT * df[hep].fillna(0.5)
    else:
        df["dili_multiplier"] = 1.0
    df["composite_score"] = (
        df["target_score"]
        * df["safety_multiplier"]
        * df["ad_confidence"]
        * df["cns_multiplier"]
        * df["dili_multiplier"]
    )
    return df


def main():
    np_lib = pd.read_csv(PROC / "np_library.csv", low_memory=False)
    print(f"library: {len(np_lib):,} unique plant compounds")

    bundles = load_models()
    if not bundles:
        raise SystemExit("no trained models found - run train.py first")
    chosen = pick_model_per_target(bundles)
    print("models in use:")
    for t, b in chosen.items():
        r = b["report"]
        sc = r["splits"].get("scaffold", {}).get(b["model_name"], {})
        print(f"   {t:6s} {b['model_name']:14s} scope={b['scope']:5s} "
              f"label={b['label_col']:10s} scaffold PR-AUC="
              f"{sc.get('pr_auc','n/a')}")

    # ---- leakage control -------------------------------------------------
    train_keys = set()
    for b in bundles.values():
        train_keys.update(b["train_inchikeys"])
    in_train = np_lib["std_inchikey"].isin(train_keys)
    print(f"\nleakage control: {int(in_train.sum()):,} library compounds appear in "
          f"ChEMBL training data and are held out as positive controls")
    controls = np_lib[in_train].copy()
    screen = np_lib[~in_train].copy()
    print(f"novel compounds to screen: {len(screen):,}")

    # ---- featurize once, score against each model ------------------------
    def score_frame(frame: pd.DataFrame) -> pd.DataFrame:
        smiles = frame["std_smiles"].tolist()
        print(f"   featurizing {len(smiles):,} compounds...", flush=True)
        X, mask = featurize_smiles(smiles)
        out = frame.copy()
        out["featurized"] = mask

        # Fingerprint the query set ONCE and reuse it for every model's
        # applicability-domain calculation. Calling max_similarity per model
        # re-fingerprints all 145k molecules each time - four passes of
        # identical RDKit work.
        print("   building query fingerprints for applicability domain...",
              flush=True)
        Q, q_idx = query_fingerprints([s if m else "" for s, m in zip(smiles, mask)])

        for t, b in chosen.items():
            probs = np.full(len(frame), np.nan, dtype=np.float32)
            if mask.any():
                probs[mask] = b["model"].predict_proba(X[mask])[:, 1]
            out[f"p_{t}"] = probs
            out[f"ad_{t}"] = max_similarity_from_fps(
                Q, q_idx, len(frame), b["ad_fps"])
            print(f"   scored {t}", flush=True)

        # Hepatic liability: three heads over DILIrank 2.0 and two Tox21 HepG2
        # assays. p_dili stays the clinical head so every existing column name
        # keeps its meaning; p_hepatic is the max over the three, which is what
        # the composite uses.
        panel_path = MODELS / "hepatic_panel.pkl"
        dili_path = MODELS / "dili.pkl"
        if panel_path.exists():
            panel = modelio.load(panel_path)
            head_cols = []
            for hname, hb in panel["heads"].items():
                col = "p_dili" if hname == "clinical" else f"p_{hname}"
                p = np.full(len(frame), np.nan, dtype=np.float32)
                if mask.any():
                    p[mask] = hb["model"].predict_proba(X[mask])[:, 1]
                out[col] = p
                head_cols.append(col)
                if hname == "clinical":
                    out["ad_dili"] = max_similarity_from_fps(
                        Q, q_idx, len(frame), hb["ad_fps"])
                print(f"   scored hepatic:{hname}", flush=True)
            out["p_hepatic"] = out[head_cols].max(axis=1)
        elif dili_path.exists():
            db = modelio.load(dili_path)
            dp = np.full(len(frame), np.nan, dtype=np.float32)
            if mask.any():
                dp[mask] = db["model"].predict_proba(X[mask])[:, 1]
            out["p_dili"] = dp
            out["p_hepatic"] = dp
            out["ad_dili"] = max_similarity_from_fps(
                Q, q_idx, len(frame), db["ad_fps"])
            print("   scored DILI (single head; no panel found)", flush=True)
        return out[out["featurized"]].drop(columns=["featurized"])

    print("\nscoring novel compounds...")
    screen = score_frame(screen)
    print("scoring positive controls...")
    controls = score_frame(controls)

    target_cols = [f"p_{t}" for t in chosen]
    ad_cols = [f"ad_{t}" for t in chosen]

    screen = add_composite_columns(screen, target_cols, ad_cols)
    controls = add_composite_columns(controls, target_cols, ad_cols)

    screen = screen.sort_values("composite_score", ascending=False)
    controls = controls.sort_values("target_score", ascending=False)

    screen.to_csv(RESULTS / "screening_all.csv.gz", index=False)
    controls.to_csv(RESULTS / "positive_controls.csv", index=False)

    # ---- shortlist -------------------------------------------------------
    shortlist = screen[
        (screen["in_domain"] == 1)
        & (screen["max_alert_severity"].fillna(0) < 3)
        & (screen["n_pains"].fillna(0) == 0)
    ].head(50)
    shortlist.to_csv(RESULTS / "shortlist.csv", index=False)

    cols = ["compound_name", "best_target", "target_score", "ad_best",
            "n_reactive_alerts", "p_hepatic", "cns_like", "composite_score",
            "source_genera", "npclass_class", "tier"]
    cols = [c for c in cols if c in screen.columns]
    print("\n=== TOP 25 SHORTLIST (in-domain, no severity-3 alert, no PAINS) ===")
    disp = shortlist[cols].head(25).copy()
    disp["compound_name"] = disp["compound_name"].astype(str).str.slice(0, 34)
    disp["source_genera"] = disp["source_genera"].astype(str).str.slice(0, 24)
    disp["npclass_class"] = disp["npclass_class"].astype(str).str.slice(0, 22)
    print(disp.to_string(index=False))

    print("\n=== POSITIVE CONTROLS (compounds the model saw in training) ===")
    print(controls[cols].head(15).to_string(index=False))

    stats = {
        "library_size": int(len(np_lib)),
        "removed_as_training_overlap": int(in_train.sum()),
        "screened": int(len(screen)),
        "in_domain": int(screen["in_domain"].sum()),
        "in_domain_rate": round(float(screen["in_domain"].mean()), 4),
        "severity3_flagged": int((screen["max_alert_severity"] == 3).sum()),
        "pains_flagged": int((screen["n_pains"] > 0).sum()),
        "shortlist_size": int(len(shortlist)),
        "weights": {
            "alert_penalty": ALERT_PENALTY,
            "ad_floor": AD_FLOOR,
            "cns_bonus": CNS_BONUS,
            "dili_weight": DILI_WEIGHT,
        },
    }
    (RESULTS / "screening_stats.json").write_text(
        json.dumps(stats, indent=2), encoding="utf-8")
    print("\n=== SCREENING STATS ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
