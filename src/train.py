"""Train per-target activity models and report honest generalisation metrics.

Evaluation protocol
-------------------
Every model is evaluated twice:

  random split    - the number most cheminformatics demos report
  scaffold split  - the number that predicts what happens when you screen a
                    natural-product library

The gap between them is reported explicitly. It is usually large, and showing
it is more informative than any single headline figure. A model that scores
0.95 random / 0.72 scaffold is not a 0.95 model; it has memorised congeneric
series, and natural products share almost no scaffolds with ChEMBL's
medicinal-chemistry space.

Metrics are PR-AUC-led because the datasets are imbalanced and we care about
the top of a ranked list, not global accuracy. Enrichment factor at 1% and 5%
is reported because that is literally the screening use case: if we hand a
chemist the top 1% of a library, how many more actives do they get than by
picking at random?

Applicability domain: for every prediction we record the maximum Tanimoto
similarity to the training set. Predictions on molecules far from anything the
model has seen are flagged, not silently trusted.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit, StratifiedKFold

from decoys import match_decoys, prepare_pool, property_overlap_report
import modelio
from featurize import featurize_smiles, scaffold, scaffold_groups

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
MODELS = ROOT / "models"
RESULTS = ROOT / "results"
MODELS.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

RANDOM_STATE = 42
MIN_CLASS = 25          # below this a target/scope combination is not modellable
TEST_FRACTION = 0.2
DECOYS_PER_ACTIVE = 1   # measured, not assumed - see note below
DECOY_WEIGHT = 0.4      # a presumed inactive is not worth a measured one

# Both numbers used to be assertions. Now they are measurements.
#
# The ratio was 3, justified in a comment as bringing a ~75% active rate down
# to roughly 20%. scripts/experiment_decoy_ratio.py tested 0, 1, 2, 3 and 5
# against the external PubChem set, where decoys cannot flatter the score:
#
#     ratio      TRPV1    CB1    FAAH
#       0        0.944   0.841   0.619      <- publication bias, unusable for CB1
#       1        0.952   0.919   0.626
#       3        0.954   0.918   0.612      <- the old setting
#       5        0.955   0.920   0.589
#
# Going from zero decoys to one is worth +0.078 on CB1. Everything past one is
# noise, and the dip at ratio 2 across all three targets confirms it. Three
# decoys per active therefore carried roughly 17,900 extra presumed-inactive
# molecules for no measurable gain, so the ratio is 1.
#
# The weight is the second half of the same problem. Every decoy used to enter
# the loss with label 0 and full weight, identical to a compound somebody had
# actually assayed and found inactive - the "presumed_inactive" tag was
# recorded in the frame and never read by anything. That asserts a certainty
# the data does not support, and it is the standard positive-unlabeled mistake.
# Decoys now carry weight 0.4, so the model is still taught what a non-ligand
# looks like but is penalised less for being unsure about a molecule nobody
# measured. Measured compounds, active or inactive, keep weight 1.0.

_morgan = GetMorganGenerator(radius=2, fpSize=2048)


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------
def enrichment_factor(y_true, y_score, frac: float) -> float:
    """How many times more actives in the top `frac` than random selection."""
    n = len(y_true)
    k = max(1, int(round(n * frac)))
    order = np.argsort(-np.asarray(y_score))
    top = np.asarray(y_true)[order[:k]]
    hit_rate_top = top.mean()
    hit_rate_all = np.asarray(y_true).mean()
    if hit_rate_all == 0:
        return float("nan")
    return float(hit_rate_top / hit_rate_all)


def evaluate(y_true, y_score, threshold: float = 0.5) -> dict:
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    y_pred = (y_score >= threshold).astype(int)
    out = {
        "n": int(len(y_true)),
        "n_active": int(y_true.sum()),
        "active_rate": round(float(y_true.mean()), 4),
        "pr_auc": round(float(average_precision_score(y_true, y_score)), 4),
        "roc_auc": round(float(roc_auc_score(y_true, y_score)), 4),
        "mcc": round(float(matthews_corrcoef(y_true, y_pred)), 4),
        "balanced_acc": round(float(balanced_accuracy_score(y_true, y_pred)), 4),
        "ef_1pct": round(enrichment_factor(y_true, y_score, 0.01), 2),
        "ef_5pct": round(enrichment_factor(y_true, y_score, 0.05), 2),
    }
    return out


# ---------------------------------------------------------------------------
# applicability domain
# ---------------------------------------------------------------------------
def build_ad_fps(smiles_list) -> np.ndarray:
    """Training-set fingerprints as a dense bit matrix.

    Stored as a numpy array rather than a list of RDKit ExplicitBitVects so
    that applicability-domain similarity can be computed as a matrix product
    (below) instead of element-by-element.
    """
    rows = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m is not None:
            rows.append(_morgan.GetFingerprintAsNumPy(m))
    if not rows:
        return np.zeros((0, 2048), dtype=np.uint8)
    # uint8, not float32: these matrices are pickled with every model, and at
    # ~20k training compounds the float32 version is 160MB per model. The
    # values are bits; they are cast to float only inside the similarity
    # matmul, one batch at a time.
    return np.asarray(rows, dtype=np.uint8)


def query_fingerprints(smiles_list) -> tuple[np.ndarray, np.ndarray]:
    """Fingerprint a query set once, returning (matrix, index-of-valid).

    Screening calls max_similarity once per model. Without this, each call
    re-parses and re-fingerprints all 145k library molecules - four times the
    RDKit work for identical results. Computing the query matrix once and
    reusing it across models is the single biggest saving in the screen.
    """
    rows, idx = [], []
    for j, s in enumerate(smiles_list):
        m = Chem.MolFromSmiles(s) if isinstance(s, str) and s else None
        if m is not None:
            rows.append(_morgan.GetFingerprintAsNumPy(m))
            idx.append(j)
    if not rows:
        return np.zeros((0, 2048), dtype=np.float32), np.asarray(idx, dtype=int)
    return np.asarray(rows, dtype=np.float32), np.asarray(idx, dtype=int)


def max_similarity_from_fps(Q: np.ndarray, q_idx: np.ndarray, n_total: int,
                            train_fps, batch: int = 4096) -> np.ndarray:
    """max_similarity for an already-fingerprinted query set."""
    sims = np.zeros(n_total, dtype=np.float32)
    T = np.asarray(train_fps, dtype=np.float32)
    if T.size == 0 or Q.size == 0:
        return sims
    t_pop = T.sum(axis=1)
    Tt = T.T
    for start in range(0, len(Q), batch):
        qb = Q[start:start + batch]
        inter = qb @ Tt
        union = qb.sum(axis=1)[:, None] + t_pop[None, :] - inter
        with np.errstate(divide="ignore", invalid="ignore"):
            tan = np.where(union > 0, inter / union, 0.0)
        sims[q_idx[start:start + batch]] = tan.max(axis=1)
    return sims


def max_similarity(smiles_list, train_fps, batch: int = 2048) -> np.ndarray:
    """Max Tanimoto similarity of each query to the training set.

    Screening 150k natural products against three models, each with thousands
    of training compounds, is billions of pairwise comparisons. Done one
    molecule at a time it dominates the runtime of the whole pipeline. Tanimoto
    on binary fingerprints is

        T = |A & B| / (|A| + |B| - |A & B|)

    and |A & B| for all pairs is exactly a matrix product of the bit matrices,
    which BLAS does orders of magnitude faster than a Python loop.
    """
    n = len(smiles_list)
    sims = np.zeros(n, dtype=np.float32)
    T = np.asarray(train_fps, dtype=np.float32)  # cast from stored uint8
    if T.size == 0:
        return sims
    t_pop = T.sum(axis=1)

    for start in range(0, n, batch):
        chunk = smiles_list[start:start + batch]
        rows, idx = [], []
        for j, s in enumerate(chunk):
            m = Chem.MolFromSmiles(s) if isinstance(s, str) and s else None
            if m is not None:
                rows.append(_morgan.GetFingerprintAsNumPy(m))
                idx.append(j)
        if not rows:
            continue
        Q = np.asarray(rows, dtype=np.float32)
        inter = Q @ T.T                      # |A & B| for every pair
        union = Q.sum(axis=1)[:, None] + t_pop[None, :] - inter
        with np.errstate(divide="ignore", invalid="ignore"):
            tan = np.where(union > 0, inter / union, 0.0)
        sims[start + np.asarray(idx)] = tan.max(axis=1)
    return sims


# ---------------------------------------------------------------------------
# training
# ---------------------------------------------------------------------------
def make_models():
    return {
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            min_samples_leaf=2,
            max_features="sqrt",
            class_weight="balanced_subsample",
            n_jobs=-1,
            random_state=RANDOM_STATE,
        ),
        "lightgbm": LGBMClassifier(
            n_estimators=600,
            learning_rate=0.05,
            num_leaves=63,
            min_child_samples=10,
            subsample=0.8,
            subsample_freq=1,
            colsample_bytree=0.5,
            class_weight="balanced",
            random_state=RANDOM_STATE,
            n_jobs=-1,
            verbose=-1,
        ),
    }


def split_indices(groups, y, kind: str):
    """Return (train_idx, test_idx) for a scaffold-grouped or random split."""
    rng = np.random.RandomState(RANDOM_STATE)
    n = len(y)
    if kind == "scaffold":
        gss = GroupShuffleSplit(
            n_splits=1, test_size=TEST_FRACTION, random_state=RANDOM_STATE
        )
        tr, te = next(gss.split(np.zeros(n), y, groups=groups))
        return tr, te
    # random split, stratified on label
    skf = StratifiedKFold(
        n_splits=int(round(1 / TEST_FRACTION)), shuffle=True,
        random_state=RANDOM_STATE
    )
    tr, te = next(skf.split(np.zeros(n), y))
    return tr, te


def train_one(target: str, scope: str, df: pd.DataFrame, label_col: str,
              decoy_pool: pd.DataFrame | None = None):
    sub = df.dropna(subset=[label_col]).copy()
    if scope == "human":
        sub = sub[sub["organism_group"] == "human"]

    # one row per compound: if a compound appears for several organisms in the
    # 'all' scope, keep the human record if present, else the most-measured one
    sub = (
        sub.sort_values(
            ["organism_group", "n_measurements"],
            key=lambda s: s.map({"human": 0, "rat": 1, "mouse": 2})
            if s.name == "organism_group" else -s,
        )
        .drop_duplicates(subset="inchikey", keep="first")
    )

    sub["evidence"] = "measured"
    y_measured = sub[label_col].astype(int).values
    if min((y_measured == 1).sum(), (y_measured == 0).sum()) < MIN_CLASS:
        print(f"  SKIP {target}/{scope}/{label_col}: "
              f"{(y_measured==1).sum()} active, {(y_measured==0).sum()} inactive")
        return None

    print(f"\n  {target}/{scope}/{label_col}: {len(sub):,} measured compounds, "
          f"{(y_measured==1).sum():,} active ({y_measured.mean():.1%})")

    decoy_report = None
    if decoy_pool is not None and len(decoy_pool):
        actives = sub[sub[label_col] == 1]
        active_scaffolds = {scaffold(s) for s in actives["std_smiles"]}
        active_scaffolds.discard("")
        dec = match_decoys(
            actives, decoy_pool, n_per_active=DECOYS_PER_ACTIVE,
            seed=RANDOM_STATE, exclude_scaffolds=active_scaffolds,
        )
        # never let a decoy duplicate a measured compound
        dec = dec[~dec["inchikey"].isin(set(sub["inchikey"]))]
        decoy_report = property_overlap_report(actives, dec)
        print(f"    + {len(dec):,} property-matched presumed-inactive decoys "
              f"(MW median gap {decoy_report['mw_median_gap']}, "
              f"cLogP gap {decoy_report['clogp_median_gap']})")
        dec_rows = pd.DataFrame({
            "inchikey": dec["inchikey"],
            "std_smiles": dec["std_smiles"],
            label_col: 0,
            "evidence": "presumed_inactive",
        })
        sub = pd.concat([sub, dec_rows], ignore_index=True)

    y = sub[label_col].astype(int).values
    smiles = sub["std_smiles"].tolist()
    n_presumed = int((sub["evidence"] == "presumed_inactive").sum())
    print(f"    final set: {len(sub):,} compounds, {(y==1).sum():,} active "
          f"({y.mean():.1%}); {n_presumed:,} negatives are presumed, "
          f"carried at weight {DECOY_WEIGHT}")

    X, mask = featurize_smiles(smiles)
    X, y = X[mask], y[mask]
    smiles = [s for s, m in zip(smiles, mask) if m]
    sub = sub[mask]
    # measured evidence counts fully; a decoy is an assumption, not a datum
    w = np.where(sub["evidence"].to_numpy() == "presumed_inactive",
                 DECOY_WEIGHT, 1.0).astype(float)
    groups = scaffold_groups(smiles)
    n_scaffolds = len(set(groups))
    print(f"    {n_scaffolds:,} distinct Bemis-Murcko scaffolds "
          f"({n_scaffolds/len(y):.1%} of compounds)")

    report = {
        "target": target,
        "scope": scope,
        "label_col": label_col,
        "n_compounds": int(len(y)),
        "n_active": int(y.sum()),
        "active_rate": round(float(y.mean()), 4),
        "n_measured": int(len(y_measured)),
        "measured_active_rate": round(float(y_measured.mean()), 4),
        "n_scaffolds": int(n_scaffolds),
        "decoy_matching": decoy_report,
        "splits": {},
    }

    best = None
    for split_kind in ("random", "scaffold"):
        tr, te = split_indices(groups, y, split_kind)
        if len(set(y[te])) < 2:
            continue
        report["splits"][split_kind] = {}
        for mname, model in make_models().items():
            model.fit(X[tr], y[tr], sample_weight=w[tr])
            score = model.predict_proba(X[te])[:, 1]
            m = evaluate(y[te], score)
            report["splits"][split_kind][mname] = m
            print(f"    {split_kind:8s} {mname:14s} "
                  f"PR-AUC {m['pr_auc']:.3f}  ROC {m['roc_auc']:.3f}  "
                  f"MCC {m['mcc']:.3f}  EF1% {m['ef_1pct']:.1f}")
            if split_kind == "scaffold":
                if best is None or m["pr_auc"] > best[2]:
                    best = (mname, model, m["pr_auc"])

    # honest gap statement
    try:
        r = report["splits"]["random"]
        s = report["splits"]["scaffold"]
        for mname in r:
            gap = r[mname]["pr_auc"] - s[mname]["pr_auc"]
            report["splits"]["scaffold"][mname]["pr_auc_gap_vs_random"] = round(gap, 4)
    except KeyError:
        pass

    # refit the winning model on everything for deployment
    if best is None:
        return report
    mname = best[0]
    final = make_models()[mname]
    final.fit(X, y, sample_weight=w)
    report["selected_model"] = mname
    report["decoys_per_active"] = DECOYS_PER_ACTIVE
    report["decoy_sample_weight"] = DECOY_WEIGHT
    report["n_presumed_inactive"] = n_presumed

    bundle = {
        "model": final,
        "model_name": mname,
        "target": target,
        "scope": scope,
        "label_col": label_col,
        "train_smiles": smiles,
        "train_inchikeys": sub["inchikey"].tolist(),
        "ad_fps": build_ad_fps(smiles),
        "report": report,
    }
    path = MODELS / f"{target}_{scope}_{label_col}.pkl"
    modelio.save(bundle, path)
    print(f"    -> saved {path.name} (model: {mname})")
    return report


def main():
    agg = pd.read_csv(PROC / "bioactivity_aggregated.csv")
    print(f"loaded {len(agg):,} aggregated bioactivity rows")

    decoy_pool = None
    pool_path = RAW / "decoy_pool.csv"
    # Standardizing 60k decoys costs ~15 minutes, and it does not change
    # between runs, so cache it rather than paying for it on every retrain.
    cache_path = ROOT / "data" / "interim" / "decoy_pool_standardized.csv"
    if cache_path.exists():
        decoy_pool = pd.read_csv(cache_path)
        decoy_pool["cell"] = [tuple(eval(c)) if isinstance(c, str) else c
                              for c in decoy_pool["cell"]]
        print(f"decoy pool: {len(decoy_pool):,} standardized (cached)")
    elif pool_path.exists():
        raw_pool = pd.read_csv(pool_path)
        print(f"decoy pool: {len(raw_pool):,} raw molecules -> standardizing",
              flush=True)
        # exclude every compound measured anywhere in the panel
        panel_keys = set(agg["inchikey"])
        decoy_pool = prepare_pool(raw_pool, panel_keys)
        decoy_pool.to_csv(cache_path, index=False)
        print(f"usable decoys after exclusion: {len(decoy_pool):,}")
    else:
        print("WARNING: no decoy pool found. Training on measured data only, "
              "which has a heavily publication-biased active rate.")

    reports = []
    for target in sorted(agg["target_name"].unique()):
        df = agg[agg["target_name"] == target]
        for scope in ("human", "all"):
            for label_col in ("label", "label_10um"):
                # the 10 uM variant is only interesting where the 1 uM cut is
                # too strict for the chemistry of interest, i.e. FAAH
                if label_col == "label_10um" and target != "FAAH":
                    continue
                rep = train_one(target, scope, df, label_col, decoy_pool)
                if rep:
                    reports.append(rep)

    out = RESULTS / "training_report.json"
    out.write_text(json.dumps(reports, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")

    # flat summary table
    rows = []
    for r in reports:
        for split, models in r["splits"].items():
            for mname, m in models.items():
                rows.append({
                    "target": r["target"], "scope": r["scope"],
                    "label": r["label_col"], "split": split, "model": mname,
                    "n": m["n"], "active_rate": m["active_rate"],
                    "pr_auc": m["pr_auc"], "roc_auc": m["roc_auc"],
                    "mcc": m["mcc"], "ef_1pct": m["ef_1pct"],
                })
    tab = pd.DataFrame(rows)
    tab.to_csv(RESULTS / "model_metrics.csv", index=False)
    print(f"wrote {RESULTS/'model_metrics.csv'}")
    print("\n=== SCAFFOLD-SPLIT PERFORMANCE (the number that matters) ===")
    sc = tab[tab["split"] == "scaffold"].sort_values("pr_auc", ascending=False)
    print(sc.to_string(index=False))


if __name__ == "__main__":
    main()
