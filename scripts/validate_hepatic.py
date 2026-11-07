"""External validation for the hepatic panel.

The activity models have a real external test: PubChem BioAssay outcomes on
structures never trained on. Until now the hepatic panel had nothing
comparable - only internal scaffold splits. Its one famous number, paracetamol
at 0.71, is not a test at all, because paracetamol is in DILIrank and
therefore in the training data. Recalling a training label proves nothing.

Two genuine external tests are available without any new data.

Test A - temporal holdout.
    DILIrank 2.0 added 300 drugs that version 1 did not contain, and the
    spreadsheet flags them. Train the clinical head on version 1 only, then
    test on the newly added drugs. The FDA assigned those classes after
    version 1 was published, so a v1-trained model has not seen them: this is
    as close to a prospective test as a retrospective project can get.
    141 of the new drugs are unambiguous (21 most-concern, 120 no-concern),
    a 14.9% positive rate - much closer to a realistic base rate than the
    39% the model trains at.

Test B - cross-dataset transfer.
    The mito and hepg2 heads are trained only on Tox21 assay outcomes and have
    never seen a clinical DILI label. So asking them to rank DILIrank
    most-concern drugs above no-concern drugs is a real external test, and it
    asks the question that actually matters: does mitochondrial disruption
    measured in a dish predict liver injury in people? Restricted to DILIrank
    drugs whose skeleton is absent from the relevant Tox21 training set, so no
    structure is scored by a model that was fitted on it.

Both tests are small. Both are reported with confidence intervals, and a wide
interval is left wide rather than described as a trend.

    python scripts/validate_hepatic.py
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import modelio  # noqa: E402
from featurize import featurize_smiles, scaffold_groups   # noqa: E402
from train import make_models                             # noqa: E402

DATA = ROOT / "data" / "processed" / "hepatic_endpoints.csv"
XLSX = ROOT / "data" / "external" / "dilirank" / "dilirank.xlsx"
MODELS = ROOT / "models"
RESULTS = ROOT / "results"
SEED = 42
N_BOOT = 2000


def stats_block(y: np.ndarray, s: np.ndarray) -> dict:
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    if len(set(y)) < 2:
        return {"n": int(len(y)), "n_positive": int(y.sum()),
                "scoreable": False, "reason": "single-class"}
    auc = float(roc_auc_score(y, s))
    rng = np.random.RandomState(SEED)
    boots = []
    for _ in range(N_BOOT):
        i = rng.randint(0, len(y), len(y))
        if len(set(y[i])) > 1:
            boots.append(roc_auc_score(y[i], s[i]))
    lo, hi = (float(np.percentile(boots, 2.5)),
              float(np.percentile(boots, 97.5))) if boots else (np.nan, np.nan)
    return {"n": int(len(y)), "n_positive": int(y.sum()),
            "n_negative": int((1 - y).sum()),
            "positive_rate": round(float(y.mean()), 4), "scoreable": True,
            "roc_auc": round(auc, 4), "ci95": [round(lo, 4), round(hi, 4)],
            "mannwhitney_p": float(sps.mannwhitneyu(
                s[y == 1], s[y == 0], alternative="greater").pvalue),
            "median_positive": round(float(np.median(s[y == 1])), 4),
            "median_negative": round(float(np.median(s[y == 0])), 4)}


def show(title: str, d: dict) -> None:
    if not d.get("scoreable"):
        print(f"  {title:34s} n={d['n']:4d}  NOT SCOREABLE ({d['reason']})")
        return
    print(f"  {title:34s} n={d['n']:4d} ({d['n_positive']} pos / "
          f"{d['n_negative']} neg, {d['positive_rate']:.1%})  "
          f"ROC-AUC {d['roc_auc']:.3f} "
          f"[{d['ci95'][0]:.3f}, {d['ci95'][1]:.3f}]  "
          f"p={d['mannwhitney_p']:.2e}")


# ---------------------------------------------------------------------------
def main() -> None:
    if not DATA.exists():
        raise SystemExit(f"missing {DATA}; run src/build_hepatic_dataset.py")
    df = pd.read_csv(DATA)
    clin = df[df["endpoint"] == "clinical"].copy()

    v2 = pd.read_excel(XLSX, sheet_name="version 2", header=1)
    v2.columns = [str(c).strip() for c in v2.columns]
    v2["name"] = v2["CompoundName"].astype(str).str.strip()
    flag = v2.set_index("name")["Comment"].astype(str).str.strip().str.lower()
    clin["is_new"] = clin["name"].map(flag).eq("new")
    clin = clin[clin["concern"].isin(["vmost-dili-concern", "vno-dili-concern"])].copy()
    clin["y"] = (clin["concern"] == "vmost-dili-concern").astype(int)
    print(f"clinical, unambiguous: {len(clin)}  "
          f"({int(clin['y'].sum())} most / {int((1 - clin['y']).sum())} none)")
    print(f"  of which newly added in DILIrank 2.0: {int(clin['is_new'].sum())}")

    summary = {}

    # ------------------------------------------------------- Test A: temporal
    print(f"\n{'=' * 78}\nTEST A - temporal holdout: train on DILIrank v1, "
          f"test on drugs added in v2\n{'=' * 78}")
    tr = clin[~clin["is_new"]].reset_index(drop=True)
    te = clin[clin["is_new"]].reset_index(drop=True)
    print(f"  train {len(tr)} ({int(tr['y'].sum())} most) | "
          f"test {len(te)} ({int(te['y'].sum())} most)")

    Xtr, mtr = featurize_smiles(tr["std_smiles"].tolist())
    Xte, mte = featurize_smiles(te["std_smiles"].tolist())
    ytr = tr["y"].to_numpy()[mtr]
    yte = te["y"].to_numpy()[mte]

    # how much of the test set is genuinely new chemistry, not just a new row
    g_tr = set(scaffold_groups([s for s, m in zip(tr["std_smiles"], mtr) if m]))
    te_smiles = [s for s, m in zip(te["std_smiles"], mte) if m]
    g_te = scaffold_groups(te_smiles)
    novel = np.array([g not in g_tr for g in g_te])
    print(f"  test compounds whose scaffold is absent from training: "
          f"{int(novel.sum())}/{len(novel)}")

    blocks = {}
    for mname, template in make_models().items():
        m = template.__class__(**template.get_params()).fit(Xtr[mtr], ytr)
        s = m.predict_proba(Xte[mte])[:, 1]
        blocks[mname] = {"all": stats_block(yte, s),
                         "novel_scaffold": stats_block(yte[novel], s[novel])}
        print(f"\n  {mname}")
        show("all newly added drugs", blocks[mname]["all"])
        show("new scaffolds only", blocks[mname]["novel_scaffold"])
    summary["test_a_temporal"] = {
        "design": ("clinical head trained on DILIrank v1 drugs, tested on the "
                   "300 drugs added in v2 (141 unambiguous)"),
        "n_train": int(len(ytr)), "n_test": int(len(yte)),
        "n_test_novel_scaffold": int(novel.sum()),
        "results": blocks}

    # -------------------------------------------- Test B: cross-dataset
    print(f"\n{'=' * 78}\nTEST B - cross-dataset: do the Tox21 heads predict "
          f"clinical DILI?\n{'=' * 78}")
    panel_path = MODELS / "hepatic_panel.pkl"
    if not panel_path.exists():
        raise SystemExit(f"missing {panel_path}; run src/train_hepatic.py")
    panel = modelio.load(panel_path)

    Xc, mc = featurize_smiles(clin["std_smiles"].tolist())
    sub = clin[mc].reset_index(drop=True)
    yc = sub["y"].to_numpy()
    sub_skel = sub["inchikey"].str[:14]

    btest = {}
    for hname in ("mito", "hepg2"):
        if hname not in panel["heads"]:
            continue
        hb = panel["heads"][hname]
        # exclude any clinical drug this head was fitted on
        trained = set(df[(df["endpoint"] == hname)]["inchikey"].str[:14])
        keep = (~sub_skel.isin(trained)).to_numpy()
        s = hb["model"].predict_proba(Xc[mc])[:, 1]
        print(f"\n  {hname} head (never saw a clinical DILI label)")
        btest[hname] = {
            "all_clinical": stats_block(yc, s),
            "excluding_tox21_trained": stats_block(yc[keep], s[keep]),
        }
        show("all DILIrank most vs none", btest[hname]["all_clinical"])
        show("structures not in its training", btest[hname]["excluding_tox21_trained"])

    # the shipped clinical head, for reference - NOT a valid external test,
    # printed so nobody mistakes its internal number for one
    ch = panel["heads"].get("clinical")
    if ch is not None:
        s = ch["model"].predict_proba(Xc[mc])[:, 1]
        d = stats_block(yc, s)
        print(f"\n  clinical head on its own training data "
              f"(NOT a test - resubstitution):")
        show("resubstitution only", d)
        btest["clinical_resubstitution_not_a_test"] = d

    summary["test_b_cross_dataset"] = {
        "design": ("Tox21-trained heads scored against DILIrank most-concern "
                   "vs no-concern; these heads never saw a clinical label"),
        "results": btest}

    (RESULTS / "hepatic_validation.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nwrote {RESULTS / 'hepatic_validation.json'}")


if __name__ == "__main__":
    main()
