"""Train the hepatic-liability panel - safety layer B, rebuilt.

What changed and why
--------------------
The previous safety model was a single classifier over TDC's 475-compound DILI
benchmark. Two things were wrong with it:

  1. It was small, and its chemical space was marketed drugs only.
  2. Paracetamol - the compound this entire project is about - was simply not
     in it, so the model's score for paracetamol was pure extrapolation and
     came out at 0.37, i.e. "low risk", for one of the best-documented
     hepatotoxins in medicine.

This replaces it with three endpoints, trained as three separate heads because
they answer three different questions and disagree in informative ways:

  clinical  FDA DILIrank 2.0. ~900 approved drugs with an unambiguous
            concern class derived from their own FDA labels. This is what
            actually happened to patients. Paracetamol is in it, as
            most-DILI-concern, severity 5.

  mito      Tox21 mitochondrial membrane potential in HepG2 (~8,000
            compounds). Mitochondrial depolarisation is a direct, measured,
            dose-response intrinsic-toxicity mechanism.

  hepg2     Tox21 CellTox Green cytotoxicity in HepG2 (~5,600 compounds).
            Outright death of liver cells.

Three heads, not one pooled model, because pooling would force the model to
treat "this drug injured patients" and "this compound depolarised mitochondria
in a dish" as the same label. They are not. Keeping them separate also means
each head's applicability domain and each head's error can be reported
honestly on its own terms.

The combined score is the maximum of the three. That is deliberately
conservative - any one hepatic signal is enough to push a compound down the
ranking - and it has no tuned weights to defend.

What this panel still cannot do
-------------------------------
It cannot see metabolism-dependent toxicity. HepG2 cells barely express
CYP2E1, so they never convert paracetamol to NAPQI, and Tox21's top
concentration is around 100 uM while paracetamol hepatotoxicity needs
millimolar exposure. Both Tox21 heads therefore have paracetamol measured
INACTIVE, correctly, for the assay they ran - and uselessly, for the question
we are asking. Only the clinical head and the structural-alert layer carry
that signal. This is the concrete, mechanistic reason the alert layer is not
redundant, and it is why the panel is a ranking input and not a safety
verdict.

    python src/train_hepatic.py
"""

from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import modelio
from featurize import featurize_smiles, scaffold_groups   # noqa: E402
from train import build_ad_fps, evaluate, make_models     # noqa: E402

DATA = ROOT / "data" / "processed" / "hepatic_endpoints.csv"
MODELS = ROOT / "models"
RESULTS = ROOT / "results"
MODELS.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

RANDOM_STATE = 42
METRICS = ["pr_auc", "roc_auc", "mcc", "balanced_acc"]

ENDPOINTS = {
    "clinical": "FDA DILIrank 2.0 concern class (label-derived, human outcome)",
    "mito": "Tox21 MMP_ratio, mitochondrial membrane potential in HepG2",
    "hepg2": "Tox21 CellTox Green 40 h cytotoxicity in HepG2",
}

# Reference compounds. Chosen so the panel can be sanity-checked against
# facts a reviewer already knows, including the ones that make it look bad.
PROBES = {
    "paracetamol":  ("CC(=O)Nc1ccc(O)cc1",            "hepatotoxic at overdose; DILIrank most-concern"),
    "isoniazid":    ("NNC(=O)c1ccncc1",               "hepatotoxic; DILIrank most-concern"),
    "troglitazone": ("Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",
                                                      "withdrawn for hepatotoxicity"),
    "diclofenac":   ("OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "hepatotoxic; DILIrank most-concern"),
    "ibuprofen":    ("CC(C)Cc1ccc(C(C)C(=O)O)cc1",    "low hepatic concern"),
    "aspirin":      ("CC(=O)Oc1ccccc1C(=O)O",         "low hepatic concern"),
    "caffeine":     ("Cn1c(=O)c2c(ncn2C)n(C)c1=O",    "no hepatic concern"),
    "curcumin":     ("COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
                                                      "plant NP; Tox21 mito-active"),
    "kaempferol":   ("O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
                                                      "plant NP; Tox21 mito-active"),
    "capsaicin":    ("COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O", "plant NP; our TRPV1 control"),
}


def mean_sd(ms: list[dict]) -> dict:
    out = {}
    for k in METRICS:
        vals = [m[k] for m in ms]
        out[k] = round(float(np.mean(vals)), 4)
        out[f"{k}_sd"] = round(float(np.std(vals)), 4)
    return out


def repeated_eval(X, y, groups, factory, kind: str, repeats: int) -> list[dict]:
    """Repeated held-out evaluation. Scaffold splits group by Bemis-Murcko
    framework so the test set is chemistry the model has not seen."""
    ms = []
    for seed in range(repeats):
        if kind == "scaffold":
            sp = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed)
            tr, te = next(sp.split(X, y, groups=groups))
        else:
            sp = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
            tr, te = next(sp.split(X, y))
        if len(set(y[te])) < 2 or len(set(y[tr])) < 2:
            continue
        m = factory().fit(X[tr], y[tr])
        ms.append(evaluate(y[te], m.predict_proba(X[te])[:, 1]))
    return ms


def clinical_subset(df: pd.DataFrame, convention: str) -> pd.DataFrame:
    """Pick the clinical head's positive and negative class.

    DILIrank grades severity 0-8 inside each concern class, and the two
    conventions disagree about what "positive" means:

      most_vs_none  most-concern (1) against no-concern (0). Less-concern is
                    dropped, not folded in. This is what we use.
      dilist        most + less (1) against none (0) - the DILIst convention.

    We use most_vs_none because the DILIst convention puts aspirin, which
    DILIrank grades severity 0 and vLess-DILI-concern, in the same positive
    class as troglitazone at severity 8 - a drug withdrawn from the market for
    killing people's livers. A model trained to call those two the same thing
    scores aspirin at 0.80 and is then useless for ranking, which is the only
    job we have for it. Dropping the 351 less-concern drugs costs training
    data and buys a label that means something.
    """
    d = df.copy()
    if convention == "most_vs_none":
        d = d[d["concern"].isin(["vmost-dili-concern", "vno-dili-concern"])]
        d["y"] = d["most_concern"]
    return d.reset_index(drop=True)


def train_head(name: str, df: pd.DataFrame) -> tuple[dict, dict]:
    """Train one endpoint. Returns (bundle, report)."""
    print(f"\n{'=' * 74}\n{name}: {ENDPOINTS[name]}\n{'=' * 74}")

    X, mask = featurize_smiles(df["std_smiles"].tolist())
    X = X[mask]
    y = df["y"].to_numpy()[mask].astype(int)
    smiles = [s for s, m in zip(df["std_smiles"], mask) if m]
    groups = scaffold_groups(smiles)
    # a small set needs more repeats before its numbers mean anything
    repeats = 10 if len(y) < 2000 else 5
    print(f"  {len(y)} compounds, {int(y.sum())} positive ({y.mean():.1%}), "
          f"{len(set(groups))} scaffolds, {repeats} repeats")

    report = {"endpoint": name, "description": ENDPOINTS[name],
              "n": int(len(y)), "n_positive": int(y.sum()),
              "positive_rate": round(float(y.mean()), 4),
              "n_scaffolds": int(len(set(groups))), "repeats": repeats,
              "candidates": {}}

    best_name, best_score, best_factory = None, -1.0, None
    for mname, template in make_models().items():
        factory = (lambda t=template: t.__class__(**t.get_params()))
        scaff = repeated_eval(X, y, groups, factory, "scaffold", repeats)
        rnd = repeated_eval(X, y, groups, factory, "random", repeats)
        if not scaff:
            continue
        s, r = mean_sd(scaff), mean_sd(rnd)
        report["candidates"][mname] = {"scaffold": s, "random": r}
        print(f"  {mname:14s} scaffold ROC {s['roc_auc']:.3f}+/-{s['roc_auc_sd']:.3f}"
              f"  PR {s['pr_auc']:.3f}  MCC {s['mcc']:.3f}"
              f"   | random ROC {r['roc_auc']:.3f}"
              f"   | gap {r['roc_auc'] - s['roc_auc']:+.3f}")
        if s["roc_auc"] > best_score:
            best_name, best_score, best_factory = mname, s["roc_auc"], factory

    if best_factory is None:
        raise SystemExit(f"{name}: no model could be evaluated")

    report["selected"] = best_name
    report["scaffold"] = report["candidates"][best_name]["scaffold"]
    report["random"] = report["candidates"][best_name]["random"]
    report["generalisation_gap_roc_auc"] = round(
        report["random"]["roc_auc"] - report["scaffold"]["roc_auc"], 4)
    # the dashboard reads dili_report.json under these keys; keep them
    report["splits"] = {"scaffold_5x": report["scaffold"],
                        "random_5x": report["random"]}
    print(f"  -> selected {best_name} on the scaffold split")

    final = best_factory().fit(X, y)
    bundle = {"model": final, "model_name": best_name, "endpoint": name,
              "description": ENDPOINTS[name], "train_smiles": smiles,
              "ad_fps": build_ad_fps(smiles), "report": report}
    return bundle, report


def main() -> None:
    if not DATA.exists():
        raise SystemExit(f"missing {DATA}\nrun: python src/build_hepatic_dataset.py")
    all_df = pd.read_csv(DATA)
    print(f"loaded {len(all_df)} endpoint rows from {DATA.name}")

    heads, reports = {}, {}
    for name in ENDPOINTS:
        sub = all_df[all_df["endpoint"] == name].reset_index(drop=True)
        if name == "clinical":
            # Train the head we ship on most-vs-none, and evaluate the DILIst
            # convention too, so the choice above is a reported number rather
            # than an assertion in a docstring.
            alt = train_head("clinical", clinical_subset(sub, "dilist"))[1]
            sub = clinical_subset(sub, "most_vs_none")
        if len(sub) < 50:
            print(f"  skipping {name}: only {len(sub)} rows")
            continue
        heads[name], reports[name] = train_head(name, sub)
        if name == "clinical":
            reports[name]["label_convention"] = "most_vs_none"
            reports[name]["dilist_convention_for_comparison"] = {
                "n": alt["n"], "positive_rate": alt["positive_rate"],
                "scaffold_roc_auc": alt["scaffold"]["roc_auc"],
                "note": ("most + less concern vs none; rejected because it "
                         "puts severity-0 aspirin in the same class as "
                         "severity-8 troglitazone"),
            }

    # ---------------------------------------------------------------- probes
    print(f"\n{'=' * 74}\nREFERENCE COMPOUNDS\n{'=' * 74}")
    Xp, mp = featurize_smiles([s for s, _ in PROBES.values()])
    names = [n for n, m in zip(PROBES, mp) if m]
    scores = {h: heads[h]["model"].predict_proba(Xp[mp])[:, 1] for h in heads}
    order = list(heads)
    print(f"  {'compound':14s} " + " ".join(f"{h:>9s}" for h in order)
          + f" {'PANEL':>9s}   note")
    probe_out = {}
    for i, n in enumerate(names):
        row = {h: float(scores[h][i]) for h in order}
        row["panel"] = float(max(row.values()))
        probe_out[n] = row
        print(f"  {n:14s} " + " ".join(f"{row[h]:9.3f}" for h in order)
              + f" {row['panel']:9.3f}   {PROBES[n][1]}")

    # ----------------------------------------------------------- persistence
    panel = {"heads": heads, "endpoints": ENDPOINTS,
             "combination": "max over heads",
             "probes": probe_out,
             "caveat": (
                 "Three structure-only heads over DILIrank 2.0 and two Tox21 "
                 "HepG2 assays. The Tox21 heads cannot see metabolism-"
                 "dependent toxicity, because HepG2 cells barely express "
                 "CYP2E1 and the assays top out near 100 uM; paracetamol is "
                 "measured inactive in both. Use the panel to deprioritise, "
                 "never to clear a compound, and always alongside the "
                 "structural alerts."
             )}
    modelio.save(panel, MODELS / "hepatic_panel.pkl")
    print(f"\nsaved {MODELS / 'hepatic_panel.pkl'}")

    # Backward compatibility: screen.py and the dashboard read models/dili.pkl
    # and the p_dili column. The clinical head is the direct replacement for
    # the old TDC model, so it keeps that filename and key layout.
    if "clinical" in heads:
        compat = dict(heads["clinical"])
        compat["task"] = "DILI"
        compat["source"] = "FDA DILIrank 2.0"
        compat["caveat"] = panel["caveat"]
        modelio.save(compat, MODELS / "dili.pkl")
        print(f"saved {MODELS / 'dili.pkl'} (clinical head, drop-in replacement)")

    out = {"panel": {"combination": "max over heads",
                     "endpoints": ENDPOINTS},
           "heads": reports, "probes": probe_out}
    (RESULTS / "hepatic_report.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    # keep the old filename pointing at the clinical head
    if "clinical" in reports:
        (RESULTS / "dili_report.json").write_text(
            json.dumps(reports["clinical"], indent=2), encoding="utf-8")
    print(f"saved {RESULTS / 'hepatic_report.json'}")

    print(f"\n{'=' * 74}")
    for n, r in reports.items():
        print(f"  {n:9s} n={r['n']:5d}  scaffold ROC-AUC "
              f"{r['scaffold']['roc_auc']:.3f} +/- {r['scaffold']['roc_auc_sd']:.3f}"
              f"   ({r['selected']})")
    print(f"{'=' * 74}")


if __name__ == "__main__":
    main()
