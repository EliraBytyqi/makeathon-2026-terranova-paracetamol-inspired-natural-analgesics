"""Consolidate the PubChem retrospective validation for the dashboard.

Writes results/pubchem_summary.json: the CB1 headline with bootstrap CIs, the
novelty-restricted subset, the score separation, the notable correct
rejections, and an explicit record of why TRPV1 and FAAH could not be scored.
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
from train import max_similarity                           # noqa: E402

RESULTS = ROOT / "results"
rng = np.random.default_rng(42)


def boot_ci(y, s, n=5000):
    y, s = np.asarray(y), np.asarray(s)
    idx, out = np.arange(len(y)), []
    for _ in range(n):
        i = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(y[i])) < 2:
            continue
        out.append(roc_auc_score(y[i], s[i]))
    return [round(float(x), 3) for x in np.percentile(out, [2.5, 97.5])]


def main() -> None:
    v = pd.read_csv(RESULTS / "pubchem_validation.csv")
    smiles = (pd.read_csv(RESULTS / "screening_all.csv",
                          usecols=["std_inchikey", "std_smiles"], low_memory=False)
              .set_index("std_inchikey")["std_smiles"])
    with open(ROOT / "models" / "CB1_human_label.pkl", "rb") as fh:
        bundle = pickle.load(fh)

    cb = v[v["target"] == "CB1"].copy()
    cb["ad"] = max_similarity(smiles.reindex(cb["std_inchikey"]).tolist(),
                              bundle["ad_fps"])
    cb["y"] = (cb["outcome"] == "Active").astype(int)

    tiers = {}
    for name, sub in [("all", cb), ("novel_ad_lt_0.5", cb[cb["ad"] < 0.5]),
                      ("novel_ad_lt_0.4", cb[cb["ad"] < 0.4])]:
        y, s = sub["y"].to_numpy(), sub["model_score"].to_numpy()
        if len(np.unique(y)) < 2:
            continue
        tiers[name] = {
            "n": int(len(y)), "n_active": int(y.sum()),
            "n_inactive": int((1 - y).sum()),
            "roc_auc": round(float(roc_auc_score(y, s)), 3),
            "ci95": boot_ci(y, s),
            "p": float(sps.mannwhitneyu(s[y == 1], s[y == 0],
                                        alternative="greater").pvalue),
        }

    # Pick ONE representative per chemical family, preferring the name a
    # chemist would recognise. The raw table is dominated by caryophyllene
    # synonyms and registry codes ("Einecs 234-173-0"), which say nothing to
    # a reader; "beta-caryophyllene" is the entire point of this panel.
    inact = cb[cb["outcome"] == "Inactive"].sort_values("model_score", ascending=False)

    # Group by the InChIKey skeleton (first 14 chars): that is the same
    # molecule regardless of how many synonyms and stereo-spellings PubChem
    # carries for it. Grouping by name fails here because the caryophyllene
    # isomers are listed under IUPAC names that never contain the word.
    inact = inact.assign(_skel=inact["std_inchikey"].str[:14])

    def readability(nm: str) -> tuple:
        """Lower sorts better: short, no registry codes, no IUPAC soup."""
        low = nm.lower()
        return (
            low.startswith("einecs") or low.startswith("unii"),
            sum(c.isdigit() for c in nm) > 4,
            nm.count("(") + nm.count("["),
            len(nm),
        )

    notable = []
    for _skel, grp in inact.groupby("_skel", sort=False):
        best = min(grp.itertuples(), key=lambda r: readability(str(r.compound_name)))
        notable.append({"compound_name": str(best.compound_name)[:60],
                        "model_score": float(best.model_score),
                        "ad": round(float(best.ad), 2),
                        "n_synonyms": int(len(grp))})
        if len(notable) >= 6:
            break
    notable.sort(key=lambda d: -d["model_score"])

    unscorable = {}
    for t in ("TRPV1", "FAAH"):
        s = v[v["target"] == t]
        unscorable[t] = {
            "n": int(len(s)),
            "n_active": int((s["outcome"] == "Active").sum()),
            "n_inactive": int((s["outcome"] == "Inactive").sum()),
            "reason": "only one outcome class present - ROC-AUC undefined",
        }

    out = {
        "source": "PubChem BioAssay, retrieved via PUG REST for all nine panel accessions",
        "design": ("library compounds with an independent measured Active/Inactive "
                   "outcome, excluded from every production training set"),
        "total_overlap": int(len(v)),
        "training_skeletons_excluded": 44516,
        "cb1": tiers,
        "separation": {
            "active_median": round(float(cb[cb["y"] == 1]["model_score"].median()), 4),
            "inactive_median": round(float(cb[cb["y"] == 0]["model_score"].median()), 4),
        },
        "notable_correct_rejections": notable,
        "unscorable": unscorable,
        "caveats": [
            "PubChem's records skew active (74.6% base rate here), so enrichment "
            "factors are meaningless; ROC-AUC is base-rate independent and is the "
            "figure to quote.",
            "Retrospective: the experiments already existed. It is a blind test "
            "because the model never saw the outcomes, not a prospective one.",
            "TRPV1 and FAAH overlaps were single-class and cannot be scored. The "
            "result is CB1 only.",
        ],
    }
    dest = RESULTS / "pubchem_summary.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    for k, t in tiers.items():
        print(f"  {k:18s} n={t['n']:3d}  ROC-AUC {t['roc_auc']:.3f} "
              f"CI{t['ci95']}  p={t['p']:.3g}")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
