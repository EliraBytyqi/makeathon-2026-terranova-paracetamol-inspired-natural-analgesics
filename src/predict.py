"""Score one molecule end to end: three targets, then hepatic liability.

This is the pipeline as a single call. Everything else in the project either
builds the models or runs them over the whole 145k library; this takes one
structure and walks it through the same steps in the same order.

    python src/predict.py "CC(=O)Nc1ccc(O)cc1"
    python src/predict.py --name paracetamol
    python src/predict.py "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O" --json

Stage 1  standardize           same standardizer as every training structure
Stage 2  activity              TRPV1, CB1, FAAH - with the applicability
                               domain for each, because a score from outside
                               the training chemistry is extrapolation
Stage 3  structural alerts     25 SMARTS patterns for reactive-metabolite
                               chemistry, plus PAINS. This layer holds the
                               veto: a severity-3 alert disqualifies a
                               compound whatever the models say
Stage 4  hepatic liability     three heads (DILIrank clinical, Tox21
                               mitochondrial, Tox21 HepG2 cytotoxicity),
                               combined as the maximum
Stage 5  composite             the same formula the library screen uses

Every number carries what is actually known about it. The confidence wording
is not decoration: it comes from the external validation in
results/pubchem_validation_expanded.json and results/hepatic_validation.json,
and it differs sharply between targets.

    CB1     externally validated, including on chemistry unlike training
            (ROC-AUC 0.93 overall, 0.69 on novel chemistry)
    TRPV1   validated near its training chemistry only (0.95); no measured
            inactives exist far from it, so performance there is unknown
    FAAH    weak (0.64) - treat its score as a hint, not a result
    hepatic weak prior (~0.65 externally); the alerts, not the model, are
            what make the safety call
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rdkit import Chem                                      # noqa: E402

import modelio                                              # noqa: E402
import safety                                               # noqa: E402
from featurize import featurize_smiles                      # noqa: E402
from standardize import standardize                         # noqa: E402
from train import max_similarity_from_fps, query_fingerprints  # noqa: E402

MODELS = ROOT / "models"
UA = "terranova-academic-research/1.0"

# must match src/screen.py
ALERT_PENALTY = 0.35
AD_FLOOR = 0.25
CNS_BONUS = 1.15
DILI_WEIGHT = 0.35

TARGETS = {"TRPV1": "TRPV1_human_label.pkl",
           "CB1": "CB1_human_label.pkl",
           "FAAH": "FAAH_human_label_10um.pkl"}

# External validation, one row per structure. See validate_pubchem_expanded.py.
VALIDATION = {
    "TRPV1": {"roc_auc": 0.918, "n": 965, "novel": None,
              "grade": "near-training only",
              "note": ("0.918 on 965 measured structures, but none of its 88 "
                       "measured inactives survive a novelty filter, so "
                       "performance on unfamiliar chemistry is untested")},
    "CB1": {"roc_auc": 0.917, "n": 2104, "novel": 0.711,
            "grade": "validated",
            "note": ("0.917 on 2,104 measured structures and 0.711 on 234 "
                     "structures unlike training - the strongest claim here")},
    "FAAH": {"roc_auc": 0.625, "n": 621, "novel": 0.622,
             "grade": "weak",
             "note": ("0.625 on 621 measured structures. Better than chance "
                      "and not much more; its internal score was inflated by "
                      "presumed-inactive decoys")},
}

HEPATIC_NOTE = (
    "Externally ~0.69: the combined panel scores 0.69 on drugs the FDA added "
    "to DILIrank after v1, and the Tox21 mitochondrial head scores 0.66 "
    "predicting clinical DILI. Treat as a weak prior. No cell-based assay "
    "can see metabolism-dependent toxicity - HepG2 barely expresses CYP2E1, "
    "which is why paracetamol is measured inactive in both Tox21 heads and "
    "why the structural alerts are not redundant."
)


def resolve_name(name: str) -> str | None:
    """Look a compound name up in PubChem and return its SMILES."""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           f"{urllib.parse.quote(name, safe='')}"
           "/property/SMILES,CanonicalSMILES,IsomericSMILES/JSON")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=30) as r:
            props = json.loads(r.read())["PropertyTable"]["Properties"][0]
    except Exception as exc:
        print(f"could not resolve {name!r} through PubChem: {exc}",
              file=sys.stderr)
        return None
    return (props.get("SMILES") or props.get("CanonicalSMILES")
            or props.get("IsomericSMILES"))


def verdict(p: float, ad: float) -> str:
    """Wording for one target. Deliberately not a yes/no."""
    if ad < AD_FLOOR:
        return "out of domain - no prediction"
    if p >= 0.90:
        return "ranks very high"
    if p >= 0.50:
        return "ranks high"
    if p >= 0.15:
        return "middling"
    return "ranks low"


def predict(smiles: str) -> dict:
    # ---------------------------------------------------- stage 1
    std = standardize(smiles)
    if not std.ok:
        return {"ok": False, "error": f"could not standardize {smiles!r}"}
    mol = Chem.MolFromSmiles(std.smiles)
    out: dict = {"ok": True, "input_smiles": smiles,
                 "std_smiles": std.smiles, "inchikey": std.inchikey}

    X, mask = featurize_smiles([std.smiles])
    if not mask[0]:
        return {"ok": False, "error": "could not featurize"}
    Q, q_idx = query_fingerprints([std.smiles])

    # ---------------------------------------------------- stage 2
    act = {}
    for target, fname in TARGETS.items():
        path = MODELS / fname
        if not path.exists():
            continue
        b = modelio.load(path)
        p = float(b["model"].predict_proba(X)[0, 1])
        ad = float(max_similarity_from_fps(Q, q_idx, 1, b["ad_fps"])[0])
        act[target] = {"p": round(p, 4), "ad": round(ad, 4),
                       "in_domain": bool(ad >= AD_FLOOR),
                       "verdict": verdict(p, ad),
                       "validation": VALIDATION[target]}
    out["activity"] = act
    best = max(act, key=lambda t: act[t]["p"]) if act else None
    out["best_target"] = best
    target_score = act[best]["p"] if best else 0.0
    ad_best = act[best]["ad"] if best else 0.0

    # ---------------------------------------------------- stage 3
    prof = safety.profile(mol)
    phys = safety.physchem(mol)
    out["safety"] = {
        "alerts": prof.alerts, "severities": prof.severities,
        "max_severity": prof.max_severity, "severity_sum": prof.severity_sum,
        "pains": prof.pains, "n_pains": len(prof.pains),
        "disqualified": bool(prof.max_severity >= 3 or prof.pains),
        "reason": ("severity-3 reactive-metabolite alert"
                   if prof.max_severity >= 3 else
                   "PAINS assay-interference pattern" if prof.pains else None),
    }
    out["physchem"] = phys

    # ---------------------------------------------------- stage 4
    hep: dict = {}
    panel_path = MODELS / "hepatic_panel.pkl"
    if panel_path.exists():
        panel = modelio.load(panel_path)
        for hname, hb in panel["heads"].items():
            hep[hname] = {
                "p": round(float(hb["model"].predict_proba(X)[0, 1]), 4),
                "ad": round(float(max_similarity_from_fps(
                    Q, q_idx, 1, hb["ad_fps"])[0]), 4),
                "measures": panel["endpoints"].get(hname, ""),
            }
        p_hep = max(v["p"] for v in hep.values())
    else:
        p_hep = 0.5
    out["hepatic"] = {"heads": hep, "p_hepatic": round(float(p_hep), 4),
                      "combination": "max over heads", "note": HEPATIC_NOTE}

    # ---------------------------------------------------- stage 5
    ad_conf = float(np.clip(ad_best / 0.5, 0.0, 1.0))
    composite = (target_score
                 * float(np.exp(-ALERT_PENALTY * prof.severity_sum))
                 * ad_conf
                 * (CNS_BONUS if phys["cns_like"] else 1.0)
                 * (1.0 - DILI_WEIGHT * p_hep))
    out["composite"] = {
        "score": round(composite, 6),
        "terms": {
            "target_score": round(target_score, 4),
            "alert_multiplier": round(float(np.exp(
                -ALERT_PENALTY * prof.severity_sum)), 4),
            "ad_confidence": round(ad_conf, 4),
            "cns_multiplier": CNS_BONUS if phys["cns_like"] else 1.0,
            "hepatic_multiplier": round(1.0 - DILI_WEIGHT * p_hep, 4),
        },
        "shortlist_eligible": bool(
            ad_best >= AD_FLOOR and prof.max_severity < 3 and not prof.pains),
    }
    return out


def render(r: dict) -> None:
    if not r.get("ok"):
        print(f"ERROR: {r.get('error')}")
        return
    w = 78
    print("=" * w)
    print(f"  {r['std_smiles']}")
    print(f"  InChIKey {r['inchikey']}")
    print("=" * w)

    print("\nSTAGE 2 - ACTIVITY ON THE THREE TARGETS")
    print(f"  {'target':7s} {'score':>7s} {'AD':>6s}  {'reading':24s} evidence")
    for t, a in r["activity"].items():
        print(f"  {t:7s} {a['p']:7.3f} {a['ad']:6.3f}  {a['verdict']:24s} "
              f"{a['validation']['grade']}")
    for t, a in r["activity"].items():
        print(f"    {t}: {a['validation']['note']}")

    s = r["safety"]
    print("\nSTAGE 3 - STRUCTURAL ALERTS (this layer holds the veto)")
    if s["alerts"]:
        for name, sev in zip(s["alerts"], s["severities"]):
            print(f"  severity {sev}  {name}")
    else:
        print("  no reactive-metabolite alerts")
    print(f"  PAINS: {s['n_pains']}" + (f"  {', '.join(s['pains'][:2])}"
                                        if s["pains"] else ""))
    if s["disqualified"]:
        print(f"  >> DISQUALIFIED from the shortlist: {s['reason']}")

    h = r["hepatic"]
    print("\nSTAGE 4 - HEPATIC LIABILITY")
    for name, v in h["heads"].items():
        print(f"  {name:9s} {v['p']:6.3f}  (AD {v['ad']:.3f})  {v['measures']}")
    print(f"  {'COMBINED':9s} {h['p_hepatic']:6.3f}  = max over heads")

    c = r["composite"]
    print("\nSTAGE 5 - COMPOSITE")
    for k, v in c["terms"].items():
        print(f"  {k:20s} {v}")
    print(f"  {'COMPOSITE SCORE':20s} {c['score']:.6f}")
    print(f"  shortlist eligible: {'yes' if c['shortlist_eligible'] else 'NO'}")

    print("\n" + "-" * w)
    print("  What this is: a ranking, not a verdict. It says where this")
    print("  molecule belongs in a queue for someone with a laboratory.")
    print(f"  Hepatic layer: {h['note'][:72]}")
    print("-" * w)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Score one molecule through the Terranova pipeline.")
    ap.add_argument("smiles", nargs="?", help="SMILES string")
    ap.add_argument("--name", help="compound name, resolved via PubChem")
    ap.add_argument("--json", action="store_true", help="emit JSON only")
    a = ap.parse_args()

    smi = a.smiles
    if not smi and a.name:
        smi = resolve_name(a.name)
        if not smi:
            raise SystemExit(2)
        if not a.json:
            print(f"resolved {a.name!r} -> {smi}")
    if not smi:
        ap.error("give a SMILES string or --name")

    r = predict(smi)
    if a.json:
        print(json.dumps(r, indent=2))
    else:
        render(r)
    if not r.get("ok"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
