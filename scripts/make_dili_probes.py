"""Precompute hepatic-panel scores for a fixed probe set.

The dashboard never loads a model to draw a table; it reads this file.

The probe set is chosen so a reviewer can check the panel against facts they
already know, including the ones that make it look bad. Paracetamol is the
important one, and the story has changed: the old TDC-trained model scored it
0.37 - "low risk" for the index compound of this whole project - because
paracetamol was absent from that benchmark. It is present in DILIrank 2.0, as
most-DILI-concern, so the clinical head now scores it correctly.

That is a closed data gap, not a validated prediction. A compound in the
training set cannot test the model. Honest validation of this panel lives in
results/hepatic_validation.json: a temporal holdout at ROC-AUC 0.68.

    python scripts/make_dili_probes.py
"""

from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import modelio  # noqa: E402
from featurize import featurize_smiles  # noqa: E402
from train import max_similarity_from_fps, query_fingerprints  # noqa: E402

PROBES = {
    "paracetamol": ("CC(=O)Nc1ccc(O)cc1",
                    "in DILIrank as most-concern; old model scored it 0.37"),
    "isoniazid": ("NNC(=O)c1ccncc1", "classic hepatotoxin; most-concern"),
    "diclofenac": ("O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "most-concern"),
    "troglitazone": ("Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",
                     "withdrawn from market for hepatotoxicity"),
    "valproate": ("CCCC(CCC)C(=O)O", "most-concern"),
    "aspirin": ("CC(=O)Oc1ccccc1C(=O)O",
                "less-concern, so excluded from training: extrapolation"),
    "ibuprofen": ("CC(C)Cc1ccc(C(C)C(=O)O)cc1", "low hepatic concern"),
    "caffeine": ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "no-concern control"),
    "curcumin": ("COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
                 "plant NP; measured mitochondria-active in Tox21"),
    "kaempferol": ("O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
                   "plant NP; measured mitochondria-active in Tox21"),
    "capsaicin": ("COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
                  "plant NP; our TRPV1 positive control"),
}


def main() -> None:
    path = ROOT / "models" / "hepatic_panel.pkl"
    if not path.exists():
        raise SystemExit(f"{path} not found - run src/train_hepatic.py first")
    panel = modelio.load(path)
    heads = panel["heads"]

    names = list(PROBES)
    smiles = [PROBES[n][0] for n in names]
    X, mask = featurize_smiles(smiles)
    Q, q_idx = query_fingerprints(smiles)

    scored = {h: heads[h]["model"].predict_proba(X[mask])[:, 1] for h in heads}
    ads = {h: max_similarity_from_fps(Q, q_idx, len(smiles), heads[h]["ad_fps"])
           for h in heads}

    out, k = [], 0
    for i, name in enumerate(names):
        if not mask[i]:
            continue
        row = {h: round(float(scored[h][k]), 4) for h in heads}
        rec = {
            "name": name, "smiles": smiles[i], "note": PROBES[name][1],
            "heads": row,
            "p_hepatic": round(float(max(row.values())), 4),
            # kept under the old key so nothing downstream breaks
            "p_dili": row.get("clinical", max(row.values())),
            "ad": round(float(ads.get("clinical", list(ads.values())[0])[i]), 4),
        }
        out.append(rec)
        k += 1

    n_train = sum(h["report"]["n"] for h in heads.values())
    payload = {
        "model": "hepatic_panel.pkl",
        "heads": {h: heads[h]["report"]["n"] for h in heads},
        "n_training_compounds": int(n_train),
        "combination": panel.get("combination", "max over heads"),
        "caveat": panel.get("caveat", ""),
        "probes": out,
    }
    dest = ROOT / "results" / "dili_probes.json"
    dest.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    cols = list(heads)
    print(f"{'compound':14s} " + " ".join(f"{c:>9s}" for c in cols)
          + f" {'PANEL':>9s}  note")
    for p in out:
        print(f"{p['name']:14s} "
              + " ".join(f"{p['heads'][c]:9.3f}" for c in cols)
              + f" {p['p_hepatic']:9.3f}  {p['note']}")
    print(f"\n{n_train:,} training compounds across {len(heads)} heads")
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
