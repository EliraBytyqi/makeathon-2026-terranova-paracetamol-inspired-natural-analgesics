"""Precompute DILI model scores for a fixed probe set.

The headline safety finding of this project is that the DILI model scores
paracetamol as *low risk* - because paracetamol is absent from the TDC
benchmark entirely. The dashboard states that, so the number behind it should
be computed rather than quoted from a comment.

This runs once and writes results/dili_probes.json. The dashboard never loads
a model itself; it reads this file.

    python scripts/make_dili_probes.py
"""

from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from featurize import featurize_smiles  # noqa: E402
from train import max_similarity_from_fps, query_fingerprints  # noqa: E402

# Chosen to make three points at once: the index compound the benchmark is
# missing (paracetamol, isoniazid), compounds the benchmark does contain
# (diclofenac, valproate, aspirin), and the natural products this project
# actually ranks (curcumin, kaempferol, capsaicin).
PROBES = {
    "paracetamol": ("CC(=O)Nc1ccc(O)cc1", "absent from benchmark"),
    "isoniazid": ("NNC(=O)c1ccncc1", "absent from benchmark"),
    "diclofenac": ("O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "benchmark positive"),
    "valproate": ("CCCC(CCC)C(=O)O", "benchmark positive"),
    "aspirin": ("CC(=O)Oc1ccccc1C(=O)O", "benchmark negative"),
    "ibuprofen": ("CC(C)Cc1ccc(C(C)C(=O)O)cc1", "reference analgesic"),
    "curcumin": ("COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
                 "natural product"),
    "kaempferol": ("O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
                   "natural product"),
    "capsaicin": ("COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O", "natural product"),
    "caffeine": ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "negative control"),
}


def main() -> None:
    path = ROOT / "models" / "dili.pkl"
    if not path.exists():
        raise SystemExit(f"{path} not found - run src/train_dili.py first")

    with open(path, "rb") as fh:
        bundle = pickle.load(fh)
    model = bundle["model"]

    names = list(PROBES)
    smiles = [PROBES[n][0] for n in names]
    X, mask = featurize_smiles(smiles)
    probs = model.predict_proba(X[mask])[:, 1]

    # applicability domain: how close is each probe to the 462 training drugs?
    Q, q_idx = query_fingerprints(smiles)
    ad = max_similarity_from_fps(Q, q_idx, len(smiles), bundle["ad_fps"])

    out = []
    k = 0
    for i, name in enumerate(names):
        if not mask[i]:
            continue
        out.append({
            "name": name,
            "smiles": smiles[i],
            "note": PROBES[name][1],
            "p_dili": round(float(probs[k]), 4),
            "ad": round(float(ad[i]), 4),
        })
        k += 1

    payload = {
        "model": "dili.pkl",
        "n_training_compounds": int(bundle["report"]["n"]),
        "caveat": bundle.get("caveat", ""),
        "probes": out,
    }
    dest = ROOT / "results" / "dili_probes.json"
    dest.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"{'compound':16s} {'p(DILI)':>8s} {'AD':>6s}  note")
    for p in out:
        print(f"{p['name']:16s} {p['p_dili']:8.3f} {p['ad']:6.3f}  {p['note']}")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
