"""Smoke test for the screening path.

Screening 150k compounds takes several minutes; a crash at the end wastes all
of it. This exercises every step of screen.py's scoring and composite logic on
a handful of molecules first.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem

from featurize import featurize_smiles
from safety import physchem, profile
from screen import add_composite_columns, load_models, pick_model_per_target
from train import max_similarity

PROBES = {
    "capsaicin": "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
    "kaempferol": "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
    "genistein": "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12",
    "curcumin": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "sucrose": "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O",
}


def main():
    bundles = load_models()
    chosen = pick_model_per_target(bundles)
    print(f"models: {[(t, b['model_name'], b['label_col']) for t, b in chosen.items()]}")

    names = list(PROBES)
    smiles = list(PROBES.values())
    df = pd.DataFrame({"compound_name": names, "std_smiles": smiles,
                       "std_inchikey": [Chem.MolToInchiKey(Chem.MolFromSmiles(s))
                                        for s in smiles]})

    # safety columns, as prepare_np.py would have produced them
    rows = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        sp, pc = profile(m), physchem(m)
        rows.append({
            "n_reactive_alerts": sp.n_alerts,
            "max_alert_severity": sp.max_severity,
            "alert_severity_sum": sp.severity_sum,
            "reactive_alerts": "|".join(sp.alerts),
            "n_pains": len(sp.pains),
            "cns_like": pc["cns_like"],
        })
    df = pd.concat([df, pd.DataFrame(rows)], axis=1)

    X, mask = featurize_smiles(smiles)
    for t, b in chosen.items():
        df[f"p_{t}"] = b["model"].predict_proba(X)[:, 1]
        df[f"ad_{t}"] = max_similarity(smiles, b["ad_fps"])

    target_cols = [f"p_{t}" for t in chosen]
    ad_cols = [f"ad_{t}" for t in chosen]
    out = add_composite_columns(df, target_cols, ad_cols)

    cols = ["compound_name", "best_target", "target_score", "ad_best",
            "n_reactive_alerts", "max_alert_severity", "cns_like",
            "composite_score"]
    print()
    print(out[cols].sort_values("composite_score", ascending=False)
          .to_string(index=False))
    print("\nper-target probabilities:")
    print(out[["compound_name"] + target_cols].to_string(index=False))
    print("\nsmoke test OK")


if __name__ == "__main__":
    main()
