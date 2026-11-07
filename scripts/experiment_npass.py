"""Close the coverage gap with real natural-product measurements from NPASS.

The diagnosis was that ChEMBL contains synthetic medicinal chemistry and the
screening library is plant chemistry, so the models are extrapolating. NPASS
(Natural Product Activity and Species Source, v3.0) is the obvious fix: it
carries ~1M measured activity records for natural products, with targets
cross-linked to UniProt - the same accessions our panel is defined by.

What it actually gives us for this panel, after joining to structures:

    TRPV1    9 compounds new to us
    CB1    280 compounds new to us
    FAAH    19 compounds new to us

Most of the new material is measured-WEAK natural products, and that is the
valuable part rather than a disappointment. Our negatives were 60k synthetic
decoys: random drug-like molecules presumed inactive. These are real
measurements on real natural products that were tested and found weak - hard
negatives drawn from exactly the chemical space the screen operates in. That
addresses the "your decoys are too easy" criticism with data rather than
argument.

LEAKAGE: every compound in the literature validation set is removed from the
NPASS import before training. Without that, adding NPASS would quietly put the
test set into training and the result would be meaningless.

    python scripts/experiment_npass.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from featurize import featurize_smiles                     # noqa: E402
from standardize import standardize                        # noqa: E402
import train as train_mod                                  # noqa: E402
from train import train_one                                # noqa: E402

# train_one() writes models/<target>_<scope>_<label>.pkl, which are the
# PRODUCTION model files. Redirect it at an experiments folder so this script
# cannot clobber the deployed pipeline.
train_mod.MODELS = ROOT / "models" / "experiments"
train_mod.MODELS.mkdir(parents=True, exist_ok=True)
EXP_MODELS = train_mod.MODELS

NPASS = ROOT / "data" / "external" / "npass"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
INTERIM = ROOT / "data" / "interim"

SEED = 42
BACKGROUND_N = 30_000
KEEP_TYPES = ["IC50", "EC50", "Ki", "Kd", "AC50"]
PANEL = {
    "TRPV1": ["Q8NER1", "O35433", "Q704Y3"],
    "CB1": ["P21554", "P20272", "P47746"],
    "FAAH": ["O00519", "P97612", "O08914"],
}
HUMAN = {"TRPV1": "Q8NER1", "CB1": "P21554", "FAAH": "O00519"}
MODEL_CUT_UM = {"TRPV1": 1.0, "CB1": 1.0, "FAAH": 10.0}
ACTIVE_CUT, INACTIVE_CUT, FAAH_CUT = 6.0, 5.0, 5.0


def build_npass_rows(exclude: set[str]) -> pd.DataFrame:
    """NPASS activities for our panel, standardized, test-set compounds removed."""
    tgt = pd.read_csv(NPASS / "NPASS3.0_target.txt", sep="\t", low_memory=False)
    act = pd.read_csv(NPASS / "NPASS3.0_activities.txt", sep="\t", low_memory=False)
    st = pd.read_csv(NPASS / "NPASS3.0_naturalproducts_structure.txt",
                     sep="\t", low_memory=False)
    smiles_of = st.set_index("np_id")["SMILES"].to_dict()

    out = []
    for target, accs in PANEL.items():
        t = tgt[tgt["uniprot_id"].isin(accs)]
        human_ids = set(t[t["uniprot_id"] == HUMAN[target]]["target_id"])
        a = act[act["target_id"].isin(set(t["target_id"]))].copy()
        a = a[a["activity_type"].isin(KEEP_TYPES) & (a["activity_units"] == "nM")]
        a["v"] = pd.to_numeric(a["activity_value"], errors="coerce")
        a = a[a["v"] > 0]
        a["pactivity"] = 9 - np.log10(a["v"])
        a = a[(a["pactivity"] >= 2) & (a["pactivity"] <= 12)]
        a["organism_group"] = np.where(a["target_id"].isin(human_ids), "human", "other")

        for np_id, g in a.groupby("np_id"):
            smi = smiles_of.get(np_id)
            if not isinstance(smi, str) or not smi.strip():
                continue
            r = standardize(smi)
            if not r.ok or r.inchikey in exclude:
                continue
            for org, gg in g.groupby("organism_group"):
                out.append({
                    "inchikey": r.inchikey,
                    "std_smiles": r.smiles,
                    "target_name": target,
                    "organism_group": org,
                    "pactivity": float(gg["pactivity"].median()),
                    "relation": "=",
                    "n_measurements": int(len(gg)),
                    "source": "NPASS",
                })
    return pd.DataFrame(out)


def label(p: float, target: str) -> float:
    cut_act = FAAH_CUT if target == "FAAH" else ACTIVE_CUT
    if p >= cut_act:
        return 1.0
    if p <= INACTIVE_CUT:
        return 0.0
    return np.nan


def main() -> None:
    lit = pd.read_csv(RESULTS / "literature_validation.csv")
    held = lit[lit["in_chembl_training"] == False].copy()
    held["fair"] = held["value_uM"] <= held["target"].map(MODEL_CUT_UM)
    exclude = set(lit["inchikey"].dropna())
    print(f"excluding {len(exclude)} validation compounds from any NPASS import")

    print("extracting and standardizing NPASS activities...", flush=True)
    npass = build_npass_rows(exclude)
    agg = pd.read_csv(PROC / "bioactivity_aggregated.csv")
    known = set(agg["inchikey"])
    npass_new = npass[~npass["inchikey"].isin(known)]
    print(f"  NPASS rows for the panel: {len(npass):,}")
    print(f"  compounds new to our training set: {npass_new['inchikey'].nunique():,}")
    for t in PANEL:
        s = npass_new[npass_new["target_name"] == t]
        if len(s):
            print(f"    {t:6s} {s['inchikey'].nunique():4d} new "
                  f"(pActivity {s['pactivity'].min():.2f}-{s['pactivity'].max():.2f})")

    # ---- merge: NPASS only contributes compounds ChEMBL does not have ----
    merged = pd.concat([agg.assign(source="ChEMBL"), npass_new], ignore_index=True)
    merged["label"] = [label(p, t) for p, t in
                       zip(merged["pactivity"], merged["target_name"])]
    merged["label_10um"] = np.where(merged["pactivity"] >= FAAH_CUT, 1.0, 0.0)
    merged.to_csv(INTERIM / "bioactivity_with_npass.csv", index=False)

    # ---- decoys ----------------------------------------------------------
    cache = INTERIM / "decoy_pool_standardized.csv"
    pool = pd.read_csv(cache)
    pool["cell"] = [tuple(eval(c)) if isinstance(c, str) else c for c in pool["cell"]]

    # ---- background, featurized once -------------------------------------
    print(f"\nsampling {BACKGROUND_N:,} background compounds...", flush=True)
    bg = pd.read_csv(RESULTS / "screening_all.csv", usecols=["std_smiles"],
                     low_memory=False).sample(n=BACKGROUND_N, random_state=SEED)
    Xbg, mbg = featurize_smiles(bg["std_smiles"].tolist())
    Xbg = Xbg[mbg]

    Xl, ml = featurize_smiles(held["std_smiles"].tolist())

    rows, summary = [], {}
    for target in PANEL:
        label_col = "label_10um" if target == "FAAH" else "label"
        df = merged[merged["target_name"] == target]
        print(f"\n{'='*66}\n{target} retrained with NPASS added\n{'='*66}", flush=True)
        rep = train_one(target, "human", df, label_col, pool)
        if rep is None:
            print(f"  {target}: not trainable")
            continue
        import pickle
        with open(EXP_MODELS / f"{target}_human_{label_col}.pkl", "rb") as fh:
            bundle = pickle.load(fh)
        model = bundle["model"]

        bg_pred = model.predict_proba(Xbg)[:, 1]
        sub = held[held["target"] == target]
        for i, (_, r) in enumerate(held.iterrows()):
            if not ml[i] or r["target"] != target:
                continue
            p = float(model.predict_proba(Xl[i:i + 1])[0, 1])
            rows.append({
                "compound_name": r["compound_name"], "target": target,
                "value_uM": r["value_uM"], "fair": bool(r["fair"]),
                "p_npass": p,
                "pct_npass": round(float((bg_pred < p).mean() * 100), 1),
                "pct_baseline": r["library_percentile"],
            })

    res = pd.DataFrame(rows)
    res.to_csv(RESULTS / "npass_experiment.csv", index=False)

    fair = res[res["fair"]]
    a = fair["pct_npass"].to_numpy(dtype=float)
    b = fair["pct_baseline"].to_numpy(dtype=float)
    w = sps.wilcoxon(a - 50, alternative="greater")
    paired = sps.wilcoxon(a - b, alternative="greater") if len(a) > 1 else None
    summary = {
        "npass_compounds_added": int(npass_new["inchikey"].nunique()),
        "n_fair": int(len(a)),
        "median_baseline": round(float(np.median(b)), 1),
        "median_with_npass": round(float(np.median(a)), 1),
        "gain": round(float(np.median(a) - np.median(b)), 1),
        "wilcoxon_vs_chance_p": round(float(w.pvalue), 4),
        "paired_vs_baseline_p": round(float(paired.pvalue), 4) if paired else None,
        "top10pct_baseline": int((b >= 90).sum()),
        "top10pct_with_npass": int((a >= 90).sum()),
        "improved": int((a > b).sum()),
    }
    (RESULTS / "npass_experiment.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print(f"\n{'='*78}\nHELD-OUT LITERATURE ACTIVES: baseline vs NPASS-augmented\n{'='*78}")
    show = res.sort_values("pct_npass", ascending=False)
    print(show[["compound_name", "target", "value_uM", "fair",
                "pct_baseline", "pct_npass"]].to_string(index=False))
    print(f"\n{'='*78}\nSUMMARY\n{'='*78}")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
