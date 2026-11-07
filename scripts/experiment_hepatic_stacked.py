"""Can mechanism features make the clinical DILI head stronger?

The bottleneck is not the algorithm, it is that DILIrank gives 506 usable
drugs and nothing else does better. A structure-only model fitted on 506
compounds reaches ROC-AUC 0.676 on the temporal holdout (train on DILIrank v1,
test on the drugs the FDA added in v2), and 0.649 on the subset of those with
scaffolds absent from training, with a confidence interval that reaches 0.50.

Tox21 has roughly 8,000 compounds measured across many liver-relevant
mechanisms. Those compounds do not have clinical DILI labels, so they cannot
be pooled with DILIrank. But their assay outcomes can be turned into features:
train one classifier per assay, then use its predicted probability as an input
to the clinical model. The clinical model then sees "this structure looks like
things that depolarise mitochondria and trigger the oxidative stress response"
instead of only the structure itself. That is the standard way to spend a
large unlabelled-for-your-task dataset on a small labelled one.

Mechanism panel, all Tox21, all in PubChem:

  mmp_ratio      mitochondrial membrane potential (HepG2)
  mmp_rhod       the same assay, second readout
  are_nrf2       Nrf2/antioxidant-response element (HepG2) - this is the
                 oxidative-stress pathway NAPQI actually triggers
  casp3          caspase-3 activation (HepG2) - apoptosis
  cytotox_flo    CellTox Green cytotoxicity (HepG2)
  viability_glo  RealTime-Glo viability (HepG2)
  pxr            PXR agonism (HepG2) - drives CYP3A4 induction
  car            CAR agonism (HepG2) - also drives CYP induction
  p53            p53 activation (HCT-116) - DNA damage response
  hek293_viab    viability in HEK293, a non-liver line. Included deliberately
                 as a control channel: a compound toxic in HepG2 but not in
                 HEK293 is liver-selective, and the model can only learn that
                 distinction if it can see both.

Leakage control. Each assay head is trained with every DILIrank structure
removed, so no mechanism feature for a clinical compound comes from a model
that was fitted on that compound. This costs a few hundred of ~8,000 assay
compounds and makes the stacked features honest.

Pre-registered comparison, so this cannot be rewritten into a success after
the fact. The number to beat is 0.676 on the temporal holdout, and 0.649 on
new scaffolds. Three variants are reported:

    structure_only     2,062 features - the current model, the baseline
    structure_plus     2,062 + 10 mechanism features
    mechanism_only     the 10 mechanism features alone

If structure_plus does not beat structure_only, that is the finding and it
gets reported as such.

    python scripts/experiment_hepatic_stacked.py
"""
from __future__ import annotations

import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from featurize import featurize_smiles, scaffold_groups   # noqa: E402
from standardize import standardize                       # noqa: E402
from train import make_models                             # noqa: E402

TOX21 = ROOT / "data" / "external" / "tox21"
SMI = ROOT / "data" / "external" / "tox21" / "smiles"
XLSX = ROOT / "data" / "external" / "dilirank" / "dilirank.xlsx"
HEP = ROOT / "data" / "processed" / "hepatic_endpoints.csv"
RESULTS = ROOT / "results"
UA = "terranova-academic-research/1.0"
SEED = 42
N_BOOT = 2000

PANEL = {
    "mmp_ratio":     (2283792, "mitochondrial membrane potential, HepG2"),
    "mmp_rhod":      (2283443, "mitochondrial membrane potential, rhodamine"),
    "are_nrf2":      (2283544, "Nrf2/ARE oxidative stress response, HepG2"),
    "casp3":         (2283989, "caspase-3 apoptosis, HepG2"),
    "cytotox_flo":   (2283888, "CellTox Green cytotoxicity, HepG2"),
    "viability_glo": (2283892, "RealTime-Glo viability, HepG2"),
    "pxr":           (2283988, "PXR agonism, HepG2 (CYP3A4 induction)"),
    "car":           (2283851, "CAR agonism, HepG2 (CYP induction)"),
    "p53":           (2283601, "p53 activation, HCT-116 (DNA damage)"),
    "hek293_viab":   (2284005, "viability in HEK293 - non-liver control"),
}


def get(url: str, dest: Path) -> bytes | None:
    if dest.exists():
        return dest.read_bytes() or None
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=180) as r:
            data = r.read()
        dest.write_bytes(data)
        time.sleep(0.3)
        return data
    except Exception as exc:
        print(f"    {dest.name}: {exc}")
        dest.write_bytes(b"")
        return None


def assay_table(short: str) -> pd.DataFrame:
    aid, _ = PANEL[short]
    path = TOX21 / f"aid_{aid}.csv"
    if not path.exists():
        print(f"    downloading AID {aid} ({short})")
        if not get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"
                   f"{aid}/concise/CSV", path):
            return pd.DataFrame()
    d = pd.read_csv(path, low_memory=False)
    d = d[d["CID"].notna() & d["Activity Outcome"].isin(["Active", "Inactive"])].copy()
    d["CID"] = d["CID"].astype(int)
    d["hit"] = (d["Activity Outcome"] == "Active").astype(int)
    return d.groupby("CID")["hit"].max().rename(short).reset_index()


def fetch_smiles(cids: list[int]) -> dict[int, str]:
    out: dict[int, str] = {}
    for i in range(0, len(cids), 150):
        block = cids[i:i + 150]
        raw = get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                  + ",".join(map(str, block))
                  + "/property/SMILES,CanonicalSMILES/JSON",
                  SMI / f"{block[0]}_{len(block)}.json")
        if not raw:
            continue
        try:
            for p in json.loads(raw)["PropertyTable"]["Properties"]:
                s = p.get("SMILES") or p.get("CanonicalSMILES")
                if s:
                    out[int(p["CID"])] = s
        except Exception:
            pass
    return out


def stats_block(y, s) -> dict:
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    if len(set(y)) < 2:
        return {"n": int(len(y)), "scoreable": False}
    auc = float(roc_auc_score(y, s))
    rng = np.random.RandomState(SEED)
    b = [roc_auc_score(y[i], s[i]) for i in
         (rng.randint(0, len(y), len(y)) for _ in range(N_BOOT))
         if len(set(y[i])) > 1]
    return {"n": int(len(y)), "n_positive": int(y.sum()), "scoreable": True,
            "roc_auc": round(auc, 4),
            "ci95": [round(float(np.percentile(b, 2.5)), 4),
                     round(float(np.percentile(b, 97.5)), 4)] if b else None,
            "mannwhitney_p": float(sps.mannwhitneyu(
                s[y == 1], s[y == 0], alternative="greater").pvalue)}


def rf():
    return make_models()["random_forest"]


def main() -> None:
    # ------------------------------------------------- 1. clinical reference
    clin = pd.read_csv(HEP)
    clin = clin[clin["endpoint"] == "clinical"].copy()
    clin = clin[clin["concern"].isin(["vmost-dili-concern", "vno-dili-concern"])]
    clin["y"] = (clin["concern"] == "vmost-dili-concern").astype(int)
    v2 = pd.read_excel(XLSX, sheet_name="version 2", header=1)
    v2.columns = [str(c).strip() for c in v2.columns]
    flag = (v2.assign(n=v2["CompoundName"].astype(str).str.strip())
              .set_index("n")["Comment"].astype(str).str.strip().str.lower())
    clin["is_new"] = clin["name"].map(flag).eq("new")
    clin = clin.reset_index(drop=True)
    clin_skel = set(clin["inchikey"].str[:14])
    print(f"clinical: {len(clin)} drugs ({int(clin['y'].sum())} most-concern), "
          f"{int(clin['is_new'].sum())} added in v2")

    # ------------------------------------------------- 2. the mechanism panel
    print(f"\n=== assembling {len(PANEL)} Tox21 endpoints ===")
    tabs, kept = [], []
    for short in PANEL:
        t = assay_table(short)
        if t.empty:
            print(f"  {short:14s} unavailable, skipped")
            continue
        print(f"  {short:14s} {len(t):5d} compounds, "
              f"{int(t[short].sum()):5d} active ({t[short].mean():.1%})")
        tabs.append(t.set_index("CID"))
        kept.append(short)
    mech = pd.concat(tabs, axis=1).reset_index()
    print(f"  union: {len(mech)} compounds over {len(kept)} endpoints")

    print("\n=== structures ===")
    smap = fetch_smiles(sorted(mech["CID"].astype(int).unique()))
    mech["raw"] = mech["CID"].astype(int).map(smap)
    mech = mech[mech["raw"].notna()].copy()
    std = [standardize(s) for s in mech["raw"]]
    mech["smiles"] = [r.smiles if r.ok else None for r in std]
    mech["skel"] = [r.inchikey[:14] if r.ok else None for r in std]
    mech = mech[mech["smiles"].notna()].drop_duplicates("skel").reset_index(drop=True)
    print(f"  usable: {len(mech)}")

    # Remove every clinical structure from the assay training pool, so no
    # mechanism feature is produced by a model that saw that compound.
    pool = mech[~mech["skel"].isin(clin_skel)].reset_index(drop=True)
    print(f"  assay training pool after removing clinical structures: {len(pool)} "
          f"({len(mech) - len(pool)} removed)")

    Xpool, mp = featurize_smiles(pool["smiles"].tolist())
    Xpool, pool = Xpool[mp], pool[mp].reset_index(drop=True)
    Xclin, mc = featurize_smiles(clin["std_smiles"].tolist())
    Xclin, clin = Xclin[mc], clin[mc].reset_index(drop=True)

    # ------------------------------------------- 3. one head per endpoint
    print(f"\n=== training {len(kept)} assay heads and scoring the clinical set ===")
    feats, used = {}, []
    for short in kept:
        yv = pool[short]
        ok = yv.notna().to_numpy()
        y = yv[ok].to_numpy().astype(int)
        if len(set(y)) < 2 or y.sum() < 30:
            print(f"  {short:14s} too few actives, skipped")
            continue
        m = rf().fit(Xpool[ok], y)
        feats[short] = m.predict_proba(Xclin)[:, 1]
        used.append(short)
        print(f"  {short:14s} trained on {int(ok.sum()):5d} "
              f"({int(y.sum()):4d} act) -> mechanism feature")
    M = np.column_stack([feats[s] for s in used])
    print(f"  mechanism feature matrix: {M.shape}")

    # ------------------------------------------- 4. the pre-registered test
    y = clin["y"].to_numpy()
    is_new = clin["is_new"].to_numpy()
    groups = scaffold_groups(clin["std_smiles"].tolist())
    variants = {
        "structure_only": Xclin,
        "structure_plus_mechanism": np.hstack([Xclin, M]),
        "mechanism_only": M,
    }

    tr, te = ~is_new, is_new
    g_tr = set(np.asarray(groups)[tr])
    novel = np.array([g not in g_tr for g in np.asarray(groups)[te]])
    print(f"\n{'=' * 78}\nTEMPORAL HOLDOUT: train on DILIrank v1, test on v2 additions"
          f"\n{'=' * 78}")
    print(f"  train {int(tr.sum())} ({int(y[tr].sum())} most-concern) | "
          f"test {int(te.sum())} ({int(y[te].sum())} most-concern), "
          f"{int(novel.sum())} on new scaffolds")

    out = {"baseline_to_beat": {"all": 0.676, "novel_scaffold": 0.649},
           "n_train": int(tr.sum()), "n_test": int(te.sum()),
           "mechanism_features": used, "temporal": {}, "scaffold_cv": {}}
    for name, X in variants.items():
        m = rf().fit(X[tr], y[tr])
        s = m.predict_proba(X[te])[:, 1]
        a = stats_block(y[te], s)
        nv = stats_block(y[te][novel], s[novel])
        out["temporal"][name] = {"all": a, "novel_scaffold": nv}
        print(f"\n  {name}")
        print(f"    all v2 additions   n={a['n']:4d}  ROC-AUC {a['roc_auc']:.3f} "
              f"[{a['ci95'][0]:.3f}, {a['ci95'][1]:.3f}]  p={a['mannwhitney_p']:.2e}")
        if nv.get("scoreable"):
            print(f"    new scaffolds only n={nv['n']:4d}  ROC-AUC {nv['roc_auc']:.3f} "
                  f"[{nv['ci95'][0]:.3f}, {nv['ci95'][1]:.3f}]  p={nv['mannwhitney_p']:.2e}")

    # ------------------------------------- 5. second read: scaffold CV
    print(f"\n{'=' * 78}\nSCAFFOLD CV (10 repeats) - a second, independent read"
          f"\n{'=' * 78}")
    for name, X in variants.items():
        aucs = []
        for seed in range(10):
            sp = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed)
            a, b = next(sp.split(X, y, groups=groups))
            if len(set(y[b])) < 2:
                continue
            m = rf().fit(X[a], y[a])
            aucs.append(roc_auc_score(y[b], m.predict_proba(X[b])[:, 1]))
        out["scaffold_cv"][name] = {"roc_auc": round(float(np.mean(aucs)), 4),
                                    "sd": round(float(np.std(aucs)), 4),
                                    "repeats": len(aucs)}
        print(f"  {name:28s} ROC-AUC {np.mean(aucs):.3f} +/- {np.std(aucs):.3f}")

    (RESULTS / "hepatic_stacked_experiment.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")

    base = out["temporal"]["structure_only"]["all"]["roc_auc"]
    plus = out["temporal"]["structure_plus_mechanism"]["all"]["roc_auc"]
    print(f"\n{'=' * 78}")
    print(f"  structure only              {base:.3f}")
    print(f"  structure + mechanism       {plus:.3f}   ({plus - base:+.3f})")
    print(f"  VERDICT: mechanism features "
          f"{'HELP' if plus > base + 0.02 else 'DO NOT HELP'} on the "
          f"pre-registered holdout")
    print(f"{'=' * 78}")
    print(f"wrote {RESULTS / 'hepatic_stacked_experiment.json'}")


if __name__ == "__main__":
    main()
