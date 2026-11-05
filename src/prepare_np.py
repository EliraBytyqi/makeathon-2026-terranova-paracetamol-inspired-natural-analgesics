"""Standardize and profile the plant-derived compound library.

Collapses LOTUS compound-organism pairs into one row per unique compound,
carrying the full list of source plants (that many-to-many relation becomes
the CONTAINS edges of the knowledge graph), then applies the identical
standardization used on the ChEMBL side, computes physicochemical properties
and runs the reactive-metabolite alert set.

Output: data/processed/np_library.csv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

from safety import physchem, profile
from standardize import standardize

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
PROC = ROOT / "data" / "processed"
PROC.mkdir(parents=True, exist_ok=True)


def main():
    src = INTERIM / "np_library_raw.csv"
    df = pd.read_csv(src, low_memory=False)
    print(f"loaded {len(df):,} compound-organism pairs")

    # ---- standardize unique structures ---------------------------------
    uniq = df["structure_smiles"].dropna().unique()
    print(f"{len(uniq):,} unique raw SMILES -> standardizing...")
    smap, reasons = {}, {}
    for i, s in enumerate(uniq):
        if i and i % 25_000 == 0:
            print(f"   {i:,}/{len(uniq):,}", flush=True)
        r = standardize(s)
        if r.ok:
            smap[s] = (r.smiles, r.inchikey)
        else:
            reasons[r.reason] = reasons.get(r.reason, 0) + 1

    print(f"standardized OK: {len(smap):,}")
    print("rejected:")
    for k, v in sorted(reasons.items(), key=lambda x: -x[1])[:12]:
        print(f"   {k:38s} {v:,}")

    df["std_smiles"] = df["structure_smiles"].map(lambda s: smap.get(s, (None,))[0])
    df["std_inchikey"] = df["structure_smiles"].map(
        lambda s: smap.get(s, (None, None))[1] if s in smap else None)
    df = df[df["std_inchikey"].notna()]
    print(f"pairs surviving standardization: {len(df):,}")

    # ---- collapse to one row per compound ------------------------------
    print("\ncollapsing to unique compounds...")

    def joinset(s, limit=None):
        vals = sorted({str(x) for x in s.dropna() if str(x).strip()})
        if limit and len(vals) > limit:
            return "|".join(vals[:limit]) + f"|...(+{len(vals)-limit})"
        return "|".join(vals)

    grouped = df.groupby("std_inchikey").agg(
        std_smiles=("std_smiles", "first"),
        compound_name=("structure_nameTraditional", "first"),
        pubchem_cid=("structure_cid", "first"),
        npclass_pathway=("structure_taxonomy_npclassifier_01pathway", "first"),
        npclass_superclass=("structure_taxonomy_npclassifier_02superclass", "first"),
        npclass_class=("structure_taxonomy_npclassifier_03class", "first"),
        classyfire_class=("structure_taxonomy_classyfire_03class", "first"),
        source_species=("organism_taxonomy_09species", lambda s: joinset(s, 12)),
        source_genera=("organism_taxonomy_08genus", lambda s: joinset(s, 12)),
        source_families=("organism_taxonomy_06family", lambda s: joinset(s, 8)),
        n_source_organisms=("organism_name", "nunique"),
        reference_dois=("reference_doi", lambda s: joinset(s, 5)),
        n_references=("reference_doi", "nunique"),
        tier=("tier", lambda s: "pain_seed" if (s == "pain_seed").any() else "plant_background"),
        seed_match_level=("seed_match_level", lambda s: joinset(s)),
        traditional_pain_use=("traditional_pain_use", lambda s: joinset(s, 3)),
        relevant_target_hypothesis=("relevant_target_hypothesis", lambda s: joinset(s, 3)),
        structure_wikidata=("structure_wikidata", "first"),
    ).reset_index()

    print(f"unique compounds: {len(grouped):,}")
    print(grouped["tier"].value_counts().to_string())

    # ---- physicochemical + safety profiling -----------------------------
    print("\nprofiling physicochemistry and reactive-metabolite alerts...")
    recs = []
    for i, (ikey, smi) in enumerate(
            zip(grouped["std_inchikey"], grouped["std_smiles"])):
        if i and i % 25_000 == 0:
            print(f"   {i:,}/{len(grouped):,}", flush=True)
        mol = Chem.MolFromSmiles(smi)
        pc = physchem(mol)
        sp = profile(mol)
        pc.update({
            "n_reactive_alerts": sp.n_alerts,
            "max_alert_severity": sp.max_severity,
            "alert_severity_sum": sp.severity_sum,
            "reactive_alerts": "|".join(sp.alerts),
            "n_pains": len(sp.pains),
            "pains_flags": "|".join(sorted(set(sp.pains))[:3]),
            "n_brenk_nih": len(sp.brenk_nih),
        })
        recs.append(pc)

    prof = pd.DataFrame(recs)
    out = pd.concat([grouped.reset_index(drop=True), prof], axis=1)
    path = PROC / "np_library.csv"
    out.to_csv(path, index=False)

    print(f"\nwrote {path}  shape={out.shape}")
    print("\n--- reactive alert prevalence in plant chemical space ---")
    print(f"compounds with >=1 alert : {(out['n_reactive_alerts']>0).mean():.1%}")
    print(f"compounds with severity 3: {(out['max_alert_severity']==3).mean():.1%}")
    print(f"PAINS-flagged            : {(out['n_pains']>0).mean():.1%}")
    print(f"CNS-like physchem        : {out['cns_like'].mean():.1%}")
    print("\nmost common alerts:")
    exploded = out["reactive_alerts"].dropna().str.split("|").explode()
    print(exploded[exploded != ""].value_counts().head(12).to_string())


if __name__ == "__main__":
    main()
