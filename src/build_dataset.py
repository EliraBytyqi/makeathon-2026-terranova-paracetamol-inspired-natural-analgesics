"""Turn raw ChEMBL activity records into modelling-ready datasets.

Every filtering decision is applied here, in one place, and counted, so the
attrition from raw records to training rows is auditable. Hidden curation is
how cheminformatics papers become irreproducible.

Key decisions and why
---------------------

ASSAY TYPES. We keep 'B' (binding) and 'F' (functional). For TRPV1 this is not
optional: it is a ligand-gated ion channel whose potency data is dominated by
calcium-flux and patch-clamp measurements registered as functional assays.
Restricting to 'B', the habitual choice for enzyme/GPCR panels, would discard
most TRPV1 pharmacology and leave a model trained almost entirely on
radioligand displacement.

ACTIVITY TYPES. IC50, EC50, Ki, Kd only. Percent-inhibition and other
single-concentration readouts are not convertible to a potency and are
dropped.

UNITS. nM only, after ChEMBL's own standardisation. pActivity = 9 - log10(nM).

CENSORED VALUES. A record of ">10000 nM" is not a measurement, but it is
strong evidence of inactivity. We use '=' records for both classes, and admit
'>' / '>=' records as inactives only. '<' records are discarded: they bound
potency from above but cannot confirm an active threshold.

AGGREGATION. Median pActivity per (compound, target, organism group). ChEMBL
contains the same compound measured many times across labs; the median is
robust to the occasional order-of-magnitude outlier.

LABELS. active <= 1 uM (pActivity >= 6); inactive >= 10 uM (pActivity <= 5).
The 5-6 band is held out rather than forced into a class - compounds there are
genuinely ambiguous and training on them teaches the model noise.

FAAH IS DIFFERENT. The natural-product FAAH inhibitors that motivate this
project (kaempferol, genistein, biochanin A, silybin, licarin A) sit at
roughly 1-15 uM. A 1 uM activity cut would classify essentially all of them as
inactive, and the model would be optimised to reject exactly the chemistry we
are looking for. For FAAH we therefore also emit a regression target on
pActivity and a 10 uM classification variant, and report both.

ORGANISM. Human is primary. Rat and mouse orthologs are kept in a separate
column and included in an 'all-species' variant, never silently pooled into
the human set. TRPV1 in particular has well documented species differences in
ligand sensitivity.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from standardize import standardize

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROC = ROOT / "data" / "processed"
PROC.mkdir(parents=True, exist_ok=True)
INTERIM.mkdir(parents=True, exist_ok=True)

KEEP_ASSAY_TYPES = {"B", "F"}
KEEP_ACT_TYPES = {"IC50", "EC50", "Ki", "Kd"}

ACTIVE_CUT = 6.0    # pActivity >= 6  => <= 1 uM
INACTIVE_CUT = 5.0  # pActivity <= 5  => >= 10 uM
FAAH_ACTIVE_CUT = 5.0  # <= 10 uM, see module docstring

AUDIT: list[tuple[str, int]] = []


def step(label: str, df: pd.DataFrame) -> pd.DataFrame:
    AUDIT.append((label, len(df)))
    print(f"  {label:48s} {len(df):>8,}")
    return df


def load_raw() -> pd.DataFrame:
    direct = RAW / "chembl_activities_raw.csv"
    if direct.exists():
        return pd.read_csv(direct, low_memory=False)
    # fall back to per-target checkpoints if the combined write never happened
    parts = sorted((RAW / "chembl_by_target").glob("*.csv"))
    if not parts:
        raise SystemExit("no ChEMBL data found - run fetch_chembl.py first")
    print(f"  (assembling from {len(parts)} checkpoint files)")
    return pd.concat([pd.read_csv(p, low_memory=False) for p in parts],
                     ignore_index=True)


def main():
    print("\n=== ATTRITION AUDIT ===")
    df = load_raw()
    step("raw activity records", df)

    df = step("has SMILES", df[df["canonical_smiles"].notna()])
    df = step("assay_type in {B,F}", df[df["assay_type"].isin(KEEP_ASSAY_TYPES)])
    df = step("standard_type in {IC50,EC50,Ki,Kd}",
              df[df["standard_type"].isin(KEEP_ACT_TYPES)])
    df = step("units == nM", df[df["standard_units"] == "nM"])
    df = step("has numeric value",
              df[pd.to_numeric(df["standard_value"], errors="coerce").notna()])

    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = step("value > 0", df[df["standard_value"] > 0])
    # ChEMBL flags records its own curation distrusts
    df = step("no data_validity_comment",
              df[df["data_validity_comment"].isna()])
    df = step("relation in {=, >, >=}",
              df[df["standard_relation"].isin(["=", ">", ">="])])

    df["pactivity"] = 9.0 - np.log10(df["standard_value"].astype(float))
    df = step("pactivity in plausible range [2, 12]",
              df[(df["pactivity"] >= 2) & (df["pactivity"] <= 12)])

    # ---- structure standardization -------------------------------------
    print("\n  standardizing unique structures...")
    uniq = df["canonical_smiles"].dropna().unique()
    print(f"  {len(uniq):,} unique raw SMILES")
    smap = {}
    reasons = {}
    for s in uniq:
        r = standardize(s)
        if r.ok:
            smap[s] = (r.smiles, r.inchikey)
        else:
            reasons[r.reason] = reasons.get(r.reason, 0) + 1
    print(f"  standardized OK: {len(smap):,}")
    if reasons:
        print("  rejected:")
        for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
            print(f"      {k:36s} {v:,}")

    df["std_smiles"] = df["canonical_smiles"].map(
        lambda s: smap.get(s, (None, None))[0])
    df["inchikey"] = df["canonical_smiles"].map(
        lambda s: smap.get(s, (None, None))[1])
    df = step("structure standardized", df[df["inchikey"].notna()])

    df.to_csv(INTERIM / "chembl_clean_records.csv", index=False)

    # ---- aggregate to one row per compound x target x organism ---------
    print("\n=== AGGREGATION ===")
    rows = []
    gcols = ["target_name", "target_organism_group", "inchikey"]
    for (target, org, ikey), g in df.groupby(gcols):
        exact = g[g["standard_relation"] == "="]
        censored = g[g["standard_relation"].isin([">", ">="])]

        if len(exact):
            p = float(exact["pactivity"].median())
            relation = "="
            n = len(exact)
            spread = float(exact["pactivity"].max() - exact["pactivity"].min())
        else:
            # only censored data: treat as an upper bound on potency
            p = float(censored["pactivity"].median())
            relation = ">"
            n = len(censored)
            spread = float(censored["pactivity"].max() - censored["pactivity"].min())

        actions = [a for a in g["action_type"].dropna().unique()]
        rows.append({
            "inchikey": ikey,
            "std_smiles": g["std_smiles"].iloc[0],
            "target_name": target,
            "organism_group": org,
            "pactivity": p,
            "relation": relation,
            "n_measurements": n,
            "pactivity_spread": spread,
            "action_type": "|".join(sorted(actions)) if actions else "",
            "assay_types": "|".join(sorted(g["assay_type"].dropna().unique())),
            "std_types": "|".join(sorted(g["standard_type"].dropna().unique())),
            "n_documents": g["document_chembl_id"].nunique(),
            "chembl_ids": "|".join(sorted(g["molecule_chembl_id"].dropna().unique())[:3]),
        })

    agg = pd.DataFrame(rows)
    print(f"  aggregated rows (compound x target x organism): {len(agg):,}")

    # ---- labelling ------------------------------------------------------
    def label_row(r):
        cut_act = ACTIVE_CUT
        cut_inact = INACTIVE_CUT
        if r["relation"] == ">":
            # censored: only usable as evidence of inactivity
            return 0 if r["pactivity"] <= cut_inact else np.nan
        if r["pactivity"] >= cut_act:
            return 1
        if r["pactivity"] <= cut_inact:
            return 0
        return np.nan  # ambiguous 1-10 uM band

    def label_row_10um(r):
        if r["relation"] == ">":
            return 0 if r["pactivity"] <= FAAH_ACTIVE_CUT else np.nan
        return 1 if r["pactivity"] >= FAAH_ACTIVE_CUT else 0

    agg["label"] = agg.apply(label_row, axis=1)
    agg["label_10um"] = agg.apply(label_row_10um, axis=1)

    agg.to_csv(PROC / "bioactivity_aggregated.csv", index=False)

    # ---- report ---------------------------------------------------------
    print("\n=== DATASET SUMMARY ===")
    summary = []
    for target in sorted(agg["target_name"].unique()):
        for org in ["human", "all"]:
            sub = agg[agg["target_name"] == target]
            if org == "human":
                sub = sub[sub["organism_group"] == "human"]
            lab = sub.dropna(subset=["label"])
            n_act = int((lab["label"] == 1).sum())
            n_inact = int((lab["label"] == 0).sum())
            n_amb = int(sub["label"].isna().sum())
            lab10 = sub.dropna(subset=["label_10um"])
            summary.append({
                "target": target,
                "scope": org,
                "compounds": sub["inchikey"].nunique(),
                "active_1uM": n_act,
                "inactive_10uM": n_inact,
                "ambiguous": n_amb,
                "active_rate": round(n_act / max(n_act + n_inact, 1), 3),
                "active_10uM_cut": int((lab10["label_10um"] == 1).sum()),
            })
    srep = pd.DataFrame(summary)
    print(srep.to_string(index=False))
    srep.to_csv(PROC / "dataset_summary.csv", index=False)

    audit = pd.DataFrame(AUDIT, columns=["step", "rows"])
    audit.to_csv(PROC / "attrition_audit.csv", index=False)
    print(f"\nwrote {PROC/'bioactivity_aggregated.csv'}")
    print(f"wrote {PROC/'dataset_summary.csv'}")
    print(f"wrote {PROC/'attrition_audit.csv'}")


if __name__ == "__main__":
    main()
