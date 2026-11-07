"""Build the experiment-ready candidate list: 50 DISTINCT molecules.

results/shortlist.csv takes the top 50 rows of the ranked screen, and that is
not the same thing as 50 compounds. The library stores stereoisomers,
tautomers and salt variants as separate rows, so the current top 50 contains
only 27 distinct chemical skeletons - 15 skeletons occupy 38 of the 50 slots.
Handed to a laboratory, that is 23 wasted assay wells.

This collapses the ranked screen to one row per InChIKey skeleton, keeping the
best-scoring representative of each, and takes the top 50 of those. It applies
the same hard requirements as the pipeline (in applicability domain, no
severity-3 reactive alert, no PAINS) and adds everything an experimentalist
needs to actually order and test the compound: the plant it comes from, the
literature DOIs, the predicted target with its applicability domain, the
hepatic flag, and the physicochemical profile.

Target diversity is reported but not enforced. The ranking is the ranking; if
one target dominates, that is a finding about the library rather than
something to correct by quota.

    python scripts/make_candidate_list.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
N = 50

COLS = [
    "rank", "compound_name", "best_target", "target_score", "ad_best",
    "p_hepatic", "max_alert_severity", "n_pains", "cns_like",
    "composite_score", "source_species", "source_genera", "source_families",
    "npclass_class", "npclass_superclass", "tier", "traditional_pain_use",
    "mw", "clogp", "tpsa", "qed", "lipinski_violations",
    "n_references", "reference_dois", "pubchem_cid", "std_inchikey",
    "std_smiles", "n_variants_collapsed", "variant_names",
]


def main() -> None:
    df = pd.read_csv(RESULTS / "screening_all.csv.gz", low_memory=False)
    print(f"ranked screen: {len(df):,} rows")

    # the pipeline's hard requirements, restated here so this file is readable
    # on its own rather than inheriting them silently
    ok = df[(df["ad_best"] >= 0.25)
            & (df["max_alert_severity"].fillna(0) < 3)
            & (df["n_pains"].fillna(0) == 0)].copy()
    print(f"  passing the hard requirements (AD>=0.25, no severity-3, no PAINS): "
          f"{len(ok):,}")

    ok["skeleton"] = ok["std_inchikey"].str[:14]
    ok = ok.sort_values("composite_score", ascending=False)

    # one row per distinct chemical skeleton, keeping the best representative
    grp = ok.groupby("skeleton", sort=False)
    best = grp.head(1).copy()
    counts = grp.size()
    names = grp["compound_name"].apply(
        lambda s: "; ".join(dict.fromkeys(str(x) for x in s))[:300])
    best["n_variants_collapsed"] = best["skeleton"].map(counts).astype(int)
    best["variant_names"] = best["skeleton"].map(names)
    print(f"  distinct skeletons: {len(best):,}")

    top = best.head(N).reset_index(drop=True)
    top.insert(0, "rank", range(1, len(top) + 1))

    cols = [c for c in COLS if c in top.columns]
    out = top[cols]
    out.to_csv(RESULTS / "candidates_top50.csv", index=False)

    # a compact JSON the dashboard and the deck can both read
    payload = {
        "n": int(len(out)),
        "selected_from": int(len(ok)),
        "distinct_skeletons_available": int(len(best)),
        "criteria": ("composite rank, one row per InChIKey skeleton, "
                     "AD >= 0.25, max alert severity < 3, zero PAINS"),
        "by_target": out["best_target"].value_counts().to_dict(),
        "compounds": json.loads(out.to_json(orient="records")),
    }
    (RESULTS / "candidates_top50.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")

    print(f"\n{'=' * 96}\nTOP {len(out)} DISTINCT CANDIDATES\n{'=' * 96}")
    print(f"  by target: {payload['by_target']}")
    n_fam = out["source_families"].astype(str).str.split("|").explode().nunique()
    print(f"  distinct plant families represented: {n_fam}")
    print(f"  with a traditional pain-use record: "
          f"{int((out['tier'] == 'pain_seed').sum())}")
    print(f"  median predicted hepatic liability: {out['p_hepatic'].median():.3f}")
    print()
    show = out.head(20).copy()
    show["compound_name"] = show["compound_name"].astype(str).str.slice(0, 36)
    show["source_genera"] = show["source_genera"].astype(str).str.slice(0, 22)
    print(show[["rank", "compound_name", "best_target", "target_score",
                "ad_best", "p_hepatic", "source_genera",
                "n_variants_collapsed"]].to_string(index=False))
    print(f"\nwrote {RESULTS / 'candidates_top50.csv'}")
    print(f"wrote {RESULTS / 'candidates_top50.json'}")


if __name__ == "__main__":
    main()
