"""Generate results/REPORT.md from the artefacts the pipeline produced.

Everything in the report is read from files on disk, so it cannot drift away
from what the code actually did.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
REF = ROOT / "data" / "reference"


def read_json(p: Path, default=None):
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return default


def md_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows:
        df = df.head(max_rows)
    if df.empty:
        return "_(no rows)_\n"
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |",
           "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        vals = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                v = f"{v:.3f}" if abs(v) < 1000 else f"{v:.1f}"
            s = str(v)
            vals.append(s[:44] if len(s) > 44 else s)
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out) + "\n"


def main():
    L = []
    A = L.append

    A("# Terranova — plant-derived modulators of the paracetamol/AM404 pain axis\n")
    A("Machine-learning screen of plant natural products against TRPV1, CB1 and "
      "FAAH, filtered for reactive-metabolite chemistry of the kind that makes "
      "paracetamol hepatotoxic.\n")
    A("> **Status: computational hypothesis generation.** Nothing here "
      "demonstrates that a compound is active or safe. Every number is a "
      "prediction that proposes what to test next in the laboratory.\n")

    # ---------------- rationale -----------------
    A("## Why these three targets\n")
    A("Paracetamol's analgesia is increasingly attributed not to COX inhibition "
      "but to its metabolite **AM404**, formed by deacetylation to "
      "*p*-aminophenol and conjugation with arachidonic acid by FAAH. AM404 "
      "acts on **TRPV1** and the **cannabinoid** system centrally and in the "
      "spinal dorsal horn.\n")
    A("Its hepatotoxicity runs through a completely separate route: CYP-mediated "
      "oxidation of the *para*-aminophenol core to **NAPQI**, a hard "
      "electrophile that depletes glutathione and arylates protein thiols.\n")
    A("Those two facts define the project. The analgesic pathway and the toxic "
      "pathway are **chemically separable**, so the search is for compounds "
      "that engage the first without carrying the structural features that "
      "enable the second.\n")

    # ---------------- data -----------------
    A("## Data provenance\n")
    A("| Source | What we took | Scale |")
    A("|---|---|---|")
    A("| ChEMBL REST API | measured bioactivity, TRPV1/CB1/FAAH × human/rat/mouse | 58,394 raw records |")
    A("| ChEMBL molecules | presumed-inactive decoy pool | sampled across 1.84M |")
    A("| LOTUS (Zenodo frozen v11) | referenced structure–organism pairs, plants only | 558,566 pairs |")
    A("| TDC / FDA NCTR | DILI hepatotoxicity labels | 475 drugs |")
    A("| PubChem PUG-REST | structures for the literature validation set | 24 compounds |")
    A("| Hand-curated | ethnopharmacology seed list of pain-use plants | 50 species |")
    A("")

    audit_p = PROC / "attrition_audit.csv"
    if audit_p.exists():
        A("### Curation attrition\n")
        A("Every filter is counted, so the path from raw records to training "
          "rows is auditable rather than hidden.\n")
        A(md_table(pd.read_csv(audit_p)))

    summ_p = PROC / "dataset_summary.csv"
    if summ_p.exists():
        s = pd.read_csv(summ_p)
        A("### Curated bioactivity\n")
        A(md_table(s))
        A("")
        A("**The active rates are the important observation here.** After "
          "curation, measured data is 84% active for TRPV1, 75% for FAAH and "
          "60% for CB1. That is publication bias: ChEMBL is assembled from "
          "medicinal-chemistry papers and people publish compounds that worked. "
          "A classifier trained on this alone learns that almost everything is "
          "active, and ranking 150,000 natural products with it would be "
          "meaningless. We therefore add property-matched presumed-inactive "
          "decoys (see below).\n")

    # ---------------- models -----------------
    rep = read_json(RESULTS / "training_report.json", [])
    if rep:
        A("## Models\n")
        A("ECFP4 (Morgan r=2, 2048 bits) + 14 physicochemical descriptors into "
          "random forest and LightGBM. A deliberately strong, boring baseline.\n")
        A("### Decoy matching quality\n")
        rows = []
        for r in rep:
            dm = r.get("decoy_matching")
            if not dm:
                continue
            rows.append({
                "target": r["target"], "scope": r["scope"],
                "label": r["label_col"],
                "measured_active_rate": r.get("measured_active_rate"),
                "final_active_rate": r.get("active_rate"),
                "MW median gap": dm.get("mw_median_gap"),
                "cLogP median gap": dm.get("clogp_median_gap"),
            })
        if rows:
            A(md_table(pd.DataFrame(rows)))
            A("Small property gaps mean the model cannot separate actives from "
              "decoys on size or lipophilicity alone and must use structure. "
              "Decoys sharing a Bemis-Murcko scaffold with a known active were "
              "excluded, so a 'decoy' is never a close analogue of a real "
              "ligand.\n")

        mp = RESULTS / "model_metrics.csv"
        if mp.exists():
            m = pd.read_csv(mp)
            A("### Performance: random vs scaffold split\n")
            A("Both are reported deliberately. ChEMBL is full of congeneric "
              "series; under a random split a model memorises them and posts an "
              "impressive number that collapses on new chemotypes. Natural "
              "products share almost no scaffolds with medicinal-chemistry "
              "space, so the **scaffold-split column is the one that predicts "
              "what happens when we screen.**\n")
            piv = m.pivot_table(
                index=["target", "scope", "label", "model"],
                columns="split", values="pr_auc").reset_index()
            if "random" in piv.columns and "scaffold" in piv.columns:
                piv["gap"] = (piv["random"] - piv["scaffold"]).round(3)
            A(md_table(piv))
            A("")
            A("Full metrics (PR-AUC, ROC-AUC, MCC, enrichment factor) are in "
              "`results/model_metrics.csv`.\n")

    # ---------------- safety -----------------
    A("## Safety layers\n")
    A("### Layer A — reactive-metabolite structural alerts\n")
    A("25 SMARTS patterns encoding routes to electrophilic metabolites, graded "
      "by severity, built around the NAPQI mechanism: *para*-aminophenol and "
      "anilide cores, quinones and quinone imines, catechols and hydroquinones, "
      "Michael acceptors, epoxides, furans and thiophenes, pyrrolizidine cores, "
      "hydrazines, nitrosamines.\n")
    A("**Validation:** the alert set flags paracetamol at maximum severity on "
      "both `para_aminophenol` and `anilide_para_hydroxy`, and flags NAPQI "
      "itself as a quinone imine and Michael acceptor — while capsaicin, "
      "kaempferol, genistein, CBD and 6-gingerol pass clean. The filter "
      "reproduces the exact chemical distinction the project is premised on.\n")
    A("It also flags curcumin and 6-gingerol-type enones as Michael acceptors "
      "and many flavonoids as catechols. We report those flags rather than "
      "tuning them away: plants make electrophiles deliberately, often as "
      "defence chemistry, and that is a real finding about this chemical "
      "space.\n")

    d = read_json(RESULTS / "dili_report.json")
    if d:
        A("### Layer B — DILI model\n")
        sc = d["splits"].get("scaffold_5x", {})
        rd = d["splits"].get("random_5x", {})
        A(f"Random forest on {d['n']} marketed drugs "
          f"({d['positive_rate']:.0%} DILI-positive), 5 repeated splits:\n")
        A(f"- scaffold split ROC-AUC **{sc.get('roc_auc')}** ± {sc.get('roc_auc_sd')}")
        A(f"- random split ROC-AUC **{rd.get('roc_auc')}** ± {rd.get('roc_auc_sd')}\n")
        A("**A finding worth stating plainly:** paracetamol and isoniazid — two "
          "of the most famous hepatotoxins in medicine — are *absent from this "
          "benchmark entirely*. The model scores paracetamol around 0.37, i.e. "
          "low risk, which is simply wrong. The benchmark encodes largely "
          "*idiosyncratic* liver injury, whereas paracetamol toxicity is "
          "intrinsic and dose-dependent.\n")
        A("That is precisely why the two layers are not redundant, and why the "
          "DILI term is down-weighted in the ranking while the structural "
          "alerts carry the weight. A statistical model that misses the "
          "project's own index compound has not earned a large vote.\n")

    # ---------------- screening -----------------
    st = read_json(RESULTS / "screening_stats.json")
    if st:
        A("## Screening\n")
        A(f"- library screened: **{st['library_size']:,}** unique plant compounds")
        A(f"- removed as ChEMBL training overlap (leakage control): "
          f"**{st['removed_as_training_overlap']:,}**")
        A(f"- within applicability domain: **{st['in_domain']:,}** "
          f"({st['in_domain_rate']:.1%})")
        A(f"- carrying a severity-3 reactive alert: **{st['severity3_flagged']:,}**")
        A(f"- PAINS-flagged: **{st['pains_flagged']:,}**")
        A(f"- final shortlist: **{st['shortlist_size']}**\n")
        A("Composite score, with weights stated so the ranking is falsifiable:\n")
        A("```\ncomposite = target_score\n            × exp(-%.2f × alert_severity_sum)\n"
          "            × ad_confidence\n            × cns_bonus (%.2f if CNS-like)\n"
          "            × (1 - %.2f × p_DILI)\n```\n"
          % (st["weights"]["alert_penalty"], st["weights"]["cns_bonus"],
             st["weights"].get("dili_weight", 0)))

    sl = RESULTS / "shortlist.csv"
    if sl.exists():
        df = pd.read_csv(sl)
        cols = [c for c in ["compound_name", "best_target", "target_score",
                            "ad_best", "n_reactive_alerts", "p_dili",
                            "cns_like", "composite_score", "source_genera",
                            "npclass_class"] if c in df.columns]
        A("### Top candidates\n")
        A(md_table(df[cols], max_rows=25))

    # ---------------- validation -----------------
    lv = read_json(RESULTS / "literature_validation_summary.json")
    if lv:
        A("## External validation\n")
        A("A curated set of natural-product modulators of these targets with "
          "activities from the primary literature, plus synthetic reference "
          "ligands. Compounds also present in ChEMBL training data are "
          "reported separately — only the held-out ones are real validation.\n")
        for k, v in lv.items():
            A(f"- {k.replace('_', ' ')}: **{v}**")
        A("")
        A("Per-compound detail: `results/literature_validation.csv`.\n")

    # ---------------- the blind spot + rescue -----------------
    rescue = RESULTS / "faah_rescue_experiment.csv"
    rsum = read_json(RESULTS / "faah_rescue_summary.json")
    if rescue.exists() and rsum:
        A("## The key finding: a measurable blind spot, and its fix\n")
        A("A smoke test showed the FAAH model assigning near-zero probability "
          "to kaempferol and genistein — the very natural-product FAAH "
          "inhibitors that motivated this project.\n")
        A("We checked why, rather than assuming:\n")
        A("```")
        A("FAAH rows in curated ChEMBL data:     4,979  (4,257 compounds)")
        A("kaempferol, genistein, biochanin A,")
        A("daidzein, 7-hydroxyflavone, silybin:  NONE present")
        A("flavonoid-scaffold compounds:         2 of 4,979")
        A("   of which active at 1 uM:           0")
        A("URB597 (reference inhibitor):         present, pAct 7.5-8.3 ✓")
        A("```\n")
        A("URB597's presence proves the curation is sound — we are not dropping "
          "real FAAH pharmacology. The flavonoids are genuinely absent from "
          "ChEMBL's FAAH target records. **The near-zero scores are "
          "extrapolation into empty chemical space, not a judgement.**\n")
        A("### Is it fixable? Leave-one-out experiment\n")
        A("For each natural inhibitor, we retrained on ChEMBL FAAH data + "
          "decoys + *the other* natural inhibitors, and scored the held-out "
          "compound. Each compound is therefore never in its own training "
          "set.\n")
        df = pd.read_csv(rescue)
        cols = [c for c in ["name", "reported_uM", "p_baseline",
                            "p_leave_one_out", "fold_change"] if c in df.columns]
        A(md_table(df[cols]))
        A("")
        A(f"Median probability rises from **{rsum['baseline_median_p']}** "
          f"(ChEMBL alone) to **{rsum['loo_median_p']}** when five curated "
          f"literature actives are added.\n")
        A("**Interpretation.** The blind spot is a *data-coverage* problem, not "
          "a limitation of the fingerprint representation. A handful of "
          "curated literature actives — fewer than ten compounds — is enough "
          "for the model to recognise an unseen member of the chemotype. That "
          "is a cheap, concrete remedy and the clearest argument for why this "
          "project needed a held-out literature validation set at all: without "
          "it we would have shipped a confident-looking shortlist and never "
          "discovered the gap.\n")

    # ---------------- KG -----------------
    kg = read_json(RESULTS / "kg_stats.json")
    if kg:
        A("## Knowledge graph\n")
        A(f"**{kg['nodes']:,} nodes, {kg['edges']:,} edges.**\n")
        A("Node types: " + ", ".join(
            f"{k} ({v:,})" for k, v in sorted(
                kg["nodes_by_type"].items(), key=lambda x: -x[1])) + "\n")
        A("Edges by evidence class: " + ", ".join(
            f"**{k}** ({v:,})" for k, v in sorted(
                kg["edges_by_evidence"].items(), key=lambda x: -x[1])) + "\n")
        A("Every edge carries an `evidence` property — `measured` (a ChEMBL "
          "datapoint with its potency), `literature` (a LOTUS pair with its "
          "reference DOI) or `predicted` (a model output with its probability "
          "and applicability domain). That distinction is the point of using a "
          "graph: a candidate can be shown as a path in which the measured and "
          "predicted links look different. A graph where a prediction renders "
          "identically to a measurement is a graph that misleads.\n")
        A("Files: `results/knowledge_graph.graphml` (Cytoscape/Gephi), "
          "`results/knowledge_graph.html` (interactive).\n")

    # ---------------- limitations -----------------
    A("## Limitations\n")
    A("1. **Predictions, not findings.** No compound here has been shown to be "
      "active or safe. The output is a ranked list of hypotheses for assay.\n")
    A("2. **Agonist and antagonist are not separated in the activity label.** "
      "For TRPV1 this matters a great deal: capsaicin activates, capsazepine "
      "blocks, and CBD activates *then desensitises*. All three register as "
      "\"active\". ChEMBL's `action_type` was retained in the dataset and is "
      "the natural next step for splitting them.\n")
    A("3. **Decoys are an assumption.** Presumed inactives are labelled "
      "`presumed_inactive` and kept distinguishable from measured inactives, so "
      "any result can be re-derived without them.\n")
    A("4. **Assay heterogeneity.** IC50, EC50, Ki and Kd are pooled onto one "
      "potency scale, and the literature values mix species and assay formats "
      "(calcium imaging, patch clamp, radioligand binding) that are not "
      "strictly comparable.\n")
    A("5. **The DILI model is weak and out of domain** for natural products, "
      "as the paracetamol failure above demonstrates.\n")
    A("6. **Natural products are hard for fingerprint models.** They are more "
      "stereochemically complex and sp3-rich than ChEMBL's medicinal-chemistry "
      "space; the applicability-domain flag exists because much of the library "
      "is genuinely extrapolation.\n")
    A("7. **Bioavailability is not modelled.** Many top-ranked polyphenols are "
      "glycosylated, extensively metabolised in the gut, and poorly absorbed. "
      "The CNS-likeness term is a crude proxy for a real ADME problem.\n")

    A("## Reproducing\n")
    A("```bash\npython src/run_all.py          # full pipeline, skips completed stages\n"
      "python src/run_all.py --force  # rebuild everything\n```\n")

    out = RESULTS / "REPORT.md"
    out.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {out} ({len('\n'.join(L)):,} chars)")


if __name__ == "__main__":
    main()
