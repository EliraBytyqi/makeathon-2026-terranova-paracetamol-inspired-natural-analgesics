# Terranova — plant-derived modulators of the paracetamol/AM404 pain axis

Machine-learning screen of plant natural products against TRPV1, CB1 and FAAH, filtered for reactive-metabolite chemistry of the kind that makes paracetamol hepatotoxic.

> **Status: computational hypothesis generation.** Nothing here demonstrates that a compound is active or safe. Every number is a prediction that proposes what to test next in the laboratory.

## Why these three targets

Paracetamol's analgesia is increasingly attributed not to COX inhibition but to its metabolite **AM404**, formed by deacetylation to *p*-aminophenol and conjugation with arachidonic acid by FAAH. AM404 acts on **TRPV1** and the **cannabinoid** system centrally and in the spinal dorsal horn.

Its hepatotoxicity runs through a completely separate route: CYP-mediated oxidation of the *para*-aminophenol core to **NAPQI**, a hard electrophile that depletes glutathione and arylates protein thiols.

Those two facts define the project. The analgesic pathway and the toxic pathway are **chemically separable**, so the search is for compounds that engage the first without carrying the structural features that enable the second.

## Data provenance

| Source | What we took | Scale |
|---|---|---|
| ChEMBL REST API | measured bioactivity, TRPV1/CB1/FAAH × human/rat/mouse | 58,394 raw records |
| ChEMBL molecules | presumed-inactive decoy pool | sampled across 1.84M |
| LOTUS (Zenodo frozen v11) | referenced structure–organism pairs, plants only | 558,566 pairs |
| TDC / FDA NCTR | DILI hepatotoxicity labels | 475 drugs |
| PubChem PUG-REST | structures for the literature validation set | 24 compounds |
| Hand-curated | ethnopharmacology seed list of pain-use plants | 50 species |

### Curation attrition

Every filter is counted, so the path from raw records to training rows is auditable rather than hidden.

| step | rows |
|---|---|
| raw activity records | 58394 |
| has SMILES | 58287 |
| assay_type in {B,F} | 57214 |
| standard_type in {IC50,EC50,Ki,Kd} | 35510 |
| units == nM | 32574 |
| has numeric value | 32386 |
| value > 0 | 32386 |
| no data_validity_comment | 31854 |
| relation in {=, >, >=} | 31510 |
| pactivity in plausible range [2, 12] | 31510 |
| structure standardized | 31271 |

### Curated bioactivity

| target | scope | compounds | active_1uM | inactive_10uM | ambiguous | active_rate | active_10uM_cut |
|---|---|---|---|---|---|---|---|
| CB1 | human | 9034 | 4311 | 2909 | 1814 | 0.597 | 5601 |
| CB1 | all | 10953 | 6067 | 3197 | 2289 | 0.655 | 7653 |
| FAAH | human | 2271 | 1352 | 448 | 471 | 0.751 | 1766 |
| FAAH | all | 4257 | 2893 | 986 | 1100 | 0.746 | 3673 |
| TRPV1 | human | 4172 | 3044 | 561 | 567 | 0.844 | 3511 |
| TRPV1 | all | 5042 | 4185 | 639 | 926 | 0.868 | 4915 |


**The active rates are the important observation here.** After curation, measured data is 84% active for TRPV1, 75% for FAAH and 60% for CB1. That is publication bias: ChEMBL is assembled from medicinal-chemistry papers and people publish compounds that worked. A classifier trained on this alone learns that almost everything is active, and ranking 150,000 natural products with it would be meaningless. We therefore add property-matched presumed-inactive decoys (see below).

## Models

ECFP4 (Morgan r=2, 2048 bits) + 14 physicochemical descriptors into random forest and LightGBM. A deliberately strong, boring baseline.

### Decoy matching quality

| target | scope | label | measured_active_rate | final_active_rate | MW median gap | cLogP median gap |
|---|---|---|---|---|---|---|
| CB1 | human | label | 0.597 | 0.377 | 2.400 | 0.180 |
| CB1 | all | label | 0.654 | 0.400 | 1.900 | 0.180 |
| FAAH | human | label | 0.751 | 0.431 | 1.100 | 0.010 |
| FAAH | human | label_10um | 0.803 | 0.450 | 1.700 | 0.080 |
| FAAH | all | label | 0.715 | 0.420 | 0.600 | 0.060 |
| FAAH | all | label_10um | 0.770 | 0.440 | 3.100 | 0.100 |
| TRPV1 | human | label | 0.844 | 0.461 | 1.100 | 0.020 |
| TRPV1 | all | label | 0.858 | 0.465 | 1.000 | 0.030 |

Small property gaps mean the model cannot separate actives from decoys on size or lipophilicity alone and must use structure. Decoys sharing a Bemis-Murcko scaffold with a known active were excluded, so a 'decoy' is never a close analogue of a real ligand.

### Performance: random vs scaffold split

Both are reported deliberately. ChEMBL is full of congeneric series; under a random split a model memorises them and posts an impressive number that collapses on new chemotypes. Natural products share almost no scaffolds with medicinal-chemistry space, so the **scaffold-split column is the one that predicts what happens when we screen.**

| target | scope | label | model | random | scaffold | gap |
|---|---|---|---|---|---|---|
| CB1 | all | label | lightgbm | 0.983 | 0.957 | 0.025 |
| CB1 | all | label | random_forest | 0.978 | 0.945 | 0.033 |
| CB1 | human | label | lightgbm | 0.978 | 0.968 | 0.009 |
| CB1 | human | label | random_forest | 0.971 | 0.962 | 0.010 |
| FAAH | all | label | lightgbm | 0.975 | 0.939 | 0.036 |
| FAAH | all | label | random_forest | 0.970 | 0.926 | 0.044 |
| FAAH | all | label_10um | lightgbm | 0.968 | 0.942 | 0.026 |
| FAAH | all | label_10um | random_forest | 0.968 | 0.931 | 0.036 |
| FAAH | human | label | lightgbm | 0.973 | 0.922 | 0.050 |
| FAAH | human | label | random_forest | 0.968 | 0.919 | 0.049 |
| FAAH | human | label_10um | lightgbm | 0.944 | 0.960 | -0.016 |
| FAAH | human | label_10um | random_forest | 0.950 | 0.952 | -0.002 |
| TRPV1 | all | label | lightgbm | 0.992 | 0.975 | 0.017 |
| TRPV1 | all | label | random_forest | 0.990 | 0.968 | 0.023 |
| TRPV1 | human | label | lightgbm | 0.990 | 0.967 | 0.023 |
| TRPV1 | human | label | random_forest | 0.988 | 0.959 | 0.029 |


Full metrics (PR-AUC, ROC-AUC, MCC, enrichment factor) are in `results/model_metrics.csv`.

## Safety layers

### Layer A — reactive-metabolite structural alerts

25 SMARTS patterns encoding routes to electrophilic metabolites, graded by severity, built around the NAPQI mechanism: *para*-aminophenol and anilide cores, quinones and quinone imines, catechols and hydroquinones, Michael acceptors, epoxides, furans and thiophenes, pyrrolizidine cores, hydrazines, nitrosamines.

**Validation:** the alert set flags paracetamol at maximum severity on both `para_aminophenol` and `anilide_para_hydroxy`, and flags NAPQI itself as a quinone imine and Michael acceptor — while capsaicin, kaempferol, genistein, CBD and 6-gingerol pass clean. The filter reproduces the exact chemical distinction the project is premised on.

It also flags curcumin and 6-gingerol-type enones as Michael acceptors and many flavonoids as catechols. We report those flags rather than tuning them away: plants make electrophiles deliberately, often as defence chemistry, and that is a real finding about this chemical space.

### Layer B — DILI model

Random forest on 506 marketed drugs (39% DILI-positive), 5 repeated splits:

- scaffold split ROC-AUC **0.7865** ± 0.0382
- random split ROC-AUC **0.8231** ± 0.0301

**A finding worth stating plainly:** paracetamol and isoniazid — two of the most famous hepatotoxins in medicine — are *absent from this benchmark entirely*. The model scores paracetamol around 0.37, i.e. low risk, which is simply wrong. The benchmark encodes largely *idiosyncratic* liver injury, whereas paracetamol toxicity is intrinsic and dose-dependent.

That is precisely why the two layers are not redundant, and why the DILI term is down-weighted in the ranking while the structural alerts carry the weight. A statistical model that misses the project's own index compound has not earned a large vote.

## Screening

- library screened: **145,380** unique plant compounds
- removed as ChEMBL training overlap (leakage control): **133**
- within applicability domain: **95,869** (66.0%)
- carrying a severity-3 reactive alert: **828**
- PAINS-flagged: **13,281**
- final distinct-candidate shortlist: **50**

Composite score, with weights stated so the ranking is falsifiable:

```
composite = target_score
            × exp(-0.35 × alert_severity_sum)
            × ad_confidence
            × cns_bonus (1.15 if CNS-like)
            × (1 - 0.35 × p_DILI)
```

Target distribution among the 50 distinct candidates: TRPV1 35 / CB1 11 / FAAH 4.

### Top distinct candidates

| compound_name | best_target | target_score | ad_best | cns_like | composite_score | source_genera | npclass_class | n_variants_collapsed |
|---|---|---|---|---|---|---|---|---|
| Hydroxy alpha sanshool | CB1 | 0.991 | 1.000 | 1 | 0.990 | Zanthoxylum | N-acyl amines | 4 |
| N-[(4-hydroxy-3-methoxyphenyl)methyl]-7-meth | TRPV1 | 1.000 | 0.935 | 1 | 0.982 | Capsicum | Capsaicins and Capsaicinoids $ Phenylalanine | 1 |
| N-(4-aminobutyl)-3-(4-hydroxy-3-methoxypheny | TRPV1 | 0.973 | 0.653 | 1 | 0.969 | Iochroma|Magnolia|Nicotiana|Solanum | Polyamines | 2 |
| N-[(4-hydroxy-3-methoxyphenyl)methyl]-6-meth | TRPV1 | 1.000 | 0.851 | 1 | 0.968 | Capsicum | Capsaicins and Capsaicinoids $ Phenylalanine | 2 |
| 3-(4-Hydroxy-3-methoxyphenyl)-1-(piperidin-1 | TRPV1 | 1.000 | 1.000 | 1 | 0.967 | Piper | Piperidine alkaloids | 1 |
| (2E,4E,8Z,10E,12E)-N-(2-hydroxypropan-2-yl)t | CB1 | 0.985 | 0.676 | 1 | 0.966 | Zanthoxylum | N-acyl amines | 1 |
| (2Z,6E,8E,10E)-N-(2-methylpropyl)dodeca-2,6, | CB1 | 0.957 | 0.879 | 1 | 0.952 | Zanthoxylum | N-acyl amines | 5 |
| 5-(1,3-benzodioxol-5-yl)-N-isobutyl-penta-2, | TRPV1 | 1.000 | 1.000 | 1 | 0.951 | Piper | nan | 3 |
| (E)-5-(4-hydroxy-3-methoxy-phenyl)-1-(1-pipe | TRPV1 | 0.990 | 0.620 | 1 | 0.951 | Piper | Piperidine alkaloids | 2 |
| (2E,4E)-5-(1,3-benzodioxol-5-yl)-N-[(2S)-but | TRPV1 | 0.989 | 0.700 | 1 | 0.949 | Piper | nan | 1 |
| 3-(4-Hydroxy-3-methoxyphenyl)prop-2-enyl 3-m | TRPV1 | 0.982 | 0.549 | 1 | 0.940 | Coreopsis|Eremanthus|Lychnophora | Cinnamic acids and derivatives | 2 |
| nan | TRPV1 | 0.996 | 0.660 | 1 | 0.936 | Asarum|Piper | nan | 4 |
| 5-(1,3-benzodioxol-5-yl)-N-(3-methylbutyl)pe | TRPV1 | 0.984 | 0.792 | 1 | 0.935 | Piper | nan | 2 |
| (2E,4E)-5-(1,3-benzodioxol-5-yl)-N-[(2R)-2-m | TRPV1 | 0.993 | 0.812 | 1 | 0.934 | Piper | nan | 1 |
| Methyl 17-ethyl-14-hydroxy-6-methoxy-3,13-di | CB1 | 0.948 | 0.570 | 1 | 0.931 | Tabernaemontana | Iboga type | 1 |
| (2E,4E)-5-(4-hydroxy-3-methoxyphenyl)-1-(pip | TRPV1 | 0.996 | 0.837 | 1 | 0.921 | Piper | Piperidine alkaloids | 2 |
| (6aS,9S,10S,10aS)-6,6,9-trimethyl-3-pentyl-7 | CB1 | 0.987 | 0.727 | 1 | 0.913 | Cannabis | Cannabinoids | 2 |
| methyl (1S,11R,13R,15S,16S,18S)-15-ethyl-19- | CB1 | 0.935 | 0.794 | 1 | 0.912 | Tabernaemontana | Iboga type | 1 |
| 3-(4-hydroxy-3-methoxyphenyl)-N-[2-(3-hydrox | TRPV1 | 0.972 | 0.627 | 1 | 0.902 | Chenopodium|Pisonia | Cinnamic acid amides | 2 |
| (E)-3-(3-hydroxy-4-methoxyphenyl)-N-[2-(4-hy | TRPV1 | 0.969 | 0.596 | 1 | 0.901 | Chenopodium | Cinnamic acid amides | 2 |
| 3-(4-hydroxy-3-methoxyphenyl)-N-(2-phenyleth | TRPV1 | 0.986 | 0.640 | 1 | 0.899 | Annona|Solanum | Cinnamic acid amides | 2 |
| Cannabitriol | CB1 | 0.983 | 0.768 | 1 | 0.895 | Cannabis | Cannabinoids | 1 |
| Ethyl palmitate | FAAH | 0.976 | 0.720 | 0 | 0.892 | Achillea|Aloe|Aristolochia|Artemisia|Azadira | Wax monoesters | 1 |
| 1-(4'-Hydroxy-3'-methoxyphenyl)-7-octen-3-on | TRPV1 | 0.924 | 0.608 | 1 | 0.887 | Zingiber | nan | 1 |
| 2-Hydroxyethyl octadeca-9,12-dienoate | TRPV1 | 0.999 | 0.581 | 0 | 0.885 | Tagetes | Wax monoesters | 2 |

## External validation

A curated set of natural-product modulators of these targets with activities from the primary literature, plus synthetic reference ligands. Compounds also present in ChEMBL training data are reported separately — only the held-out ones are real validation.

- n claims: **26**
- n held out: **16**
- n seen in training: **10**
- held out median probability: **0.0006**
- held out frac above 0.5: **0.0**
- held out median library percentile: **75.0**
- held out frac top10pct of library: **0.1875**
- seen median probability: **0.9999**

Per-compound detail: `results/literature_validation.csv`.

## The key finding: a measurable blind spot, and its fix

A smoke test showed the FAAH model assigning near-zero probability to kaempferol and genistein — the very natural-product FAAH inhibitors that motivated this project.

We checked why, rather than assuming:

```
FAAH rows in curated ChEMBL data:     4,979  (4,257 compounds)
kaempferol, genistein, biochanin A,
daidzein, 7-hydroxyflavone, silybin:  NONE present
flavonoid-scaffold compounds:         2 of 4,979
   of which active at 1 uM:           0
URB597 (reference inhibitor):         present, pAct 7.5-8.3 ✓
```

URB597's presence proves the curation is sound — we are not dropping real FAAH pharmacology. The flavonoids are genuinely absent from ChEMBL's FAAH target records. **The near-zero scores are extrapolation into empty chemical space, not a judgement.**

### Is it fixable? Leave-one-out experiment

For each natural inhibitor, we retrained on ChEMBL FAAH data + decoys + *the other* natural inhibitors, and scored the held-out compound. Each compound is therefore never in its own training set.

| name | reported_uM | p_baseline | p_leave_one_out | fold_change |
|---|---|---|---|---|
| kaempferol | 5.000 | 0.000 | 0.316 | 8538.7 |
| genistein | 1.300 | 0.000 | 0.999 | 7032.4 |
| biochanin A | 2.400 | 0.000 | 0.344 | 1292.9 |
| daidzein | 14.000 | 0.000 | 0.415 | 12509.3 |
| 7-hydroxyflavone | 2.040 | 0.000 | 0.003 | 141.859 |
| silybin | 5.080 | 0.000 | 0.000 | 0.893 |


Median probability rises from **3.509e-05** (ChEMBL alone) to **0.329796** when five curated literature actives are added.

**Interpretation.** The blind spot is a *data-coverage* problem, not a limitation of the fingerprint representation. A handful of curated literature actives — fewer than ten compounds — is enough for the model to recognise an unseen member of the chemotype. That is a cheap, concrete remedy and the clearest argument for why this project needed a held-out literature validation set at all: without it we would have shipped a confident-looking shortlist and never discovered the gap.

## Knowledge graph

**9,090 nodes, 9,846 edges.**

Node types: Compound (8,880), Plant (178), ChemicalClass (20), ? (4), Target (3), Pathway (3), Disease (1), SafetyFlag (1)

Edges by evidence class: **measured** (8,707), **literature** (578), **predicted** (561)

Every edge carries an `evidence` property — `measured` (a ChEMBL datapoint with its potency), `literature` (a LOTUS pair with its reference DOI) or `predicted` (a model output with its probability and applicability domain). That distinction is the point of using a graph: a candidate can be shown as a path in which the measured and predicted links look different. A graph where a prediction renders identically to a measurement is a graph that misleads.

Files: `results/knowledge_graph.graphml` (Cytoscape/Gephi), `results/knowledge_graph.html` (interactive).

## Limitations

1. **Predictions, not findings.** No compound here has been shown to be active or safe. The output is a ranked list of hypotheses for assay.

2. **Agonist and antagonist are not separated in the activity label.** For TRPV1 this matters a great deal: capsaicin activates, capsazepine blocks, and CBD activates *then desensitises*. All three register as "active". ChEMBL's `action_type` was retained in the dataset and is the natural next step for splitting them.

3. **Decoys are an assumption.** Presumed inactives are labelled `presumed_inactive` and kept distinguishable from measured inactives, so any result can be re-derived without them.

4. **Assay heterogeneity.** IC50, EC50, Ki and Kd are pooled onto one potency scale, and the literature values mix species and assay formats (calcium imaging, patch clamp, radioligand binding) that are not strictly comparable.

5. **The DILI model is weak and out of domain** for natural products, as the paracetamol failure above demonstrates.

6. **Natural products are hard for fingerprint models.** They are more stereochemically complex and sp3-rich than ChEMBL's medicinal-chemistry space; the applicability-domain flag exists because much of the library is genuinely extrapolation.

7. **Bioavailability is not modelled.** Many top-ranked polyphenols are glycosylated, extensively metabolised in the gut, and poorly absorbed. The CNS-likeness term is a crude proxy for a real ADME problem.

## Reproducing

```bash
python src/run_all.py          # full pipeline, skips completed stages
python src/run_all.py --force  # rebuild everything
```
