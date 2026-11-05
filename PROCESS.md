# Terranova — the entire process, explained

**Track 02: Natural Products for Drug Discovery**

This document explains what we built, why each decision was made, and what the
results do and do not mean. Read it top to bottom and you can defend any part
of the pipeline.

> **One-line summary:** we find plant-derived compounds that may act on the
> three targets behind paracetamol's *analgesic* mechanism (TRPV1, CB1, FAAH),
> rank them with machine learning, and filter out the chemistry behind its
> *hepatotoxic* one.

---

## Part 1 — The scientific idea

### Why paracetamol is the right starting point

Paracetamol (acetaminophen) is the most-used analgesic in the world and we
still argue about how it works. The modern account is that its analgesia is
not mainly COX inhibition at all. Instead:

1. Paracetamol is deacetylated to ***p*-aminophenol**.
2. **FAAH** (fatty acid amide hydrolase) conjugates *p*-aminophenol with
   arachidonic acid to form **AM404**.
3. **AM404** acts centrally and in the spinal dorsal horn, engaging **TRPV1**
   and the **cannabinoid (CB1)** system.

Its liver toxicity runs through a *completely different* route:

1. A minor fraction is oxidised by CYP2E1/CYP3A4 to **NAPQI**
   (*N*-acetyl-*p*-benzoquinone imine).
2. NAPQI is a hard electrophile. At therapeutic dose, hepatic **glutathione**
   conjugates it harmlessly.
3. In overdose the glutathione pool is depleted, NAPQI arylates protein
   thiols, and you get oxidative stress, mitochondrial dysfunction and
   centrilobular necrosis.

### The key insight the project is built on

**The analgesic pathway and the toxic pathway are chemically separable.**

The analgesia needs a molecule that engages TRPV1 / CB1 / FAAH. The toxicity
needs a *para*-aminophenol core that can be oxidised to a quinone imine. There
is no law of nature saying a molecule must have both.

So the search is well posed: find plant compounds that hit the first set of
targets **without** carrying the structural features that enable the second.

That gives us two independent questions, which is why the project has two
independent model families:

| Question | What answers it |
|---|---|
| Can this compound plausibly engage a pain target? | Three per-target ML classifiers |
| Does it carry liver-risk chemistry? | Structural alert set + DILI model |

---

## Part 2 — Where the data came from

Nothing here is synthetic or hand-waved. Every row traces to a public source.

| Source | What we took | Scale |
|---|---|---|
| **ChEMBL REST API** | measured bioactivity for TRPV1/CB1/FAAH × human/rat/mouse | **58,394 raw records** |
| **ChEMBL molecules** | presumed-inactive decoy pool sampled across 1.84M | **59,923 molecules** |
| **LOTUS** (Zenodo frozen v11, Apr 2026) | referenced compound–organism pairs, plants only | **558,566 pairs** |
| **TDC / FDA NCTR** | DILI (liver injury) labels for marketed drugs | **475 drugs** |
| **PubChem PUG-REST** | structures for the literature validation set | **24 compounds** |
| **Hand-curated** | ethnopharmacology seed list of pain-use plants | **50 species** |

### How we resolved the targets

We did *not* hardcode ChEMBL IDs from memory. We resolved them from **UniProt
accessions**, which is the stable identifier, and kept the rodent orthologs
separately:

| Target | Human | Rat | Mouse |
|---|---|---|---|
| TRPV1 | CHEMBL4794 (Q8NER1) | CHEMBL5102 | CHEMBL1781864 |
| CB1 | CHEMBL218 (P21554) | CHEMBL3571 | CHEMBL3037 |
| FAAH | CHEMBL2243 (O00519) | CHEMBL3229 | CHEMBL3455 |

Rodent data is kept in its own column and never silently pooled with human.
TRPV1 in particular has well-documented species differences in ligand
sensitivity — rat and human TRPV1 respond differently to several antagonists.

### Why LOTUS rather than just scraping compound names

LOTUS gives **referenced structure–organism pairs**: every (compound, plant)
link carries a literature DOI, plus full Open Tree of Life taxonomy down to
species, plus NPClassifier chemical classification. That provenance is what
lets the knowledge graph be defensible rather than decorative — we can point at
the paper behind every "this plant contains this compound" edge.

We filtered to **Archaeplastida** (land plants), discarding fungal, bacterial
and marine sources, because the project premise is *plant*-derived analgesics.
That took 709,740 pairs down to 558,566.

### The plant seed list

We hand-curated 50 species traditionally used for pain, each with its family,
its traditional use, and a target hypothesis — *Capsicum annuum* (capsaicin,
TRPV1), *Curcuma longa* (curcumin), *Zingiber officinale* (gingerols),
*Syzygium aromaticum* (eugenol), *Cannabis sativa*, *Silybum marianum*
(silybin, FAAH), *Glycine max* (genistein/daidzein), *Corydalis yanhusuo*,
*Harpagophytum procumbens*, and so on.

Compounds are tagged in two tiers:

- **`pain_seed`** — from a seed-list plant (16,729 compounds)
- **`plant_background`** — from any other plant (128,651 compounds)

The background tier is not padding. Without it, "top ranked" would be
meaningless — we would only ever be ranking compounds we already believed in.

---

## Part 3 — The engineering problems we hit (and how we solved them)

These are worth telling because they are real and they show judgment.

### Problem 1: The ChEMBL API was failing 70–80% of requests

Mid-build, EBI's list endpoints started returning HTTP 500 and read timeouts
on most requests. We measured it rather than guessing: success rate was
~20–30%, **independent of page size**, and the failures were **random** — a
request could fail nine times then succeed, while an adjacent offset succeeded
first try.

That diagnosis dictated the fix. Random failures are *not* congestion, so
exponential backoff is exactly wrong — it just multiplies dead waiting time.
What works instead:

- **constant short retries** (2s) with many attempts, not exponential backoff
- a **short 30s timeout** — a healthy response arrives in 4–15s, so anything
  still open at 30s is a hung connection, and waiting 60s on it is pure waste
- **concurrency across offsets**, since pages are independent
- **per-target checkpointing**, so an hour-long harvest survives interruption

Result: a harvest that was stalling indefinitely completed cleanly.

### Problem 2: The applicability-domain calculation was going to take hours

Scoring 150k natural products against 3 models, each with thousands of
training compounds, is **billions of pairwise Tanimoto comparisons**. Done
molecule-by-molecule it would dominate the entire runtime.

Tanimoto on binary fingerprints is

```
T(A,B) = |A ∩ B| / (|A| + |B| - |A ∩ B|)
```

and `|A ∩ B|` for *all pairs at once* is exactly a **matrix product** of the
two bit matrices. We rewrote it as BLAS matrix multiplication and verified it
against RDKit's own `BulkTanimotoSimilarity`: **max absolute difference 0.0**.
Same answer, orders of magnitude faster.

### Problem 3 (the scientific one): the active rate was a lie

This is the most important thing we found.

After curation, the measured data was:

| Target | Active rate |
|---|---|
| TRPV1 | **84.4%** |
| FAAH | **75.1%** |
| CB1 | **59.7%** |

That is not biology. **That is publication bias.** ChEMBL is assembled from
medicinal-chemistry papers, and people publish compounds that worked. Genuine
inactives are massively under-reported.

A classifier trained on 84% actives learns that nearly everything is active.
Point it at 150,000 natural products and it returns high probabilities almost
everywhere — the resulting "ranking" is noise. **Class weighting does not fix
this**, because the problem is not the loss function; it is that the model has
never been shown what a non-ligand looks like.

**The fix: property-matched decoys.** Real screening hit rates are well under
1%, so a random drug-like molecule is overwhelmingly likely to be inactive
against TRPV1. We sample presumed-inactives from ChEMBL at 3 per active, which
brings the active rate down to ~21%.

Two safeguards make the decoys honest:

1. **Property matching (DUD-E style).** Decoys are drawn from the same
   (molecular weight, cLogP) cell as the active they match. Without matching,
   the model separates classes on gross physical properties — "heavy and
   greasy means active" — which inflates benchmark scores and transfers
   nothing. Measured matching quality for CB1/human: **MW median gap 12.1 Da,
   cLogP median gap 0.44**.
2. **Scaffold exclusion.** Any decoy sharing a Bemis-Murcko scaffold with a
   known active for that target is dropped (282 dropped for CB1/human). Without
   this, a "decoy" can be a close analogue of a real ligand that simply was
   not measured — teaching the model that real actives are inactive.

Decoys are labelled **`presumed_inactive`**, kept distinguishable from measured
inactives, so any result can be re-derived without the assumption. A presumed
inactive is an assumption; a measured inactive is data.

---

## Part 4 — How the data was curated

Every filter is counted and written to `data/processed/attrition_audit.csv`.
Hidden curation is how cheminformatics results become irreproducible.

```
raw activity records                      58,394
has SMILES                                58,287
assay_type in {B,F}                       57,214
standard_type in {IC50,EC50,Ki,Kd}        35,510
units == nM                               32,574
has numeric value                         32,386
no data_validity_comment                  31,854
relation in {=, >, >=}                     31,510
pactivity in plausible range [2,12]       31,510
structure standardized                     31,271
   → aggregated to 22,282 compound×target×organism rows
```

### Decisions worth defending

**We kept functional assays, not just binding.** TRPV1 is a ligand-gated ion
channel; most of its potency data is calcium-flux and patch-clamp, registered
as `assay_type = F`. The habitual choice for GPCR/enzyme panels is binding
only (`B`) — that would have discarded most TRPV1 pharmacology.

**Censored values are used correctly.** A record of ">10,000 nM" is not a
measurement, but it is strong evidence of *inactivity*. We admit `>` and `>=`
records as **inactives only**, never as actives. `<` records are discarded
entirely — they bound potency from above and cannot confirm an active.

**Median aggregation.** The same compound is measured many times across labs.
The median is robust to the occasional order-of-magnitude outlier; the mean is
not.

**The ambiguous band is held out, not forced.** Active ≤ 1 µM, inactive
≥ 10 µM, and the 1–10 µM band is dropped from training rather than shoved into
a class. Compounds there are genuinely ambiguous, and training on them teaches
the model noise.

**FAAH is labelled differently on purpose.** The natural-product FAAH
inhibitors that motivate this project sit at roughly 1–15 µM:

| Compound | Reported FAAH activity |
|---|---|
| Kaempferol | Ki 5 µM (review table: IC50 1.06 µM) |
| Genistein | IC50 1.3 µM / Ki 8 µM |
| Biochanin A | IC50 2.1–2.4 µM |
| 7-Hydroxyflavone | IC50 2.04 µM |
| Silybin | IC50 5.08 µM |
| Licarin A | IC50 7.02 µM |
| N-benzyl-linoleamide | IC50 7.2 µM |
| Daidzein | IC50 14 µM |

A 1 µM cut classifies **nearly all of them as inactive**, optimising the model
to reject exactly the chemistry we are hunting. So FAAH also gets a **10 µM
classification variant**, and both are reported.

### Structure standardization

Both datasets pass through the **same** normalization — otherwise a compound
could appear as a salt in ChEMBL and a free base in LOTUS and never be
recognised as the same molecule. That failure is silent and it is the most
common way an NP screening pipeline reports a "novel hit" that was in its own
training set.

Pipeline: parse → Cleanup → largest-fragment (strip salts/solvates) →
uncharge → reionize → canonical isomeric SMILES + InChIKey. Plus element
allow-list, heavy-atom bounds (6–70) and an MW ceiling.

Verified: paracetamol and its HCl salt collapse to the same InChIKey
(`RZVAJINKPMORJF-UHFFFAOYSA-N`).

---

## Part 5 — The models

**Features:** ECFP4 (Morgan radius 2, 2048 bits) + 14 physicochemical
descriptors. A deliberately boring, strong baseline. A graph neural network on
~10k compounds per target would not reliably beat it and would cost the whole
hackathon.

**Algorithms:** Random Forest and LightGBM, both class-balanced. The better
scaffold-split model per target is selected and refit on everything.

### The single most important methodological choice: scaffold splitting

ChEMBL is built from medicinal-chemistry papers, each contributing a
**congeneric series** of close analogues. Under a **random split** those
analogues straddle train and test, the model memorises the series, and you
report an ROC-AUC near 0.95 that **collapses on anything genuinely new**.

Screening a natural-product library is *exactly* the "anything genuinely new"
case — natural products share almost no scaffolds with medicinal-chemistry
space.

So we group compounds by **Bemis-Murcko scaffold** and assign whole scaffold
groups to train or test. We report **both splits and the gap between them**.
Showing the gap is more informative than any single headline number.

> Worked example of why scaffolds matter: kaempferol (a flavonol,
> 2-phenylchromen-4-one) and genistein (an isoflavone, 3-phenylchromen-4-one)
> receive **different** Murcko scaffolds, because the phenyl sits at a different
> position. They are genuinely distinct chemotypes, and the split correctly
> refuses to let the model learn one from the other.

### Metrics, and why these ones

- **PR-AUC leads** — the data is imbalanced and we care about the top of a
  ranked list, not global accuracy.
- **ROC-AUC, MCC, balanced accuracy** for completeness.
- **Enrichment factor @1% and @5%** — literally the screening question: if we
  hand a chemist the top 1% of a library, how many more actives do they get
  than by picking at random?
  - *Caveat we report honestly:* at a ~21% active rate the theoretical maximum
    EF is only ~4.6, and we are measuring ~5.0. **EF is saturated and no longer
    informative** at this ratio — PR-AUC and MCC are the meaningful numbers.
- **Applicability domain** — every prediction carries its max Tanimoto
  similarity to the training set. Below ~0.25 it is flagged as extrapolation,
  not silently trusted.

### Results so far

All eight model configurations trained; LightGBM won every one. Scaffold-split
(the honest number):

| Target | Scope | Label | Active rate | Scaffold PR-AUC | ROC-AUC | MCC |
|---|---|---|---|---|---|---|
| TRPV1 | all | 1 µM | 22.7% | **0.983** | 0.995 | 0.901 |
| CB1 | all | 1 µM | 23.0% | **0.981** | 0.993 | 0.915 |
| TRPV1 | human | 1 µM | 26.2% | **0.973** | 0.992 | 0.905 |
| FAAH | all | 1 µM | 16.5% | **0.962** | 0.987 | 0.901 |
| CB1 | human | 1 µM | 19.9% | **0.961** | 0.984 | 0.890 |
| FAAH | human | 1 µM | 22.7% | **0.956** | 0.986 | 0.904 |
| FAAH | all | 10 µM | 21.3% | **0.930** | 0.977 | 0.854 |
| FAAH | human | 10 µM | 27.5% | **0.895** | 0.931 | 0.801 |

The small random-vs-scaffold gap (0.015 for CB1/human) is the encouraging part:
the model is not merely memorising congeneric series. The 10 µM FAAH variant is
hardest, as expected — a looser cut admits weaker, more structurally diverse
chemistry.

### The most important result in the project: a measurable blind spot

A smoke test on known compounds produced this:

| Compound | p(FAAH) | Known FAAH activity |
|---|---|---|
| Kaempferol | **0.000006** | Ki 5 µM / IC50 1.06 µM |
| Genistein | **0.000065** | IC50 1.3 µM |
| Capsaicin | 0.032 | (TRPV1 compound; p(TRPV1) = **0.99999**) |

The model assigns **essentially zero probability to the exact natural-product
FAAH inhibitors that motivated this entire project.**

Before calling that a finding we checked which of three things was happening
(`src/check_faah_flavonoids.py`): are they in training and mislabelled, in
training and scored wrong, or simply absent?

**The answer is unambiguous:**

```
FAAH rows in curated ChEMBL data:        4,979  (4,257 compounds)
kaempferol, genistein, biochanin A,
daidzein, 7-hydroxyflavone, silybin:     NONE present
flavonoid-scaffold compounds present:    2 of 4,979
   of which active at 1 µM:              0
URB597 (reference inhibitor):            present, pAct 7.5–8.3, active ✓
```

URB597's presence proves our curation is sound — we are not dropping real FAAH
pharmacology. The flavonoids are genuinely **not in ChEMBL's FAAH target
records**. The natural-product FAAH literature (Thors et al. 2008; Zada et al.
2022) simply has not been curated into that target.

So the near-zero scores are **extrapolation into an unoccupied region of
chemical space, not a considered negative judgement.** The model has literally
never seen a flavonoid engage FAAH.

**Why this matters more than a pretty shortlist.** The honest headline of this
project is not "here are 50 promising compounds". It is:

> A conventional ChEMBL-trained ML screen has a specific, measurable blind spot
> for exactly the natural-product chemistry it is being asked to find — and we
> can show it, quantify it, and say what to do about it.

What to do about it:
1. **Trust the applicability-domain flag.** It exists precisely for this; most
   flavonoids fall below the similarity floor against the FAAH model.
2. **Augment training with curated literature actives.** Our
   `literature_actives.csv` is the seed. A leave-one-out experiment (train on
   7 of the 8 known flavonoid inhibitors, test on the 8th) would measure
   whether even a handful of natural actives rescues generalisation.
3. **Complement ML with similarity/pharmacophore search** anchored on the known
   natural inhibitors, which does not require the target to be populated.

This is also the clearest possible argument for *why the project needed a
held-out literature validation set at all*. Without it, we would have shipped a
confident-looking shortlist and never discovered that the FAAH model cannot see
flavonoids.

**An honest caveat about those high numbers:** adding decoys changes the task.
The model is now answering *"is this a CB1-like molecule at all?"*, which is
genuinely easier than *"is this potent?"*. That is the right question for
library screening, but it means the AUC is not comparable to a potency-ranking
benchmark.

---

## Part 6 — The safety layers

### Layer A — reactive-metabolite structural alerts

25 hand-written SMARTS patterns, graded by severity, built around the NAPQI
mechanism rather than copied from a generic list:

| Severity | Examples |
|---|---|
| **3** (direct bioactivation precedent) | *para*-aminophenol, *para*-hydroxy anilide, quinone imine, *para*/*ortho*-quinone, acyl halide, nitrosamine, pyrrolizidine core |
| **2** (established reactive route) | catechol, hydroquinone, primary aniline, nitroaromatic, Michael acceptor enone, α-methylene lactone, epoxide, furan, thiophene, hydrazine, thiourea |
| **1** (context-dependent) | aldehyde, terminal alkyne, acrylate ester, azo |

Plus RDKit's **PAINS / Brenk / NIH** catalogs, because an NP hit list that has
not been checked for assay-interference chemistry is not credible — many
polyphenols are notorious frequent hitters.

**This layer validates beautifully.** Run against known compounds:

| Compound | Alerts | Max severity |
|---|---|---|
| **Paracetamol** | para_aminophenol, anilide_para_hydroxy | **3** |
| **NAPQI** | quinone_imine, michael_acceptor + PAINS hit | **3** |
| Curcumin | michael_acceptor_enone | 2 |
| 6-Gingerol | — | 0 |
| Kaempferol | — | 0 |
| Genistein | — | 0 |
| Capsaicin | — | 0 |
| CBD | — | 0 |

**The filter flags paracetamol at maximum severity on the exact mechanism that
makes it hepatotoxic, while the natural analgesics we care about pass clean.**
That is the headline slide.

It also flags curcumin as a Michael acceptor and many flavonoids as catechols.
**We report those flags rather than tuning them away.** Plants make
electrophiles deliberately, often as defence chemistry — that is a real finding
about this chemical space, not a bug.

### Layer B — the DILI model

Random forest on the TDC drug-induced-liver-injury benchmark (475 marketed
drugs, ~50% positive). Scaffold-split ROC-AUC **0.889 ± 0.024**, random
**0.903 ± 0.044**.

**But here is what we found, and it matters more than the AUC:**

> **Paracetamol and isoniazid — two of the most famous hepatotoxins in
> medicine — are absent from this benchmark entirely.**

We verified it directly (`src/check_dili_labels.py`). Diclofenac, troglitazone
and valproate are present as positives; aspirin as a negative; paracetamol and
isoniazid are simply not there. The model scores paracetamol **~0.37 — "low
risk"**, which is flatly wrong.

Why? The benchmark largely encodes ***idiosyncratic*** DILI: unpredictable,
host-dependent, not dose-proportional. Paracetamol toxicity is the opposite —
**intrinsic, dose-dependent, mechanistically understood.**

**This is the strongest argument for the two-layer design.** The statistical
model cannot see the project's own index compound; the structural alert set
flags it at maximum severity. So in the ranking, the DILI term is deliberately
**down-weighted (0.30)** while the alerts carry the weight. A model that misses
paracetamol has not earned a large vote.

### Layer C — CNS-likeness

The entire premise is *central* analgesia via the AM404 axis. A compound that
cannot plausibly cross the blood-brain barrier is a weaker candidate for *this
hypothesis* regardless of its target score. We compute a simplified
Wager-style CNS-MPO check (MW ≤ 400, 1 ≤ cLogP ≤ 4, TPSA ≤ 90, HBD ≤ 3,
RotB ≤ 8) alongside Lipinski.

---

## Part 7 — Screening and ranking

```
152,414 plant compounds
  → standardize (same pipeline as ChEMBL)
  → REMOVE anything present in training data        ← leakage control
  → score against TRPV1 / CB1 / FAAH models
  → applicability-domain similarity per model
  → reactive alerts + PAINS + DILI + CNS-likeness
  → composite rank
```

### Leakage control is not a formality

A natural-product hit list that silently contains the model's own training
compounds will rank **capsaicin first on TRPV1** and look spectacular while
demonstrating nothing. We remove exact InChIKey matches, report how many, and
re-add them separately as **positive controls** so we can check the model
recovers known actives.

### The composite score, stated so it is falsifiable

```
composite = target_score
          × exp(-0.35 × alert_severity_sum)
          × ad_confidence
          × 1.15 if CNS-like
          × (1 - 0.30 × p_DILI)
```

Weights are printed in `results/screening_stats.json`. An opaque ranking
formula is unfalsifiable; this one can be argued with, and the top-N can be
re-derived under different weights.

The shortlist additionally requires: **in applicability domain**, **no
severity-3 alert**, **no PAINS flag**.

---

## Part 8 — The knowledge graph

**Nodes:** Plant, Compound, Target, Pathway, Disease(Pain), SafetyFlag,
ChemicalClass

**Edges:**

| Relation | From → To | Evidence |
|---|---|---|
| `CONTAINS` | Plant → Compound | literature (+ reference DOI) |
| `ACTS_ON` | Compound → Target | **measured** (+ pActivity, assay type) |
| `PREDICTED_ACTIVE` | Compound → Target | **predicted** (+ probability, AD) |
| `HAS_FLAG` | Compound → SafetyFlag | predicted |
| `MEMBER_OF` | Compound → ChemicalClass | literature (NPClassifier) |
| `INVOLVED_IN` | Target → Pathway | literature |
| `MODULATES` | Pathway → Pain | literature |
| `TRADITIONALLY_USED_FOR` | Plant → Pain | literature (ethnopharmacology) |

**Every edge carries an `evidence` property: `measured`, `literature` or
`predicted`** — and they render differently (green solid / blue solid / yellow
dashed).

That distinction is the entire reason to build a graph rather than a
spreadsheet. It lets a candidate be presented as a *path* — "this compound is
reported in this plant by this paper, the plant is traditionally used for
pain, and our model predicts it engages FAAH" — with the measured and
predicted links **visibly different**. A graph where a prediction looks
identical to a measurement is a graph that misleads.

Outputs: `knowledge_graph.graphml` (Cytoscape/Gephi) and an interactive
`knowledge_graph.html`.

---

## Part 9 — External validation

A curated set of 24 compounds with activities from the primary literature —
natural-product modulators plus synthetic reference ligands (capsaicin,
resiniferatoxin, capsazepine, URB597, PF-3845, rimonabant, THC, anandamide,
and **AM404 itself**). All resolved to structures via PubChem.

**It is deliberately not training data.** If these compounds were folded into
training, the headline claim — "the model recovers known natural analgesic
chemistry" — would be circular. Compounds that *do* appear in ChEMBL training
data are identified and reported separately; only held-out ones count.

The question is **not** "what is the AUC" — with 24 compounds that number is
meaningless. The question is: **does the model rank known modulators above the
bulk of the plant library?** We answer it by reporting each compound's
**percentile against the full 150k background**, which is exactly the screening
decision the pipeline exists to make.

The literature table also records assay format, species, and value type
(Ki vs IC50 vs EC50) per entry, because — as we confirmed while researching —
published values mix rat/human/mouse enzymes and calcium imaging vs patch
clamp vs radioligand binding, and **are not directly comparable**. Kaempferol
alone is reported at both Ki 5 µM and IC50 1.06 µM in different sources.

---

## Part 10 — Status

### Done

- [x] Target resolution from UniProt → 9 ChEMBL targets
- [x] ChEMBL harvest — **58,394 activity records**
- [x] Decoy pool — **59,923 molecules** sampled across 1.84M
- [x] Curation + attrition audit → **22,282 compound×target rows**
- [x] LOTUS extraction → **558,566 plant pairs / 152,414 unique compounds**
- [x] Shared standardization (verified on paracetamol salt forms)
- [x] Reactive-metabolite alert set (validated on paracetamol/NAPQI)
- [x] DILI model + the paracetamol-absence finding
- [x] Literature validation set resolved (24/24)
- [x] Fast Tanimoto AD (verified exact vs RDKit)
- [x] CB1 models trained — scaffold PR-AUC 0.961 (human), 0.981 (all species)

### In progress

- [ ] TRPV1 and FAAH models (CB1 done; ~20 min)
- [ ] Plant library safety profiling (~95k of 145k done)

### Remaining

- [ ] Screen 152k compounds → ranked shortlist
- [ ] Literature validation run
- [ ] Knowledge graph build + interactive export
- [ ] Figures + auto-generated `results/REPORT.md`

---

## Part 11 — Limitations (say these out loud)

1. **These are predictions, not findings.** No compound here has been shown to
   be active or safe. The output is a ranked list of hypotheses for assay.
2. **Agonist vs antagonist is not separated in the activity label.** For TRPV1
   this matters enormously: capsaicin *activates*, capsazepine *blocks*, and
   CBD *activates then desensitises*. All three register as "active". We
   retained ChEMBL's `action_type` field in the dataset — splitting on it is
   the obvious next step.
3. **Decoys are an assumption**, labelled as such and separable from measured
   inactives.
4. **Assay heterogeneity.** IC50, EC50, Ki and Kd are pooled onto one potency
   scale.
5. **The DILI model is weak and out of domain** for natural products — as the
   paracetamol failure proves.
6. **Natural products are hard for fingerprint models** — more stereochemically
   complex and sp3-rich than ChEMBL space. The AD flag exists because much of
   the library is genuine extrapolation.
7. **Bioavailability is not modelled.** Many top-ranked polyphenols are
   glycosylated, heavily gut-metabolised and poorly absorbed. CNS-likeness is a
   crude proxy for a real ADME problem.
8. **Traditional use is not evidence of mechanism.** The ethnopharmacology seed
   list prioritises where to look; it does not support any pharmacological
   claim on its own.

---

## Part 12 — Reproducing

```bash
py -3.13 -m venv .venv
.venv\Scripts\python.exe -m pip install rdkit pandas numpy scikit-learn \
    lightgbm requests networkx pyvis matplotlib seaborn tqdm

python src/run_all.py          # full pipeline; skips completed stages
python src/run_all.py --force  # rebuild everything
```

Stage order and the artefact each produces is defined in `src/run_all.py`.
`results/REPORT.md` is generated entirely from files on disk, so it cannot
drift from what the code actually did.
