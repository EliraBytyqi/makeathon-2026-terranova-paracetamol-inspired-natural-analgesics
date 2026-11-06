# Terranova

**Track 02 — Natural Products for Drug Discovery**

Finding plant-derived compounds that may act on the pain targets behind
paracetamol's *analgesic* mechanism (TRPV1, CB1, FAAH), while screening out the
chemistry behind its *hepatotoxic* one.

> Everything this project outputs is a computational prediction. Nothing here
> establishes that a compound is active, or that it is safe.

---

## The idea in one paragraph

Paracetamol's pain relief is increasingly attributed not to COX inhibition but
to its metabolite **AM404**, formed by deacetylation to *p*-aminophenol and
conjugation with arachidonic acid by **FAAH**. AM404 acts centrally and in the
spinal dorsal horn through **TRPV1** and the **cannabinoid** system. Its liver
toxicity runs through a *different* route entirely: CYP oxidation of the
*para*-aminophenol core to **NAPQI**, a hard electrophile that depletes
glutathione and arylates protein thiols. Those two pathways are chemically
separable — so we look for plant compounds that engage the first without
carrying the structural features that enable the second.

## Pipeline

```
ChEMBL (TRPV1/CB1/FAAH, 3 species)  ──┐
ChEMBL decoy pool (60k molecules)   ──┤──> curate ──> scaffold split ──> models
                                       │
LOTUS plant natural products ─────────┴──> standardize ──> screen ──> rank
                                                              │
                              reactive-metabolite alerts ─────┤
                              DILI model ─────────────────────┤
                              applicability domain ───────────┘
                                                              │
                                                  knowledge graph + shortlist
```

## Quick start

```bash
py -3.13 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt

python src/run_all.py          # full pipeline; skips stages already done
python src/run_all.py --force  # rebuild everything
```

Read `results/REPORT.md` when it finishes.

## Dashboard

A local full-stack app presents the finished results — the hypothesis, data
provenance, model evaluation, external validation, safety layers, an explorer
over all 145,164 screened compounds, and a glossary for readers without a
chemistry background.

### Running it

From the **repository root**:

```bash
.venv\Scripts\python.exe -m app.run
```

It starts on <http://127.0.0.1:8000> and opens your browser. `Ctrl+C` stops it.

| Flag | Effect |
|---|---|
| `--port 9000` | Use a different port (if 8000 is taken) |
| `--no-browser` | Don't open a browser window |
| `--reload` | Auto-reload on code changes, for development |

### Running it from a fresh clone

You do **not** need the modelling stack. The dashboard loads no model and
imports no RDKit — it reads precomputed files — so three packages are enough:

```bash
git clone https://github.com/EliraBytyqi/makeathon-2026-terranova-paracetamol-inspired-natural-analgesics.git
cd makeathon-2026-terranova-paracetamol-inspired-natural-analgesics

git lfs install
git lfs pull                       # REQUIRED - see below

py -3.13 -m venv .venv
.venv\Scripts\python.exe -m pip install pandas fastapi "uvicorn[standard]"
.venv\Scripts\python.exe -m app.run
```

> **`git lfs pull` is not optional.** `results/screening_all.csv` (104 MB) and
> `results/positive_controls.csv` are stored in Git LFS. Without it you get
> small pointer files instead of data, and the screen explorer will fail to
> load while the rest of the app appears fine.

### If something is wrong

| Symptom | Cause and fix |
|---|---|
| Amber banner listing missing files | Those pipeline outputs aren't on disk. Run `python src/run_all.py` (needs the full `requirements.txt`). |
| Screen explorer empty, other views fine | LFS files weren't fetched. Run `git lfs pull`. |
| `Address already in use` | The app is already running, or something else holds the port. Use it, or pass `--port 9000`. |
| `No module named app` | You're not in the repository root. `cd` there first. |
| Knowledge graph is blank | Run `python src/run_all.py kg.py --force` to regenerate it and its local assets. |

It serves precomputed artefacts only, so it starts in about a second and nothing
loads a model at request time. No build step and no CDN — it runs with no
network access at all. See [`app/README.md`](app/README.md) for the API and the
design decisions.

## Layout

| Path | What it is |
|---|---|
| `src/chembl_client.py` | Resilient ChEMBL REST client (see note below) |
| `src/resolve_targets.py` | UniProt accession → ChEMBL target IDs |
| `src/fetch_chembl.py` | Bioactivity harvest, checkpointed per target |
| `src/fetch_decoys.py` | Presumed-inactive molecule pool |
| `src/build_np_library.py` | LOTUS → plant compound library, two tiers |
| `src/standardize.py` | Shared structure normalization |
| `src/featurize.py` | ECFP4 + descriptors, Bemis-Murcko scaffolds |
| `src/decoys.py` | DUD-E-style property matching |
| `src/safety.py` | Reactive-metabolite alerts, PAINS, CNS-likeness |
| `src/build_dataset.py` | Curation with a full attrition audit |
| `src/train.py` | Target models, random vs scaffold evaluation |
| `src/train_dili.py` | Hepatotoxicity model |
| `src/screen.py` | Library screening and composite ranking |
| `src/validate_literature.py` | External validation set |
| `src/kg.py` | Knowledge graph build and export |
| `src/make_report.py` | Generates `results/REPORT.md` from artefacts |
| `scripts/make_dili_probes.py` | Precomputes DILI probe scores for the dashboard |
| `app/` | FastAPI backend + static frontend (see `app/README.md`) |
| `data/reference/` | Hand-curated
| `results/` | Metrics, shortlist, graph, report |

## Four decisions that carry the project

**1. Scaffold splitting, not random.** ChEMBL is built from medicinal-chemistry
papers, each contributing a congeneric series. Under a random split those
analogues straddle train and test, the model memorises the series, and you
report an ROC-AUC near 0.95 that collapses on anything new. Since natural
products share almost no scaffolds with ChEMBL, that is exactly the wrong
measurement. We report both splits and the gap between them.

**2. Decoys, because the active rate is a lie.** After curation, measured data
is 84% active for TRPV1 and 75% for FAAH. That is publication bias, not
biology. A model trained on it calls everything active. We add
property-matched (MW, cLogP) presumed-inactive decoys, excluding any that
share a scaffold with a known active, and keep them labelled
`presumed_inactive` so results can be re-derived without the assumption.

**3. Two safety layers, because one is not enough.** The structural alert set
flags paracetamol at maximum severity on the exact NAPQI mechanism. The DILI
model, trained on a public benchmark, scores paracetamol as *low risk* —
because paracetamol and isoniazid are **absent from that benchmark entirely**.
We kept both layers, and weighted the alerts above the model.

**4. Evidence is typed throughout.** Every knowledge-graph edge is `measured`,
`literature` or `predicted`, and renders differently. A graph where a
prediction looks like a measurement is a graph that misleads.

## Note on the ChEMBL API

During development the EBI list endpoints were intermittently returning HTTP
500 and read timeouts on roughly 70–80% of requests, with failures distributed
randomly rather than by quota. Exponential backoff is the wrong response to
that failure mode — it multiplies dead time. `chembl_client.py` instead uses
many cheap constant-interval retries, a short timeout so hung connections fail
fast, and concurrency across independent offsets, with per-target
checkpointing so a long harvest survives interruption.

## FAAH is labelled differently, on purpose

The natural-product FAAH inhibitors motivating this work — kaempferol,
genistein, biochanin A, silybin, licarin A — sit at roughly 1–15 µM. A 1 µM
activity cut classifies nearly all of them as inactive, training the model to
reject exactly the chemistry we are hunting. FAAH therefore also gets a 10 µM
classification variant, and both are reported.
