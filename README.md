# Terranova

Screening plant natural products for the analgesic pharmacology of paracetamol,
while filtering out the chemistry that makes it hepatotoxic.

**Makeathon 2026 · Track 02 · Natural Products**
Elira Bytyqi · Shkurta Jashari · Fjolla Nebihi

---

## AI tools

We used **Claude** as a coding assistant during development. It helped us write,
review and debug parts of the project code. We reviewed and tested the code and
remain responsible for the implementation and results.

## Project description

### The chosen challenge

Track 02 — Natural Products. Use computational methods to find plant-derived
compounds with therapeutic potential.

### The problem being solved

Paracetamol is one of the most widely used drugs in the world and the leading
cause of acute liver failure in the West. The two facts come from the same
molecule taking two different metabolic routes:

- **Analgesia.** Deacetylation gives p-aminophenol, which crosses into the CNS.
  There **FAAH** conjugates it with arachidonic acid to form **AM404** — a
  potent **TRPV1** agonist, an anandamide reuptake inhibitor and a weak **CB1**
  ligand.
- **Hepatotoxicity.** **CYP2E1** and CYP3A4 oxidise the same parent molecule to
  **NAPQI**, an electrophilic quinone imine. Once glutathione is depleted it
  arylates mitochondrial protein thiols and kills hepatocytes.

Those two routes need different chemistry. The analgesic route needs an amide
FAAH will accept as a substrate; the toxic route needs a *para*-aminophenol that
CYP can oxidise. **In principle a molecule can offer the first without offering
the second** — and that gap is what this project screens for.

The pipeline scores **145,247 plant compounds** against TRPV1, CB1 and FAAH,
applies a reactive-metabolite and hepatic-liability safety layer, and returns a
ranked, provenance-tracked candidate list. It is a hypothesis list for a
laboratory, not an experimental result: **no compound was assayed by this team.**

### What it produces

| Output | File |
|---|---|
| 50 distinct candidate molecules, with plant source and literature DOIs | `results/candidates_top50.csv` |
| Full scored screen | `results/screening_all.csv` |
| External validation against measured outcomes | `results/pubchem_validation_expanded.json` |
| Knowledge graph (plant → compound → target → pain) | `results/knowledge_graph.html` |
| Written report | `results/REPORT.md` |

Headline result: on **3,691 structures measured by other groups and absent from
all training data**, the CB1 model reaches **ROC-AUC 0.917**, and **0.711** on
chemistry unlike anything it was trained on. TRPV1 reaches 0.918 but has no
measured inactives in novel chemical space, so it is validated only near its
training set. FAAH reaches 0.625 — real, and weak. Each is reported with that
grading attached.

---

## Technologies

### Programming language

- **Python 3.13**

### Libraries and frameworks

| Area | Packages |
|---|---|
| Cheminformatics | `rdkit` 2026.9.1 — structure standardisation, ECFP4 fingerprints, descriptors, SMARTS alerts, PAINS/Brenk/NIH catalogues |
| Machine learning | `scikit-learn` 1.9.1 (Random Forest, grouped splits, metrics), `lightgbm` 4.7.0 |
| Data | `pandas` 3.1.1, `numpy` 2.5.4, `scipy` 1.18.1, `openpyxl` 3.1.5 |
| Web backend | `fastapi` 0.162.0, `uvicorn` 0.54.0 |
| Web frontend | Plain ES2019, hand-rolled SVG charts — **no framework, no bundler, no build step** |
| Knowledge graph | `networkx` 3.7, `pyvis` 0.3.2 (vis-network vendored locally, runs offline) |
| Documents | `reportlab` 5.1.0, `pypdf` 7.0.0 |
| HTTP | `requests` 2.34.2 |

### Data sources

| Source | Role |
|---|---|
| **LOTUS** | Screening library — 145,380 referenced plant structures, each with a DOI |
| **ChEMBL** | Training labels (58,394 raw activities, 9 UniProt accessions) and the presumed-inactive pool |
| **PubChem BioAssay** | External test set — the only source that records *Inactive* outcomes |
| **FDA DILIrank 2.0** | Clinical hepatotoxicity labels |
| **Tox21** | HepG2 mitochondrial-potential and cytotoxicity assays |

COCONUT and NPASS were evaluated and **not used**: COCONUT carries structures
with no bioactivity data, and NPASS augmentation measurably reduced external
performance.

---

## Execution instructions

### Installation

```bash
git clone <repository-url>
cd terranova

python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate

pip install -r requirements.txt
```

Python 3.11 or newer is required (developed on 3.13).

### Run the application

The dashboard reads precomputed results from `results/`, so it starts
immediately after cloning:

```bash
python -m app.run
```

Then open **http://127.0.0.1:8000**.

During development, `--reload` picks up code changes automatically:

```bash
python -m uvicorn app.backend.main:app --reload
```

### Score a single molecule

```bash
python src/predict.py "CC(=O)Nc1ccc(O)cc1"      # SMILES
python src/predict.py --name paracetamol         # resolved via PubChem
python src/predict.py "<smiles>" --json          # machine-readable
```

### Rebuild the pipeline from source data

```bash
python src/run_all.py            # skips any stage whose output already exists
python src/run_all.py --force    # rebuild everything, including network fetches
python src/run_all.py screen.py  # run one named stage
```

Stages run in order: resolve targets → fetch ChEMBL → fetch decoys → build
dataset → build natural-product library → standardise → resolve literature →
build hepatic dataset → train hepatic panel → train activity models → screen →
validate → knowledge graph → report.

A full `--force` run re-downloads tens of thousands of records and takes
several hours. The default incremental run takes about 15 minutes.

---

## API management

**No API keys, accounts or secrets are required.** Every data source used here
is a public, unauthenticated endpoint:

| Service | Endpoint | Auth |
|---|---|---|
| ChEMBL | `https://www.ebi.ac.uk/chembl/api/data` | none |
| PubChem PUG REST | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | none |
| NCBI E-utilities | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils` | none |
| FDA DILIrank | static file download | none |
| LOTUS | static dataset download | none |

Requests are rate-limited in code and identified by a descriptive
`User-Agent`. All responses are cached under `data/external/`, so a rebuild
does not re-fetch what is already on disk. There is no `.env` file and nothing
to configure.

The project's own HTTP API is served by the dashboard and is also unauthenticated:

| Route | Returns |
|---|---|
| `GET /api/health` | status and which result files are present |
| `GET /api/overview` | headline screening statistics |
| `GET /api/candidates` | the 50 distinct candidate molecules |
| `GET /api/compounds/{source}` | paged, filterable compound tables |
| `GET /api/predict?smiles=...` | **live scoring of any molecule** |
| `GET /api/validation` | external validation results |
| `GET /kg` | interactive knowledge graph |

---

## Repository layout

```
src/            pipeline stages, in run order (run_all.py drives them)
app/            FastAPI backend + dependency-free frontend
scripts/        validation experiments and document generators
data/           raw, interim, processed and cached external downloads
models/         trained models (.pkl)
results/        every output: screen, candidates, validation, graph, report
docs/           presentation documents
```

---
