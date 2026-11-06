# Terranova dashboard

A local full-stack app that presents the finished pipeline: the hypothesis, the
data provenance, the model evaluation, the safety layers, and an explorer over
all 145,164 screened compounds.

```bash
.venv\Scripts\python.exe -m app.run
```

Then open <http://127.0.0.1:8000> (it opens by itself unless you pass
`--no-browser`).

## What it is, and what it deliberately is not

**It serves precomputed artefacts only.** No model is loaded in the web process.
Every number on screen was produced by `src/run_all.py` and read off disk, which
means the app starts in about a second, uses ~200 MB, and cannot fail on stage
because a 200 MB pickle was slow to unpickle.

**There is no build step and no CDN.** Plain HTML, CSS and ES2019, with the
charts hand-drawn as SVG. It renders on a venue network that resolves nothing.

## Layout

| Path | What it is |
|---|---|
| `app/run.py` | Launcher (`--host`, `--port`, `--reload`, `--no-browser`) |
| `app/backend/main.py` | FastAPI routes + static mount |
| `app/backend/data.py` | Artefact loaders, cached; lazy on the 101 MB screen file |
| `app/frontend/index.html` | All six views |
| `app/frontend/styles.css` | Design tokens, light + dark |
| `app/frontend/js/charts.js` | SVG chart primitives |
| `app/frontend/js/app.js` | Views, filtering, the explorer |

## API

| Endpoint | Returns |
|---|---|
| `GET /api/health` | Which artefacts are present on disk |
| `GET /api/overview` | Headline counts, scaffold-split bests, shortlist breakdown |
| `GET /api/datasets` | Attrition audit, per-target summary, seed plants |
| `GET /api/models` | All metrics, training report, DILI report + probes, FAAH rescue |
| `GET /api/safety` | Alert catalogue, flag counts, composite weights |
| `GET /api/compounds/{shortlist\|screen\|controls}` | Paginated, filterable, sortable |
| `GET /api/compound/{inchikey}` | One compound + its LOTUS DOIs |
| `GET /api/distribution?field=…` | Histogram over the full screen |

Query parameters on `/api/compounds/…`: `q`, `target`, `tier`, `in_domain`,
`no_alerts`, `cns_only`, `sort`, `desc`, `page`, `page_size`.

## Design decisions worth knowing

**The alert catalogue is parsed, not copied.** `data.py` reads the `ALERTS`
table out of `src/safety.py` with `ast`, so the Safety view cannot drift from
the real alert set. Importing the module directly would drag in RDKit and build
the PAINS catalogs for data we only need to read.

**Colour follows the target, never the rank.** TRPV1 is always blue, CB1 always
orange, FAAH always aqua. Filtering never repaints a survivor. The three hues
clear the colourblind and normal-vision separation floors in both light and dark
mode, and every chart also ships a table view, so nothing is encoded by colour
alone.

**DOIs are streamed, not held.** The reference-DOI column roughly doubles the
footprint of the screening table and is only ever needed one row at a time, so
`compound_references()` chunks the CSV instead of keeping it in memory.

## If a view is empty

The app reports missing artefacts in a banner at the top. Regenerate with:

```bash
python src/run_all.py              # full pipeline
python scripts/make_dili_probes.py # DILI probe scores for the Models view
```

> Everything the dashboard displays is a computational prediction. Nothing in it
> establishes that a compound is active, or that it is safe.
