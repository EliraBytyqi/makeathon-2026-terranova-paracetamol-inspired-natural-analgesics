"""Terranova dashboard API.

Serves the precomputed pipeline outputs to the frontend. No model is loaded
anywhere in this process - every number returned here was produced by
``src/run_all.py`` and read off disk.

Run with:
    python -m app.run
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from . import data

FRONTEND = Path(__file__).resolve().parents[1] / "frontend"

app = FastAPI(
    title="Terranova",
    description="Plant-derived analgesic candidates on the paracetamol/AM404 axis.",
    version="1.0.0",
)


# ---------------------------------------------------------------------------
# meta
# ---------------------------------------------------------------------------
@app.get("/api/health")
def health() -> dict:
    return {"status": "ok", "artefacts": data.availability()}


@app.get("/api/overview")
def overview() -> dict:
    """Headline numbers for the landing view."""
    stats = data.screening_stats()
    summary = data.dataset_summary()
    metrics = data.model_metrics()

    measured_compounds = (
        int(summary[summary["scope"] == "all"]["compounds"].sum())
        if len(summary) else 0
    )
    scaffold_best = {}
    if len(metrics):
        sc = metrics[metrics["split"] == "scaffold"]
        for target in data.TARGETS:
            rows = sc[sc["target"] == target]
            if len(rows):
                best = rows.loc[rows["pr_auc"].idxmax()]
                scaffold_best[target] = {
                    "pr_auc": float(best["pr_auc"]),
                    "roc_auc": float(best["roc_auc"]),
                    "model": str(best["model"]),
                    "scope": str(best["scope"]),
                    "label": str(best["label"]),
                }

    attr = data.attrition()
    short = data.shortlist()

    return {
        "targets": list(data.TARGETS),
        "screening": stats,
        "measured_compounds": measured_compounds,
        "raw_activity_records": int(attr["rows"].iloc[0]) if len(attr) else 0,
        "curated_activity_records": int(attr["rows"].iloc[-1]) if len(attr) else 0,
        "scaffold_best": scaffold_best,
        "shortlist_size": int(len(short)),
        "shortlist_by_target": (
            short["best_target"].value_counts().to_dict() if len(short) else {}
        ),
        "shortlist_by_tier": (
            short["tier"].value_counts().to_dict() if len(short) else {}
        ),
        "n_models": len(list((data.ROOT / "models").glob("*.pkl"))),
    }


# ---------------------------------------------------------------------------
# data provenance
# ---------------------------------------------------------------------------
@app.get("/api/datasets")
def datasets() -> dict:
    return {
        "attrition": data.records(data.attrition()),
        "summary": data.records(data.dataset_summary()),
        "seed_plants": data.seed_plants(),
        "literature_actives": data.literature_actives(),
    }


@app.get("/api/models")
def models() -> dict:
    metrics = data.model_metrics()
    return {
        "metrics": data.records(metrics),
        "training_report": data.training_report(),
        "dili": data.dili_report(),
        "dili_probes": data.dili_probes(),
        "faah_rescue": data.faah_rescue(),
    }


@app.get("/api/validation")
def validation() -> dict:
    """External validation plus knowledge-graph statistics."""
    out = data.literature_validation()
    out["kg"] = data.kg_stats()
    out["corrected"] = data.validation_corrected()
    return out


@app.get("/kg")
def knowledge_graph() -> FileResponse:
    """The interactive pyvis knowledge graph, served as its own page."""
    path = data.RESULTS / "knowledge_graph.html"
    if not path.exists():
        raise HTTPException(
            status_code=404,
            detail="knowledge_graph.html not found — run python src/run_all.py kg.py",
        )
    return FileResponse(path, media_type="text/html")


@app.get("/api/interpretation")
def interpretation() -> dict:
    """Computed evidence behind the claims on the ML-decoder view."""
    return data.interpretation()


@app.get("/api/safety")
def safety() -> dict:
    stats = data.screening_stats()
    return {
        "alerts": data.alert_catalog(),
        "flagged": {
            "severity3": stats.get("severity3_flagged"),
            "pains": stats.get("pains_flagged"),
            "screened": stats.get("screened"),
        },
        "weights": stats.get("weights", {}),
    }


# ---------------------------------------------------------------------------
# compound tables
# ---------------------------------------------------------------------------
SORTABLE = {
    "rank", "composite_score", "target_score", "ad_best", "mw", "clogp",
    "qed", "p_dili", "p_TRPV1", "p_CB1", "p_FAAH", "n_reactive_alerts",
    "n_references", "tpsa",
}


def _frame(source: str) -> pd.DataFrame:
    if source == "shortlist":
        return data.shortlist()
    if source == "controls":
        return data.positive_controls()
    if source == "screen":
        return data.screen()
    raise HTTPException(status_code=404, detail=f"unknown source '{source}'")


@app.get("/api/compounds/{source}")
def compounds(
    source: str,
    q: str = Query("", description="name / genus / class substring"),
    target: str = Query("", description="TRPV1 | CB1 | FAAH"),
    tier: str = Query("", description="pain_seed | plant_background"),
    in_domain: bool = Query(False, description="applicability domain only"),
    no_alerts: bool = Query(False, description="exclude severity-3 and PAINS"),
    cns_only: bool = Query(False),
    sort: str = Query("rank"),
    desc: bool = Query(False),
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
) -> dict:
    df = _frame(source)
    if not len(df):
        return {"total": 0, "page": page, "page_size": page_size, "rows": []}

    if q:
        needle = q.strip().lower()
        cols = [c for c in ("compound_name", "source_genera", "npclass_class",
                            "source_species", "traditional_pain_use",
                            "std_inchikey")
                if c in df.columns]
        mask = pd.Series(False, index=df.index)
        for col in cols:
            mask |= df[col].astype(str).str.lower().str.contains(needle, na=False)
        df = df[mask]
    if target:
        df = df[df["best_target"] == target]
    if tier:
        df = df[df["tier"] == tier]
    if in_domain and "in_domain" in df.columns:
        df = df[df["in_domain"] == 1]
    if no_alerts:
        df = df[(df["max_alert_severity"].fillna(0) < 3)
                & (df["n_pains"].fillna(0) == 0)]
    if cns_only and "cns_like" in df.columns:
        df = df[df["cns_like"] == 1]

    total = int(len(df))
    key = sort if sort in SORTABLE and sort in df.columns else "rank"
    df = df.sort_values(key, ascending=not desc, kind="stable")

    start = (page - 1) * page_size
    page_df = df.iloc[start:start + page_size]
    return {
        "total": total,
        "page": page,
        "page_size": page_size,
        "sort": key,
        "desc": desc,
        "rows": data.records(page_df),
    }


@app.get("/api/compound/{inchikey}")
def compound(inchikey: str) -> dict:
    for source in ("screen", "shortlist", "controls"):
        df = _frame(source)
        if not len(df) or "std_inchikey" not in df.columns:
            continue
        hit = df[df["std_inchikey"] == inchikey]
        if len(hit):
            row = data.records(hit.head(1))[0]
            row["_source"] = source
            row["reference_dois"] = data.compound_references(inchikey)
            return row
    raise HTTPException(status_code=404, detail="compound not found")


@app.get("/api/distribution")
def distribution(
    field: str = Query("composite_score"),
    bins: int = Query(40, ge=5, le=100),
) -> dict:
    """Histogram of one numeric field across the full screen."""
    df = data.screen()
    if not len(df) or field not in df.columns:
        raise HTTPException(status_code=404, detail=f"no field '{field}'")
    series = pd.to_numeric(df[field], errors="coerce").dropna()
    if not len(series):
        return {"field": field, "bins": []}
    counts, edges = pd.cut(series, bins=bins, retbins=True)
    hist = counts.value_counts().sort_index()
    return {
        "field": field,
        "n": int(len(series)),
        "bins": [
            {"lo": float(edges[i]), "hi": float(edges[i + 1]), "count": int(c)}
            for i, c in enumerate(hist.to_list())
        ],
    }


# ---------------------------------------------------------------------------
# frontend
# ---------------------------------------------------------------------------
@app.get("/")
def index() -> FileResponse:
    target = FRONTEND / "index.html"
    if not target.exists():
        return JSONResponse(
            {"error": "frontend not built", "expected": str(target)},
            status_code=500,
        )
    return FileResponse(target)


# The knowledge graph references lib/vis-9.1.2/... relative to itself, so those
# vendored assets have to be reachable at /lib for /kg to render. Mounted before
# the catch-all frontend mount, which would otherwise swallow the path.
_KG_LIB = data.RESULTS / "lib"
if _KG_LIB.exists():
    app.mount("/lib", StaticFiles(directory=_KG_LIB), name="kg-assets")

if FRONTEND.exists():
    app.mount("/", StaticFiles(directory=FRONTEND, html=True), name="frontend")
