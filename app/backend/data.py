"""Load every precomputed artefact the dashboard serves.

Design notes
------------
Nothing here loads a model. The pipeline has already run; this module reads its
outputs. That keeps the demo fast to start and impossible to break on stage.

The screening table is 145k rows / ~101 MB on disk, so it is loaded lazily on
first request rather than at import, and the app degrades gracefully to the
shortlist alone if it is absent.

The reactive-alert catalogue is parsed out of ``src/safety.py`` with ``ast``
rather than copied. Copying it would drift the moment someone edits the real
alert set, and importing the module directly would drag in RDKit and build the
PAINS catalogs (several seconds) for data we could read statically.
"""

from __future__ import annotations

import ast
import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
PROC = ROOT / "data" / "processed"
REFERENCE = ROOT / "data" / "reference"
SRC = ROOT / "src"

TARGETS = ("TRPV1", "CB1", "FAAH")


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _clean(value: Any) -> Any:
    """Make a pandas value JSON-safe (NaN/NaT -> None, numpy scalar -> python)."""
    if value is None:
        return None
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (ValueError, AttributeError):
            return str(value)
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def records(df: pd.DataFrame) -> list[dict]:
    """DataFrame -> list of JSON-safe dicts."""
    return [{k: _clean(v) for k, v in row.items()}
            for row in df.to_dict(orient="records")]


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _read_csv(path: Path, **kw) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, **kw)


# ---------------------------------------------------------------------------
# static artefacts - small, cached at first touch
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def screening_stats() -> dict:
    return _read_json(RESULTS / "screening_stats.json")


@lru_cache(maxsize=1)
def dili_report() -> dict:
    return _read_json(RESULTS / "dili_report.json")


@lru_cache(maxsize=1)
def dili_probes() -> list[dict]:
    """Probe scores written by scripts/make_dili_probes.py (optional)."""
    data = _read_json(RESULTS / "dili_probes.json")
    return data.get("probes", [])


@lru_cache(maxsize=1)
def training_report() -> list[dict]:
    path = RESULTS / "training_report.json"
    if not path.exists():
        return []
    return json.loads(path.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def model_metrics() -> pd.DataFrame:
    return _read_csv(RESULTS / "model_metrics.csv")


@lru_cache(maxsize=1)
def attrition() -> pd.DataFrame:
    return _read_csv(PROC / "attrition_audit.csv")


@lru_cache(maxsize=1)
def dataset_summary() -> pd.DataFrame:
    return _read_csv(PROC / "dataset_summary.csv")


@lru_cache(maxsize=1)
def faah_rescue() -> dict:
    return {
        "summary": _read_json(RESULTS / "faah_rescue_summary.json"),
        "compounds": records(_read_csv(RESULTS / "faah_rescue_experiment.csv")),
    }


@lru_cache(maxsize=1)
def seed_plants() -> list[dict]:
    return records(_read_csv(REFERENCE / "pain_plants_seed.csv"))


@lru_cache(maxsize=1)
def literature_actives() -> list[dict]:
    return records(_read_csv(REFERENCE / "literature_actives.csv"))


@lru_cache(maxsize=1)
def alert_catalog() -> list[dict]:
    """Parse the ALERTS table out of src/safety.py without importing RDKit."""
    path = SRC / "safety.py"
    if not path.exists():
        return []
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:
        targets = getattr(node, "targets", []) or (
            [node.target] if hasattr(node, "target") else []
        )
        names = [t.id for t in targets if isinstance(t, ast.Name)]
        if "ALERTS" not in names:
            continue
        raw = ast.literal_eval(node.value)
        return [
            {"name": n, "smarts": s, "severity": sev, "note": note}
            for n, s, sev, note in raw
        ]
    return []


# ---------------------------------------------------------------------------
# compound tables
# ---------------------------------------------------------------------------
# Columns the explorer needs. reference_dois is excluded here and fetched on
# demand for a single compound: it is a pipe-joined DOI list, by far the widest
# column in the file, and nothing in the table view renders it.
EXPLORER_COLS = [
    "std_inchikey", "std_smiles", "compound_name", "pubchem_cid",
    "npclass_pathway", "npclass_superclass", "npclass_class",
    "source_species", "source_genera", "source_families",
    "n_source_organisms", "n_references", "tier", "seed_match_level",
    "traditional_pain_use", "relevant_target_hypothesis",
    "mw", "clogp", "tpsa", "qed", "fsp3", "heavy_atoms",
    "lipinski_violations", "cns_like",
    "n_reactive_alerts", "max_alert_severity", "alert_severity_sum",
    "reactive_alerts", "n_pains", "pains_flags",
    "p_CB1", "ad_CB1", "p_FAAH", "ad_FAAH", "p_TRPV1", "ad_TRPV1",
    "p_dili", "best_target", "target_score", "ad_best", "ad_confidence",
    "in_domain", "safety_multiplier", "cns_multiplier", "dili_multiplier",
    "composite_score",
]


def _load_compound_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=EXPLORER_COLS)
    header = pd.read_csv(path, nrows=0).columns
    usecols = [c for c in EXPLORER_COLS if c in header]
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df


@lru_cache(maxsize=1)
def shortlist() -> pd.DataFrame:
    return _load_compound_csv(RESULTS / "shortlist.csv")


@lru_cache(maxsize=1)
def candidates() -> pd.DataFrame:
    """The experiment-ready list: 50 DISTINCT molecules.

    shortlist.csv is the top 50 ROWS of the ranked screen, which is not the
    same thing as 50 compounds - the library stores stereoisomers, tautomers
    and salt variants separately, so those 50 rows hold only 27 distinct
    chemical skeletons. Handed to a laboratory that is 23 wasted assay wells.

    This file collapses to one row per InChIKey skeleton, keeps the
    best-scoring representative, and carries the provenance an
    experimentalist actually needs: source species, literature DOIs, the
    predicted target with its applicability domain, and the hepatic flag.
    Written by scripts/make_candidate_list.py.
    """
    return _load_compound_csv(RESULTS / "candidates_top50.csv")


@lru_cache(maxsize=1)
def candidates_summary() -> dict:
    return _read_json(RESULTS / "candidates_top50.json")


@lru_cache(maxsize=1)
def positive_controls() -> pd.DataFrame:
    return _load_compound_csv(RESULTS / "positive_controls.csv")


@lru_cache(maxsize=1)
def screen() -> pd.DataFrame:
    """The full 145k screen. Lazy: only paid for if someone opens the explorer."""
    return _load_compound_csv(RESULTS / "screening_all.csv.gz")


def compound_references(inchikey: str) -> list[str]:
    """Pull the DOI list for one compound by streaming the screen file.

    Chunked rather than held in memory: the DOI column roughly doubles the
    footprint of the screening table and is only ever needed one row at a time.
    """
    out: list[str] = []
    for path in (RESULTS / "screening_all.csv.gz", RESULTS / "positive_controls.csv"):
        if not path.exists():
            continue
        header = pd.read_csv(path, nrows=0).columns
        if "reference_dois" not in header:
            continue
        for chunk in pd.read_csv(
            path, usecols=["std_inchikey", "reference_dois"],
            chunksize=50_000, low_memory=False,
        ):
            hit = chunk[chunk["std_inchikey"] == inchikey]
            if len(hit):
                raw = hit.iloc[0]["reference_dois"]
                if isinstance(raw, str) and raw.strip():
                    out = [d for d in raw.split("|") if d.strip()]
                return out
    return out


@lru_cache(maxsize=1)
def literature_validation() -> dict:
    """External validation: literature actives scored by models that never saw them.

    This is the only genuinely held-out evaluation in the project. The
    percentile column is the one that matters — it asks the actual screening
    question ("would this compound surface above the library background?")
    rather than an AUC over 16 rows, which would mean nothing.
    """
    df = _read_csv(RESULTS / "literature_validation.csv")
    summary = _read_json(RESULTS / "literature_validation_summary.json")
    cols = [c for c in [
        "compound_name", "target", "interaction", "value_uM",
        "p_claimed_target", "library_percentile", "ad_claimed_target",
        "in_chembl_training", "n_reactive_alerts", "p_dili", "confidence",
    ] if c in df.columns]
    if len(df):
        df = df[cols].sort_values("library_percentile", ascending=False)
    return {
        "summary": summary,
        "held_out": records(df[df["in_chembl_training"] == False]) if len(df) else [],
        "seen": records(df[df["in_chembl_training"] == True]) if len(df) else [],
    }


@lru_cache(maxsize=1)
def validation_corrected() -> dict:
    """Corrected reading of the external validation.

    The raw summary counts all 16 held-out claims against a 0.5 probability
    threshold. Seven of those compounds are weaker than the cut their model was
    trained on, so a low score is correct rather than a miss, and 0.5 sits near
    the 99th percentile of the library rather than at any natural boundary.
    scripts/make_validation_summary.py recomputes both.
    """
    return _read_json(RESULTS / "validation_corrected.json")


@lru_cache(maxsize=1)
def pubchem_validation() -> dict:
    """Retrospective validation against PubChem BioAssay.

    The only analysis in the project with a real denominator: PubChem records
    Inactive outcomes as well as Active ones, so compounds that were screened
    and did nothing are in the data. Written by scripts/validate_pubchem.py
    then summarised by scripts/make_pubchem_summary.py.
    """
    return _read_json(RESULTS / "pubchem_summary.json")


@lru_cache(maxsize=1)
def pubchem_validation_expanded() -> dict:
    """The validation that made all three targets scoreable.

    Supersedes pubchem_validation() on two counts, both of which changed the
    numbers. First, outcomes from PubChem are keyed on the InChIKey skeleton,
    which collapses stereoisomers and salts, and the plant library holds
    several entries per skeleton - so the earlier set counted one measured
    compound up to several times. The CB1 overlap reported as n=71 is 31
    distinct structures. Everything here is deduplicated before any statistic
    is computed.

    Second, the earlier set required the test compound to be in the plant
    library, which is what left TRPV1 and FAAH single-class. Dropping that
    restriction - keeping "measured" and "never trained on" - turns 11
    scoreable TRPV1 structures into 966 and 2 FAAH into 621.

    Written by scripts/validate_pubchem_expanded.py.
    """
    return _read_json(RESULTS / "pubchem_validation_expanded.json")


@lru_cache(maxsize=1)
def hepatic_report() -> dict:
    """Training report for the three-head hepatic panel.

    Replaces the single TDC model: DILIrank 2.0 plus two Tox21 HepG2 assays,
    12,125 compounds against the old 462.
    """
    return _read_json(RESULTS / "hepatic_report.json")


@lru_cache(maxsize=1)
def hepatic_validation() -> dict:
    """External validation of the hepatic panel - the honest numbers.

    Two designs. A temporal holdout trains the clinical head on DILIrank v1
    and tests it on the drugs the FDA added in v2. A cross-dataset test scores
    the Tox21-trained heads, which have never seen a clinical label, against
    DILIrank. Both land near 0.65, well below the internal scaffold splits,
    and the file also records the 0.998 resubstitution score that shows why
    the paracetamol probe is not evidence of anything.

    Written by scripts/validate_hepatic.py.
    """
    return _read_json(RESULTS / "hepatic_validation.json")


@lru_cache(maxsize=1)
def hepatic_stacked() -> dict:
    """A pre-registered experiment that failed, kept because it failed.

    Ten Tox21 mechanism endpoints were turned into features for the clinical
    head. The target to beat was fixed in advance at 0.676. Adding them scored
    0.656. Reported rather than buried.

    Written by scripts/experiment_hepatic_stacked.py.
    """
    return _read_json(RESULTS / "hepatic_stacked_experiment.json")


@lru_cache(maxsize=1)
def kg_stats() -> dict:
    return _read_json(RESULTS / "kg_stats.json")


@lru_cache(maxsize=1)
def interpretation() -> dict:
    """Evidence behind the claims on the ML-decoder view.

    Computed rather than written down, so the page cannot drift from the data.

    The subtlety worth knowing: the 216 "positive controls" are library compounds
    that appeared *anywhere* in training, which includes compounds that were
    measured and found INACTIVE. Scoring them as one group understates the
    models badly, so they are split by what the curated table actually says
    about them.
    """
    out: dict[str, Any] = {}
    df = screen()
    sl = shortlist()
    pc = positive_controls()

    # --- applicability domain ------------------------------------------
    if len(df):
        out["domain"] = {
            "total": int(len(df)),
            "in_domain": int((df["in_domain"] == 1).sum()),
            "extrapolating": int((df["in_domain"] == 0).sum()),
            "in_domain_rate": round(float((df["in_domain"] == 1).mean()), 4),
            "median_ad": round(float(df["ad_best"].median()), 3),
        }

        # --- is the model discriminating, or calling everything active? ---
        disc = {}
        for t in TARGETS:
            col = f"p_{t}"
            if col not in df.columns:
                continue
            v = pd.to_numeric(df[col], errors="coerce").dropna()
            if not len(v):
                continue
            disc[t] = {
                "median": float(v.median()),
                "p99": round(float(v.quantile(0.99)), 4),
                "frac_above_half": round(float((v > 0.5).mean()), 4),
            }
        out["discrimination"] = disc

        # --- tier enrichment, against the pool the shortlist was drawn from ---
        eligible = df[(df["in_domain"] == 1)
                      & (df["max_alert_severity"].fillna(0) < 3)
                      & (df["n_pains"].fillna(0) == 0)]
        base = float((df["tier"] == "pain_seed").mean()) if len(df) else 0.0
        elig_base = float((eligible["tier"] == "pain_seed").mean()) if len(eligible) else 0.0
        short = float((sl["tier"] == "pain_seed").mean()) if len(sl) else 0.0
        out["enrichment"] = {
            "library_base_rate": round(base, 4),
            "eligible_pool": int(len(eligible)),
            "eligible_base_rate": round(elig_base, 4),
            "shortlist_rate": round(short, 4),
            "fold": round(short / elig_base, 2) if elig_base else None,
        }

    # --- positive controls, split by their TRAINING label ----------------
    agg_path = PROC / "bioactivity_aggregated.csv"
    if len(pc) and agg_path.exists():
        agg = pd.read_csv(agg_path, usecols=["inchikey", "label"])
        keys = set(pc["std_inchikey"])
        m = agg[agg["inchikey"].isin(keys)]
        per = m.groupby("inchikey")["label"].max()
        actives = set(per[per == 1].index)
        inactives = set(per[per == 0].index)

        def group(keyset: set) -> dict:
            sub = pc[pc["std_inchikey"].isin(keyset)]
            if not len(sub):
                return {"n": 0}
            s = pd.to_numeric(sub["target_score"], errors="coerce").dropna()
            return {
                "n": int(len(sub)),
                "median_score": round(float(s.median()), 4),
                "above_90": int((s > 0.9).sum()),
                "above_50": int((s > 0.5).sum()),
                "frac_above_90": round(float((s > 0.9).mean()), 4),
                "frac_above_50": round(float((s > 0.5).mean()), 4),
            }

        out["controls"] = {
            "total_held_out": int(len(pc)),
            "matched_in_curated_table": int(per.shape[0]),
            "measured_active": group(actives),
            "measured_inactive": group(inactives),
        }

    return out


def availability() -> dict:
    """Which artefacts are on disk - drives the warning banner in the UI."""
    checks = {
        "shortlist": RESULTS / "shortlist.csv",
        "candidates_top50": RESULTS / "candidates_top50.csv",
        "screening_all": RESULTS / "screening_all.csv.gz",
        "positive_controls": RESULTS / "positive_controls.csv",
        "screening_stats": RESULTS / "screening_stats.json",
        "model_metrics": RESULTS / "model_metrics.csv",
        "training_report": RESULTS / "training_report.json",
        "dili_report": RESULTS / "dili_report.json",
        "hepatic_report": RESULTS / "hepatic_report.json",
        "hepatic_validation": RESULTS / "hepatic_validation.json",
        "pubchem_validation_expanded": RESULTS / "pubchem_validation_expanded.json",
        "attrition_audit": PROC / "attrition_audit.csv",
        "dataset_summary": PROC / "dataset_summary.csv",
        "literature_validation": RESULTS / "literature_validation.csv",
        "validation_corrected": RESULTS / "validation_corrected.json",
        "pubchem_validation": RESULTS / "pubchem_summary.json",
        "knowledge_graph": RESULTS / "knowledge_graph.html",
        "report": RESULTS / "REPORT.md",
    }
    return {name: path.exists() for name, path in checks.items()}
