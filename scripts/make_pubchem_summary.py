"""Condense the external validation into the summary the UI and PDFs read.

results/pubchem_validation_expanded.json is the full record: every target, at
four applicability-domain cut-offs, with bootstrap intervals. The dashboard
and both documents need a smaller, stabler shape than that, so it is derived
here rather than each consumer reaching into the full file.

Run after scripts/validate_pubchem_expanded.py.

    python scripts/make_pubchem_summary.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

GRADE = {
    "CB1": "validated, including on chemistry unlike training",
    "TRPV1": ("validated near training chemistry only; no measured inactives "
              "survive an AD < 0.4 filter"),
    "FAAH": "measured and weak - better than chance, not much more",
}
SCOPES = {
    "all": "external_all",
    "novel_ad_lt_0.5": "external_novel_ad_lt_0.5",
    "novel_ad_lt_0.4": "external_novel_ad_lt_0.4",
    "in_plant_library": "in_library",
}


def block(d: dict | None) -> dict | None:
    """Keep only the fields a reader needs; drop a single-class scope."""
    if not d or not d.get("scoreable"):
        return None
    return {"n": d["n"], "n_active": d["n_active"], "n_inactive": d["n_inactive"],
            "roc_auc": d["roc_auc"], "ci95": d.get("ci95"),
            "p": d.get("mannwhitney_p")}


def main() -> None:
    src = RESULTS / "pubchem_validation_expanded.json"
    if not src.exists():
        raise SystemExit(f"missing {src}\n"
                         "run: python scripts/validate_pubchem_expanded.py")
    targets = json.loads(src.read_text(encoding="utf-8")).get("targets", {})
    if not targets:
        raise SystemExit(f"{src.name} has no target results")

    out: dict = {
        "source": "PubChem BioAssay, via PUG REST over all nine panel accessions",
        "design": ("structures with an independent measured Active/Inactive "
                   "outcome, absent from every production training set, "
                   "deduplicated to one row per InChIKey skeleton"),
        "total_scored": sum(t["external_all"]["n"] for t in targets.values()
                            if t.get("external_all", {}).get("scoreable")),
        "targets": {},
    }

    for name in ("CB1", "TRPV1", "FAAH"):
        t = targets.get(name)
        if not t:
            continue
        out["targets"][name] = {"grade": GRADE.get(name, "")} | {
            key: block(t.get(src_key)) for key, src_key in SCOPES.items()
        }

    # legacy key: older readers expect a top-level cb1 object
    if "CB1" in out["targets"]:
        out["cb1"] = {k: v for k, v in out["targets"]["CB1"].items()
                      if k != "grade"}

    dest = RESULTS / "pubchem_summary.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print(f"{'target':8s} {'n':>6s} {'act/inact':>12s} {'ROC-AUC':>8s}  grade")
    for name, t in out["targets"].items():
        a = t.get("all")
        if not a:
            print(f"{name:8s} {'-':>6s} {'-':>12s} {'-':>8s}  not scoreable")
            continue
        print(f"{name:8s} {a['n']:6,d} {a['n_active']:5d}/{a['n_inactive']:<6d} "
              f"{a['roc_auc']:8.3f}  {t['grade']}")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
