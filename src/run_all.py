"""Run the whole pipeline end to end, in dependency order.

Each stage is skipped if its output already exists, so this is safe to re-run.
Pass --force to rebuild everything.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
PY = str(ROOT / ".venv" / "Scripts" / "python.exe")

# (script, output that proves it ran)
STAGES = [
    ("resolve_targets.py", ROOT / "data/raw/targets.json"),
    ("fetch_chembl.py", ROOT / "data/raw/chembl_activities_raw.csv"),
    ("fetch_decoys.py", ROOT / "data/raw/decoy_pool.csv"),
    ("build_dataset.py", ROOT / "data/processed/bioactivity_aggregated.csv"),
    ("build_np_library.py", ROOT / "data/interim/np_library_raw.csv"),
    ("prepare_np.py", ROOT / "data/processed/np_library.csv"),
    ("resolve_literature.py", ROOT / "data/interim/literature_actives_resolved.csv"),
    ("build_hepatic_dataset.py", ROOT / "data/processed/hepatic_endpoints.csv"),
    ("train_hepatic.py", ROOT / "models/hepatic_panel.pkl"),
    ("train.py", ROOT / "results/training_report.json"),
    ("screen.py", ROOT / "results/screening_all.csv.gz"),
    ("../scripts/make_candidate_list.py", ROOT / "results/candidates_top50.csv"),
    ("validate_literature.py", ROOT / "results/literature_validation.csv"),
    ("kg.py", ROOT / "results/knowledge_graph.graphml"),
    ("make_report.py", ROOT / "results/REPORT.md"),
]


def main():
    force = "--force" in sys.argv
    only = [a for a in sys.argv[1:] if not a.startswith("--")]

    for script, output in STAGES:
        if only and script not in only:
            continue
        if output.exists() and not force:
            print(f"[skip] {script}  ({output.relative_to(ROOT)} exists)")
            continue
        print(f"\n{'='*70}\n[run ] {script}\n{'='*70}", flush=True)
        r = subprocess.run([PY, script], cwd=SRC)
        if r.returncode != 0:
            print(f"\n[FAIL] {script} exited {r.returncode}")
            sys.exit(r.returncode)
    print("\nall stages complete")


if __name__ == "__main__":
    main()
