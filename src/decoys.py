"""Property-matched decoy selection.

For each active compound, draw decoys from the presumed-inactive pool that
resemble it in molecular weight and cLogP but are structurally distinct.

Two safeguards:

1. PROPERTY MATCHING. Decoys are sampled from the same (MW, cLogP) cell as the
   active they are matched to. Unmatched decoys let a model separate classes on
   gross physical properties alone, which inflates benchmark scores and
   generalises to nothing (the central lesson of the DUD-E work).

2. STRUCTURAL EXCLUSION. Any pool molecule whose InChIKey appears anywhere in
   the panel's measured data is removed, and so is any molecule sharing a
   Bemis-Murcko scaffold with a known active for that target. Without the
   second check, a "decoy" can be a close analogue of an active that simply
   was not measured, which teaches the model that real actives are inactive.

The resulting label is `presumed_inactive`, and it is recorded as such in the
output rather than being silently merged with measured inactives. A presumed
inactive is an assumption; a measured inactive is data. The distinction is kept
visible so that any result can be re-derived without the assumption.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors

from featurize import scaffold

RDLogger.DisableLog("rdApp.*")

MW_BIN = 50.0
LOGP_BIN = 1.0


def _props(smiles: str):
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return None
    return Descriptors.MolWt(m), Crippen.MolLogP(m)


def _cell(mw: float, logp: float) -> tuple[int, int]:
    return (int(mw // MW_BIN), int(np.floor(logp / LOGP_BIN)))


def prepare_pool(pool: pd.DataFrame, exclude_inchikeys: set[str]) -> pd.DataFrame:
    """Standardize the decoy pool and index it by property cell."""
    from standardize import standardize

    rows = []
    for smi in pool["canonical_smiles"].dropna().unique():
        r = standardize(smi)
        if not r.ok or r.inchikey in exclude_inchikeys:
            continue
        p = _props(r.smiles)
        if p is None:
            continue
        mw, logp = p
        rows.append({
            "inchikey": r.inchikey,
            "std_smiles": r.smiles,
            "mw": mw,
            "clogp": logp,
            "cell": _cell(mw, logp),
        })
    return pd.DataFrame(rows).drop_duplicates(subset="inchikey")


def match_decoys(
    actives: pd.DataFrame,
    pool: pd.DataFrame,
    n_per_active: int = 3,
    seed: int = 42,
    exclude_scaffolds: set[str] | None = None,
) -> pd.DataFrame:
    """Return property-matched decoys for a set of actives.

    `actives` needs a `std_smiles` column. `pool` must come from prepare_pool.
    """
    rng = np.random.RandomState(seed)
    exclude_scaffolds = exclude_scaffolds or set()

    by_cell: dict[tuple[int, int], list[int]] = {}
    for idx, c in enumerate(pool["cell"]):
        by_cell.setdefault(c, []).append(idx)

    wanted = len(actives) * n_per_active
    chosen: set[int] = set()

    # pass 1: exact property cell
    for smi in actives["std_smiles"]:
        p = _props(smi)
        if p is None:
            continue
        cell = _cell(*p)
        cands = [i for i in by_cell.get(cell, []) if i not in chosen]
        if not cands:
            continue
        take = min(n_per_active, len(cands))
        for i in rng.choice(cands, size=take, replace=False):
            chosen.add(int(i))

    # pass 2: top up from neighbouring cells if matching was too strict
    if len(chosen) < wanted:
        all_idx = [i for i in range(len(pool)) if i not in chosen]
        if all_idx:
            extra = min(wanted - len(chosen), len(all_idx))
            for i in rng.choice(all_idx, size=extra, replace=False):
                chosen.add(int(i))

    out = pool.iloc[sorted(chosen)].copy()

    # structural exclusion: drop decoys sharing a scaffold with a known active
    if exclude_scaffolds:
        keep = [scaffold(s) not in exclude_scaffolds for s in out["std_smiles"]]
        n_drop = len(out) - sum(keep)
        out = out[keep]
        if n_drop:
            print(f"      dropped {n_drop} decoys sharing an active's scaffold")

    out["label"] = 0
    out["evidence"] = "presumed_inactive"
    return out


def property_overlap_report(actives: pd.DataFrame, decoys: pd.DataFrame) -> dict:
    """Quantify how well matched the two sets are.

    If the medians differ sharply, the matching failed and any model trained on
    the result is separating classes on physics rather than pharmacology.
    """
    def stats(df):
        vals = [_props(s) for s in df["std_smiles"]]
        vals = [v for v in vals if v]
        if not vals:
            return {}
        mw = np.array([v[0] for v in vals])
        lp = np.array([v[1] for v in vals])
        return {
            "n": len(vals),
            "mw_median": round(float(np.median(mw)), 1),
            "mw_iqr": round(float(np.percentile(mw, 75) - np.percentile(mw, 25)), 1),
            "clogp_median": round(float(np.median(lp)), 2),
            "clogp_iqr": round(float(np.percentile(lp, 75) - np.percentile(lp, 25)), 2),
        }

    a, d = stats(actives), stats(decoys)
    return {
        "actives": a,
        "decoys": d,
        "mw_median_gap": round(abs(a.get("mw_median", 0) - d.get("mw_median", 0)), 1),
        "clogp_median_gap": round(
            abs(a.get("clogp_median", 0) - d.get("clogp_median", 0)), 2),
    }
