"""Assemble the hepatic-liability training data: DILIrank 2.0 + Tox21.

Replaces the 475-compound TDC DILI benchmark that the old safety model was
trained on. Three endpoints, all with measured human or cell-based outcomes:

  clinical  FDA DILIrank 2.0 - 1,336 approved drugs ranked by the evidence in
            their own FDA labels. 982 of them carry an unambiguous concern
            class (most / less / none); the remaining 354 are explicitly
            "ambiguous" and are dropped rather than guessed at.
            https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb

  mito      Tox21 AID 2283792 (TOX21_MMP_ratio) - mitochondrial membrane
            potential in HepG2, ~8,100 compounds, Active/Inactive only.
            Mitochondrial depolarisation is a direct intrinsic-toxicity
            mechanism, measured rather than inferred from a drug label.

  hepg2     Tox21 AID 2283888 (CellTox Green, 40 h) - outright death of HepG2
            liver cells, ~5,600 scoreable compounds.

Why all three. The clinical set says what actually happened to patients but is
small and covers only marketed drugs. The two Tox21 sets are an order of
magnitude larger and are real dose-response measurements in liver cells, but
they see only what a hepatoma line in a dish can show. Neither alone is
enough; the disagreements between them are informative in themselves.

Structures: DILIrank ships drug names, Tox21 ships CIDs, so both are resolved
through PubChem and then pushed through the project's own standardizer, so
every structure here is directly comparable to the screening library.

    python src/build_hepatic_dataset.py
"""

from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from standardize import standardize  # noqa: E402

DILIRANK_XLSX = ROOT / "data" / "external" / "dilirank" / "dilirank.xlsx"
DILIRANK_URL = "https://www.fda.gov/media/113052/download?attachment"
TOX21 = ROOT / "data" / "external" / "tox21"
CACHE = ROOT / "data" / "external" / "hepatic_cache"
OUT = ROOT / "data" / "processed"
UA = "terranova-academic-research/1.0"

# AID -> (short name, what it measures)
ASSAYS = {
    "mito": (2283792, "TOX21_MMP_ratio: mitochondrial membrane potential, HepG2"),
    "hepg2": (2283888, "TOX21_RT_HEPG2_FLO_40hr: CellTox Green cytotoxicity, HepG2"),
}

# salt and formulation words to strip if the full name does not resolve
SALT_WORDS = [
    "sulfate", "sulphate", "hydrochloride", "hydrobromide", "sodium",
    "potassium", "calcium", "magnesium", "maleate", "mesylate", "tartrate",
    "citrate", "acetate", "phosphate", "succinate", "fumarate", "besylate",
    "tosylate", "oxalate", "nitrate", "bromide", "chloride", "iodide",
    "dihydrate", "monohydrate", "hydrate", "anhydrous", "disodium",
    "dipotassium", "palmitate", "stearate", "valerate", "propionate",
    "dipropionate", "furoate", "pivalate", "lactate", "gluconate",
    "bitartrate", "pamoate", "napsylate", "edisylate", "camsylate",
]


def get(url: str, dest: Path, pause: float = 0.22) -> bytes | None:
    """Fetch with an on-disk cache. Returns None on failure (incl. 404)."""
    if dest.exists():
        if dest.stat().st_size == 0:
            return None                      # cached miss
        return dest.read_bytes()
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
        dest.write_bytes(data)
        time.sleep(pause)
        return data
    except Exception:
        dest.write_bytes(b"")                # cache the miss, don't retry it
        time.sleep(pause)
        return None


def smiles_from_payload(raw: bytes | None) -> dict[str, str]:
    """CID/key -> SMILES out of a PUG REST property payload."""
    if not raw:
        return {}
    try:
        props = json.loads(raw)["PropertyTable"]["Properties"]
    except Exception:
        return {}
    out = {}
    for p in props:
        s = p.get("SMILES") or p.get("CanonicalSMILES") or p.get("IsomericSMILES")
        if s:
            out[str(p.get("CID"))] = s
    return out


# ---------------------------------------------------------------------------
# 1. clinical endpoint: FDA DILIrank 2.0
# ---------------------------------------------------------------------------
def name_to_smiles(name: str) -> str | None:
    """Resolve a drug name through PubChem, retrying without its salt word."""
    tries = [name]
    low = name.lower()
    for w in SALT_WORDS:
        if low.endswith(" " + w):
            tries.append(name[: -(len(w) + 1)].strip())
            break
    for attempt in tries:
        slug = urllib.parse.quote(attempt, safe="")
        dest = CACHE / "names" / (
            "".join(c if c.isalnum() else "_" for c in attempt.lower())[:80] + ".json")
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
               f"{slug}/property/SMILES,CanonicalSMILES,IsomericSMILES/JSON")
        got = smiles_from_payload(get(url, dest))
        if got:
            return next(iter(got.values()))
    return None


def clinical() -> pd.DataFrame:
    if not DILIRANK_XLSX.exists():
        print(f"  downloading DILIrank -> {DILIRANK_XLSX}")
        DILIRANK_XLSX.parent.mkdir(parents=True, exist_ok=True)
        data = get(DILIRANK_URL, DILIRANK_XLSX)
        if not data:
            raise SystemExit(f"could not download DILIrank from {DILIRANK_URL}")

    d = pd.read_excel(DILIRANK_XLSX, sheet_name="version 2", header=1)
    d.columns = [str(c).strip() for c in d.columns]
    d = d.rename(columns={"CompoundName": "name", "SeverityClass": "severity",
                          "vDILI-Concern": "concern"})
    # the sheet has case drift in the concern labels ("vMOST" vs "vMost")
    d["concern"] = d["concern"].astype(str).str.strip().str.lower()
    d["name"] = d["name"].astype(str).str.strip()
    print(f"  DILIrank 2.0: {len(d)} drugs")
    print("   " + d["concern"].value_counts().to_string().replace("\n", "\n   "))

    # All three unambiguous classes are kept in the CSV. Which of them become
    # the positive and negative class is a modelling decision, made in
    # train_hepatic.py, not here - see the note there on why less-concern is
    # excluded rather than folded into the positives.
    keep = {"vmost-dili-concern", "vless-dili-concern", "vno-dili-concern"}
    d = d[d["concern"].isin(keep)].copy()
    d["most_concern"] = (d["concern"] == "vmost-dili-concern").astype(int)
    # y: the DILIst convention (most + less vs none), kept for comparison
    d["y"] = (d["concern"] != "vno-dili-concern").astype(int)
    d["severity"] = pd.to_numeric(d["severity"], errors="coerce").fillna(0)
    print(f"  unambiguous: {len(d)}  "
          f"most={int(d['most_concern'].sum())} "
          f"less={int((d['concern'] == 'vless-dili-concern').sum())} "
          f"none={int((d['concern'] == 'vno-dili-concern').sum())}")

    print("  resolving names through PubChem (cached; first run is slow)...")
    smi, miss = [], 0
    for i, n in enumerate(d["name"], 1):
        s = name_to_smiles(n)
        smi.append(s)
        if s is None:
            miss += 1
        if i % 200 == 0:
            print(f"    {i}/{len(d)} resolved, {miss} unresolved", flush=True)
    d["raw_smiles"] = smi
    before = len(d)
    d = d[d["raw_smiles"].notna()].copy()
    print(f"  structures found for {len(d)}/{before} "
          f"({before - len(d)} unresolved - mostly biologics and mixtures)")
    d["endpoint"] = "clinical"
    return d[["name", "raw_smiles", "y", "most_concern", "severity",
              "concern", "endpoint"]]


# ---------------------------------------------------------------------------
# 2+3. in-vitro endpoints: Tox21
# ---------------------------------------------------------------------------
def assay(short: str) -> pd.DataFrame:
    aid, label = ASSAYS[short]
    path = TOX21 / f"aid_{aid}.csv"
    if not path.exists():
        print(f"  downloading AID {aid}")
        if not get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"
                   f"{aid}/concise/CSV", path):
            raise SystemExit(f"could not download AID {aid}")

    d = pd.read_csv(path, low_memory=False)
    d = d[d["CID"].notna()].copy()
    d["CID"] = d["CID"].astype(int)
    print(f"  {short} (AID {aid}): {len(d)} rows, "
          f"{d['Activity Outcome'].value_counts().to_dict()}")

    # Tox21 screens most compounds several times. Inconclusive means the
    # dose-response curve could not be called - that is not a negative, so it
    # is dropped rather than folded into either class.
    d = d[d["Activity Outcome"].isin(["Active", "Inactive"])].copy()
    d["hit"] = (d["Activity Outcome"] == "Active").astype(int)

    # A compound counts as active if any replicate called it active: a real
    # curve in one run is better evidence than silence in another.
    g = d.groupby("CID").agg(y=("hit", "max"), n_runs=("hit", "size"),
                             ac50=("Activity Value [uM]", "median")).reset_index()
    print(f"    {len(g)} unique compounds, {int(g['y'].sum())} active "
          f"({g['y'].mean():.1%})")

    print(f"    fetching structures for {len(g)} CIDs...")
    cids = g["CID"].tolist()
    smap: dict[str, str] = {}
    for i in range(0, len(cids), 180):
        block = cids[i:i + 180]
        dest = CACHE / short / f"{i:06d}.json"
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
               + ",".join(map(str, block))
               + "/property/SMILES,CanonicalSMILES/JSON")
        smap.update(smiles_from_payload(get(url, dest)))
    g["raw_smiles"] = g["CID"].astype(str).map(smap)
    before = len(g)
    g = g[g["raw_smiles"].notna()].copy()
    print(f"    structures for {len(g)}/{before}")

    g["name"] = "CID" + g["CID"].astype(str)
    g["most_concern"] = g["y"]
    g["severity"] = 0.0
    g["concern"] = "n/a"
    g["endpoint"] = short
    return g[["name", "raw_smiles", "y", "most_concern", "severity",
              "concern", "endpoint"]]


# ---------------------------------------------------------------------------
def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    print("=== clinical endpoint: FDA DILIrank 2.0 ===")
    frames = [clinical()]
    for short in ASSAYS:
        print(f"\n=== in-vitro endpoint: {short} ===")
        frames.append(assay(short))

    df = pd.concat(frames, ignore_index=True)

    print(f"\n=== standardizing {len(df)} structures ===")
    res = [standardize(s) for s in df["raw_smiles"]]
    df["std_smiles"] = [r.smiles if r.ok else None for r in res]
    df["inchikey"] = [r.inchikey if r.ok else None for r in res]
    ok = df["std_smiles"].notna()
    print(f"  standardized OK: {int(ok.sum())}/{len(df)}")
    df = df[ok].copy()

    # one row per (endpoint, structure); keep the worse label on a collision
    df = (df.sort_values("y", ascending=False)
            .drop_duplicates(subset=["endpoint", "inchikey"], keep="first")
            .reset_index(drop=True))

    path = OUT / "hepatic_endpoints.csv"
    df.to_csv(path, index=False)
    print(f"\n=== final ===")
    for ep, s in df.groupby("endpoint"):
        print(f"  {ep:9s} n={len(s):5d}  positive={int(s['y'].sum()):5d} "
              f"({s['y'].mean():.1%})")
    print(f"  wrote {path}")

    # how much of the plant library's chemistry these endpoints actually cover
    # is answered in train_hepatic.py via the applicability domain; here we
    # only check the compound this whole project is about.
    para = df[df["inchikey"].str.startswith("RZVAJINKPMORJF")]
    print("\n  paracetamol in the training data:")
    if not len(para):
        print("    absent from all three endpoints")
    for _, r in para.iterrows():
        print(f"    {r['endpoint']:9s} label={int(r['y'])} "
              f"severity={r['severity']:.0f}  ({r['name']})")


if __name__ == "__main__":
    main()
