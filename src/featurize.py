"""Molecular featurization and scaffold assignment.

Features: ECFP4 (Morgan radius 2, 2048 bits) concatenated with a compact
physicochemical descriptor block. This is a deliberately boring, strong
baseline. A graph neural network on ~10k compounds per target would not beat
it reliably and would cost the whole hackathon.

Scaffolds: Bemis-Murcko frameworks, used for grouped splitting. This is the
single most important methodological choice in the project. ChEMBL is built
from medicinal-chemistry papers, each contributing a congeneric series of
close analogues. Under a random split those analogues straddle train and test,
the model memorises the series, and you report an ROC-AUC near 0.95 that
collapses on anything genuinely new - which is exactly what screening a
natural-product library asks it to do. Grouped scaffold splitting measures the
quantity we actually care about: generalisation to unseen chemotypes.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

FP_BITS = 2048
FP_RADIUS = 2
_morgan = GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)

DESCRIPTORS = [
    ("mw", Descriptors.MolWt),
    ("clogp", Crippen.MolLogP),
    ("tpsa", Descriptors.TPSA),
    ("hbd", Descriptors.NumHDonors),
    ("hba", Descriptors.NumHAcceptors),
    ("rotb", Descriptors.NumRotatableBonds),
    ("aromatic_rings", Descriptors.NumAromaticRings),
    ("rings", Descriptors.RingCount),
    ("fsp3", Descriptors.FractionCSP3),
    ("heavy_atoms", lambda m: float(m.GetNumHeavyAtoms())),
    ("qed", Descriptors.qed),
    ("molmr", Crippen.MolMR),
    ("heteroatoms", Descriptors.NumHeteroatoms),
    ("stereocenters", lambda m: float(len(Chem.FindMolChiralCenters(
        m, includeUnassigned=True, useLegacyImplementation=False)))),
]
DESC_NAMES = [n for n, _ in DESCRIPTORS]
FEATURE_NAMES = [f"fp_{i}" for i in range(FP_BITS)] + DESC_NAMES


def fingerprint(mol: Chem.Mol) -> np.ndarray:
    return np.array(_morgan.GetFingerprintAsNumPy(mol), dtype=np.uint8)


def descriptors(mol: Chem.Mol) -> np.ndarray:
    out = []
    for _, fn in DESCRIPTORS:
        try:
            v = float(fn(mol))
            if not np.isfinite(v):
                v = 0.0
        except Exception:
            v = 0.0
        out.append(v)
    return np.array(out, dtype=np.float32)


def featurize_smiles(smiles_list, with_descriptors: bool = True):
    """Return (X, valid_mask). Invalid molecules get a zero row and mask=False."""
    n = len(smiles_list)
    width = FP_BITS + (len(DESCRIPTORS) if with_descriptors else 0)
    X = np.zeros((n, width), dtype=np.float32)
    mask = np.zeros(n, dtype=bool)
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if mol is None:
            continue
        row = fingerprint(mol).astype(np.float32)
        if with_descriptors:
            row = np.concatenate([row, descriptors(mol)])
        X[i] = row
        mask[i] = True
    return X, mask


def scaffold(smiles: str, generic: bool = False) -> str:
    """Bemis-Murcko scaffold SMILES. Acyclic molecules return "" (own group)."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return ""
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        if generic:
            core = MurckoScaffold.MakeScaffoldGeneric(core)
        return Chem.MolToSmiles(core)
    except Exception:
        return ""


def scaffold_groups(smiles_list) -> np.ndarray:
    """Group id per molecule. Acyclic/failed molecules each get a unique group
    so they are never used to leak information between train and test."""
    groups = np.empty(len(smiles_list), dtype=object)
    lookup: dict[str, int] = {}
    nxt = 0
    for i, smi in enumerate(smiles_list):
        s = scaffold(smi)
        if not s:
            groups[i] = f"__singleton_{i}"
            continue
        if s not in lookup:
            lookup[s] = nxt
            nxt += 1
        groups[i] = s
    return groups


if __name__ == "__main__":
    probes = [
        ("capsaicin", "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O"),
        ("kaempferol", "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
        ("genistein", "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12"),
        ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ]
    X, mask = featurize_smiles([s for _, s in probes])
    print("feature matrix:", X.shape, "valid:", mask.sum())
    print(f"\n{'compound':14s} {'bits_on':>7s}  scaffold")
    for (name, smi), row in zip(probes, X):
        print(f"{name:14s} {int(row[:FP_BITS].sum()):7d}  {scaffold(smi)}")
    print("\nkaempferol vs genistein share a flavone scaffold ->",
          scaffold(probes[1][1]) == scaffold(probes[2][1]))
