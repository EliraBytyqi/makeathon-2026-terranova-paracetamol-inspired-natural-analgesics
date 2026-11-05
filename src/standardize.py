"""Structure standardization, shared by the ChEMBL and natural-product sides.

Both datasets must pass through the *same* normalization, otherwise the
train/screen comparison is corrupted: a compound could appear as a salt in
ChEMBL and a free base in LOTUS and never be recognised as the same molecule.
That failure mode is silent and it is the single most common way an NP
screening pipeline reports a "novel hit" that was in its own training set.

Pipeline per molecule:
  parse -> Cleanup (sanitize, normalize functional groups, disconnect metals)
        -> FragmentParent (strip salts/solvates, keep largest organic fragment)
        -> Uncharger (neutralize where chemically sensible)
        -> reionize
        -> canonical isomeric SMILES + InChIKey

Rejections are counted and reported rather than dropped silently.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

_cleanup = rdMolStandardize.CleanupParameters()
_lfc = rdMolStandardize.LargestFragmentChooser()
_uncharger = rdMolStandardize.Uncharger()

# Elements we allow. Anything else (metals, boron, heavy halogens in odd
# contexts) is out of domain for a Morgan-fingerprint model trained on
# drug-like organics.
ALLOWED = {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "Se"}

MIN_HEAVY = 6     # below this is a fragment/solvent, not a lead
MAX_HEAVY = 70    # above this is a saponin/polysaccharide; fingerprints degrade
MAX_MW = 1000.0


@dataclass
class StdResult:
    ok: bool
    smiles: str | None = None
    inchikey: str | None = None
    reason: str | None = None


def standardize(smiles: str) -> StdResult:
    if not smiles or not isinstance(smiles, str):
        return StdResult(False, reason="empty_smiles")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return StdResult(False, reason="unparseable")

    try:
        mol = rdMolStandardize.Cleanup(mol, _cleanup)
        mol = _lfc.choose(mol)
        mol = _uncharger.uncharge(mol)
        mol = rdMolStandardize.Reionize(mol)
        Chem.SanitizeMol(mol)
    except Exception:
        return StdResult(False, reason="standardize_failed")

    if mol.GetNumAtoms() == 0:
        return StdResult(False, reason="empty_after_cleanup")

    syms = {a.GetSymbol() for a in mol.GetAtoms()}
    bad = syms - ALLOWED
    if bad:
        return StdResult(False, reason=f"disallowed_element:{'|'.join(sorted(bad))}")

    heavy = mol.GetNumHeavyAtoms()
    if heavy < MIN_HEAVY:
        return StdResult(False, reason="too_small")
    if heavy > MAX_HEAVY:
        return StdResult(False, reason="too_large")

    mw = Descriptors.MolWt(mol)
    if mw > MAX_MW:
        return StdResult(False, reason="mw_too_high")

    # Mixtures that survived fragment selection as a single component are fine;
    # anything still containing a dot is a genuine multi-component record.
    can = Chem.MolToSmiles(mol, isomericSmiles=True)
    if "." in can:
        return StdResult(False, reason="multicomponent")

    try:
        ikey = Chem.MolToInchiKey(mol)
    except Exception:
        return StdResult(False, reason="inchikey_failed")
    if not ikey:
        return StdResult(False, reason="inchikey_empty")

    return StdResult(True, smiles=can, inchikey=ikey)


def standardize_series(smiles_iter):
    """Vectorised-ish helper returning (smiles, inchikey, reason) tuples."""
    out = []
    for s in smiles_iter:
        r = standardize(s)
        out.append((r.smiles, r.inchikey, r.reason))
    return out


if __name__ == "__main__":
    # sanity checks, including the molecule that motivates the whole project
    tests = {
        "paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "paracetamol HCl salt form": "CC(=O)Nc1ccc(O)cc1.Cl",
        "NAPQI": "CC(=O)N=C1C=CC(=O)C=C1",
        "capsaicin": "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
        "sodium salt": "[Na+].CC(=O)[O-]",
        "bad": "not_a_smiles",
        "tiny": "CCO",
    }
    for name, smi in tests.items():
        r = standardize(smi)
        print(f"{name:28s} ok={r.ok!s:5s} {r.inchikey or r.reason}")
