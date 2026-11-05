"""Diagnostic: does the TDC DILI benchmark actually contain paracetamol,
and how is it labelled?

This matters for how the DILI layer may be interpreted. "DILI" in the clinical
literature often means *idiosyncratic* drug-induced liver injury - unpredictable,
host-dependent, not dose-proportional. Paracetamol hepatotoxicity is the
opposite: intrinsic, dose-dependent and mechanistically understood (NAPQI).
If the benchmark excludes or negatively labels paracetamol, then the DILI model
cannot see the very mechanism this project is built around, and the structural
alert layer is doing that work instead.
"""

from pathlib import Path

import pandas as pd

from standardize import standardize

ROOT = Path(__file__).resolve().parents[1]

PROBES = {
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "isoniazid": "NNC(=O)c1ccncc1",
    "ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "diclofenac": "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
    "troglitazone": "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2",
    "valproic acid": "CCCC(CCC)C(=O)O",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
}


def main():
    df = pd.read_csv(ROOT / "data" / "external" / "dili" / "tdc_dili.csv")
    df["inchikey"] = [(standardize(s).inchikey or "") for s in df["SMILES"]]

    print(f"{'compound':16s} {'in dataset':>11s} {'label':>7s}")
    for name, smi in PROBES.items():
        ik = standardize(smi).inchikey
        hit = df[df["inchikey"] == ik]
        if len(hit):
            labels = hit["Y"].tolist()
            print(f"{name:16s} {'yes':>11s} {str(labels):>7s}")
        else:
            print(f"{name:16s} {'NO':>11s} {'-':>7s}")

    print(f"\ndataset size: {len(df)}, positives: {int(df['Y'].sum())}")


if __name__ == "__main__":
    main()
