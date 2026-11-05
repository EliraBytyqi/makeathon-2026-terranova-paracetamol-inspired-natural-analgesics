"""Diagnostic: why does the FAAH model score known flavonoid inhibitors at ~0?

The smoke test showed the FAAH model assigning kaempferol p=6e-6 and genistein
p=6.5e-5, despite both being reported FAAH inhibitors in the low micromolar
range. Before treating that as a finding we need to know which of three things
is happening:

  (a) they are in the FAAH training set labelled ACTIVE and the model is simply
      wrong  -> model failure, needs fixing
  (b) they are in the FAAH training set labelled INACTIVE (because their
      potency sits above our cut)  -> labelling choice, explainable
  (c) they are not in the FAAH training set at all  -> extrapolation, and the
      question becomes whether ChEMBL's FAAH chemical space contains any
      flavonoids to learn from

This script answers that directly from the curated data.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from standardize import standardize

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"

PROBES = {
    "kaempferol": "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
    "genistein": "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12",
    "biochanin A": "COc1ccc(-c2coc3cc(O)cc(O)c3c2=O)cc1",
    "daidzein": "O=c1c(-c2ccc(O)cc2)coc2cc(O)ccc12",
    "7-hydroxyflavone": "O=c1cc(-c2ccccc2)oc2cc(O)ccc12",
    "silybin": "COc1cc(C2Oc3ccc(C4Oc5cc(O)cc(O)c5C(=O)C4O)cc3OC2CO)ccc1O",
    # URB597 = cyclohexylcarbamic acid 3'-carbamoylbiphenyl-3-yl ester.
    # The reference FAAH inhibitor, included as a positive control: if this is
    # missing too, the gap is in our curation, not in ChEMBL's coverage.
    "URB597": "NC(=O)c1cccc(-c2cccc(OC(=O)NC3CCCCC3)c2)c1",
    "quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
}


def main():
    agg = pd.read_csv(PROC / "bioactivity_aggregated.csv")
    faah = agg[agg["target_name"] == "FAAH"]
    print(f"FAAH rows in curated data: {len(faah):,} "
          f"({faah['inchikey'].nunique():,} compounds)\n")

    print(f"{'compound':18s} {'in FAAH data':>13s} {'pAct':>6s} "
          f"{'label':>6s} {'10uM':>5s} {'organism':>8s}")
    for name, smi in PROBES.items():
        ik = standardize(smi).inchikey
        hit = faah[faah["inchikey"] == ik]
        if len(hit):
            for _, r in hit.iterrows():
                print(f"{name:18s} {'yes':>13s} {r['pactivity']:6.2f} "
                      f"{str(r['label']):>6s} {str(r['label_10um']):>5s} "
                      f"{r['organism_group']:>8s}")
        else:
            print(f"{name:18s} {'NO':>13s} {'-':>6s} {'-':>6s} {'-':>5s} {'-':>8s}")

    # How much flavonoid chemistry does ChEMBL's FAAH space contain at all?
    # Crude but informative: count compounds matching a flavone/chromen-4-one core.
    from rdkit import Chem
    core = Chem.MolFromSmarts("O=c1cc(-c2ccccc2)oc2ccccc12")   # flavone
    core2 = Chem.MolFromSmarts("O=c1c(-c2ccccc2)coc2ccccc12")  # isoflavone
    core3 = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc21")   # flavanone

    n_flav = 0
    n_flav_active = 0
    for _, r in faah.iterrows():
        m = Chem.MolFromSmiles(r["std_smiles"])
        if m is None:
            continue
        if (m.HasSubstructMatch(core) or m.HasSubstructMatch(core2)
                or m.HasSubstructMatch(core3)):
            n_flav += 1
            if r["label"] == 1:
                n_flav_active += 1

    print(f"\nflavonoid-scaffold compounds in ChEMBL FAAH data: {n_flav} "
          f"(of {len(faah)} rows)")
    print(f"   of which labelled active at 1 uM: {n_flav_active}")
    print("\nIf this count is near zero, the model has essentially never seen a")
    print("flavonoid engage FAAH, and scoring kaempferol low is extrapolation")
    print("rather than a considered negative judgement.")


if __name__ == "__main__":
    main()
