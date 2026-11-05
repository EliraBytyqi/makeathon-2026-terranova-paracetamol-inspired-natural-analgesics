"""Hepatic-liability screening: structural alerts for reactive-metabolite formation.

Motivation
----------
Paracetamol's analgesia is not what hurts the liver. The damage comes from a
minor CYP2E1/CYP3A4-mediated oxidation of the para-aminophenol core to NAPQI
(N-acetyl-p-benzoquinone imine), a hard electrophile. At therapeutic dose,
hepatic glutathione conjugates it. In overdose the glutathione pool is
depleted, NAPQI arylates cellular protein thiols, and you get oxidative stress,
mitochondrial dysfunction and centrilobular necrosis.

So the structural question is specific and answerable: *can this scaffold be
oxidised to a quinone/quinone-imine or otherwise act as a soft or hard
electrophile?* That is what the alert set below encodes.

Two honest caveats, which belong in the writeup:

1. Structural alerts are a hazard flag, not a risk assessment. Dose, exposure,
   clearance route and detoxification capacity all matter and none are encoded
   here. A flagged compound is not "toxic"; it is "worth a reactivity check".
2. The alert set flags several compounds this project likes. Curcumin and
   6-gingerol are Michael acceptors. Many flavonoids are catechols. That is a
   genuine finding about the chemical space of plant analgesics, not a bug to
   be tuned away - plants make electrophiles on purpose, often as defence
   chemistry. We report the flags and keep the compounds, tagged.

Also included: a PAINS/Brenk/NIH filter-catalog pass, because an NP screening
hit list that has not been checked for assay-interference chemistry
(catechols, quinones, rhodanines) is not credible. Many polyphenols are
frequent hitters.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# Reactive-metabolite / hepatic liability alerts
# severity: 3 = direct NAPQI-type bioactivation precedent
#           2 = established reactive-metabolite route
#           1 = context-dependent / weaker concern
# ---------------------------------------------------------------------------
ALERTS: list[tuple[str, str, int, str]] = [
    # --- the paracetamol axis itself -------------------------------------
    (
        "para_aminophenol",
        "[OX2H]c1ccc([NX3;!$(N=*)])cc1",
        3,
        "direct NAPQI-type quinone-imine precursor; the paracetamol motif",
    ),
    (
        "anilide_para_hydroxy",
        "[OX2H]c1ccc(NC(=O)[#6])cc1",
        3,
        "acetaminophen-like acylated aminophenol; oxidises to quinone imine",
    ),
    (
        "quinone_imine",
        "O=C1C=CC(=N)C=C1",
        3,
        "pre-formed quinone imine - the reactive species itself",
    ),
    # --- quinone / polyphenol oxidation routes ----------------------------
    (
        "para_quinone",
        "O=C1C=CC(=O)C=C1",
        3,
        "quinone; redox cycling and thiol arylation",
    ),
    (
        "ortho_quinone",
        "O=C1C(=O)C=CC=C1",
        3,
        "ortho-quinone; arylates protein thiols",
    ),
    (
        "catechol",
        "[OX2H]c1ccccc1[OX2H]",
        2,
        "oxidises to ortho-quinone; also a frequent assay interferer",
    ),
    (
        "hydroquinone",
        "[OX2H]c1ccc([OX2H])cc1",
        2,
        "oxidises to para-quinone",
    ),
    (
        "aniline_primary",
        "[NX3;H2][c]",
        2,
        "aromatic amine; nitrenium/quinone-imine bioactivation",
    ),
    (
        "nitroaromatic",
        "[$([NX3](=O)=O),$([NX3+](=O)[O-])][c]",
        2,
        "nitroreduction to nitroso/hydroxylamine",
    ),
    # --- direct electrophiles ---------------------------------------------
    (
        "michael_acceptor_enone",
        "[CX3]=[CX3][CX3](=O)[#6,#1]",
        2,
        "alpha,beta-unsaturated carbonyl; direct thiol Michael addition",
    ),
    (
        "michael_acceptor_ester",
        "[CX3]=[CX3][CX3](=O)[OX2][#6]",
        1,
        "acrylate-type Michael acceptor",
    ),
    (
        "alpha_methylene_lactone",
        "C=C1CCOC1=O",
        2,
        "sesquiterpene lactone warhead; potent thiol alkylator",
    ),
    (
        "epoxide",
        "C1OC1",
        2,
        "strained electrophile; alkylates nucleophiles directly",
    ),
    (
        "aldehyde",
        "[CX3H1](=O)[#6]",
        1,
        "forms Schiff bases with lysine; context dependent",
    ),
    (
        "alpha_halo_carbonyl",
        "[CX3](=O)[CX4][F,Cl,Br,I]",
        2,
        "direct alkylating agent",
    ),
    (
        "acyl_halide",
        "[CX3](=O)[F,Cl,Br,I]",
        3,
        "acylating agent; hydrolytically and chemically reactive",
    ),
    (
        "isocyanate_isothiocyanate",
        "[NX2]=[CX2]=[OX1,SX1]",
        2,
        "carbamoylates protein nucleophiles",
    ),
    # --- CYP bioactivation of small heteroaromatics -----------------------
    (
        "furan",
        "c1ccoc1",
        2,
        "CYP epoxidation to reactive cis-enedione (furan hepatotoxicity route)",
    ),
    (
        "thiophene",
        "c1ccsc1",
        2,
        "CYP oxidation to S-oxide / thiophene epoxide",
    ),
    (
        "terminal_alkyne",
        "[CX2;H1]#[CX2]",
        1,
        "mechanism-based CYP inactivation via ketene",
    ),
    # --- classic hepatotoxicophores ---------------------------------------
    (
        "hydrazine",
        "[NX3][NX3]",
        2,
        "hydrazine/hydrazide; isoniazid-type hepatotoxicity",
    ),
    (
        "pyrrolizidine_core",
        "C1CN2CCC=C2C1",
        3,
        "pyrrolizidine alkaloid core; well documented plant hepatotoxin class",
    ),
    (
        "thiourea",
        "[NX3][CX3](=[SX1])[NX3]",
        2,
        "oxidative desulfuration to reactive sulfenic/sulfinic species",
    ),
    (
        "nitrosamine",
        "[NX3][NX2]=[OX1]",
        3,
        "nitrosamine; metabolic activation to alkylating diazonium",
    ),
    (
        "azo",
        "[#6][NX2]=[NX2][#6]",
        1,
        "azoreduction to aromatic amines",
    ),
]

_COMPILED = [
    (name, Chem.MolFromSmarts(sma), sev, note) for name, sma, sev, note in ALERTS
]
_BAD = [n for n, p, _, _ in _COMPILED if p is None]
if _BAD:
    raise RuntimeError(f"invalid SMARTS in ALERTS: {_BAD}")


def _build_catalog(cats) -> FilterCatalog:
    params = FilterCatalogParams()
    for cat in cats:
        params.AddCatalog(cat)
    return FilterCatalog(params)


# Kept as two catalogs rather than one, so a hit can be attributed to
# assay-interference (PAINS) versus general medicinal-chemistry undesirability
# (Brenk/NIH). RDKit's FilterCatalogEntry does not expose the originating
# catalog, so separating them at construction is the reliable way to tell.
_PAINS_CATALOG = _build_catalog(
    (
        FilterCatalogParams.FilterCatalogs.PAINS_A,
        FilterCatalogParams.FilterCatalogs.PAINS_B,
        FilterCatalogParams.FilterCatalogs.PAINS_C,
    )
)
_BRENK_NIH_CATALOG = _build_catalog(
    (
        FilterCatalogParams.FilterCatalogs.BRENK,
        FilterCatalogParams.FilterCatalogs.NIH,
    )
)


@dataclass
class SafetyProfile:
    alerts: list[str] = field(default_factory=list)
    severities: list[int] = field(default_factory=list)
    pains: list[str] = field(default_factory=list)
    brenk_nih: list[str] = field(default_factory=list)

    @property
    def n_alerts(self) -> int:
        return len(self.alerts)

    @property
    def max_severity(self) -> int:
        return max(self.severities, default=0)

    @property
    def severity_sum(self) -> int:
        return sum(self.severities)

    @property
    def is_pains(self) -> bool:
        return bool(self.pains)


def profile(mol: Chem.Mol) -> SafetyProfile:
    """Return reactive-metabolite alerts and interference flags for a molecule."""
    p = SafetyProfile()
    if mol is None:
        return p

    for name, patt, sev, _note in _COMPILED:
        if mol.HasSubstructMatch(patt):
            p.alerts.append(name)
            p.severities.append(sev)

    p.pains = [e.GetDescription() for e in _PAINS_CATALOG.GetMatches(mol)]
    p.brenk_nih = [e.GetDescription() for e in _BRENK_NIH_CATALOG.GetMatches(mol)]

    return p


def physchem(mol: Chem.Mol) -> dict:
    """Descriptors relevant to a *central* analgesic: CNS exposure matters here.

    The whole premise is central analgesia via the AM404 axis, so a compound
    that cannot plausibly reach the CNS is a weaker candidate for this
    hypothesis regardless of how well it scores on target. We compute a
    CNS-MPO-flavoured check rather than only Lipinski.
    """
    if mol is None:
        return {}
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    arom = Descriptors.NumAromaticRings(mol)

    lipinski_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    # Pragmatic CNS-penetration heuristic (Wager-style, simplified):
    cns_like = (
        mw <= 400 and 1 <= logp <= 4 and tpsa <= 90 and hbd <= 3 and rotb <= 8
    )

    return {
        "mw": round(mw, 2),
        "clogp": round(logp, 3),
        "tpsa": round(tpsa, 2),
        "hbd": hbd,
        "hba": hba,
        "rotb": rotb,
        "aromatic_rings": arom,
        "fsp3": round(Descriptors.FractionCSP3(mol), 3),
        "qed": round(Descriptors.qed(mol), 4),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "lipinski_violations": lipinski_viol,
        "cns_like": int(cns_like),
    }


if __name__ == "__main__":
    probes = {
        "paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "NAPQI": "CC(=O)N=C1C=CC(=O)C=C1",
        "curcumin": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
        "6-gingerol": "CCCCC[C@@H](O)CC(=O)CCc1ccc(O)c(OC)c1",
        "kaempferol": "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12",
        "genistein": "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12",
        "capsaicin": "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
        "CBD": "CCCCCc1cc(O)c([C@@H]2C=C(C)CC[C@H]2C(=C)C)c(O)c1",
        "ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    }
    print(f"{'compound':14s} {'nAlert':>6s} {'maxSev':>6s} {'PAINS':>5s} "
          f"{'cns':>3s}  alerts")
    for name, smi in probes.items():
        m = Chem.MolFromSmiles(smi)
        p = profile(m)
        pc = physchem(m)
        print(f"{name:14s} {p.n_alerts:6d} {p.max_severity:6d} "
              f"{len(p.pains):5d} {pc['cns_like']:3d}  {','.join(p.alerts)}")
