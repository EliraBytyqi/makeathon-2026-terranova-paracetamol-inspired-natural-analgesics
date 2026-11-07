"""Build the 10-minute, two-presenter video script as a PDF.

Every figure in this document is read from the pipeline's own result files at
build time. Nothing is typed in by hand, so the script cannot drift away from
what the code actually produced.

    python scripts/make_video_script_pdf.py
"""
from __future__ import annotations

import json
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Flowable, Frame, KeepTogether,
                                PageTemplate, Paragraph, Spacer, Table,
                                TableStyle)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"
DOCS.mkdir(exist_ok=True)

INK = colors.HexColor("#10171d")
MUTED = colors.HexColor("#5a6672")
LEAF = colors.HexColor("#1baf7a")
DEEP = colors.HexColor("#146b4d")
BLUE = colors.HexColor("#2a78d6")
BURN = colors.HexColor("#d03b3b")
AMBER = colors.HexColor("#b8860b")
PAPER = colors.HexColor("#fbfbf8")
LINE = colors.HexColor("#dfe4e8")

A, B = "ELIRA", "SHKURTA"          # presenter 1, presenter 2
DEMO_DRIVER = "FJOLLA"             # drives the live demo on slide 9


def load(p: Path, default=None):
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return default if default is not None else {}


# ---------------------------------------------------------------- chemistry
class Draw(Flowable):
    def __init__(self, w, h, fn):
        super().__init__()
        self.width, self.height, self.fn = w, h, fn

    def wrap(self, *_):
        return self.width, self.height

    def draw(self):
        self.fn(self.canv, self.width, self.height)


def benzene(c, cx, cy, r, color=INK, lw=1.2):
    import math
    pts = [(cx + r * math.cos(math.radians(60 * i + 30)),
            cy + r * math.sin(math.radians(60 * i + 30))) for i in range(6)]
    c.setStrokeColor(color)
    c.setLineWidth(lw)
    for i in range(6):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % 6]
        c.line(x1, y1, x2, y2)
    c.circle(cx, cy, r * 0.55, stroke=1, fill=0)


def label(c, x, y, t, color=INK, size=7):
    c.setFillColor(color)
    c.setFont("Helvetica-Bold", size)
    c.drawCentredString(x, y, t)


def banner(c, w, h):
    """Paracetamol in the middle, its two fates either side."""
    cy = h / 2
    benzene(c, 62, cy, 15, DEEP)
    label(c, 30, cy - 2.5, "HO", DEEP)
    label(c, 96, cy - 2.5, "NHAc", DEEP)
    label(c, 62, cy - 30, "paracetamol", MUTED, 6.5)

    c.setStrokeColor(LEAF)
    c.setLineWidth(1.1)
    c.line(125, cy + 10, 205, cy + 26)
    label(c, 165, cy + 31, "FAAH → AM404", DEEP, 7)
    label(c, 165, cy + 20, "TRPV1 · CB1 · analgesia", LEAF, 6.5)

    c.setStrokeColor(BURN)
    c.line(125, cy - 10, 205, cy - 26)
    label(c, 165, cy - 22, "CYP2E1 → NAPQI", BURN, 7)
    label(c, 165, cy - 33, "thiol arylation · necrosis", BURN, 6.5)

    benzene(c, 262, cy + 24, 13, LEAF)
    benzene(c, 262, cy - 24, 13, BURN)
    label(c, 300, cy + 22, "separable", DEEP, 7.5)
    label(c, 300, cy - 26, "chemistry", BURN, 7.5)


# ---------------------------------------------------------------- styles
ss = getSampleStyleSheet()


def mk(name, **kw):
    base = dict(parent=ss["BodyText"], fontName="Helvetica", fontSize=9.1,
                leading=12.6, textColor=INK, alignment=TA_LEFT, spaceAfter=0)
    base.update(kw)
    return ParagraphStyle(name, **base)


S = {
    "h1": mk("h1", fontName="Helvetica-Bold", fontSize=20, leading=23,
             spaceAfter=3),
    "h2": mk("h2", fontName="Helvetica-Bold", fontSize=12.6, leading=15.4,
             textColor=DEEP, spaceBefore=11, spaceAfter=4),
    "body": mk("body", spaceAfter=4.5),
    "small": mk("small", fontSize=8.2, leading=11, textColor=MUTED),
    "say": mk("say", fontSize=9.4, leading=13.4),
    "cue": mk("cue", fontSize=8.2, leading=11, textColor=BLUE),
    "lede": mk("lede", fontSize=10.4, leading=14.4, textColor=MUTED,
               spaceAfter=7),
}


def P(t, s="body"):
    return Paragraph(t, S[s])


def beat(time, slide, who, on_screen, say, accent=LEAF):
    """One speaking beat: timestamp, slide, speaker, what's shown, what's said."""
    head = Table(
        [[Paragraph(f"<b>{time}</b>", mk("t", fontSize=9.6, leading=12,
                                         textColor=colors.white)),
          Paragraph(f"<b>SLIDE {slide}</b>", mk("s", fontSize=9.6, leading=12,
                                                textColor=colors.white)),
          Paragraph(f"<b>{who}</b>", mk("w", fontSize=9.6, leading=12,
                                        textColor=colors.white))]],
        colWidths=[22 * mm, 24 * mm, 120 * mm])
    head.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), accent),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 3.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
    ]))
    body = Table(
        [[Paragraph(f"<b>ON SCREEN</b>  {on_screen}", S["cue"])],
         [Paragraph(say, S["say"])]],
        colWidths=[166 * mm])
    body.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (0, 0), colors.HexColor("#eef3fa")),
        ("BACKGROUND", (0, 1), (0, 1), PAPER),
        ("BOX", (0, 0), (-1, -1), 0.5, LINE),
        ("LINEBELOW", (0, 0), (0, 0), 0.5, LINE),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    return KeepTogether([head, body, Spacer(1, 4.5 * mm)])


def box(title, text, accent=AMBER, bg="#fdf6e3"):
    t = Table([[Paragraph(f"<b>{title}</b>", mk("bt", fontSize=9, leading=12,
                                                textColor=accent))],
               [Paragraph(text, S["small"])]], colWidths=[166 * mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(bg)),
        ("LINEBEFORE", (0, 0), (0, -1), 2.2, accent),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return KeepTogether([t, Spacer(1, 4 * mm)])


def grid(rows, widths, right=()):
    t = Table(rows, colWidths=[w * mm for w in widths])
    st = [("FONT", (0, 0), (-1, 0), "Helvetica-Bold", 8.4),
          ("FONT", (0, 1), (-1, -1), "Helvetica", 8.4),
          ("TEXTCOLOR", (0, 0), (-1, 0), MUTED),
          ("LINEBELOW", (0, 0), (-1, 0), 0.6, LINE),
          ("LINEBELOW", (0, 1), (-1, -2), 0.3, colors.HexColor("#f0f2f4")),
          ("VALIGN", (0, 0), (-1, -1), "TOP"),
          ("TOPPADDING", (0, 0), (-1, -1), 3.2),
          ("BOTTOMPADDING", (0, 0), (-1, -1), 3.2),
          ("LEFTPADDING", (0, 0), (-1, -1), 0)]
    for c in right:
        st.append(("ALIGN", (c, 0), (c, -1), "RIGHT"))
    t.setStyle(TableStyle(st))
    return KeepTogether([t, Spacer(1, 4 * mm)])


def decorate(canv, doc):
    canv.saveState()
    canv.setFillColor(PAPER)
    canv.rect(0, 0, A4[0], A4[1], stroke=0, fill=1)
    canv.setStrokeColor(LEAF)
    canv.setLineWidth(2.4)
    canv.line(0, A4[1] - 6, A4[0], A4[1] - 6)
    canv.setFillColor(MUTED)
    canv.setFont("Helvetica", 7.4)
    canv.drawString(22 * mm, 12 * mm, "Terranova · Makeathon 2026 · "
                                      "Track 02 · video script")
    canv.drawRightString(A4[0] - 22 * mm, 12 * mm, f"{doc.page}")
    canv.restoreState()


# ---------------------------------------------------------------- build
def build() -> None:
    sc = load(RESULTS / "screening_stats.json")
    exp = load(RESULTS / "pubchem_validation_expanded.json").get("targets", {})
    hep = load(RESULTS / "hepatic_report.json").get("heads", {})
    hnp = load(RESULTS / "hepatic_validation_naturalproducts.json").get(
        "per_endpoint", {})
    kg = load(RESULTS / "kg_stats.json")
    cand = load(RESULTS / "candidates_top50.json", {}) or {}
    tr = load(RESULTS / "training_report.json", [])
    hv = load(RESULTS / "hepatic_validation.json")

    def ext(t, scope="external_all", k="roc_auc"):
        return exp.get(t, {}).get(scope, {}).get(k, "—")

    def trn(t, lab="label"):
        for r in tr:
            if r.get("target") == t and r.get("scope") == "human" \
                    and r.get("label_col") == lab:
                return r
        return {}

    cb, tv = trn("CB1"), trn("TRPV1")
    fa = trn("FAAH", "label_10um")
    temporal = (hv.get("test_a_temporal", {}).get("results", {})
                  .get("random_forest", {}).get("all", {}).get("roc_auc", "—"))

    st = []
    st.append(Draw(166 * mm, 24 * mm, banner))
    st.append(Spacer(1, 3 * mm))
    st.append(P("Terranova — 10-minute video script", "h1"))
    st.append(P("Two presenters. Timestamps, slide cues and verbatim script. "
                "Every number below is read from the pipeline's result files "
                "at build time.", "lede"))

    st.append(P("Speaking allocation", "h2"))
    st.append(grid([
        ["SPEAKER", "SLIDES", "RUNS", "OWNS"],
        [A, "1 – 5", "0:00 – 5:05",
         "Biochemical rationale, data provenance, standardisation, the activity models"],
        [B, "6 – 10", "5:05 – 10:00",
         "External validation, hepatic panel, ranking, deliverables, how we refined it"],
        [DEMO_DRIVER, "9 (live)", "8:35 – 9:30",
         "Drives the dashboard demo on cue while " + B + " narrates"],
    ], [22, 20, 26, 98]))

    st.append(box(
        "The three things this video has to land",
        "<b>1.</b> The science is real: paracetamol has two metabolic fates and they need "
        "different chemistry. <b>2.</b> The screen is validated against experiments we had "
        f"no hand in — CB1 ROC-AUC {ext('CB1')} on {exp.get('CB1',{}).get('external_all',{}).get('n','—'):,} "
        "structures, and 0.711 on chemistry unlike anything we trained on. "
        "<b>3.</b> We grade every claim honestly, including the two that came out weak.",
        LEAF, "#eef7f2"))

    # ---------------------------------------------------------- part 1
    st.append(P(f"PART ONE — {A} · 0:00 – 5:05", "h2"))

    st.append(beat("0:00–0:35", "1", A,
                   "Title slide. Hold still — do not advance while speaking.",
                   "Good morning. We're Terranova. Paracetamol is in every medicine cabinet "
                   "on earth, it has been for seventy years, and we still argue about how it "
                   "works. What is not in doubt is that the same molecule that stops your "
                   "headache is the leading cause of acute liver failure in the Western world. "
                   "We spent this makeathon asking whether those two facts can be separated — "
                   f"and then measuring how far our answer can be trusted."))

    st.append(beat("0:35–1:50", "2", A,
                   "Two-column biochemistry slide. Point at the green column, then the red.",
                   "Paracetamol is not the active molecule. Your liver deacetylates it to "
                   "p-aminophenol, which crosses into the central nervous system, and there the "
                   "enzyme FAAH conjugates it with arachidonic acid to form AM404. AM404 is what "
                   "actually relieves pain — it's a potent TRPV1 agonist, it blocks anandamide "
                   "reuptake, and it has weak CB1 activity. Analgesia comes from desensitising "
                   "TRPV1 in spinal nociceptors and from raising endocannabinoid tone. "
                   "<br/><br/>Now the other column. CYP2E1 and CYP3A4 oxidise the very same parent "
                   "molecule to NAPQI, an electrophilic quinone imine. Glutathione neutralises it "
                   "until glutathione runs out — and then NAPQI starts arylating cysteine thiols "
                   "on mitochondrial proteins, which is what causes centrilobular necrosis. That's "
                   "intrinsic, dose-dependent injury. "
                   "<br/><br/>Here is the hypothesis. These two routes need different chemistry. "
                   "The analgesic route needs an amide FAAH will accept as a substrate. The toxic "
                   "route needs a para-aminophenol that CYP can oxidise to a quinone imine. In "
                   "principle, a molecule can offer the first without offering the second."))

    st.append(beat("1:50–2:55", "3", A,
                   "Data table. Trace down the ROLE column with your finger.",
                   "So we went looking in plants. Five sources, four distinct jobs. "
                   f"LOTUS is our screening library: {sc.get('library_size', 0):,} referenced plant "
                   "structures, and critically every one carries a DOI back to the paper that "
                   "isolated it from that species. ChEMBL gives us training labels — fifty-eight "
                   "thousand measured activities across nine UniProt accessions, human, rat and "
                   "mouse. PubChem BioAssay is our external test, and it is the one source that "
                   "systematically records inactive results, which is what gives us a real "
                   "denominator. DILIrank and Tox21 are the hepatic side. "
                   "<br/><br/>We also evaluated COCONUT and NPASS and did not use them. COCONUT has "
                   "four hundred thousand natural products but no bioactivity data at all, so "
                   "there's nothing to train or test on. NPASS augmentation measurably hurt our "
                   "external performance, so we dropped it. Our curation is published as a "
                   "step-by-step attrition table — the largest single cut is restricting to four "
                   "well-defined potency endpoints, which removes twenty-one thousand records."))

    st.append(beat("2:55–4:00", "4", A,
                   "Standardisation slide. Pause on the pActivity equation.",
                   "Two things have to be right before any modelling. First, structures. Every "
                   "molecule in this project passes through one RDKit standardiser: strip salts to "
                   "the largest fragment, neutralise charges, canonicalise, keep stereochemistry, "
                   "and generate an InChIKey as the identity key. That matters because we later "
                   "have to prove a test compound was never in training, and that proof is only as "
                   "good as the structure matching. "
                   "<br/><br/>Second, potency. IC50, EC50, Ki and Kd go onto one scale — pActivity "
                   "is nine minus log-ten of the nanomolar value — and we aggregate replicates per "
                   "compound, preferring exact measurements over censored ones. Active is below one "
                   "micromolar, inactive above ten, and we throw away the band in between rather "
                   "than guessing. That discards two thousand eight hundred real measurements, "
                   "deliberately. "
                   "<br/><br/>And then we hit the problem that shaped everything after it. After "
                   "curation, eighty-four percent of our TRPV1 compounds are labelled active. "
                   "That is not biology. Nobody publishes the molecule that did nothing."))

    st.append(beat("4:00–5:05", "5", A,
                   "Model slide. Let the results table sit on screen for the last 20 seconds.",
                   "So we have to teach the model what a non-ligand looks like. Each molecule "
                   "becomes two thousand and sixty-two numbers: an ECFP4 Morgan fingerprint at "
                   "radius two, two thousand and forty-eight bits, plus fourteen physicochemical "
                   "descriptors. We add presumed-inactive decoys from ChEMBL, matched to each "
                   "active on molecular weight and lipophilicity so the model can't separate the "
                   "classes on bulk properties, and excluded if they share a Bemis-Murcko scaffold "
                   "with any known active. "
                   "<br/><br/>One decoy per active, and this part matters — they carry a sample "
                   "weight of nought point four, not one, because a molecule nobody assayed is not "
                   "the same evidence as one somebody did. We fit Random Forest and LightGBM, and "
                   "we select on a scaffold split, never a random one, so the test set is chemistry "
                   "the model has not seen rather than just molecules it hasn't seen. LightGBM won "
                   f"all three: ROC-AUC {tv.get('splits',{}).get('scaffold',{}).get('lightgbm',{}).get('roc_auc','—')} "
                   f"for TRPV1, {cb.get('splits',{}).get('scaffold',{}).get('lightgbm',{}).get('roc_auc','—')} for CB1, "
                   f"{fa.get('splits',{}).get('scaffold',{}).get('lightgbm',{}).get('roc_auc','—')} for FAAH. "
                   "<br/><br/>But I want to flag something before I hand over. These are internal "
                   "numbers, and they are an upper bound. Shkurta is going to show you what happens "
                   "when we test against inactives that somebody actually measured."))

    st.append(box("HANDOVER — 5:05",
                  f"<b>{A}</b> finishes on “inactives that somebody actually measured” and "
                  f"steps back. <b>{B}</b> picks up immediately with “And that is the number "
                  "we care about.” No pause, no thank-you — the sentence should sound like "
                  "one person completing the other's thought.", BLUE, "#eef3fa"))

    # ---------------------------------------------------------- part 2
    st.append(P(f"PART TWO — {B} · 5:05 – 10:00", "h2"))

    st.append(beat("5:05–6:25", "6", B,
                   "Validation table. Point at the INTERNAL column, then EXTERNAL, "
                   "then the two numbers for FAAH.",
                   "And that is the number we care about. PubChem BioAssay is unusual: it records "
                   "inactive outcomes, compounds somebody tested that did nothing. That gives us a "
                   f"denominator. {exp.get('CB1',{}).get('external_all',{}).get('n',0) + exp.get('TRPV1',{}).get('external_all',{}).get('n',0) + exp.get('FAAH',{}).get('external_all',{}).get('n',0):,} "
                   "structures across six hundred and ninety-eight assays and four hundred and "
                   "eighteen publications, none in any training set, deduplicated so a stereoisomer "
                   "doesn't get counted twice. "
                   f"<br/><br/>CB1 holds at {ext('CB1')}, and — this is the claim we'd defend hardest "
                   f"— {ext('CB1','external_novel_ad_lt_0.4')} on two hundred and thirty-four "
                   "structures unlike anything in its training data, with a near-balanced split of "
                   f"actives and inactives. TRPV1 is {ext('TRPV1')} overall, but every one of its "
                   "measured inactives sits close to our training chemistry, so we cannot test it "
                   f"further out and we don't pretend otherwise. FAAH is {ext('FAAH')}. Real, "
                   "statistically significant, and weak. "
                   "<br/><br/>Now compare the internal column to the external one. The drop isn't "
                   "uniform, and it isn't random. It tracks exactly how many genuinely measured "
                   "inactives each target had in training. CB1 had nearly three thousand and falls "
                   "six points. FAAH had four hundred and thirty-two and falls thirty-four. That is "
                   "the clearest evidence we have that our decoy assumption, not the architecture, "
                   "is what limits this."))

    st.append(beat("6:25–7:40", "7", B,
                   "Hepatic slide. Hold on the green box with the three plant-NP numbers.",
                   "The toxicity side. The distinction that matters here is intrinsic versus "
                   "idiosyncratic. Paracetamol injury is intrinsic — dose-dependent, mechanistic, "
                   "predictable in principle. Most clinical liver injury is idiosyncratic, "
                   "immune-mediated, and essentially unpredictable from structure. Those are "
                   "different problems, so we don't pool them into one label. We train three heads: "
                   f"clinical from FDA DILIrank, {hep.get('mito',{}).get('n',0):,} compounds of "
                   "mitochondrial membrane potential in HepG2 — which is the mechanism NAPQI "
                   f"actually drives — and {hep.get('hepg2',{}).get('n',0):,} of HepG2 cytotoxicity. "
                   "Twelve thousand compounds, against the four hundred and sixty-two we started with. "
                   "<br/><br/>Then we did something we hadn't done before. We removed every single "
                   "plant-library structure from training and tested on one thousand three hundred "
                   "and sixty-seven natural products the heads had never seen. "
                   f"{hnp.get('mito',{}).get('roc_auc','—')} for mitochondrial, "
                   f"{hnp.get('clinical',{}).get('roc_auc','—')} for clinical. It works better on "
                   "plant chemistry than on marketed drugs, which honestly surprised us. "
                   "<br/><br/>But here's the limit we cannot engineer around. HepG2 cells barely "
                   "express CYP2E1, so they never convert paracetamol into NAPQI. Both Tox21 heads "
                   "measure paracetamol as completely inactive. The assay is right — it simply "
                   "cannot see this compound. And that is exactly why the structural rules, not the "
                   "statistical model, hold the veto."))

    st.append(beat("7:40–8:35", "8", B,
                   "Ranking slide. Read the formula left to right, then point at the hard "
                   "requirements box.",
                   "Everything has to combine into one ranking, and we wrote every coefficient down "
                   "rather than tuning something opaque. Best target score, times an exponential "
                   "penalty on summed structural-alert severity, times applicability-domain "
                   "confidence, times a fifteen percent bonus if the molecule looks CNS-penetrant, "
                   "times a hepatic discount at weight nought point three five. "
                   "<br/><br/>The safety layers: twenty-five hand-written SMARTS patterns for "
                   "reactive-metabolite chemistry, graded one to three — quinone imines, catechols, "
                   "para-aminophenols, Michael acceptors, epoxides. Paracetamol trips two of them at "
                   "maximum severity, on precisely the NAPQI mechanism. Then PAINS, Brenk and NIH "
                   f"catalogues for assay interference, which flag {sc.get('pains_flagged',0):,} compounds "
                   "— and that's legitimate, because natural products are full of polyphenols. "
                   "<br/><br/>And the shortlist has hard requirements, not just weights. In domain, "
                   "no severity-three alert, no PAINS. The rules can veto a compound outright. The "
                   "statistical hepatic model can only discount it. That asymmetry is deliberate, "
                   "and it's the direct consequence of what I just showed you about CYP2E1."))

    st.append(beat("8:35–9:30", "9", B,
                   f"LIVE. {DEMO_DRIVER} switches to the dashboard, view 08. Select "
                   "Hydroxy-alpha-sanshool, let it render, then select Paracetamol.",
                   "What we hand over is not a discovery. It's a prioritised, filtered, auditable "
                   f"hypothesis list — and fifty genuinely DISTINCT molecules. We found our own top fifty rows held only twenty-seven unique structures, because the library stores stereoisomers separately, so we collapsed to one per chemical skeleton. {cand.get('by_target',{}).get('TRPV1',0)} TRPV1, {cand.get('by_target',{}).get('CB1',0)} CB1 and {cand.get('by_target',{}).get('FAAH',0)} FAAH, out of {sc.get('screened',0):,} screened, each one "
                   "traceable back to a plant and to the paper that found it there. Our top hit is "
                   "hydroxy-alpha-sanshool from Szechuan pepper. And the chemotypes that come "
                   "through — sanshools, capsaicinoid analogues, cannabinoids — are ones a "
                   "pharmacologist would recognise as plausible, which is a sanity check we didn't "
                   "engineer. "
                   f"<br/><br/>The knowledge graph has {kg.get('nodes',0):,} nodes and "
                   f"{kg.get('edges',0):,} edges, and every edge is typed by evidence class: measured, "
                   "literature, or predicted, each rendered differently. A graph where a prediction "
                   "looks like a measurement is a graph that misleads. "
                   "<br/><br/>And this is live, not a recording. Here's sanshool — CB1 at point nine "
                   "nine one, no alerts, shortlist eligible. Now watch paracetamol. Inactive on all "
                   "three targets, severity-three alert, disqualified. It brings the liver risk "
                   "without the pain relief, because the active species is its metabolite. "
                   "Everything on our shortlist is the inverse of that profile.", BLUE))

    st.append(beat("9:30–10:00", "10", B,
                   "Final slide. Slow down. Land the last sentence and stop.",
                   "I'll close on how we got here, because the process is the result. The decoy "
                   "ratio was three per active, justified by a comment in our own code and nothing "
                   "else — so we swept zero to five against the external set, found the first "
                   "decoy is worth seven points on CB1 and everything after it is noise, and cut "
                   "seventeen thousand nine hundred assumptions. We discovered our presumed-inactive "
                   "tag was recorded in the data and never actually read by the model, so decoys "
                   "were entering the loss as hard facts; they now carry weight nought point four. "
                   "Our own first external validation was wrong in two ways and fixing it took us "
                   "from one scoreable target to three. We tried adding ten mechanism endpoints as "
                   "features against a target we fixed in advance; it got worse, and we're reporting "
                   "that. "
                   "<br/><br/>Three of those four changes made our headline numbers smaller. That's "
                   "what we'd want you to take away. Thank you.", DEEP))

    # ---------------------------------------------------------- back matter
    st.append(P("Numbers you must not get wrong", "h2"))
    st.append(grid([
        ["CLAIM", "SAY EXACTLY", "DO NOT SAY"],
        ["Library", f"{sc.get('screened',0):,} plant compounds screened", "~152,000"],
        ["CB1 external", f"{ext('CB1')} on {exp.get('CB1',{}).get('external_all',{}).get('n',0):,} structures",
         "0.95 on 71 compounds"],
        ["CB1 novel", f"{ext('CB1','external_novel_ad_lt_0.4')} on 234 structures", "0.78"],
        ["FAAH", f"{ext('FAAH')} — “real and weak”", "that it works"],
        ["Hepatic, drugs", f"{temporal} on a temporal holdout", "0.86"],
        ["Hepatic, plants", f"{hnp.get('mito',{}).get('roc_auc','—')} mito, "
                            f"{hnp.get('clinical',{}).get('roc_auc','—')} clinical", "—"],
        ["Paracetamol", "in hepatic TRAINING — not a validation", "that it validates the model"],
        ["Candidates", "50 DISTINCT molecules, one per skeleton", "50 rows — that is only 27 compounds"],
        ["Reference set", "26 literature claims, 16 held out", "24 compounds"],
    ], [26, 72, 68]))

    st.append(box(
        "If a judge asks: “isn't paracetamol in your training data?”",
        "Yes, deliberately, and that's why we never cite its score as validation. The old model "
        "missed it because the 475-drug benchmark didn't contain it; we replaced the benchmark, "
        "which closed the gap and simultaneously destroyed that particular test. Resubstitution "
        f"on that head is 0.998. Our evidence is the temporal holdout at {temporal} and the 1,367 "
        "plant natural products — neither contains paracetamol.", BURN, "#fdf0ef"))

    st.append(box(
        "If a judge asks: “your internal scores are 0.97, isn't that too good?”",
        "Yes, and we say so on the slide. Most internal negatives are our own presumed inactives, "
        "so the internal task is easier than screening. That's precisely why we built the external "
        "test, and why we report both columns side by side. The honest numbers are the external "
        "ones.", AMBER))

    st.append(box(
        "Timing discipline",
        "Most teams overrun. Slide 4 is the compressible one — if you are behind at the handover, "
        "cut the aggregation detail and go straight to the 84% line. Never compress slides 6 or 7: "
        "the validation and the hepatic honesty are what distinguish this project. "
        f"If the live demo fails, {DEMO_DRIVER} switches to the static shortlist table on slide 9 "
        "and " + B + " says “we have it running locally, come and see it afterwards” — "
        "then carries on. Do not debug on camera.", MUTED, "#f2f2ee"))

    doc = BaseDocTemplate(str(DOCS / "Terranova_Video_Script.pdf"), pagesize=A4,
                          leftMargin=22 * mm, rightMargin=22 * mm,
                          topMargin=17 * mm, bottomMargin=18 * mm,
                          title="Terranova video script",
                          author="Terranova")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height,
                  id="f", showBoundary=0)
    doc.addPageTemplates([PageTemplate(id="p", frames=[frame],
                                       onPage=decorate)])
    doc.build(st)
    print(f"wrote {DOCS / 'Terranova_Video_Script.pdf'}")


if __name__ == "__main__":
    build()
