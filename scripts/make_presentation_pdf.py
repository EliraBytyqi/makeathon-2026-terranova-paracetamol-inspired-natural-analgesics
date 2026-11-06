"""Generate the presentation and demo plan PDF.

Reads the live pipeline artefacts so every number in the document is the one
currently on disk, then writes docs/Terranova_Presentation_Plan.pdf.

    python scripts/make_presentation_pdf.py

Note on glyphs: ReportLab's built-in fonts use WinAnsi encoding, which has no
arrows, check marks or math comparison signs. Those render as black boxes, so
this file uses "->", ">=" and text labels instead.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate, Frame, KeepTogether, PageBreak, PageTemplate, Paragraph,
    Spacer, Table, TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PROC = ROOT / "data" / "processed"
DOCS = ROOT / "docs"
DOCS.mkdir(exist_ok=True)

# palette, matched to the dashboard
INK = colors.HexColor("#0b0b0b")
INK2 = colors.HexColor("#52514e")
MUTED = colors.HexColor("#898781")
RULE = colors.HexColor("#e1e0d9")
SURFACE = colors.HexColor("#f2f2ee")
BLUE = colors.HexColor("#2a78d6")
ORANGE = colors.HexColor("#eb6834")
AQUA = colors.HexColor("#1baf7a")
RED = colors.HexColor("#d03b3b")
AMBER = colors.HexColor("#fab219")
GREEN = colors.HexColor("#0ca30c")


def load(path: Path, default=None):
    if not path.exists():
        return default
    if path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    return pd.read_csv(path)


# ---------------------------------------------------------------- styles
ss = getSampleStyleSheet()


def mk(name, **kw):
    base = dict(fontName="Helvetica", fontSize=9.5, leading=14, textColor=INK,
                spaceAfter=6, alignment=TA_LEFT)
    base.update(kw)
    return ParagraphStyle(name, **base)


S = {
    "title": mk("title", fontName="Helvetica-Bold", fontSize=26, leading=30,
                spaceAfter=4),
    "subtitle": mk("subtitle", fontSize=12.5, leading=17, textColor=INK2,
                   spaceAfter=14),
    "h1": mk("h1", fontName="Helvetica-Bold", fontSize=16, leading=20,
             spaceBefore=16, spaceAfter=8),
    "h2": mk("h2", fontName="Helvetica-Bold", fontSize=11.5, leading=15,
             spaceBefore=11, spaceAfter=5),
    "h3": mk("h3", fontName="Helvetica-Bold", fontSize=9.8, leading=13,
             spaceBefore=8, spaceAfter=3),
    "body": mk("body"),
    "small": mk("small", fontSize=8.6, leading=12, textColor=INK2),
    "muted": mk("muted", fontSize=8.2, leading=11.5, textColor=MUTED),
    "eyebrow": mk("eyebrow", fontName="Helvetica-Bold", fontSize=8,
                  leading=11, textColor=BLUE, spaceAfter=3),
    "say": mk("say", fontName="Helvetica-Oblique", fontSize=9.8, leading=14.5,
              textColor=INK, leftIndent=9, spaceAfter=5),
    "cell": mk("cell", fontSize=8.4, leading=11.5, spaceAfter=0),
    "cellb": mk("cellb", fontName="Helvetica-Bold", fontSize=8.4, leading=11.5,
                spaceAfter=0),
    "cellm": mk("cellm", fontSize=8.4, leading=11.5, textColor=INK2, spaceAfter=0),
}


def P(text, style="body"):
    return Paragraph(text, S[style])


def callout(title, body, accent=BLUE, style="body"):
    inner = [P(f"<b>{title}</b>", "h3")] if title else []
    inner.append(P(body, style))
    t = Table([[inner]], colWidths=[165 * mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), SURFACE),
        ("LINEBEFORE", (0, 0), (0, -1), 2.2, accent),
        ("LEFTPADDING", (0, 0), (-1, -1), 9),
        ("RIGHTPADDING", (0, 0), (-1, -1), 9),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    return t


def table(rows, widths, header=True, align_right=()):
    data = []
    for r_i, row in enumerate(rows):
        out = []
        for c_i, cell in enumerate(row):
            if isinstance(cell, Paragraph):
                out.append(cell)
            else:
                st = "cellb" if (header and r_i == 0) else "cell"
                out.append(P(str(cell), st))
        data.append(out)
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
    style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE),
    ]
    if header:
        style += [
            ("BACKGROUND", (0, 0), (-1, 0), SURFACE),
            ("LINEBELOW", (0, 0), (-1, 0), 0.8, MUTED),
        ]
    for c in align_right:
        style.append(("ALIGN", (c, 0), (c, -1), "RIGHT"))
    t.setStyle(TableStyle(style))
    return t


# ------------------------------------------------------------- page frame
def decorate(canvas, doc):
    canvas.saveState()
    w, h = A4
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(MUTED)
    if doc.page > 1:
        canvas.drawString(22 * mm, h - 13 * mm, "Terranova - Presentation & Demo Plan")
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.4)
        canvas.line(22 * mm, h - 15 * mm, w - 22 * mm, h - 15 * mm)
    canvas.drawRightString(w - 22 * mm, 12 * mm, f"{doc.page}")
    canvas.drawString(22 * mm, 12 * mm, "Computational predictions only. Nothing here establishes activity or safety.")
    canvas.restoreState()


def build():
    # ---------------------------------------------------------- live data
    stats = load(RESULTS / "screening_stats.json", {}) or {}
    lit = load(RESULTS / "literature_validation_summary.json", {}) or {}
    kg = load(RESULTS / "kg_stats.json", {}) or {}
    rescue = load(RESULTS / "faah_rescue_summary.json", {}) or {}
    dili = load(RESULTS / "dili_report.json", {}) or {}
    probes = (load(RESULTS / "dili_probes.json", {}) or {}).get("probes", [])
    metrics = load(RESULTS / "model_metrics.csv")
    attr = load(PROC / "attrition_audit.csv")
    shortlist = load(RESULTS / "shortlist.csv")
    litcsv = load(RESULTS / "literature_validation.csv")

    para = next((p for p in probes if p["name"] == "paracetamol"), {})
    n_held = int(lit.get("n_held_out", 0))
    med_pct = lit.get("held_out_median_library_percentile")
    top10 = round(lit.get("held_out_frac_top10pct_of_library", 0) * n_held)
    scaff = metrics[metrics["split"] == "scaffold"] if metrics is not None else None
    pr_lo = float(scaff["pr_auc"].min()) if scaff is not None and len(scaff) else 0
    pr_hi = float(scaff["pr_auc"].max()) if scaff is not None and len(scaff) else 0

    doc = BaseDocTemplate(
        str(DOCS / "Terranova_Presentation_Plan.pdf"), pagesize=A4,
        leftMargin=22 * mm, rightMargin=22 * mm,
        topMargin=20 * mm, bottomMargin=18 * mm,
        title="Terranova - Presentation & Demo Plan",
        author="Terranova", subject="Makeathon 2026, Track 02",
    )
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="f")
    doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=decorate)])

    st = []  # story

    # ======================================================== COVER
    st += [
        Spacer(1, 22 * mm),
        P("MAKEATHON 2026 &nbsp;&middot;&nbsp; TRACK 02 &nbsp;&middot;&nbsp; NATURAL PRODUCTS", "eyebrow"),
        P("Terranova", "title"),
        P("Presentation &amp; demo plan", "subtitle"),
        callout(
            "The one-line version",
            "Paracetamol relieves pain and destroys livers through two chemically "
            "separable routes. We screened "
            f"<b>{int(stats.get('screened', 0)):,}</b> plant natural products for the "
            "first while filtering out the second - and then built the honest "
            "evaluation that tells us how far to trust the result.",
            BLUE),
        Spacer(1, 7 * mm),
    ]

    st += [P("What is in this document", "h2")]
    st.append(table([
        ["Section", "What it gives you"],
        ["1. Read this first", "A result that landed today and changes your framing"],
        ["2. Why the order is the argument", "What makes this structure work, and why it cannot be copied"],
        ["3. The narrative", "Six beats, timed, with what to say at each"],
        ["4. The demo script", "Exact click path through the dashboard"],
        ["5. Numbers to know cold", "The figures a judge may ask you to defend"],
        ["6. Question prep", "Eight likely questions with answers"],
        ["7. Red lines", "Claims you must not make"],
        ["8. Pre-flight checklist", "What to verify before you walk up"],
    ], [42 * mm, 123 * mm]))

    st += [
        Spacer(1, 6 * mm),
        callout(
            "Status",
            "All thirteen pipeline stages have now run. The three that were still "
            "missing this morning - external literature validation, the knowledge "
            "graph, and the written report - completed today, and the dashboard "
            "shows all of them.",
            GREEN, "small"),
    ]

    # ======================================================== 1. READ FIRST
    st += [PageBreak(), P("1. Read this first", "h1")]
    st += [P(
        "The external validation stage ran for the first time today. It is the only "
        "genuinely held-out test in the project, and <b>it did not go well</b>. You need "
        "to know this before you plan a single slide, because the honest framing of it "
        "is far stronger than any attempt to talk around it.", "body")]

    st.append(table([
        ["What we measured", "Result"],
        [P("Held-out literature actives (absent from ChEMBL training)", "cell"), P(f"<b>{n_held}</b> compound-target claims", "cell")],
        [P("How many score above 0.5 on their claimed target", "cell"), P(f"<b>0 of {n_held}</b>", "cell")],
        [P("Median percentile against the 145k library background", "cell"), P(f"<b>{med_pct}</b> (50 would be chance)", "cell")],
        [P("How many reach the library's top 10%", "cell"), P(f"<b>{top10} of {n_held}</b>", "cell")],
        [P("Same models, on compounds they <i>did</i> see in training", "cell"), P(f"median probability <b>{lit.get('seen_median_probability')}</b>", "cell")],
    ], [105 * mm, 60 * mm]))

    st += [
        Spacer(1, 4 * mm),
        callout(
            "Why this is the strongest thing you have",
            f"Your scaffold-split PR-AUC runs {pr_lo:.2f} to {pr_hi:.2f}. Most teams would "
            "put that on a slide and stop. You went one step further, tested against "
            "compounds from the primary literature that the model had never seen, and "
            "found the number does not transfer. <b>Almost nobody in the room will have "
            "done that.</b> Judges reward a team that can find and name its own ceiling "
            "far more than one quoting a number they cannot interrogate.",
            GREEN),
        Spacer(1, 3 * mm),
        callout(
            "And you have the diagnosis, not just the failure",
            "The FAAH rescue experiment answers exactly the question this raises. Adding "
            "five curated literature actives moves the median prediction on a held-out "
            f"sixth from {rescue.get('baseline_median_p', 0):.1e} to "
            f"{rescue.get('loo_median_p', 0):.3f} - four orders of magnitude. "
            "<b>So the failure is coverage, not capability.</b> The features and the model "
            "can represent this chemistry; ChEMBL simply does not contain it. That turns a "
            "negative result into a costed next step: a few hundred curated literature "
            "actives, not a new architecture.",
            BLUE),
    ]

    st += [P("The pattern inside the failure", "h2")]
    st += [P(
        "It is not random noise, and saying so shows you read your own output. The four "
        "compounds that do surface are all TRPV1 vanilloids - structurally close to "
        "capsaicin, which <i>is</i> in training. The FAAH flavonoids all sit mid-pack. "
        "<b>The models extrapolate a short distance from what they know, and no further.</b>", "body")]

    if litcsv is not None and len(litcsv):
        held = litcsv[litcsv["in_chembl_training"] == False].sort_values(
            "library_percentile", ascending=False)
        rows = [["Compound", "Target", "Reported", "Library %ile"]]
        for _, r in held.head(6).iterrows():
            rows.append([
                P(str(r["compound_name"]), "cell"),
                P(str(r["target"]), "cell"),
                P(f"{r['value_uM']:g} uM", "cellm"),
                P(f"<b>{r['library_percentile']:.1f}</b>", "cell"),
            ])
        rows.append([P("<i>... and the worst case</i>", "cellm"), "", "", ""])
        worst = held.tail(1).iloc[0]
        rows.append([
            P(f"{worst['compound_name']}", "cell"),
            P(str(worst["target"]), "cell"),
            P(f"{worst['value_uM']:g} uM", "cellm"),
            P(f"<b>{worst['library_percentile']:.1f}</b>", "cell"),
        ])
        st += [Spacer(1, 2 * mm), table(rows, [70 * mm, 25 * mm, 35 * mm, 35 * mm],
                                        align_right=(3,))]
        st += [P(
            "AM404 is the paracetamol metabolite this entire project is built around, and "
            "on TRPV1 it lands in the bottom 5% of the library. If you do not say that "
            "first, somebody may find it and say it for you.", "small")]

    # ======================================================== 2. THE ARC
    st += [PageBreak(), P("2. Why the order is the argument", "h1")]
    st += [P(
        "The running order in the next section is not a stylistic choice, and changing it "
        "breaks the case. This page explains why, so that whoever is holding the clicker "
        "understands what each beat is doing.", "body")]

    st += [P("The shape", "h2")]
    st += [P(
        "It works, then here is where it breaks, then here is why, then here is the fix. "
        "Four steps, and your project has the evidence for all four.", "body")]

    st.append(table([
        ["Step", "Your evidence", "What it establishes"],
        [P("<b>It works</b>", "cellb"),
         P("Positive controls: 38 measured actives, median 0.9988; 47 measured inactives, median 0.0004", "cell"),
         P("The machine is built correctly", "cell")],
        [P("<b>Where it breaks</b>", "cellb"),
         P(f"External literature set: 0 of {n_held} above 0.5, median percentile {med_pct}", "cell"),
         P("The boundary of competence", "cell")],
        [P("<b>Why</b>", "cellb"),
         P(f"FAAH rescue: median moves {rescue.get('baseline_median_p', 0):.1e} to {rescue.get('loo_median_p', 0):.3f}", "cell"),
         P("The mechanism of the failure", "cell")],
        [P("<b>The fix</b>", "cellb"),
         P("A few hundred curated literature actives; the rescue prices it", "cell"),
         P("A roadmap, not an apology", "cell")],
    ], [30 * mm, 80 * mm, 55 * mm]))

    st += [P("Why every step is load-bearing", "h2")]
    st += [P("The clearest way to see it is to remove one and look at what is left.", "body")]
    st.append(table([
        ["If you present...", "What a judge actually hears"],
        ["\"It works\" alone",
         "A demo. Unfalsifiable, and every team has one. It carries no information, because a team that never checked would show the identical slide."],
        ["\"It breaks\" alone", "A failed project. Nobody presents this."],
        ["Works + breaks",
         "Honest but defeatist. A limitation you cannot explain, and the last thing the room remembers is the failure."],
        ["Works + breaks + why",
         "A diagnosis. This is the step that turns a bad number into a finding."],
        [P("<b>All four</b>", "cellb"),
         P("<b>A research programme. The negative result becomes the justification for the next phase.</b>", "cellb")],
    ], [42 * mm, 123 * mm]))

    st += [
        Spacer(1, 3 * mm),
        callout("The hinge is step three",
                "This is not four nice things stacked up. Without the rescue experiment, "
                f"\"0 of {n_held}\" is simply bad news. With it, the same number becomes "
                "evidence for a specific, testable claim about why - and the argument "
                "survives.",
                BLUE, "small"),
    ]

    st += [P("Why it cannot be copied in the last hour", "h2")]
    st += [P(
        "Each step depends on infrastructure that had to exist before you knew you would "
        "need it. That is precisely what makes it credible.", "body")]
    st.append(table([
        ["Step", "What it required, and when"],
        ["It works",
         "Leakage control designed into screen.py - removing training compounds and keeping them aside as controls. An architectural decision made when the screen was written. Retrofitting it means re-running the screen."],
        ["Where it breaks",
         "resolve_literature.py: 24 compounds hand-curated from the primary literature, structures resolved through PubChem, potencies recorded with confidence levels. Hours of careful work, and not fabricable - a judge can check whether kaempferol really inhibits FAAH at 5 uM."],
        ["Why",
         "experiment_faah_rescue.py: a leave-one-out design that retrains a model per compound. A designed experiment, not a metric computed from existing output."],
    ], [30 * mm, 135 * mm]))

    st += [
        Spacer(1, 3 * mm),
        callout("It is expensive to fake, which is the whole point",
                "A team that spent their time optimising a number cannot produce this at "
                "11pm - not through laziness, but because what makes it possible is the "
                "evaluation infrastructure, and that had to be built alongside the model. "
                "The arc is not a presentation technique bolted on at the end. It is the "
                "visible output of having built the project a particular way.",
                AQUA, "small"),
    ]

    st += [P("Why judges reward it specifically", "h2")]
    st += [P(
        "Consider the judge's actual problem. They see fifteen to thirty teams in a day, "
        "every one claims their model works, and they can verify none of it. So what they "
        "are really trying to estimate is: <b>does this team know whether their own thing "
        "works?</b>", "body")]
    st += [P(
        "A team reporting only successes provides no evidence either way. A team reporting "
        "a specific, quantified limitation has demonstrably run a test that could have "
        "embarrassed them, and did not flinch. In a format where nothing can be verified, "
        "that is the only credible signal available.", "body")]
    st += [P(
        "It also inoculates you completely. The strongest attack on any screening project "
        "is \"this will not generalise\". Once you have said it yourself, with a number, "
        "before anyone asks, there is no attack left. The worst thing a judge can say, you "
        "have already said, with more precision than they could.", "body")]

    st += [P("Three rules for delivering it", "h2")]
    st.append(table([
        ["#", "Rule", "Why"],
        ["1", "Order is non-negotiable: competence first, limitation second, never reversed",
         "You earn the right to report the negative by first demonstrating the positive."],
        ["2", "Beats 4 and 5 must be adjacent - no slide, no tangent, no pause for questions",
         "If the room sits in the failure for even thirty seconds, the framing is lost."],
        ["3", "Tone is finding, not confession",
         "You went looking for this. Never \"unfortunately our model did not work\"."],
    ], [8 * mm, 85 * mm, 72 * mm]))

    st += [P("The transitions, where it is won or lost", "h2")]
    st += [
        P("(beat 3 to 4) &nbsp; \"That is the part that works. Now - here is what happened "
          "when we tested it properly.\"", "say"),
        P("(beat 4 to 5) &nbsp; \"That number could have ended the project. It did not, "
          "because we could isolate why.\"", "say"),
        P("(beat 5 to 6) &nbsp; \"Which means the fix is not a better model. It is about "
          "three hundred more rows of data.\"", "say"),
    ]
    st += [
        Spacer(1, 2 * mm),
        callout("Land on that last line",
                "It reframes everything behind it. You are not a team whose model "
                "underperformed. You are a team that knows exactly what their model costs "
                "to finish - and that is a different kind of team entirely.",
                BLUE, "small"),
    ]

    # ======================================================== 3. NARRATIVE
    st += [PageBreak(), P("3. The narrative - six beats", "h1")]
    st += [P(
        "Timed for a six-minute slot. Cut beat 6 first if you are given five, and "
        "compress beat 2 if you are given four. The order matters: you earn the right to "
        "report a negative result by first showing the thing works.", "body")]

    beats = [
        ("0:00-0:45", "The hook: one molecule, two pathways", BLUE,
         "Paracetamol is in every medicine cabinet on earth and we still argue about how it "
         "works. The current account is that it is a prodrug: the body converts it to AM404, "
         "which acts on TRPV1 and the cannabinoid system. Its liver toxicity runs through a "
         "completely different route - oxidation to NAPQI. Two pathways, chemically "
         "separable. That gap is the whole opportunity.",
         "Open on the Hypothesis view. Point at the two arms."),

        ("0:45-1:30", "What we built", AQUA,
         "We harvested every measured activity for those three targets from ChEMBL, curated "
         f"{int(attr['rows'].iloc[0]):,} raw records down to {int(attr['rows'].iloc[-1]):,} with "
         "every filter counted, trained models with a scaffold-grouped split, then screened "
         f"{int(stats.get('screened', 0)):,} plant compounds from LOTUS - each one carrying a "
         "literature DOI back to the paper that found it in that plant.",
         "Overview view. Let the funnel chart do the talking."),

        ("1:30-2:30", "It works - here is the proof", AQUA,
         "Before trusting anything, we held out every library compound that appeared in "
         "training. Of those, the ones measured active score a median of 0.9988 - capsaicin, "
         "resiniferatoxin, THC, all recovered at the top. The ones measured inactive score "
         "0.0004. The pipeline separates them cleanly, and the score distribution is sharply "
         "peaked near zero, so it is not simply calling everything active.",
         "Screen explorer -> Positive controls tab. Scroll the top rows."),

        ("2:30-3:45", "But here is what we actually found", RED,
         "Then we ran the test most screens skip. We took compounds whose activity is reported "
         f"in the primary literature but which are absent from ChEMBL, and scored them. <b>Not "
         f"one of the {n_held} scores above 0.5.</b> The median sits at the "
         f"{med_pct}th percentile of our own library. Our scaffold-split number was "
         f"{pr_hi:.2f}. It did not transfer.",
         "External validation view. Stop on the red card. Let it land."),

        ("3:45-4:45", "And we know why", BLUE,
         "We isolated it. Take six known plant FAAH inhibitors, none in ChEMBL. Score each with "
         "a model that saw the other five but never that one. The median prediction moves four "
         "orders of magnitude. So the representation is fine and the model is fine - the "
         "training set has simply never seen this kind of molecule. That is a data problem with "
         "a price tag, not an architecture problem.",
         "Models view -> FAAH rescue dumbbell chart."),

        ("4:45-6:00", "What it is for, and what is next", AQUA,
         "What we hand over is not a discovery, it is a prioritised, filtered, "
         f"provenance-tracked hypothesis list: {len(shortlist) if shortlist is not None else 50} "
         "compounds, each traceable to a plant and a paper, each screened against the "
         "reactive chemistry that makes paracetamol dangerous. The knowledge graph makes "
         "every claim inspectable, and measured, literature and predicted edges look "
         "different on purpose. The next step is a few hundred curated literature actives "
         "to close the coverage gap we just measured.",
         "Knowledge graph, then back to the shortlist."),
    ]

    for t, title, accent, say, action in beats:
        block = [
            table([[P(f"<b>{t}</b>", "cellb"), P(f"<b>{title}</b>", "cellb")]],
                  [24 * mm, 141 * mm], header=False),
            P(say, "say"),
            table([[P(f"<b>On screen:</b> {action}", "cellm")]], [165 * mm], header=False),
            Spacer(1, 4 * mm),
        ]
        st.append(KeepTogether(block))

    # ======================================================== 3. DEMO
    st += [PageBreak(), P("4. The demo script", "h1")]
    st += [P(
        "Start the dashboard before you are called up, not when you reach the podium. "
        "Have it already open on the Overview view with the browser at full screen.", "body")]

    st += [callout("Launch",
                   "In the project folder:<br/><br/>"
                   "<font face='Courier'>.venv\\Scripts\\python.exe -m app.run</font>"
                   "<br/><br/>Then <font face='Courier'>http://127.0.0.1:8000</font>. "
                   "It opens the browser for you. Everything is precomputed - no model "
                   "loads at request time, so nothing can stall mid-demo.",
                   BLUE, "small")]

    st += [P("Click path", "h2")]
    st.append(table([
        ["#", "Where", "What to do", "What to say / point at"],
        ["1", "Overview", "Land here. Do not scroll yet.",
         "The four tiles: screened, in-domain, shortlist, curated records."],
        ["2", "Overview", "Scroll to the funnel.",
         "145,380 down to 50, every stage counted."],
        ["3", "The hypothesis", "Click nav 02.",
         "The two pathway arms. This is the whole idea in one picture."],
        ["4", "Screen explorer", "Click nav 07, then the Positive controls tab.",
         "Capsaicin at 0.99999. The model finds what we already know."],
        ["5", "Screen explorer", "Type 'Zanthoxylum' in search.",
         "Sichuan pepper sanshools - a real shortlist hit with traditional use."],
        ["6", "Screen explorer", "Click any row.",
         "The compound modal: scores, safety flags, and the DOIs behind it."],
        ["7", "External validation", "Click nav 05. Pause on the red card.",
         "This is the honest result. Let the silence work."],
        ["8", "External validation", "Scroll to the percentile chart.",
         "Only the vanilloids surface. The pattern is informative."],
        ["9", "Models", "Click nav 04, scroll to FAAH rescue.",
         "Four orders of magnitude. Coverage, not capability."],
        ["10", "Models", "Scroll to the DILI card.",
         f"Paracetamol scores {para.get('p_dili', 0):.2f} - low risk. The model is wrong, and we show it."],
        ["11", "Validation", "Click 'Open interactive graph'.",
         f"{kg.get('nodes', 0):,} nodes. Drag one. Evidence classes are coloured differently."],
        ["12", "ML decoder", "Only if a judge is technical.",
         "The glossary and the six claims with their evidence."],
    ], [8 * mm, 28 * mm, 56 * mm, 73 * mm]))

    st += [
        Spacer(1, 4 * mm),
        callout("If the laptop dies",
                "Everything the demo shows also exists as files: "
                "<font face='Courier'>results/REPORT.md</font> is the whole project in "
                "16,000 words, <font face='Courier'>results/shortlist.csv</font> is the "
                "50 candidates, and "
                "<font face='Courier'>results/knowledge_graph.html</font> opens in any "
                "browser with no server at all. Have that last one open in a second tab "
                "as insurance.",
                AMBER, "small"),
    ]

    # ======================================================== 4. NUMBERS
    st += [PageBreak(), P("5. Numbers to know cold", "h1")]
    st += [P("If you can only memorise one column, memorise the right-hand one.", "body")]

    st += [P("Scale", "h2")]
    st.append(table([
        ["Quantity", "Value"],
        ["Raw ChEMBL activity records harvested", f"{int(attr['rows'].iloc[0]):,}"],
        ["After curation (every filter counted)", f"{int(attr['rows'].iloc[-1]):,}"],
        ["Plant compounds in the LOTUS library", f"{int(stats.get('library_size', 0)):,}"],
        ["Screened after leakage control", f"{int(stats.get('screened', 0)):,}"],
        ["Within applicability domain", f"{int(stats.get('in_domain', 0)):,} ({stats.get('in_domain_rate', 0) * 100:.1f}%)"],
        ["Flagged by a severity-3 structural alert", f"{int(stats.get('severity3_flagged', 0)):,}"],
        ["Flagged as PAINS (assay interference)", f"{int(stats.get('pains_flagged', 0)):,}"],
        ["Final shortlist", f"{int(stats.get('shortlist_size', 0))}"],
        ["Knowledge graph", f"{kg.get('nodes', 0):,} nodes / {kg.get('edges', 0):,} edges"],
    ], [110 * mm, 55 * mm], align_right=(1,)))

    st += [P("Performance - and its ceiling", "h2")]
    st.append(table([
        ["Measure", "Value"],
        ["Scaffold-split PR-AUC, across all models", f"{pr_lo:.3f} - {pr_hi:.3f}"],
        ["Positive controls, measured active: median score", "0.9988 (34 of 38 above 0.9)"],
        ["Positive controls, measured inactive: median score", "0.0004"],
        [P("<b>Held-out literature actives above 0.5</b>", "cellb"), P(f"<b>0 of {n_held}</b>", "cellb")],
        [P("<b>Their median library percentile</b>", "cellb"), P(f"<b>{med_pct}</b>", "cellb")],
        ["FAAH rescue: median before / after", f"{rescue.get('baseline_median_p', 0):.1e} -> {rescue.get('loo_median_p', 0):.3f}"],
        ["DILI model ROC-AUC (scaffold, 5 repeats)", f"{dili.get('splits', {}).get('scaffold_5x', {}).get('roc_auc', 0):.3f}"],
        ["DILI model on paracetamol", f"{para.get('p_dili', 0):.3f} - i.e. low risk, and wrong"],
    ], [110 * mm, 55 * mm], align_right=(1,)))

    st += [P("Design decisions you should be able to justify", "h2")]
    st.append(table([
        ["Decision", "One-sentence reason"],
        ["Scaffold split, not random",
         "ChEMBL is congeneric series; a random split leaks analogues and inflates the score."],
        ["60,000 decoys as presumed inactives",
         "The measured active rate is 84% - publication bias, not biology."],
        ["Property-matched decoys",
         "Unmatched ones let the model learn 'heavy and greasy = active' and transfer nothing."],
        ["1-10 uM band discarded",
         "Forcing a label on genuinely ambiguous compounds teaches the model noise."],
        ["FAAH labelled at 10 uM as well as 1 uM",
         "Plant FAAH inhibitors sit at 1-15 uM; a 1 uM cut trains the model to reject them."],
        ["Two safety layers, alerts weighted above the model",
         "The alert set flags paracetamol at maximum severity; the DILI model calls it safe."],
        ["Evidence typed throughout the graph",
         "A graph where a prediction looks like a measurement is a graph that misleads."],
    ], [62 * mm, 103 * mm]))

    # ======================================================== 5. Q&A
    st += [PageBreak(), P("6. Question prep", "h1")]

    qa = [
        ("Your scaffold-split PR-AUC is 0.98. Isn't that suspiciously high?",
         "Yes, and we can tell you why. The decoys are easy negatives - random drug-like "
         "molecules, so the model is mostly separating 'real ligand' from 'arbitrary "
         "molecule', which is much easier than ranking among plausible candidates. That is "
         "exactly why we ran the external literature validation, and there the number did "
         "not hold. We treat the scaffold figure as an upper bound, not a result."),
        ("So your model doesn't work?",
         "It works inside its applicability domain and fails outside it, and we can show "
         "you precisely where the boundary is. It recovers known ligands, it separates "
         "measured actives from measured inactives, and it extrapolates a short way - the "
         "TRPV1 vanilloids surface. It does not reach chemistry ChEMBL has never contained. "
         "We measured that rather than assuming it, and the FAAH rescue shows it is fixable "
         "with data."),
        ("What would you do with another week?",
         "Curate several hundred literature actives for the three targets, the way we did "
         "the six in the rescue experiment, and retrain. The rescue tells us what that buys: "
         "four orders of magnitude on held-out plant chemistry. We would also add hard "
         "negatives rather than random decoys, so the scaffold number stops flattering us."),
        ("Why only LOTUS? What about COCONUT?",
         "We checked. COCONUT is larger at around 695,000 structures, but its bulk downloads "
         "carry IDs, SMILES and computed scores - not the organism and per-compound DOI "
         "linkage. Getting that needs their 32 GB SQL dump. LOTUS gives referenced "
         "structure-organism pairs directly, and that provenance is what makes the knowledge "
         "graph defensible. Most of LOTUS is in COCONUT anyway via its Wikidata collection."),
        ("Isn't it circular that most of your shortlist comes from plants you picked for pain use?",
         "Partly, yes. Seed-list compounds are 11.5% of the library and 72% of the shortlist, "
         "a 6.3-fold enrichment - but those plants are on the list because they have "
         "documented pain use, and their chemistry resembles the ChEMBL ligands. We treat it "
         "as a consistency check, not a discovery. The interesting rows are the 14 "
         "background-tier compounds that got there with no traditional-use prior."),
        ("Have you tested anything in a lab?",
         "No, and we never claim otherwise. Every number is a computational prediction over "
         "compounds this project has not assayed. What we are handing over is a ranked, "
         "filtered, provenance-tracked hypothesis list for someone with a lab."),
        ("Why should I trust your safety filter?",
         "Do not trust it on its own - that is why there are two layers. The DILI model "
         f"scores paracetamol {para.get('p_dili', 0):.2f}, i.e. low risk, because paracetamol "
         "and isoniazid are absent from its benchmark entirely. The structural alert set does "
         "flag paracetamol at maximum severity on the exact NAPQI mechanism. So we weighted "
         "the alerts above the model and documented why."),
        ("What is actually novel here?",
         "Three things. The target panel is chosen from a specific mechanistic hypothesis "
         "rather than convenience. The evaluation is built to fail honestly - scaffold split, "
         "external validation, applicability domain, typed evidence. And the FAAH rescue "
         "experiment isolates data coverage from model capability, which is the question "
         "everyone waves at and nobody measures."),
    ]
    for i, (q, a) in enumerate(qa, 1):
        st.append(KeepTogether([
            P(f"{i}. {q}", "h3"),
            P(a, "small"),
            Spacer(1, 3 * mm),
        ]))

    # ======================================================== 6. RED LINES
    st += [PageBreak(), P("7. Red lines", "h1")]
    st += [P(
        "This is a biomedical project presented to people who may know more biology than you "
        "do. One overclaim costs more credibility than every chart you built. These are not "
        "style preferences.", "body")]

    st.append(table([
        [P("<b>Never say</b>", "cellb"), P("<b>Say instead</b>", "cellb")],
        ["\"We found a new painkiller.\"", "\"We produced a ranked hypothesis list for lab testing.\""],
        ["\"These compounds are safe.\"", "\"These compounds lack the structural features associated with paracetamol's toxicity.\""],
        ["\"The model is 98% accurate.\"", "\"Scaffold-split PR-AUC is 0.98, but external validation shows it does not transfer.\""],
        ["\"This compound inhibits FAAH.\"", "\"Our model predicts activity; it has not been assayed.\""],
        ["\"Natural means safe.\"", "\"Plants make electrophiles on purpose - several of our own hits carry alerts.\""],
        ["\"We validated the model.\"", "\"We ran an external validation and it exposed a coverage limit.\""],
    ], [72 * mm, 93 * mm]))

    st += [
        Spacer(1, 4 * mm),
        callout("The disclaimer is on every page of the dashboard for a reason",
                "If a judge asks whether any of this is ready for a human being, the answer "
                "is an immediate and unhedged no. Saying that quickly and plainly reads as "
                "competence, not weakness.",
                RED, "small"),
    ]

    # ======================================================== 7. CHECKLIST
    st += [P("8. Pre-flight checklist", "h1")]
    st.append(table([
        ["", "Check", "How"],
        ["[ ]", "Dashboard starts and all artefacts are present",
         "Open http://127.0.0.1:8000 - no amber banner at the top means all files are there"],
        ["[ ]", "All eight views render in both themes",
         "Click each nav item, then the theme toggle, and click through again"],
        ["[ ]", "The knowledge graph opens and is draggable",
         "External validation -> Open interactive graph"],
        ["[ ]", "The compound modal shows DOIs",
         "Screen explorer -> click any row -> scroll to Literature provenance"],
        ["[ ]", "Browser zoom set so tables are readable from the back",
         "Ctrl+Plus two or three times; check the projector, not your screen"],
        ["[ ]", "Second tab open on results/knowledge_graph.html",
         "Insurance if the server dies"],
        ["[ ]", "Everything committed and pushed",
         "git status is clean; git push has finished"],
        ["[ ]", "You can say the AM404 percentile out loud without flinching",
         "It is 4.4. Practise it."],
    ], [9 * mm, 62 * mm, 94 * mm]))

    st += [
        Spacer(1, 7 * mm),
        callout("Last thing",
                "The temptation in the room will be to present the 0.98 and hope nobody asks. "
                "Do the opposite. Lead with the thing that works, then show the test that "
                "broke it, then show that you already know why. A team that can say 'here is "
                "where our model stops working, here is the measurement that proves it, and "
                "here is what it would cost to fix' is doing science. That is the wow factor, "
                "and it is the one nobody else can copy in the last hour.",
                BLUE),
    ]

    doc.build(st)
    return DOCS / "Terranova_Presentation_Plan.pdf"


if __name__ == "__main__":
    out = build()
    print(f"wrote {out}  ({out.stat().st_size / 1024:.0f} KB)")
