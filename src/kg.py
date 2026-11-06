"""Build the plant -> compound -> target -> pathway -> pain knowledge graph.

Every edge carries an `evidence` property with one of three values:

    measured    a ChEMBL bioactivity measurement, with its potency
    literature  a LOTUS structure-organism pair, with its reference DOI
    predicted   a model output, with its probability and applicability domain

That distinction is the entire point of building a graph rather than a
spreadsheet. It lets a candidate be presented as a path - "this compound is
reported in this plant by this paper, the plant is used traditionally for
pain, and our model predicts the compound engages FAAH" - with the
measured and predicted links visibly different. A graph where a prediction
looks the same as a measurement is a graph that misleads.

Outputs:
    results/knowledge_graph.graphml   full graph, loadable in Cytoscape/Gephi
    results/knowledge_graph.html      interactive pyvis view of the top subgraph
    results/kg_stats.json
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import networkx as nx
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
REF = ROOT / "data" / "reference"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

DISEASE = "Disease::Pain"
PATHWAYS = {
    "TRPV1": "Pathway::Vanilloid nociceptor signalling",
    "CB1": "Pathway::Endocannabinoid signalling",
    "FAAH": "Pathway::Endocannabinoid degradation",
}

NODE_COLORS = {
    "Plant": "#4C9F70",
    "Compound": "#3C78D8",
    "Target": "#E07B39",
    "Pathway": "#9B59B6",
    "Disease": "#C0392B",
    "SafetyFlag": "#7F8C8D",
    "ChemicalClass": "#1ABC9C",
}

TOP_N_FOR_VIEW = 30


def add_node(G, nid, ntype, **attrs):
    if nid not in G:
        G.add_node(nid, node_type=ntype, label=nid.split("::", 1)[-1],
                   color=NODE_COLORS.get(ntype, "#999999"), **attrs)
    else:
        G.nodes[nid].update({k: v for k, v in attrs.items() if v is not None})


def build(shortlist_n: int = 200) -> nx.MultiDiGraph:
    G = nx.MultiDiGraph()

    # ---- disease and pathway backbone ----------------------------------
    add_node(G, DISEASE, "Disease")
    for target, pathway in PATHWAYS.items():
        tnode = f"Target::{target}"
        add_node(G, tnode, "Target")
        add_node(G, pathway, "Pathway")
        G.add_edge(tnode, pathway, relation="INVOLVED_IN", evidence="literature",
                   source="curated from AM404 mechanism literature")
        G.add_edge(pathway, DISEASE, relation="MODULATES", evidence="literature",
                   source="curated from AM404 mechanism literature")

    # ---- measured activity from ChEMBL ---------------------------------
    agg_path = PROC / "bioactivity_aggregated.csv"
    n_measured = 0
    if agg_path.exists():
        agg = pd.read_csv(agg_path)
        # only the confirmed actives become graph edges; inactives are training
        # signal, not knowledge worth asserting
        act = agg[(agg["label"] == 1) & (agg["organism_group"] == "human")]
        for _, r in act.iterrows():
            cnode = f"Compound::{r['inchikey']}"
            add_node(G, cnode, "Compound", smiles=r["std_smiles"],
                     inchikey=r["inchikey"], origin="chembl")
            G.add_edge(
                cnode, f"Target::{r['target_name']}",
                relation="ACTS_ON", evidence="measured",
                pactivity=float(r["pactivity"]),
                n_measurements=int(r["n_measurements"]),
                assay_types=str(r["assay_types"]),
                action_type=str(r.get("action_type", "")),
                source="ChEMBL",
            )
            n_measured += 1

    # ---- predicted activity + plant provenance -------------------------
    screen_path = RESULTS / "screening_all.csv"
    n_predicted = 0
    n_contains = 0
    n_flags = 0
    if screen_path.exists():
        screen = pd.read_csv(screen_path, low_memory=False)
        top = screen.head(shortlist_n)
        targets = [c[2:] for c in screen.columns if c.startswith("p_")]

        for _, r in top.iterrows():
            cnode = f"Compound::{r['std_inchikey']}"
            add_node(
                G, cnode, "Compound",
                smiles=r["std_smiles"],
                inchikey=r["std_inchikey"],
                compound_name=str(r.get("compound_name", "")),
                origin="plant",
                composite_score=float(r.get("composite_score", 0) or 0),
                tier=str(r.get("tier", "")),
            )
            if pd.notna(r.get("compound_name")):
                G.nodes[cnode]["label"] = str(r["compound_name"])[:40]

            # predicted target edges
            for t in targets:
                p = r.get(f"p_{t}")
                if pd.isna(p) or p < 0.5:
                    continue
                G.add_edge(
                    cnode, f"Target::{t}",
                    relation="PREDICTED_ACTIVE", evidence="predicted",
                    probability=round(float(p), 4),
                    applicability_domain=round(float(r.get(f"ad_{t}", 0) or 0), 4),
                    source="terranova model",
                )
                n_predicted += 1

            # plant provenance, with the DOI that supports it
            species = str(r.get("source_species", "") or "")
            dois = str(r.get("reference_dois", "") or "")
            first_doi = dois.split("|")[0] if dois else ""
            for sp in [s for s in species.split("|") if s and not s.startswith("...")][:6]:
                pnode = f"Plant::{sp}"
                add_node(G, pnode, "Plant",
                         family=str(r.get("source_families", "")).split("|")[0])
                G.add_edge(pnode, cnode, relation="CONTAINS",
                           evidence="literature", reference_doi=first_doi,
                           source="LOTUS")
                n_contains += 1

            # chemical class
            klass = r.get("npclass_class")
            if pd.notna(klass) and str(klass).strip():
                knode = f"ChemicalClass::{klass}"
                add_node(G, knode, "ChemicalClass")
                G.add_edge(cnode, knode, relation="MEMBER_OF",
                           evidence="literature", source="NPClassifier")

            # safety flags
            alerts = str(r.get("reactive_alerts", "") or "")
            for a in [x for x in alerts.split("|") if x]:
                snode = f"SafetyFlag::{a}"
                add_node(G, snode, "SafetyFlag")
                G.add_edge(cnode, snode, relation="HAS_FLAG",
                           evidence="predicted",
                           source="reactive-metabolite SMARTS alert set")
                n_flags += 1

    # ---- traditional-use edges from the ethnopharmacology seed ----------
    seed = pd.read_csv(REF / "pain_plants_seed.csv")
    n_trad = 0
    for _, r in seed.iterrows():
        pnode = f"Plant::{r['species']}"
        if pnode not in G:
            continue
        G.nodes[pnode]["common_name"] = str(r["common_name"])
        G.add_edge(pnode, DISEASE, relation="TRADITIONALLY_USED_FOR",
                   evidence="literature", use=str(r["traditional_pain_use"]),
                   source="curated ethnopharmacology seed list")
        n_trad += 1

    G.graph["n_measured_edges"] = n_measured
    G.graph["n_predicted_edges"] = n_predicted
    G.graph["n_contains_edges"] = n_contains
    G.graph["n_flag_edges"] = n_flags
    G.graph["n_traditional_edges"] = n_trad
    return G


def export_html(G: nx.MultiDiGraph, path: Path, top_n: int = TOP_N_FOR_VIEW):
    """Interactive view of the highest-scoring compounds and their context."""
    from pyvis.network import Network

    compounds = [
        (n, d.get("composite_score", 0))
        for n, d in G.nodes(data=True)
        if d.get("node_type") == "Compound" and d.get("origin") == "plant"
    ]
    compounds.sort(key=lambda x: -x[1])
    keep = {n for n, _ in compounds[:top_n]}
    for n in list(keep):
        keep.update(G.predecessors(n))
        keep.update(G.successors(n))
    keep.add(DISEASE)
    keep.update(PATHWAYS.values())

    sub = G.subgraph(keep)
    net = Network(height="820px", width="100%", directed=True,
                  bgcolor="#111317", font_color="#e8e8e8")
    net.barnes_hut(gravity=-9000, spring_length=190, spring_strength=0.012)

    for n, d in sub.nodes(data=True):
        ntype = d.get("node_type", "")
        size = {"Disease": 46, "Pathway": 34, "Target": 40,
                "Compound": 22, "Plant": 24, "SafetyFlag": 16,
                "ChemicalClass": 18}.get(ntype, 16)
        title = f"{ntype}: {d.get('label', n)}"
        if ntype == "Compound":
            title += (f"\\nscore={d.get('composite_score', 0):.3f}"
                      f"\\ntier={d.get('tier','')}"
                      f"\\n{d.get('smiles','')[:60]}")
        net.add_node(n, label=d.get("label", n)[:30], color=d.get("color"),
                     size=size, title=title, shape="dot")

    for u, v, d in sub.edges(data=True):
        ev = d.get("evidence", "")
        color = {"measured": "#2ECC71", "predicted": "#F1C40F",
                 "literature": "#5DADE2"}.get(ev, "#777777")
        dashes = ev == "predicted"
        title = f"{d.get('relation','')} [{ev}]"
        if "probability" in d:
            title += f" p={d['probability']}"
        if "pactivity" in d:
            title += f" pAct={d['pactivity']:.2f}"
        if d.get("reference_doi"):
            title += f" doi:{d['reference_doi']}"
        net.add_edge(u, v, color=color, title=title, dashes=dashes, width=2)

    net.set_options("""
    {"physics": {"stabilization": {"iterations": 250}},
     "interaction": {"hover": true, "tooltipDelay": 120}}
    """)
    net.write_html(str(path), notebook=False)
    _localize_assets(path)


def _localize_assets(path: Path):
    """Make the exported graph work without a network connection.

    pyvis writes its `lib/` folder relative to the *current working directory*,
    not relative to the output file, and links vis-network from a CDN. Run from
    src/, that leaves results/knowledge_graph.html pointing at a lib/ directory
    that is not there and at cdnjs for the library itself - so the graph breaks
    on any machine that is offline or behind a captive portal. Since the real
    library is already vendored, point at it and drop the dead references.
    """
    import shutil

    out_dir = path.parent
    for candidate in (Path.cwd() / "lib", ROOT / "src" / "lib", ROOT / "lib"):
        if candidate.is_dir() and candidate.resolve() != (out_dir / "lib").resolve():
            shutil.copytree(candidate, out_dir / "lib", dirs_exist_ok=True)
            shutil.rmtree(candidate, ignore_errors=True)
            break

    html = path.read_text(encoding="utf-8")
    replacements = {
        "https://cdnjs.cloudflare.com/ajax/libs/vis-network/9.1.2/dist/dist/vis-network.min.css":
            "lib/vis-9.1.2/vis-network.css",
        "https://cdnjs.cloudflare.com/ajax/libs/vis-network/9.1.2/dist/vis-network.min.js":
            "lib/vis-9.1.2/vis-network.min.js",
    }
    for remote, local in replacements.items():
        html = html.replace(remote, local)

    # Bootstrap is only used for pyvis's optional filter widget, which this
    # export does not enable. Its tags span several lines with the attributes
    # on their own rows, so they need a pattern rather than a literal match.
    # (The ../node_modules references pyvis also emits sit inside an HTML
    # comment and never generate a request, so they are left alone.)
    html = re.sub(
        r"\s*<link[^>]*?cdn\.jsdelivr\.net/npm/bootstrap[^>]*?/?>", "", html,
        flags=re.DOTALL)
    html = re.sub(
        r"\s*<script[^>]*?cdn\.jsdelivr\.net/npm/bootstrap[^>]*?>\s*</script>", "",
        html, flags=re.DOTALL)

    path.write_text(html, encoding="utf-8")

    if re.search(r'(?<!<!--)\s(?:src|href)="https?://', html):
        print("  WARNING: graph still references a remote asset")

    missing = [
        ref for ref in ("lib/vis-9.1.2/vis-network.min.js", "lib/vis-9.1.2/vis-network.css")
        if ref in html and not (out_dir / ref).exists()
    ]
    if missing:
        print(f"  WARNING: graph references missing local assets: {missing}")
    else:
        print("  assets localised - graph renders with no network access")


def main():
    G = build()
    print(f"nodes: {G.number_of_nodes():,}   edges: {G.number_of_edges():,}")

    by_type = {}
    for _, d in G.nodes(data=True):
        by_type[d.get("node_type", "?")] = by_type.get(d.get("node_type", "?"), 0) + 1
    print("\nnodes by type:")
    for k, v in sorted(by_type.items(), key=lambda x: -x[1]):
        print(f"   {k:16s} {v:,}")

    by_rel = {}
    by_ev = {}
    for _, _, d in G.edges(data=True):
        by_rel[d.get("relation", "?")] = by_rel.get(d.get("relation", "?"), 0) + 1
        by_ev[d.get("evidence", "?")] = by_ev.get(d.get("evidence", "?"), 0) + 1
    print("\nedges by relation:")
    for k, v in sorted(by_rel.items(), key=lambda x: -x[1]):
        print(f"   {k:24s} {v:,}")
    print("\nedges by evidence class:")
    for k, v in sorted(by_ev.items(), key=lambda x: -x[1]):
        print(f"   {k:16s} {v:,}")

    # GraphML cannot hold list/None attributes; coerce everything to primitives
    H = nx.MultiDiGraph()
    H.add_nodes_from(
        (n, {k: ("" if v is None else v if isinstance(v, (str, int, float, bool)) else str(v))
             for k, v in d.items()})
        for n, d in G.nodes(data=True)
    )
    for u, v, d in G.edges(data=True):
        H.add_edge(u, v, **{
            k: ("" if val is None else val if isinstance(val, (str, int, float, bool)) else str(val))
            for k, val in d.items()
        })
    gpath = RESULTS / "knowledge_graph.graphml"
    nx.write_graphml(H, gpath)
    print(f"\nwrote {gpath}")

    hpath = RESULTS / "knowledge_graph.html"
    try:
        export_html(G, hpath)
        print(f"wrote {hpath}")
    except Exception as exc:
        print(f"pyvis export skipped: {exc}")

    stats = {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "nodes_by_type": by_type,
        "edges_by_relation": by_rel,
        "edges_by_evidence": by_ev,
    }
    (RESULTS / "kg_stats.json").write_text(json.dumps(stats, indent=2),
                                           encoding="utf-8")


if __name__ == "__main__":
    main()
