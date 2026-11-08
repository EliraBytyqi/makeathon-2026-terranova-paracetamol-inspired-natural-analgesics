/* Terranova dashboard — views, data binding and the screen explorer.
 *
 * Plain ES2019, no framework and no bundler: the demo has to run from a clean
 * checkout with one command, offline.
 */
(function () {
  "use strict";

  const C = window.Charts;
  const $ = (sel, root) => (root || document).querySelector(sel);
  const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));

  const TARGET_COLORS = { TRPV1: "--series-1", CB1: "--series-2", FAAH: "--series-3" };
  const state = {
    source: "candidates",
    page: 1,
    pageSize: 50,
    sort: "rank",
    desc: false,
    filters: { q: "", target: "", tier: "", in_domain: false, no_alerts: false, cns_only: false },
    cache: {},
    rendered: {},
  };

  /* ------------------------------------------------------------- utility */
  const esc = C.esc;
  const fmt = C.fmt;

  function num(n, digits) {
    if (n === null || n === undefined || n === "" || Number.isNaN(Number(n))) return "—";
    return Number(n).toFixed(digits === undefined ? 3 : digits);
  }
  function int(n) {
    if (n === null || n === undefined || Number.isNaN(Number(n))) return "—";
    return Math.round(Number(n)).toLocaleString();
  }
  function titleCase(s) {
    return String(s || "").replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase());
  }
  function targetBadge(t) {
    if (!t) return '<span class="dim">—</span>';
    return `<span class="badge badge-${esc(t)}"><span class="badge-dot"></span>${esc(t)}</span>`;
  }
  function compoundName(row) {
    return row.compound_name || row.std_inchikey || "(unnamed)";
  }

  async function api(path) {
    if (state.cache[path]) return state.cache[path];
    const res = await fetch(path);
    if (!res.ok) throw new Error(`${res.status} ${res.statusText} — ${path}`);
    const json = await res.json();
    state.cache[path] = json;
    return json;
  }

  function stat(label, value, sub, accent) {
    return `<div class="stat">
      <div class="stat-label">${esc(label)}</div>
      <div class="stat-value ${accent ? "stat-accent-" + accent : ""}">${value}</div>
      ${sub ? `<div class="stat-sub">${sub}</div>` : ""}
    </div>`;
  }

  function errorBox(host, err) {
    host.innerHTML = `<div class="empty">Could not load this section.<br><code>${esc(err.message)}</code></div>`;
  }

  /* --------------------------------------------------------------- theme */
  function initTheme() {
    let saved = null;
    try { saved = localStorage.getItem("tn-theme"); } catch (e) { /* private mode */ }
    if (saved) document.documentElement.setAttribute("data-theme", saved);
    paintThemeButton();
    $("#theme-toggle").addEventListener("click", () => {
      const cur = document.documentElement.getAttribute("data-theme");
      const sysDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
      const next = cur ? (cur === "dark" ? "light" : "dark") : sysDark ? "light" : "dark";
      document.documentElement.setAttribute("data-theme", next);
      try { localStorage.setItem("tn-theme", next); } catch (e) { /* ignore */ }
      paintThemeButton();
      // charts read CSS variables at draw time, so they must be redrawn
      state.rendered = {};
      renderView(currentView, true);
    });
  }
  function paintThemeButton() {
    const cur = document.documentElement.getAttribute("data-theme");
    const dark = cur ? cur === "dark" : window.matchMedia("(prefers-color-scheme: dark)").matches;
    $("#theme-icon").textContent = dark ? "☀" : "☾";
    $("#theme-label").textContent = dark ? "Light" : "Dark";
  }

  /* ----------------------------------------------------------- navigation */
  let currentView = "overview";
  function go(view) {
    currentView = view;
    $$(".nav-btn").forEach((b) => b.setAttribute("aria-current", String(b.dataset.view === view)));
    $$(".view").forEach((v) => v.classList.toggle("active", v.id === "view-" + view));
    if (location.hash.slice(1) !== view) history.replaceState(null, "", "#" + view);
    window.scrollTo({ top: 0, behavior: "instant" in window ? "instant" : "auto" });
    renderView(view);
  }

  function renderView(view, force) {
    if (state.rendered[view] && !force) return;
    state.rendered[view] = true;
    const fn = {
      overview: renderOverview,
      science: renderScience,
      data: renderData,
      models: renderModels,
      validation: renderValidation,
      safety: renderSafety,
      screen: renderScreen,
      demo: renderDemo,
      glossary: renderGlossaryPage,
    }[view];
    if (fn) fn();
  }

  /* ------------------------------------------------------------ 01 overview */
  async function renderOverview() {
    try {
      const o = await api("/api/overview");
      const s = o.screening || {};

      $("#lede-lib").textContent = int(s.library_size);

      $("#overview-stats").innerHTML = [
        stat("Compounds screened", int(s.screened),
             `${int(s.removed_as_training_overlap)} held out as controls`, 1),
        stat("Within applicability domain", int(s.in_domain),
             `${((s.in_domain_rate || 0) * 100).toFixed(1)}% of the library`, 3),
        stat("Shortlist", int(o.shortlist_size),
             "in-domain · no severity-3 · no PAINS", 2),
        stat("Curated activity records", int(o.curated_activity_records),
             `from ${int(o.raw_activity_records)} raw ChEMBL rows`),
      ].join("");

      C.funnel($("#pipeline-funnel"), {
        title: "Pipeline attrition",
        ariaLabel: "Compounds surviving each pipeline stage",
        labelWidth: 250,
        stages: [
          { label: "LOTUS plant library", value: s.library_size, note: "unique plant compounds" },
          { label: "After leakage control", value: s.screened, note: "training overlap removed" },
          { label: "Within applicability domain", value: s.in_domain, note: "max Tanimoto ≥ 0.25" },
          { label: "No PAINS flag", value: Math.max(0, (s.in_domain || 0) - (s.pains_flagged || 0)), note: "assay-interference filter" },
          { label: "Final shortlist", value: s.shortlist_size, note: "top 50 by composite score" },
        ],
      });

      const byTarget = o.shortlist_by_target || {};
      C.hbar($("#shortlist-by-target"), {
        title: "Shortlist by best target",
        labelWidth: 70,
        labelHeader: "Target",
        valueHeader: "Compounds",
        data: ["TRPV1", "CB1", "FAAH"]
          .filter((t) => byTarget[t])
          .map((t) => ({
            label: t,
            value: byTarget[t],
            color: C.cssVar(TARGET_COLORS[t]),
            tip: [["compounds", byTarget[t]],
                  ["share", ((byTarget[t] / o.shortlist_size) * 100).toFixed(0) + "%"]],
          })),
        legend: ["TRPV1", "CB1", "FAAH"].filter((t) => byTarget[t])
          .map((t) => ({ label: t, color: C.cssVar(TARGET_COLORS[t]) })),
      });

      const ctl = await api("/api/compounds/controls?page_size=8&sort=target_score&desc=true");
      $("#controls-preview").innerHTML = controlsTable(ctl.rows);

      const sl = await api("/api/compounds/shortlist?page_size=10");
      $("#overview-shortlist").innerHTML = shortlistTable(sl.rows);
      wireRowClicks($("#overview-shortlist"));
    } catch (e) { errorBox($("#overview-stats"), e); }
  }

  function controlsTable(rows) {
    if (!rows || !rows.length) return '<div class="empty">No positive controls found.</div>';
    return `<div class="table-wrap"><table>
      <caption>Known ligands the models never saw in the screen, ranked by target score.</caption>
      <thead><tr><th scope="col">Compound</th><th scope="col">Target</th><th scope="col" class="num">Score</th></tr></thead>
      <tbody>${rows.map((r) => `<tr>
        <th scope="row" style="font-weight:500">${esc(compoundName(r))}</th>
        <td>${targetBadge(r.best_target)}</td>
        <td class="num">${num(r.target_score, 4)}</td>
      </tr>`).join("")}</tbody></table></div>`;
  }

  function shortlistTable(rows) {
    if (!rows || !rows.length) return '<div class="empty">Shortlist not found. Run <code>python src/run_all.py</code>.</div>';
    return `<div class="table-wrap"><table>
      <thead><tr>
        <th scope="col" class="num">#</th><th scope="col">Compound</th><th scope="col">Target</th>
        <th scope="col" class="num">Target score</th><th scope="col" class="num">AD</th>
        <th scope="col">Source genus</th><th scope="col">Tier</th>
      </tr></thead>
      <tbody>${rows.map((r) => `<tr class="clickable" tabindex="0" data-key="${esc(r.std_inchikey)}">
        <td class="num dim">${r.rank}</td>
        <th scope="row" style="font-weight:500"><span class="truncate">${esc(compoundName(r))}</span></th>
        <td>${targetBadge(r.best_target)}</td>
        <td class="num">${num(r.target_score, 4)}</td>
        <td class="num">${num(r.ad_best, 2)}</td>
        <td><span class="truncate" style="max-width:170px">${esc((r.source_genera || "").split("|")[0] || "—")}</span></td>
        <td><span class="badge">${esc(titleCase(r.tier))}</span></td>
      </tr>`).join("")}</tbody></table></div>`;
  }

  /* ------------------------------------------------------------- 02 science */
  async function renderScience() {
    try {
      const d = await api("/api/datasets");
      const seeds = d.seed_plants || [];
      $("#seed-count").textContent = seeds.length;
      $("#seed-table").innerHTML = `<div class="table-wrap" style="max-height:440px;overflow-y:auto">
        <table><thead><tr>
          <th scope="col">Species</th><th scope="col">Common name</th>
          <th scope="col">Traditional pain use</th><th scope="col">Hypothesised target</th>
        </tr></thead><tbody>
        ${seeds.map((s) => `<tr>
          <th scope="row" style="font-weight:500"><em>${esc(s.species)}</em></th>
          <td>${esc(s.common_name)}</td>
          <td class="dim">${esc(s.traditional_pain_use)}</td>
          <td>${String(s.relevant_target_hypothesis || "")
                .split(";")
                .map((t) => (TARGET_COLORS[t] ? targetBadge(t) : `<span class="badge">${esc(t)}</span>`))
                .join(" ")}</td>
        </tr>`).join("")}
        </tbody></table></div>`;
    } catch (e) { errorBox($("#seed-table"), e); }
  }

  /* ---------------------------------------------------------------- 03 data */
  const SOURCES = [
    {
      name: "ChEMBL bioactivity",
      role: "Training labels",
      accent: 1,
      detail: "Measured potencies for TRPV1, CB1 and FAAH across human, rat and mouse. Harvested from the ChEMBL REST API after resolving nine UniProt accessions to ChEMBL target IDs.",
      meta: ["58,394 raw records", "9 targets", "IC50 / EC50 / Ki / Kd"],
    },
    {
      name: "ChEMBL decoy pool",
      role: "Presumed inactives",
      accent: 2,
      detail: "60,000 drug-like molecules (MW 150–700) sampled across ChEMBL's whole ID range — not the first 60k, which would over-represent pre-2000 medicinal chemistry. Property-matched to each target's actives.",
      meta: ["60,000 molecules", "DUD-E style matching", "scaffold-excluded"],
    },
    {
      name: "LOTUS natural products",
      role: "The screening library",
      accent: 3,
      detail: "Referenced structure–organism pairs: every compound–plant edge carries a literature DOI. Filtered to Archaeplastida (land plants); fungal, bacterial and marine sources excluded.",
      meta: ["145,380 plant compounds", "DOI-backed provenance", "two source tiers"],
    },
    {
      name: "FDA DILIrank 2.0",
      role: "Clinical hepatotoxicity labels",
      accent: 4,
      detail: "1,336 approved drugs ranked by the evidence in their own FDA labels. 982 carry an unambiguous concern class; the 354 marked ambiguous are dropped rather than guessed at. Replaces the 475-drug TDC benchmark, which did not contain paracetamol.",
      meta: ["506 modelable (most vs none)", "217 most-concern", "contains paracetamol"],
    },
    {
      name: "Tox21 \u2014 HepG2 assay panel",
      role: "Measured liver-cell mechanism",
      accent: 5,
      detail: "Mitochondrial membrane potential and CellTox Green cytotoxicity, both run in HepG2 liver cells as dose\u2013response screens. Measured outcomes rather than inferences from a drug label, and an order of magnitude more of them.",
      meta: ["11,619 measured compounds", "2 endpoints", "Active / Inactive calls"],
    },
  ];

  async function renderData() {
    try {
      const d = await api("/api/datasets");
      const o = await api("/api/overview");

      $("#source-cards").innerHTML = SOURCES.map((s) => `<div class="card" style="margin:0">
        <div class="stat-label stat-accent-${s.accent}">${esc(s.role)}</div>
        <h3 style="margin:4px 0 8px">${esc(s.name)}</h3>
        <p class="card-note">${esc(s.detail)}</p>
        <div class="pill-row" style="margin-top:11px">
          ${s.meta.map((m) => `<span class="badge">${esc(m)}</span>`).join("")}
        </div>
      </div>`).join("");

      const attr = d.attrition || [];
      C.hbar($("#attrition-chart"), {
        title: "ChEMBL curation attrition",
        ariaLabel: "Records surviving each curation step",
        ordinal: true,
        labelWidth: 260,
        labelHeader: "Curation step",
        valueHeader: "Records",
        data: attr.map((a, i) => ({
          label: a.step,
          value: a.rows,
          display: int(a.rows),
          emphasis: i === attr.length - 1,
          tip: [
            ["records", int(a.rows)],
            ["dropped here", i ? int(attr[i - 1].rows - a.rows) : "—"],
            ["of raw", ((a.rows / attr[0].rows) * 100).toFixed(1) + "%"],
          ],
        })),
      });

      const sum = (d.summary || []).filter((r) => r.scope === "human" || r.scope === "all");
      C.stacked($("#dataset-chart"), {
        title: "Labelled compounds per target and scope",
        ariaLabel: "Active, inactive and ambiguous compound counts per target",
        labelWidth: 120,
        labelHeader: "Target / scope",
        data: sum.map((r) => ({
          label: `${r.target} · ${r.scope}`,
          active: r.active_1uM,
          inactive: r.inactive_10uM,
          ambiguous: r.ambiguous,
        })),
        series: [
          { key: "active", label: "Active (≤ 1 µM)", colorVar: "--series-1" },
          { key: "inactive", label: "Inactive (≥ 10 µM)", colorVar: "--series-2" },
          { key: "ambiguous", label: "Ambiguous (held out)", colorVar: "--series-3" },
        ],
      });

      const human = (d.summary || []).filter((r) => r.scope === "human");
      $("#bias-stats").innerHTML = human.map((r) =>
        stat(`${r.target} measured active rate`,
             ((r.active_rate || 0) * 100).toFixed(1) + "%",
             `${int(r.active_1uM)} active vs ${int(r.inactive_10uM)} inactive`,
             r.target === "TRPV1" ? 1 : r.target === "CB1" ? 2 : 3)).join("");
    } catch (e) { errorBox($("#source-cards"), e); }
  }

  /* -------------------------------------------------------------- 04 models */
  let metricsSplit = "scaffold";

  async function renderModels() {
    try {
      const m = await api("/api/models");
      const metrics = m.metrics || [];

      // dumbbell: random vs scaffold PR-AUC, one row per model configuration
      const byConfig = {};
      metrics.forEach((r) => {
        const key = `${r.target}/${r.scope}/${r.label}/${r.model}`;
        byConfig[key] = byConfig[key] || { label: key, target: r.target };
        byConfig[key][r.split] = r.pr_auc;
      });
      const pairs = Object.values(byConfig)
        .filter((c) => c.random != null && c.scaffold != null)
        .map((c) => ({
          label: c.label.replace(/\//g, " · ").replace("random_forest", "RF").replace("lightgbm", "LGBM"),
          a: c.random,
          b: c.scaffold,
          extra: [["target", c.target]],
        }))
        .sort((x, y) => x.b - y.b);

      C.dumbbell($("#split-chart"), {
        title: "PR-AUC: random split vs scaffold split",
        ariaLabel: "Comparison of random and scaffold split PR-AUC for every model",
        labelWidth: 250,
        labelHeader: "Target · scope · label · model",
        aLabel: "Random split",
        bLabel: "Scaffold split",
        aColorVar: "--series-neutral",
        bColorVar: "--series-1",
        min: 0.8,
        max: 1.0,
        ticks: [0.8, 0.85, 0.9, 0.95, 1.0],
        data: pairs,
      });

      renderMetricsTable(metrics);
      $$("#view-models .seg [data-split]").forEach((b) => {
        b.addEventListener("click", () => {
          metricsSplit = b.dataset.split;
          $$("#view-models .seg [data-split]").forEach((x) =>
            x.setAttribute("aria-pressed", String(x === b)));
          renderMetricsTable(metrics);
        });
      });

      // FAAH rescue
      const rescue = m.faah_rescue || {};
      const comps = rescue.compounds || [];
      if (comps.length) {
        C.dumbbell($("#rescue-chart"), {
          title: "FAAH rescue: predicted probability before and after",
          ariaLabel: "Predicted FAAH activity probability, ChEMBL-only model versus leave-one-out model",
          labelWidth: 150,
          labelHeader: "Compound",
          aLabel: "ChEMBL-only model",
          bLabel: "Leave-one-out model",
          aColorVar: "--series-neutral",
          bColorVar: "--series-3",
          log: true,
          min: 1e-6,
          max: 1,
          ticks: [1e-6, 1e-4, 1e-2, 1],
          tickFmt: (t) => (t >= 1 ? "1" : t.toExponential(0).replace("e-", "e−")),
          valueFmt: (v) => (v < 0.001 ? Number(v).toExponential(1) : Number(v).toFixed(3)),
          showDelta: false,
          data: comps.map((c) => ({
            label: c.name,
            a: c.p_baseline,
            b: c.p_leave_one_out,
            extra: [["reported", c.reported_uM + " µM"],
                    ["fold change", fmt(c.fold_change)]],
          })),
        });
        const s = rescue.summary || {};
        $("#rescue-note").innerHTML = `<strong>Median predicted probability moves from
          ${Number(s.baseline_median_p).toExponential(1)} to ${num(s.loo_median_p, 3)}</strong>
          once the model has seen just five curated literature actives — none of which
          overlap ChEMBL's FAAH set. The blind spot is a <em>data-coverage</em> problem,
          not a limitation of the representation. Silybin is the exception and stays
          low: it is a large flavonolignan far outside the training distribution.`;
      }

    } catch (e) { errorBox($("#split-chart"), e); }
  }

  function renderMetricsTable(metrics) {
    const rows = metrics.filter((r) => !metricsSplit || r.split === metricsSplit)
      .slice()
      .sort((a, b) => b.pr_auc - a.pr_auc);
    $("#metrics-table").innerHTML = `<div class="table-wrap" style="max-height:460px;overflow-y:auto">
      <table><thead><tr>
        <th scope="col">Target</th><th scope="col">Scope</th><th scope="col">Label</th>
        <th scope="col">Model</th><th scope="col">Split</th>
        <th scope="col" class="num">n</th><th scope="col" class="num">PR-AUC</th>
        <th scope="col" class="num">ROC-AUC</th><th scope="col" class="num">MCC</th>
        <th scope="col" class="num">EF 1%</th>
      </tr></thead><tbody>
      ${rows.map((r) => `<tr>
        <th scope="row">${targetBadge(r.target)}</th>
        <td>${esc(r.scope)}</td>
        <td class="mono dim">${esc(r.label)}</td>
        <td>${esc(r.model === "random_forest" ? "Random forest" : "LightGBM")}</td>
        <td>${r.split === "scaffold"
              ? '<span class="badge badge-good">scaffold</span>'
              : '<span class="badge">random</span>'}</td>
        <td class="num">${int(r.n)}</td>
        <td class="num"><strong>${num(r.pr_auc)}</strong></td>
        <td class="num">${num(r.roc_auc)}</td>
        <td class="num">${num(r.mcc)}</td>
        <td class="num">${num(r.ef_1pct, 2)}</td>
      </tr>`).join("")}
      </tbody></table></div>`;
  }

  /* ---------------------------------------------------------- 05 validation */
  let valTab = "held";
  let valBound = false;

  async function renderValidation() {
    try {
      const v = await api("/api/validation");
      const c = v.corrected || {};
      const held = v.held_out || [];
      const kg = v.kg || {};
      const t = c.tests || {};
      const thr = c.threshold_experiment || null;

      renderPubChem(v.pubchem || {}, v.expanded || {});

      $("#validation-stats").innerHTML = [
        stat("Fair denominator", int(c.fair_denominator),
             `${int(c.raw_claims)} claims, ${int(c.excluded_weaker_than_cut)} below their model&rsquo;s cut-off`, 1),
        stat("Median percentile", num(c.median_percentile, 1),
             "against 50 for chance", 2),
        stat("Not separable from chance", `p = ${num(t.wilcoxon_p, 2)}`,
             `Wilcoxon; sign test p = ${num(t.sign_test_p, 2)}`, 2),
        stat("Reaching the top 10%", `${int(c.in_top_10pct)} of ${int(c.fair_denominator)}`,
             "of the 145k library", 3),
      ].join("");

      const ex = c.excluded || [];

      const half = c.p_half_is_percentile || {};


      const cutByTarget = { TRPV1: 1, CB1: 1, FAAH: 10 };
      C.hbar($("#validation-chart"), {
        title: "Held-out literature actives: percentile against the 145k library",
        ariaLabel: "Library percentile for each held-out literature active",
        labelWidth: 250,
        labelHeader: "Compound (claimed target)",
        valueHeader: "Library percentile",
        max: 100,
        data: held.map((r) => {
          const fair = r.value_uM <= (cutByTarget[r.target] || 1);
          return {
            label: `${(r.compound_name || "?").slice(0, 20)} · ${r.target}${fair ? "" : " (below cut)"}`,
            value: r.library_percentile == null ? 0 : r.library_percentile,
            display: r.library_percentile == null ? "—" : num(r.library_percentile, 1),
            color: C.cssVar(TARGET_COLORS[r.target] || "--series-neutral"),
            muted: !fair,
            emphasis: fair && (r.library_percentile || 0) >= 90,
            tip: [
              ["claimed target", r.target],
              ["reported potency", r.value_uM == null ? "—" : r.value_uM + " µM"],
              ["model's cut-off", (cutByTarget[r.target] || 1) + " µM"],
              ["counted in the test?", fair ? "yes" : "no — weaker than the cut"],
              ["model probability", r.p_claimed_target == null ? "—" : Number(r.p_claimed_target).toExponential(2)],
              ["library percentile", num(r.library_percentile, 1)],
              ["applicability domain", num(r.ad_claimed_target, 2)],
            ],
          };
        }),
        legend: ["TRPV1", "CB1", "FAAH"].map((x) => ({ label: x, color: C.cssVar(TARGET_COLORS[x]) }))
          .concat([{ label: "faded = below its model's cut-off, excluded", color: C.cssVar("--series-neutral") }]),
      });

      if (thr) {
        $("#threshold-result").innerHTML = `
          <div class="grid grid-4" style="margin-bottom:14px">
            ${stat("Median percentile", `${num(thr.median_pct_1um, 1)} → ${num(thr.median_pct_10um, 1)}`,
                   `${thr.n_improved} of ${thr.n} compounds improved`, 1)}
            ${stat("True actives", `${num(thr.true_actives_median_before, 1)} → ${num(thr.true_actives_median_after, 1)}`,
                   "the ones it should rank highly", 3)}
            ${stat("Too-weak compounds", `${num(thr.too_weak_median_before, 1)} → ${num(thr.too_weak_median_after, 1)}`,
                   "should rank LOW — this got worse", 2)}
            ${stat("Reaching top 10%", `${int(thr.top10pct_before)} → ${int(thr.top10pct_after)}`,
                   "unchanged", 2)}
          </div>
          <div class="callout callout-warning">
            <h4>Verdict: a real gain, but not the bottleneck</h4>
            <p class="card-note">${esc(thr.verdict)}</p>
          </div>`;
      } else {
        $("#threshold-result").innerHTML = `<div class="empty">
          Not computed for this run.</div>`;
      }

      renderValidationTable(v);
      if (!valBound) {
        valBound = true;
        $$("#view-validation [data-val]").forEach((b) => {
          b.addEventListener("click", () => {
            valTab = b.dataset.val;
            $$("#view-validation [data-val]").forEach((x) =>
              x.setAttribute("aria-pressed", String(x === b)));
            renderValidationTable(v);
          });
        });
      }

      const byType = kg.nodes_by_type || {};
      const byEv = kg.edges_by_evidence || {};
      $("#kg-stats").innerHTML = [
        stat("Nodes", int(kg.nodes), `${int(byType.Compound)} compounds · ${int(byType.Plant)} plants`, 1),
        stat("Edges", int(kg.edges), "every one typed by evidence", 3),
        stat("Measured edges", int(byEv.measured), "ChEMBL bioactivity", 1),
        stat("Predicted edges", int(byEv.predicted), "model output, rendered differently", 2),
      ].join("");

      C.stacked($("#kg-evidence"), {
        title: "Knowledge-graph edges by evidence class",
        ariaLabel: "Edge counts by evidence class",
        labelWidth: 110,
        labelHeader: "Graph",
        data: [{
          label: "All edges",
          measured: byEv.measured || 0,
          literature: byEv.literature || 0,
          predicted: byEv.predicted || 0,
        }],
        series: [
          { key: "measured", label: "Measured (ChEMBL)", colorVar: "--series-1" },
          { key: "literature", label: "Literature (LOTUS DOI)", colorVar: "--series-2" },
          { key: "predicted", label: "Predicted (our model)", colorVar: "--series-3" },
        ],
      });
    } catch (e) { errorBox($("#validation-stats"), e); }
  }

  /* The external validation, rewritten after two problems were found in the
   * first version. Outcomes from PubChem key on the InChIKey skeleton, which
   * collapses stereoisomers, and the plant library holds several entries per
   * skeleton - so the set reported as n=71 for CB1 was 31 distinct structures
   * and its p-value was overstated by about five orders of magnitude. And
   * requiring the test compound to be in the plant library is what left TRPV1
   * and FAAH single-class. Both are fixed here, and both numbers are shown. */

  const TARGET_VERDICT = {
    CB1: ["Validated", "The strongest claim in the project: it holds on chemistry " +
          "genuinely unlike anything it trained on."],
    TRPV1: ["Validated near training chemistry only", "Of 88 measured inactives, " +
            "none survive an AD&nbsp;&lt;&nbsp;0.4 filter, so performance on novel " +
            "chemistry is untested. Publication bias, not a model fault: nobody " +
            "publishes a plant compound that failed against the capsaicin receptor."],
    FAAH: ["Measured, and weak", "0.64 is better than chance and not much more. " +
           "Its internal score of 0.93 was inflated by decoy negatives; this is " +
           "the first test that could show that."],
  };

  function ciCell(d) {
    if (!d || !d.ci95) return "&mdash;";
    return "[" + num(d.ci95[0], 3) + ", " + num(d.ci95[1], 3) + "]";
  }
  function pCell(d) {
    if (!d || d.mannwhitney_p == null) return "&mdash;";
    return d.mannwhitney_p < 0.01
      ? Number(d.mannwhitney_p).toExponential(1)
      : Number(d.mannwhitney_p).toFixed(3);
  }

  function renderPubChem(pc, exp) {
    const host = $("#pubchem-card");
    const targets = (exp || {}).targets || {};
    if (!Object.keys(targets).length) {
      // fall back to the old single-target summary if the expanded run is absent
      host.hidden = !(pc && pc.cb1);
      return;
    }
    host.hidden = false;

    const SCOPES = [
      ["external_all", "All untrained structures"],
      ["external_novel_ad_lt_0.5", "Novel chemistry (AD &lt; 0.5)"],
      ["external_novel_ad_lt_0.4", "Novel chemistry (AD &lt; 0.4)"],
      ["in_library", "Plant library only"],
    ];

    const cards = ["CB1", "TRPV1", "FAAH"].filter(function (k) { return targets[k]; })
      .map(function (tn) {
        const blocks = targets[tn];
        const a = blocks.external_all || {};
        const v = TARGET_VERDICT[tn] || ["", ""];
        const rows = SCOPES.map(function (s) {
          const d = blocks[s[0]];
          if (!d) return "";
          if (!d.scoreable) {
            return "<tr><td>" + s[1] + "</td><td class='num'>" + int(d.n) + "</td>" +
              "<td class='num dim'>" + int(d.n_active) + " / " + int(d.n_inactive) + "</td>" +
              "<td colspan='3' class='dim'>single-class &mdash; ROC-AUC undefined</td></tr>";
          }
          return "<tr><td>" + s[1] + "</td><td class='num'>" + int(d.n) + "</td>" +
            "<td class='num dim'>" + int(d.n_active) + " / " + int(d.n_inactive) + "</td>" +
            "<td class='num'><strong>" + num(d.roc_auc, 3) + "</strong></td>" +
            "<td>" + ciCell(d) + "</td><td class='num'>" + pCell(d) + "</td></tr>";
        }).join("");

        return '<div class="card" style="margin-top:14px">' +
          '<div class="card-head"><div><h3>' + targetBadge(tn) + " &nbsp;" +
          esc(v[0]) + "</h3>" +
          '<p class="card-note">' + v[1] + "</p></div>" +
          '<div class="stat-value" style="color:var(' + TARGET_COLORS[tn] + ')">' +
          (a.scoreable ? num(a.roc_auc, 3) : "&mdash;") + "</div></div>" +
          '<div class="table-wrap" style="margin-top:10px"><table><thead><tr>' +
          '<th scope="col">Test set</th><th scope="col" class="num">n</th>' +
          '<th scope="col" class="num">act / inact</th>' +
          '<th scope="col" class="num">ROC-AUC</th><th scope="col">95% CI</th>' +
          '<th scope="col" class="num">p</th></tr></thead><tbody>' + rows +
          "</tbody></table></div></div>";
      }).join("");

    host.innerHTML =
      '<div class="card-head"><div>' +
      "<h2>All three targets, tested against measured outcomes</h2>" +
      '<p class="card-note">PubChem BioAssay records <em>Inactive</em> as well as ' +
      "<em>Active</em>, so compounds that were screened and did nothing are in the " +
      "data. That gives a real denominator. Every structure here was measured by " +
      "somebody else and appears in no training set, and the table is deduplicated " +
      "to one row per chemical skeleton.</p></div></div>" +
      cards +
      '<div class="callout callout-warning" style="margin-top:14px">' +
      "<h4>Two corrections we made to our own result</h4>" +
      "<p><strong>We were counting molecules twice.</strong> Outcomes key on the " +
      "InChIKey skeleton, which collapses stereoisomers and salts, and the library " +
      "holds several entries per skeleton. The CB1 set we first reported as n=71 is " +
      "31 distinct structures. The effect size barely moved (0.952 &rarr; 0.946) but " +
      "the p-value fell about five orders of magnitude, from 5&times;10<sup>&minus;9</sup> " +
      "to 2.5&times;10<sup>&minus;4</sup>. Everything above is deduplicated.</p>" +
      "<p style='margin-top:10px'><strong>We were over-constraining the test.</strong> " +
      "Requiring test compounds to be plant natural products is what left TRPV1 and " +
      "FAAH single-class. The question &ldquo;does this rank real actives above real " +
      "inactives&rdquo; never needed that. Dropping it took TRPV1 from 11 scoreable " +
      "structures to 966 and FAAH from 2 to 621.</p>" +
      "<p style='margin-top:10px'>The external set is drug-like, not plant-like, so it " +
      "validates <em>the models</em> rather than the library. Both scopes are in the " +
      "tables. And this is retrospective throughout: a blind test, not a prospective " +
      "one. No compound was assayed by this team.</p></div>";
  }

  function renderValidationTable(v) {
    const rows = valTab === "held" ? (v.held_out || []) : (v.seen || []);
    $("#validation-table").innerHTML = `<div class="table-wrap" style="max-height:480px;overflow-y:auto">
      <table><thead><tr>
        <th scope="col">Compound</th><th scope="col">Target</th>
        <th scope="col">Interaction</th><th scope="col" class="num">Reported µM</th>
        <th scope="col" class="num">Model p</th><th scope="col" class="num">Library %ile</th>
        <th scope="col" class="num">AD</th><th scope="col" class="num">p(DILI)</th>
      </tr></thead><tbody>
      ${rows.map((r) => {
        const pct = r.library_percentile;
        const cls = pct >= 90 ? "badge-good" : pct >= 50 ? "badge-warning" : "badge-critical";
        return `<tr>
          <th scope="row" style="font-weight:500">${esc(r.compound_name)}</th>
          <td>${targetBadge(r.target)}</td>
          <td class="dim">${esc(r.interaction || "—")}</td>
          <td class="num">${r.value_uM == null ? "—" : num(r.value_uM, 3)}</td>
          <td class="num">${r.p_claimed_target == null ? "—" : Number(r.p_claimed_target).toExponential(1)}</td>
          <td class="num"><span class="badge ${cls}">${num(pct, 1)}</span></td>
          <td class="num">${num(r.ad_claimed_target, 2)}</td>
          <td class="num">${num(r.p_dili, 2)}</td>
        </tr>`;
      }).join("")}
      </tbody></table></div>`;
  }

  /* -------------------------------------------------------------- 06 safety */
  let sevFilter = "";
  async function renderSafety() {
    try {
      const s = await api("/api/safety");
      const f = s.flagged || {};
      const w = s.weights || {};

      $("#safety-stats").innerHTML = [
        stat("Severity-3 flagged", int(f.severity3),
             `${(((f.severity3 || 0) / (f.screened || 1)) * 100).toFixed(1)}% of the screen — direct NAPQI-type precedent`, 2),
        stat("PAINS flagged", int(f.pains),
             `${(((f.pains || 0) / (f.screened || 1)) * 100).toFixed(1)}% — assay-interference chemistry`, 1),
        stat("Alert patterns", (s.alerts || []).length,
             "SMARTS, read live from src/safety.py", 3),
      ].join("");

      $("#weight-stats").innerHTML = [
        stat("Alert penalty", w.alert_penalty, "per unit of summed severity"),
        stat("AD floor", w.ad_floor, "below this is extrapolation"),
        stat("CNS bonus", w.cns_bonus, "central analgesia is the hypothesis"),
        stat("Hepatic weight", w.dili_weight, "raised from 0.30; capped by a 0.65 external score"),
      ].join("");

      $("#alert-count").textContent = (s.alerts || []).length;
      renderAlerts(s.alerts || []);
      $$("#view-safety .seg [data-sev]").forEach((b) => {
        b.addEventListener("click", () => {
          sevFilter = b.dataset.sev;
          $$("#view-safety .seg [data-sev]").forEach((x) =>
            x.setAttribute("aria-pressed", String(x === b)));
          renderAlerts(s.alerts || []);
        });
      });
    } catch (e) { errorBox($("#safety-stats"), e); }
  }

  function renderAlerts(alerts) {
    const rows = alerts.filter((a) => !sevFilter || String(a.severity) === sevFilter);
    const sevBadge = (n) =>
      n === 3 ? '<span class="badge badge-critical">3 · direct precedent</span>'
      : n === 2 ? '<span class="badge badge-serious">2 · established route</span>'
      : '<span class="badge badge-warning">1 · context-dependent</span>';
    $("#alerts-table").innerHTML = `<div class="table-wrap" style="max-height:520px;overflow-y:auto">
      <table><thead><tr>
        <th scope="col">Alert</th><th scope="col">Severity</th>
        <th scope="col">SMARTS</th><th scope="col">Mechanism</th>
      </tr></thead><tbody>
      ${rows.map((a) => `<tr>
        <th scope="row" class="mono" style="font-weight:500">${esc(a.name)}</th>
        <td>${sevBadge(a.severity)}</td>
        <td class="mono dim"><span class="truncate" style="max-width:230px">${esc(a.smarts)}</span></td>
        <td class="card-note">${esc(a.note)}</td>
      </tr>`).join("")}
      </tbody></table></div>`;
  }

  /* -------------------------------------------------------------- 06 screen */
  async function renderScreen() {
    bindScreenControls();
    await loadScreenPage();
    try {
      const d = await api("/api/distribution?field=composite_score&bins=46");
      C.histogram($("#dist-chart"), {
        title: "Composite score distribution across the screen",
        ariaLabel: "Histogram of composite score across all screened compounds",
        bins: d.bins || [],
      });
    } catch (e) { /* the explorer still works without the histogram */ }
  }

  let controlsBound = false;
  function bindScreenControls() {
    if (controlsBound) return;
    controlsBound = true;

    $$("#view-screen .seg [data-source]").forEach((b) => {
      b.addEventListener("click", () => {
        state.source = b.dataset.source;
        state.page = 1;
        state.sort = state.source === "controls" ? "target_score" : "rank";
        state.desc = state.source === "controls";
        $$("#view-screen .seg [data-source]").forEach((x) =>
          x.setAttribute("aria-pressed", String(x === b)));
        loadScreenPage();
      });
    });

    let debounce;
    $("#f-q").addEventListener("input", (e) => {
      clearTimeout(debounce);
      debounce = setTimeout(() => {
        state.filters.q = e.target.value.trim();
        state.page = 1;
        loadScreenPage();
      }, 260);
    });

    [["f-target", "target"], ["f-tier", "tier"]].forEach(([id, key]) => {
      $("#" + id).addEventListener("change", (e) => {
        state.filters[key] = e.target.value;
        state.page = 1;
        loadScreenPage();
      });
    });
    [["f-domain", "in_domain"], ["f-alerts", "no_alerts"], ["f-cns", "cns_only"]].forEach(([id, key]) => {
      $("#" + id).addEventListener("change", (e) => {
        state.filters[key] = e.target.checked;
        state.page = 1;
        loadScreenPage();
      });
    });

    $("#f-reset").addEventListener("click", () => {
      state.filters = { q: "", target: "", tier: "", in_domain: false, no_alerts: false, cns_only: false };
      state.page = 1;
      $("#f-q").value = "";
      $("#f-target").value = "";
      $("#f-tier").value = "";
      ["f-domain", "f-alerts", "f-cns"].forEach((id) => { $("#" + id).checked = false; });
      loadScreenPage();
    });

    $("#page-prev").addEventListener("click", () => {
      if (state.page > 1) { state.page--; loadScreenPage(); }
    });
    $("#page-next").addEventListener("click", () => { state.page++; loadScreenPage(); });
  }

  function screenQuery() {
    const f = state.filters;
    const p = new URLSearchParams();
    if (f.q) p.set("q", f.q);
    if (f.target) p.set("target", f.target);
    if (f.tier) p.set("tier", f.tier);
    if (f.in_domain) p.set("in_domain", "true");
    if (f.no_alerts) p.set("no_alerts", "true");
    if (f.cns_only) p.set("cns_only", "true");
    p.set("sort", state.sort);
    p.set("desc", String(state.desc));
    p.set("page", String(state.page));
    p.set("page_size", String(state.pageSize));
    return `/api/compounds/${state.source}?${p.toString()}`;
  }

  const COLUMNS = [
    { key: "rank", label: "#", type: "rank", sortable: true },
    { key: "compound_name", label: "Compound", type: "name" },
    { key: "best_target", label: "Target", type: "target" },
    { key: "target_score", label: "Target score", type: "num4", sortable: true },
    { key: "composite_score", label: "Composite", type: "num3", sortable: true },
    { key: "ad_best", label: "AD", type: "num2", sortable: true },
    { key: "p_dili", label: "p(DILI)", type: "num2", sortable: true },
    { key: "n_reactive_alerts", label: "Alerts", type: "alerts", sortable: true },
    { key: "mw", label: "MW", type: "num0", sortable: true },
    { key: "clogp", label: "cLogP", type: "num1", sortable: true },
    { key: "source_genera", label: "Genus", type: "genus" },
    { key: "tier", label: "Tier", type: "tier" },
  ];

  async function loadScreenPage() {
    const host = $("#screen-table");
    host.classList.add("loading"); // hold the previous render, no skeleton flash
    try {
      const url = screenQuery();
      const res = await fetch(url);
      if (!res.ok) throw new Error(`${res.status} ${res.statusText}`);
      const data = await res.json();

      const start = (data.page - 1) * data.page_size;
      const end = Math.min(start + data.rows.length, data.total);
      $("#screen-info").textContent = data.total
        ? `${int(start + 1)}–${int(end)} of ${int(data.total)} compounds`
        : "No compounds match these filters";
      $("#page-prev").disabled = data.page <= 1;
      $("#page-next").disabled = end >= data.total;

      if (!data.rows.length) {
        host.innerHTML = '<div class="empty">No compounds match these filters.</div>';
        return;
      }

      host.innerHTML = `<div class="table-wrap" style="border:0;max-height:620px;overflow-y:auto">
        <table><thead><tr>
        ${COLUMNS.map((c) => {
          const sorted = state.sort === c.key;
          return `<th scope="col" class="${c.sortable ? "sortable" : ""} ${c.type.indexOf("num") === 0 || c.type === "rank" ? "num" : ""}"
            ${sorted ? `aria-sort="${state.desc ? "descending" : "ascending"}"` : ""}
            ${c.sortable ? `data-sort="${c.key}" tabindex="0" role="button"` : ""}>
            ${esc(c.label)}${c.sortable ? ` <span class="sort-caret">${sorted ? (state.desc ? "▾" : "▴") : "⇅"}</span>` : ""}
          </th>`;
        }).join("")}
        </tr></thead><tbody>
        ${data.rows.map((r) => `<tr class="clickable" tabindex="0" data-key="${esc(r.std_inchikey)}">
          ${COLUMNS.map((c) => cell(r, c)).join("")}
        </tr>`).join("")}
        </tbody></table></div>`;

      $$("th[data-sort]", host).forEach((th) => {
        const act = () => {
          const key = th.dataset.sort;
          if (state.sort === key) state.desc = !state.desc;
          else { state.sort = key; state.desc = key !== "rank"; }
          state.page = 1;
          loadScreenPage();
        };
        th.addEventListener("click", act);
        th.addEventListener("keydown", (e) => {
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); act(); }
        });
      });
      wireRowClicks(host);
    } catch (e) {
      errorBox(host, e);
    } finally {
      host.classList.remove("loading");
    }
  }

  function cell(r, c) {
    const v = r[c.key];
    switch (c.type) {
      case "rank": return `<td class="num dim">${r.rank}</td>`;
      case "name": return `<th scope="row" class="col-name" style="font-weight:500"><span class="truncate">${esc(compoundName(r))}</span></th>`;
      case "target": return `<td>${targetBadge(v)}</td>`;
      case "num4": return `<td class="num">${num(v, 4)}</td>`;
      case "num3": return `<td class="num">${num(v, 3)}</td>`;
      case "num2": return `<td class="num">${num(v, 2)}</td>`;
      case "num1": return `<td class="num">${num(v, 1)}</td>`;
      case "num0": return `<td class="num">${v == null ? "—" : Math.round(v)}</td>`;
      case "alerts": {
        const n = v || 0;
        const sev = r.max_alert_severity || 0;
        if (!n) return '<td class="num dim">0</td>';
        const cls = sev >= 3 ? "badge-critical" : sev === 2 ? "badge-serious" : "badge-warning";
        return `<td class="num"><span class="badge ${cls}">${n}</span></td>`;
      }
      case "genus": return `<td><span class="truncate" style="max-width:150px">${esc((v || "").split("|")[0] || "—")}</span></td>`;
      case "tier": return `<td><span class="badge">${esc(v === "pain_seed" ? "Seed" : "Background")}</span></td>`;
      default: return `<td>${esc(v)}</td>`;
    }
  }

  /* ------------------------------------------------------------ 09 glossary */
  let glossaryCat = "";
  let glossaryQ = "";

  let glossaryBound = false;
  /* ------------------------------------------------------------- 08 demo */
  /* The only view that runs a model. Everything else reads a file the
   * pipeline already wrote; this one posts a structure and waits. The first
   * call loads the models server-side, so it is slow once and fast after. */

  const DEMO_GRADE = {
    "validated": "validated externally, including on novel chemistry",
    "near-training only": "validated near training chemistry only",
    "weak": "weak - treat as a hint, not a result",
  };

  function demoBar(p, tone) {
    var pct = Math.max(0, Math.min(100, Number(p) * 100));
    return '<div style="height:10px;border-radius:999px;background:var(--surface-3);overflow:hidden;margin-top:8px">'
      + '<div style="height:100%;width:' + pct.toFixed(1) + '%;border-radius:999px;background:var(' + tone + ')"></div>'
      + '</div>';
  }

  function demoReading(p, ad) {
    if (ad < 0.25) return ["OUT OF DOMAIN", "bad"];
    if (p >= 0.5) return ["ACTIVE", "ok"];
    if (p >= 0.15) return ["BORDERLINE", "warn"];
    return ["INACTIVE", "muted"];
  }

  function renderDemoResult(r) {
    var host = $("#demo-out");
    var acts = r.activity || {};
    var order = ["TRPV1", "CB1", "FAAH"].filter(function (k) { return acts[k]; });

    var cards = order.map(function (tname) {
      var a = acts[tname];
      var rd = demoReading(a.p, a.ad);
      var grade = (a.validation || {}).grade || "";
      return '<div class="card" style="flex:1;min-width:232px">'
        + '<div class="stat-label">' + esc(tname) + '</div>'
        + '<div class="stat-value" style="color:var(' + TARGET_COLORS[tname] + ')">' + num(a.p, 3) + '</div>'
        + demoBar(a.p, TARGET_COLORS[tname])
        + '<p style="margin:10px 0 4px;font-weight:680;letter-spacing:0.02em">' + esc(rd[0]) + '</p>'
        + '<p class="card-note" style="margin:0">applicability domain ' + num(a.ad, 2) + '</p>'
        + '<p class="card-note" style="margin:6px 0 0">' + esc(DEMO_GRADE[grade] || grade) + '</p>'
        + '</div>';
    }).join("");

    var s = r.safety || {};
    var alerts = s.alerts || [];
    var sev = s.severities || [];
    var alertRows = alerts.length
      ? alerts.map(function (nm, i) {
          return '<li><b>severity ' + sev[i] + '</b> &mdash; ' + esc(String(nm).replace(/_/g, " ")) + '</li>';
        }).join("")
      : '<li>no reactive-metabolite alerts</li>';

    var h = r.hepatic || {};
    var heads = h.heads || {};
    var headRows = Object.keys(heads).map(function (k) {
      return '<tr><td>' + esc(k) + '</td><td class="num">' + num(heads[k].p, 3) + '</td>'
        + '<td class="card-note">' + esc(heads[k].measures || "") + '</td></tr>';
    }).join("");
    var hepHigh = Number(h.p_hepatic) >= 0.6;

    var c = r.composite || {};
    var terms = c.terms || {};
    var termRows = Object.keys(terms).map(function (k) {
      return '<tr><td>' + esc(titleCase(k)) + '</td><td class="num">' + num(terms[k], 4) + '</td></tr>';
    }).join("");

    host.innerHTML =
      '<div class="card" style="margin-top:18px">'
      + '<h3>Stage 1 &mdash; the standardised structure</h3>'
      + '<p style="word-break:break-all;margin:6px 0"><code>' + esc(r.std_smiles) + '</code></p>'
      + '<p class="card-note">InChIKey ' + esc(r.inchikey) + '</p></div>'

      + '<h3 style="margin:22px 0 10px">Stage 2 &mdash; active or inactive, on each target</h3>'
      + '<div style="display:flex;gap:14px;flex-wrap:wrap">' + cards + '</div>'

      + '<div class="card" style="margin-top:18px">'
      + '<h3>Stage 3 &mdash; structural alerts <span class="card-note">(this layer holds the veto)</span></h3>'
      + '<ul style="margin:8px 0 0 18px">' + alertRows + '</ul>'
      + '<p class="card-note" style="margin-top:8px">PAINS patterns: ' + int(s.n_pains) + '</p>'
      + (s.disqualified
          ? '<p style="font-weight:680;margin-top:10px">DISQUALIFIED from the shortlist &mdash; ' + esc(s.reason || "") + '</p>'
          : '<p style="font-weight:680;margin-top:10px">No disqualifying alert.</p>')
      + '</div>'

      + '<div class="card" style="margin-top:18px">'
      + '<h3>Stage 4 &mdash; hepatic liability</h3>'
      + '<div class="stat-value">' + num(h.p_hepatic, 3) + '</div>'
      + demoBar(h.p_hepatic, hepHigh ? "--series-4" : "--series-3")
      + '<table class="table" style="margin-top:12px"><thead><tr><th>Head</th>'
      + '<th class="num">Score</th><th>What it measures</th></tr></thead><tbody>'
      + headRows + '</tbody></table>'
      + '<p class="card-note" style="margin-top:10px">Combined as the maximum over the three heads. '
      + esc(h.note || "") + '</p></div>'

      + '<div class="card" style="margin-top:18px">'
      + '<h3>Stage 5 &mdash; composite rank</h3>'
      + '<div class="stat-value">' + num(c.score, 5) + '</div>'
      + '<table class="table" style="margin-top:12px"><tbody>' + termRows + '</tbody></table>'
      + '<p style="margin-top:12px;font-weight:680">'
      + (c.shortlist_eligible ? 'Shortlist eligible.' : 'Not shortlist eligible.') + '</p>'
      + '<p class="card-note">A ranking, not a verdict: it says where this molecule belongs '
      + 'in a queue for someone with a laboratory.</p></div>';
  }

  async function runDemo(smiles, note) {
    var out = $("#demo-out");
    $("#demo-note").textContent = "Running the pipeline on this structure...";
    out.innerHTML = '<div class="card" style="margin-top:18px">Scoring. The first run loads '
      + 'the models and can take a few seconds; later runs are immediate.</div>';
    try {
      var r = await api("/api/predict?smiles=" + encodeURIComponent(smiles));
      renderDemoResult(r);
      $("#demo-note").textContent = note || "";
    } catch (e) {
      out.innerHTML = '<div class="card empty" style="margin-top:18px">Could not score that structure.<br>'
        + esc(String((e && e.message) || e)) + '</div>';
      $("#demo-note").textContent = "";
    }
  }

  async function renderDemo() {
    var sel = $("#demo-pick");
    var mols = [];
    try {
      mols = await api("/api/demo_molecules");
    } catch (e) {
      mols = [];
    }
    state.demoMols = mols;
    sel.innerHTML = mols.length
      ? mols.map(function (m, i) {
          return '<option value="' + i + '">' + esc(m.label) + '</option>';
        }).join("")
      : '<option value="">none available</option>';

    function pick() {
      var m = state.demoMols[Number(sel.value)];
      if (!m) return;
      $("#demo-smiles").value = m.smiles;
      runDemo(m.smiles, m.note);
    }
    sel.onchange = pick;
    $("#demo-run").onclick = function () {
      var raw = String($("#demo-smiles").value || "").trim();
      if (raw) runDemo(raw, "Custom structure.");
    };
    $("#demo-smiles").onkeydown = function (ev) {
      if (ev.key === "Enter") $("#demo-run").click();
    };
    if (mols.length) {
      sel.value = "0";
      pick();
    }
  }

  function renderGlossaryPage() {
    renderGlossary();
    if (glossaryBound) return;
    glossaryBound = true;
    $("#g-q").addEventListener("input", (e) => {
      glossaryQ = e.target.value.trim().toLowerCase();
      renderGlossary();
    });
    $$("#view-glossary [data-gcat]").forEach((b) => {
      b.addEventListener("click", () => {
        glossaryCat = b.dataset.gcat;
        $$("#view-glossary [data-gcat]").forEach((x) =>
          x.setAttribute("aria-pressed", String(x === b)));
        renderGlossary();
      });
    });
  }
  function renderGlossary() {
    const all = window.GLOSSARY || [];
    const rows = all.filter((g) => {
      if (glossaryCat && g.cat !== glossaryCat) return false;
      if (!glossaryQ) return true;
      return (g.term + " " + g.plain + " " + g.ml).toLowerCase().indexOf(glossaryQ) >= 0;
    });
    if (!rows.length) {
      $("#glossary-table").innerHTML = '<div class="empty">No terms match that search.</div>';
      return;
    }
    $("#glossary-table").innerHTML = `<div class="table-wrap">
      <table><thead><tr>
        <th scope="col" style="min-width:150px">Term</th>
        <th scope="col" style="min-width:250px">What it means</th>
        <th scope="col" style="min-width:260px">The ML equivalent</th>
      </tr></thead><tbody>
      ${rows.map((g) => `<tr>
        <th scope="row" style="font-weight:600">${esc(g.term)}
          <div><span class="badge" style="margin-top:4px">${esc(g.cat)}</span></div></th>
        <td class="card-note">${g.plain}</td>
        <td class="card-note">${g.ml}</td>
      </tr>`).join("")}
      </tbody></table></div>
      <p class="card-note" style="margin-top:9px">${rows.length} of ${all.length} terms.</p>`;
  }


  /* --------------------------------------------------------------- modal */
  function wireRowClicks(root) {
    $$("tr.clickable", root).forEach((tr) => {
      const open = () => openCompound(tr.dataset.key);
      tr.addEventListener("click", open);
      tr.addEventListener("keydown", (e) => {
        if (e.key === "Enter") { e.preventDefault(); open(); }
      });
    });
  }

  let lastFocus = null;
  async function openCompound(key) {
    if (!key) return;
    lastFocus = document.activeElement;
    const modal = $("#modal");
    modal.hidden = false;
    $("#modal-title").textContent = "Loading…";
    $("#modal-body").innerHTML = '<div class="empty"><span class="spinner"></span></div>';
    $("#modal-close").focus();
    C.hideTip();

    try {
      const r = await api(`/api/compound/${encodeURIComponent(key)}`);
      $("#modal-title").textContent = compoundName(r);
      const dois = r.reference_dois || [];
      const alerts = (r.reactive_alerts || "").split("|").filter(Boolean);
      const pains = (r.pains_flags || "").split("|").filter(Boolean);

      $("#modal-body").innerHTML = `
        <div class="pill-row" style="margin-bottom:14px">
          ${targetBadge(r.best_target)}
          <span class="badge">${esc(titleCase(r.tier))}</span>
          ${r.in_domain ? '<span class="badge badge-good">in domain</span>'
                        : '<span class="badge badge-warning">extrapolation</span>'}
          ${r.cns_like ? '<span class="badge badge-good">CNS-like</span>' : ""}
        </div>

        <div class="grid grid-2" style="gap:14px">
          <div>
            <h4>Model scores</h4>
            <dl class="kv" style="margin-top:7px">
              <dt>TRPV1</dt><dd>${num(r.p_TRPV1, 4)} <span class="dim">(AD ${num(r.ad_TRPV1, 2)})</span></dd>
              <dt>CB1</dt><dd>${num(r.p_CB1, 4)} <span class="dim">(AD ${num(r.ad_CB1, 2)})</span></dd>
              <dt>FAAH</dt><dd>${num(r.p_FAAH, 4)} <span class="dim">(AD ${num(r.ad_FAAH, 2)})</span></dd>
              <dt>p(DILI)</dt><dd>${num(r.p_dili, 3)}</dd>
              <dt>Composite</dt><dd><strong>${num(r.composite_score, 4)}</strong></dd>
            </dl>
          </div>
          <div>
            <h4>Properties</h4>
            <dl class="kv" style="margin-top:7px">
              <dt>MW</dt><dd>${num(r.mw, 1)}</dd>
              <dt>cLogP</dt><dd>${num(r.clogp, 2)}</dd>
              <dt>TPSA</dt><dd>${num(r.tpsa, 1)}</dd>
              <dt>QED</dt><dd>${num(r.qed, 3)}</dd>
              <dt>Lipinski violations</dt><dd>${r.lipinski_violations == null ? "—" : r.lipinski_violations}</dd>
            </dl>
          </div>
        </div>

        <h4 style="margin-top:16px">Structure</h4>
        <div class="smiles" style="margin-top:6px">${esc(r.std_smiles)}</div>
        <div class="dim" style="font-size:0.76rem;margin-top:5px">InChIKey ${esc(r.std_inchikey)}</div>

        <h4 style="margin-top:16px">Source organisms</h4>
        <p class="card-note">
          ${esc((r.source_genera || "").split("|").join(", ") || "—")}
          ${r.n_source_organisms ? `<span class="dim"> · ${r.n_source_organisms} organism(s)</span>` : ""}
        </p>
        ${r.traditional_pain_use ? `<p class="card-note"><strong>Traditional use:</strong> ${esc(r.traditional_pain_use)}</p>` : ""}
        ${r.npclass_class ? `<p class="card-note"><strong>Class:</strong> ${esc(String(r.npclass_class).split("$").join(" · "))}</p>` : ""}

        <h4 style="margin-top:16px">Safety flags</h4>
        ${alerts.length || pains.length
          ? `<div class="pill-row" style="margin-top:6px">
               ${alerts.map((a) => `<span class="badge ${r.max_alert_severity >= 3 ? "badge-critical" : "badge-serious"}">${esc(a)}</span>`).join("")}
               ${pains.map((p) => `<span class="badge badge-warning">PAINS: ${esc(p)}</span>`).join("")}
             </div>`
          : '<p class="card-note">No reactive-metabolite alert and no PAINS flag.</p>'}

        <h4 style="margin-top:16px">Literature provenance</h4>
        ${dois.length
          ? `<p class="card-note">${dois.length} referenced structure–organism pair(s) in LOTUS.</p>
             <div class="pill-row" style="margin-top:6px">
               ${dois.slice(0, 8).map((d) => `<a class="badge mono" href="https://doi.org/${encodeURIComponent(d.replace(/^https?:\/\/doi\.org\//, ""))}" target="_blank" rel="noopener">${esc(d.replace(/^https?:\/\/doi\.org\//, ""))}</a>`).join("")}
               ${dois.length > 8 ? `<span class="badge dim">+${dois.length - 8} more</span>` : ""}
             </div>`
          : '<p class="card-note dim">No DOI recorded for this compound.</p>'}

        <div class="callout callout-warning" style="margin-top:18px;font-size:0.82rem">
          These are model predictions on a compound that has not been assayed here.
          Nothing on this page establishes activity or safety.
        </div>`;
    } catch (e) {
      $("#modal-title").textContent = "Not found";
      $("#modal-body").innerHTML = `<div class="empty">${esc(e.message)}</div>`;
    }
  }

  function closeModal() {
    $("#modal").hidden = true;
    if (lastFocus && lastFocus.focus) lastFocus.focus();
  }

  /* ----------------------------------------------------------------- init */
  async function checkArtefacts() {
    try {
      const h = await api("/api/health");
      const missing = Object.keys(h.artefacts || {}).filter((k) => !h.artefacts[k]);
      if (missing.length) {
        const b = $("#artefact-banner");
        b.hidden = false;
        b.innerHTML = `<strong>Missing pipeline output:</strong>
          ${missing.map((m) => `<code>${esc(m)}</code>`).join(", ")}.
          Run <code>python src/run_all.py</code> to regenerate. Sections that depend on
          these will be empty.`;
      }
    } catch (e) { /* health is advisory only */ }
  }

  function init() {
    initTheme();
    $$(".nav-btn").forEach((b) => b.addEventListener("click", () => go(b.dataset.view)));
    $$("[data-goto]").forEach((b) => b.addEventListener("click", () => go(b.dataset.goto)));
    $("#modal-close").addEventListener("click", closeModal);
    $("#modal").addEventListener("click", (e) => { if (e.target.id === "modal") closeModal(); });
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape" && !$("#modal").hidden) closeModal();
    });
    window.addEventListener("hashchange", () => {
      const v = location.hash.slice(1);
      if (v && v !== currentView) go(v);
    });

    let resizeTimer;
    window.addEventListener("resize", () => {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(() => {
        state.rendered = {};
        renderView(currentView, true);
      }, 220);
    });

    checkArtefacts();
    const initial = location.hash.slice(1);
    go(["overview", "science", "data", "models", "validation", "safety", "screen", "glossary"]
       .indexOf(initial) >= 0 ? initial : "overview");
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
