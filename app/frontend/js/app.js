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
    source: "shortlist",
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
      decoder: renderDecoder,
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
      name: "TDC DILI benchmark",
      role: "Hepatotoxicity labels",
      accent: 4,
      detail: "475 marketed drugs aggregated by the FDA's National Center for Toxicological Research, derived from Xu et al. 2015. Small, and entirely outside natural-product chemical space.",
      meta: ["462 usable compounds", "50.7% positive", "331 scaffolds"],
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

      // DILI
      const dili = m.dili || {};
      const sc = (dili.splits || {}).scaffold_5x || {};
      const rnd = (dili.splits || {}).random_5x || {};
      $("#dili-n").textContent = int(dili.n);
      $("#dili-metrics").innerHTML = [
        ["Compounds", int(dili.n)],
        ["Scaffolds", int(dili.n_scaffolds)],
        ["Positive rate", ((dili.positive_rate || 0) * 100).toFixed(1) + "%"],
        ["ROC-AUC (scaffold, 5×)", `${num(sc.roc_auc)} ± ${num(sc.roc_auc_sd)}`],
        ["ROC-AUC (random, 5×)", `${num(rnd.roc_auc)} ± ${num(rnd.roc_auc_sd)}`],
        ["MCC (scaffold, 5×)", `${num(sc.mcc)} ± ${num(sc.mcc_sd)}`],
      ].map(([k, v]) => `<dt>${esc(k)}</dt><dd>${v}</dd>`).join("");

      const probes = m.dili_probes || [];
      if (probes.length) {
        const para = probes.find((p) => p.name === "paracetamol");
        if (para) $("#para-score").textContent = num(para.p_dili, 2);
        $("#dili-probe-host").innerHTML = `<h4>Probe scores</h4>
          <p class="card-note">Predicted probability of hepatic liability, with the
          applicability domain beside it.</p>
          <div class="table-wrap" style="margin-top:10px">
          <table><thead><tr>
            <th scope="col">Compound</th><th scope="col" class="num">p(DILI)</th>
            <th scope="col" class="num">AD</th><th scope="col">In benchmark?</th>
          </tr></thead><tbody>
          ${probes.map((p) => {
            const high = p.p_dili >= 0.5;
            const outOfDomain = p.ad < 0.4;
            return `<tr>
              <th scope="row" style="font-weight:500">${esc(p.name)}</th>
              <td class="num">${num(p.p_dili, 3)}</td>
              <td class="num ${outOfDomain ? "dim" : ""}">${num(p.ad, 2)}</td>
              <td><span class="badge ${p.note.indexOf("absent") === 0 ? "badge-critical" : ""}">${esc(p.note)}</span></td>
            </tr>`;
          }).join("")}
          </tbody></table></div>`;
      } else {
        $("#dili-probe-host").innerHTML = `<div class="empty">
          Probe scores not computed.<br>Run <code>python scripts/make_dili_probes.py</code>.</div>`;
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

      renderPubChem(v.pubchem || {});

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
      $("#excluded-table").innerHTML = `<div class="table-wrap"><table>
        <thead><tr><th scope="col">Compound</th><th scope="col" class="num">Reported</th>
        <th scope="col" class="num">Model cut-off</th></tr></thead><tbody>
        ${ex.map((r) => `<tr>
          <th scope="row" style="font-weight:500">${esc(r.compound_name)} <span class="dim">${esc(r.target)}</span></th>
          <td class="num"><strong>${num(r.value_uM, 2)}</strong></td>
          <td class="num dim">${num(r.model_cut_uM, 0)}</td>
        </tr>`).join("")}
        </tbody></table></div>
        <p class="card-note" style="margin-top:8px">Values in µM. Lower means more potent.</p>`;

      const half = c.p_half_is_percentile || {};
      $("#half-table").innerHTML = `<div class="table-wrap"><table>
        <thead><tr><th scope="col">Target</th><th scope="col" class="num">A score of 0.5 is the…</th></tr></thead>
        <tbody>${Object.keys(half).map((k) => `<tr>
          <th scope="row">${targetBadge(k)}</th>
          <td class="num"><strong>${num(half[k], 2)}th</strong> percentile</td>
        </tr>`).join("")}</tbody></table></div>`;

      $("#corrected-verdict").innerHTML = `
        <div class="callout callout-warning">
          <h4>Inconclusive — not a failure, and not a success</h4>
          <p class="card-note">
            On the ${int(c.fair_denominator)} compounds the models could fairly be asked
            about, the median library percentile is
            <strong>${num(c.median_percentile, 1)}</strong> against <strong>50</strong>
            for chance. A Wilcoxon signed-rank test gives
            <strong>p = ${num(t.wilcoxon_p, 3)}</strong>, a sign test
            <strong>p = ${num(t.sign_test_p, 2)}</strong>
            (${int(t.n_above_chance)} of ${int(c.fair_denominator)} above chance).
            <br><br>
            <strong>We cannot reject the null, and we cannot show signal either.</strong>
            Nine compounds is too few to settle it. The honest conclusion is that we have
            no evidence the screen generalises to plant chemistry — which is different
            from evidence that it does not.
          </p>
        </div>`;

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
          Not computed. Run <code>python scripts/experiment_threshold.py</code>.</div>`;
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

  function renderPubChem(pc) {
    const host = $("#pubchem-card");
    const cb = pc.cb1 || {};
    const all = cb.all, novel = cb["novel_ad_lt_0.5"];
    if (!all) {
      host.hidden = true;
      return;
    }
    const sep = pc.separation || {};
    const fmtp = (x) => (x < 0.001 ? Number(x).toExponential(1) : Number(x).toFixed(4));

    $("#pubchem-stats").innerHTML = [
      stat("CB1 · all overlap", num(all.roc_auc, 3),
           `ROC-AUC · n=${all.n} (${all.n_active} active / ${all.n_inactive} inactive) · p=${fmtp(all.p)}`, 3),
      novel ? stat("CB1 · novel chemistry", num(novel.roc_auc, 3),
           `ROC-AUC · n=${novel.n} (${novel.n_active}/${novel.n_inactive}) · p=${fmtp(novel.p)}`, 1) : "",
      stat("Score separation", `${num(sep.active_median, 3)} vs ${Number(sep.inactive_median).toExponential(0)}`,
           "median predicted p(CB1): measured active vs measured inactive", 2),
    ].join("");

    const rej = (pc.notable_correct_rejections || []).map((r) =>
      `<tr><th scope="row" style="font-weight:500">${esc(r.compound_name)}</th>
       <td class="num">${Number(r.model_score).toExponential(1)}</td>
       <td class="num dim">${num(r.ad, 2)}</td></tr>`).join("");
    const un = pc.unscorable || {};
    const unrows = Object.keys(un).map((k) =>
      `<tr><th scope="row">${targetBadge(k)}</th>
       <td class="num">${int(un[k].n)}</td>
       <td class="num">${int(un[k].n_active)} / ${int(un[k].n_inactive)}</td>
       <td class="card-note">${esc(un[k].reason)}</td></tr>`).join("");

    $("#pubchem-detail").innerHTML = `
      <div class="grid grid-2">
        <div>
          <h4>Correctly rejected — the convincing part</h4>
          <p class="card-note">Measured <strong>inactive</strong> at CB1, and the model
          said so. β-caryophyllene is the trap: widely called a &#8220;dietary
          cannabinoid&#8221; because it hits CB2, but it is not a CB1 ligand.</p>
          <div class="table-wrap" style="margin-top:9px"><table>
            <thead><tr><th scope="col">Compound</th><th scope="col" class="num">p(CB1)</th>
            <th scope="col" class="num">AD</th></tr></thead><tbody>${rej}</tbody></table></div>
        </div>
        <div>
          <h4>What could not be scored</h4>
          <p class="card-note">Both overlaps were single-class, so ROC-AUC is undefined.
          <strong>The result is CB1 only</strong> — do not generalise it.</p>
          <div class="table-wrap" style="margin-top:9px"><table>
            <thead><tr><th scope="col">Target</th><th scope="col" class="num">n</th>
            <th scope="col" class="num">act/inact</th><th scope="col">Why</th></tr></thead>
            <tbody>${unrows}</tbody></table></div>
          <div class="callout callout-warning" style="margin-top:12px;font-size:0.82rem">
            PubChem skews active (74.6% base rate here), so enrichment factors are
            meaningless. ROC-AUC is base-rate independent — quote that. And this is
            retrospective: a blind test, not a prospective one.
          </div>
        </div>
      </div>`;
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
        stat("DILI weight", w.dili_weight, "deliberately small"),
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

  /* ------------------------------------------------------------ 07 decoder */
  let glossaryCat = "";
  let glossaryQ = "";

  const MODEL_PURPOSE = [
    { file: "TRPV1_human_label.pkl", task: "P(acts on TRPV1)", rows: "Human TRPV1 measurements", thr: "1 µM" },
    { file: "TRPV1_all_label.pkl", task: "P(acts on TRPV1)", rows: "+ rat and mouse", thr: "1 µM" },
    { file: "CB1_human_label.pkl", task: "P(acts on CB1)", rows: "Human CB1 measurements", thr: "1 µM" },
    { file: "CB1_all_label.pkl", task: "P(acts on CB1)", rows: "+ rat and mouse", thr: "1 µM" },
    { file: "FAAH_human_label.pkl", task: "P(acts on FAAH)", rows: "Human FAAH measurements", thr: "1 µM" },
    { file: "FAAH_all_label.pkl", task: "P(acts on FAAH)", rows: "+ rat and mouse", thr: "1 µM" },
    { file: "FAAH_human_label_10um.pkl", task: "P(acts on FAAH)", rows: "Human FAAH measurements", thr: "10 µM" },
    { file: "FAAH_all_label_10um.pkl", task: "P(acts on FAAH)", rows: "+ rat and mouse", thr: "10 µM" },
    { file: "dili.pkl", task: "P(causes liver injury)", rows: "462 marketed drugs (TDC)", thr: "n/a — given label" },
  ];

  let decoderBound = false;
  async function renderDecoder() {
    renderGlossary();
    if (!decoderBound) {
      decoderBound = true;
      $("#g-q").addEventListener("input", (e) => {
        glossaryQ = e.target.value.trim().toLowerCase();
        renderGlossary();
      });
      $$("#view-decoder [data-gcat]").forEach((b) => {
        b.addEventListener("click", () => {
          glossaryCat = b.dataset.gcat;
          $$("#view-decoder [data-gcat]").forEach((x) =>
            x.setAttribute("aria-pressed", String(x === b)));
          renderGlossary();
        });
      });
    }

    $("#model-purpose-table").innerHTML = `<div class="table-wrap"><table>
      <thead><tr>
        <th scope="col">File</th><th scope="col">Predicts</th>
        <th scope="col">Trained on</th><th scope="col">Positive threshold</th>
      </tr></thead><tbody>
      ${MODEL_PURPOSE.map((m) => `<tr>
        <th scope="row" class="mono" style="font-weight:500">${esc(m.file)}</th>
        <td>${esc(m.task)}</td>
        <td class="card-note">${esc(m.rows)} <span class="dim">+ 3 decoys per active</span></td>
        <td>${esc(m.thr)}</td>
      </tr>`).join("")}
      </tbody></table></div>`;

    try {
      const v = await api("/api/interpretation");
      const m = await api("/api/models");
      renderConclusions(v, m);
      renderAttacks(v, m);
    } catch (e) { errorBox($("#conclusions"), e); }
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

  function claim(n, verdict, title, body, evidence) {
    const cls = verdict === "supported" ? "callout-good"
      : verdict === "partial" ? "callout-warning" : "callout-critical";
    const badge = verdict === "supported" ? '<span class="badge badge-good">supported</span>'
      : verdict === "partial" ? '<span class="badge badge-warning">supported with a caveat</span>'
      : '<span class="badge badge-critical">not supported</span>';
    return `<div class="callout ${cls}" style="margin-bottom:14px">
      <h4>${n}. ${esc(title)} ${badge}</h4>
      <p class="card-note" style="margin-bottom:8px">${body}</p>
      ${evidence ? `<dl class="kv" style="margin-top:8px">${evidence}</dl>` : ""}
    </div>`;
  }
  const ev = (pairs) => pairs.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${v}</dd>`).join("");

  function renderConclusions(v, m) {
    const c = v.controls || {};
    const a = c.measured_active || {};
    const i = c.measured_inactive || {};
    const d = v.domain || {};
    const en = v.enrichment || {};
    const disc = v.discrimination || {};
    const rescue = (m.faah_rescue || {}).summary || {};
    const probes = m.dili_probes || [];
    const para = probes.find((p) => p.name === "paracetamol") || {};
    const kaem = probes.find((p) => p.name === "kaempferol") || {};

    $("#conclusions").innerHTML = [
      claim(1, "supported",
        "The models separate known actives from known inactives",
        `Of the ${int(c.total_held_out)} library compounds removed for leakage control,
         ${int(c.matched_in_curated_table)} appear in the curated activity table — and
         <strong>most of those were measured <em>inactive</em></strong>, not active. Split
         by what the data actually says, the separation is clean.
         <br><br><strong>What it does not establish:</strong> these compounds were in
         training. This is a memorisation check that confirms the pipeline is wired up
         correctly — not evidence of generalisation.`,
        ev([
          ["Measured active", `n = ${int(a.n)} · median score <strong>${num(a.median_score, 4)}</strong> · ${a.above_90 || 0}/${int(a.n)} above 0.9`],
          ["Measured inactive", `n = ${int(i.n)} · median score <strong>${num(i.median_score, 4)}</strong> · ${i.above_50 || 0}/${int(i.n)} above 0.5`],
        ])),

      claim(2, "supported",
        "The decoys fixed the everything-is-active failure",
        `This was the specific risk of training on 84%-positive data: a model that
         returns high probabilities across the whole library, producing a ranking driven
         by noise. It did not happen. The score distributions are sharply peaked near
         zero with a thin tail, which is what a discriminating ranker looks like.`,
        ev(Object.keys(disc).map((t) => [t,
          `median <strong>${Number(disc[t].median).toExponential(1)}</strong> ·
           99th pct ${num(disc[t].p99, 3)} ·
           ${(disc[t].frac_above_half * 100).toFixed(2)}% above 0.5`]))),

      claim(3, "partial",
        "The shortlist is enriched in traditional pain plants",
        `Compounds from the hand-curated ethnopharmacology seed list make up
         ${((en.library_base_rate || 0) * 100).toFixed(1)}% of the library but
         <strong>${((en.shortlist_rate || 0) * 100).toFixed(0)}% of the shortlist</strong> —
         a ${en.fold}× enrichment. The base rate is identical in the filtered pool the
         shortlist was drawn from, so this is not an artefact of the safety filters.
         <br><br><strong>The caveat, and it is a real one:</strong> this is partly
         circular. Those plants were put on the seed list <em>because</em> they have
         documented pain use, and their compounds (capsaicinoids, cannabinoids,
         sanshools) are chemically close to the known ligands in ChEMBL. The models are
         rediscovering a relationship the seed list already encodes. Treat it as a
         consistency check, not as independent validation.`,
        ev([
          ["Library base rate", `${((en.library_base_rate || 0) * 100).toFixed(2)}%`],
          ["Eligible pool", `${int(en.eligible_pool)} compounds, ${((en.eligible_base_rate || 0) * 100).toFixed(2)}% seed`],
          ["Shortlist", `${((en.shortlist_rate || 0) * 100).toFixed(0)}% seed`],
          ["Enrichment", `<strong>${en.fold}×</strong>`],
        ])),

      claim(4, "supported",
        "The FAAH blind spot is a data-coverage problem, not a representation problem",
        `This is the strongest methodological result in the project, and the one most
         worth presenting. Six known plant FAAH inhibitors, none in ChEMBL's FAAH set.
         Each is scored by a model that never saw it. Adding just five curated literature
         actives moves the median prediction by four orders of magnitude.
         <br><br><strong>Why it matters:</strong> it cleanly separates two hypotheses
         that usually get confused — "the features cannot express this chemistry" versus
         "the training set has never seen this chemistry". It is the latter, and that is
         fixable with data rather than with a new architecture.`,
        ev([
          ["Median p, ChEMBL-only model", `<strong>${Number(rescue.baseline_median_p).toExponential(1)}</strong>`],
          ["Median p, +5 literature actives", `<strong>${num(rescue.loo_median_p, 3)}</strong>`],
          ["Overlap with ChEMBL FAAH set", `${rescue.overlap_with_chembl} compounds`],
        ])),

      claim(5, "not",
        "The DILI model cannot be used as a safety filter",
        `It misses the project's own index compound. Paracetamol and isoniazid — two of
         the most famous hepatotoxins in medicine — are absent from the training
         benchmark entirely, and the model calls paracetamol low-risk. Worse, it is
         <em>confidently wrong in the other direction</em> on natural products: kaempferol
         scores ${num(kaem.p_dili, 2)} at an applicability domain of only
         ${num(kaem.ad, 2)}, meaning that score is extrapolation dressed up as a
         prediction.
         <br><br><strong>The project handles this correctly</strong> by weighting it at
         0.30 against the structural alerts, which do flag paracetamol at maximum
         severity. Keeping a model you have shown to be unreliable, at a reduced weight,
         with the reason documented, is a better answer than deleting it or trusting it.`,
        ev([
          ["paracetamol", `p = ${num(para.p_dili, 3)} <span class="dim">(absent from benchmark)</span>`],
          ["kaempferol", `p = ${num(kaem.p_dili, 3)} at AD ${num(kaem.ad, 2)} <span class="dim">(extrapolation)</span>`],
          ["Weight in composite", "0.30, versus an exponential penalty on alerts"],
        ])),

      claim(6, "partial",
        "One fifth of the library is outside the models' competence",
        `The applicability domain is doing real work rather than decorating the output.
         A median max-similarity of ${num(d.median_ad, 2)} means the typical plant
         compound is only loosely related to anything in the training set — which is the
         expected consequence of training on synthetic drug chemistry and deploying on
         plant chemistry.
         <br><br><strong>Read this as a scope limit:</strong> for ${int(d.extrapolating)}
         compounds the model is guessing, and the shortlist correctly excludes them.`,
        ev([
          ["In domain", `${int(d.in_domain)} (${((d.in_domain_rate || 0) * 100).toFixed(1)}%)`],
          ["Extrapolating", `${int(d.extrapolating)} (${(100 - (d.in_domain_rate || 0) * 100).toFixed(1)}%)`],
          ["Median max-similarity", num(d.median_ad, 3)],
        ])),
    ].join("");
  }

  function renderAttacks(v, m) {
    const metrics = (m.metrics || []).filter((r) => r.split === "scaffold");
    const weakest = metrics.length
      ? metrics.reduce((a, b) => (a.pr_auc < b.pr_auc ? a : b)) : null;
    const en = v.enrichment || {};

    const items = [
      {
        q: "Your scaffold-split PR-AUC is 0.88–0.98. That is too good. What is inflating it?",
        a: `<strong>The decoys are easy negatives.</strong> They are random drug-like
            molecules, so the model is mostly separating "real ligand" from "arbitrary
            molecule" — a much easier task than ranking among plausible candidates, which
            is what the screen actually faces. The honest answer is that these numbers are
            an upper bound on screening performance.
            ${weakest ? `The most realistic row in the table is
            <code>${esc(weakest.target)}/${esc(weakest.scope)}/${esc(weakest.label)}</code>
            at PR-AUC ${num(weakest.pr_auc)}, because the 10 µM threshold forces the
            hardest discrimination.` : ""}`,
      },
      {
        q: "Isn't the seed-list enrichment circular?",
        a: `<strong>Partly, yes — concede it.</strong> Those plants are on the list because
            they have documented pain use, and their chemistry resembles the known ligands
            in training. The ${en.fold}× enrichment shows the pipeline is internally
            consistent, not that it discovered anything. The genuinely interesting rows are
            the ${Math.round((1 - (en.shortlist_rate || 0)) * 50)} background-tier compounds
            that reached the shortlist without any traditional-use prior.`,
      },
      {
        q: "You label 60,000 molecules inactive without measuring them. Isn't that made up?",
        a: `It is an assumption, stated as one. Real screening hit rates are well under 1%,
            so a random drug-like molecule is overwhelmingly likely to be inactive against
            any given target. They are labelled <code>presumed_inactive</code> rather than
            <code>inactive</code>, so every result can be re-derived without them. The
            alternative — training on 84%-positive data — is demonstrably worse.`,
      },
      {
        q: "Have you tested any of this in a lab?",
        a: `<strong>No, and nothing here claims otherwise.</strong> There is no prospective
            validation. Every number is a computational prediction over compounds that have
            not been assayed by this project. The deliverable is a ranked, filtered,
            provenance-tracked hypothesis list for someone with a lab — not a result.`,
      },
      {
        q: "Your positive controls were in the training data. What does recovering them prove?",
        a: `That the pipeline is wired up correctly — featurisation, model loading, scoring
            and ranking all work end to end. It is a smoke test, not a performance claim,
            and it is reported separately from the screen for exactly that reason. The
            generalisation claim rests on the scaffold split and the FAAH rescue, not on
            these.`,
      },
    ];

    $("#attacks").innerHTML = items.map((it, n) => `<div class="callout" style="margin-bottom:12px">
      <h4>${n + 1}. ${esc(it.q)}</h4>
      <p class="card-note">${it.a}</p>
    </div>`).join("");
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
    go(["overview", "science", "data", "models", "validation", "safety", "screen", "decoder"]
       .indexOf(initial) >= 0 ? initial : "overview");
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
