/* Hand-rolled SVG charts.
 *
 * No chart library and no CDN: the whole dashboard has to render on a venue
 * network that may not resolve anything, so every mark here is drawn from the
 * data with plain SVG.
 *
 * House rules enforced throughout (data-viz method):
 *   - thin marks, 4px rounded data-ends anchored to the baseline
 *   - a 2px surface gap between adjacent fills, never a border
 *   - hairline, solid gridlines one shade off the surface (never dashed)
 *   - a hover/focus tooltip on every chart, with hit targets >= 24px
 *   - a legend whenever >= 2 series are present; direct labels used selectively
 *   - every chart ships a table-view twin, so no value is colour-gated
 *   - colour follows the entity (the target), never the row's rank
 */
(function (global) {
  "use strict";

  const NS = "http://www.w3.org/2000/svg";

  /* ------------------------------------------------------------ tooltip */
  let tipEl = null;
  function tip() {
    if (!tipEl) {
      tipEl = document.createElement("div");
      tipEl.className = "tooltip";
      tipEl.setAttribute("role", "status");
      document.body.appendChild(tipEl);
    }
    return tipEl;
  }
  function showTip(html, evt) {
    const t = tip();
    t.innerHTML = html;
    t.classList.add("show");
    moveTip(evt);
  }
  function moveTip(evt) {
    const t = tip();
    const pad = 14;
    const r = t.getBoundingClientRect();
    let x = evt.clientX + pad;
    let y = evt.clientY + pad;
    if (x + r.width > window.innerWidth - 8) x = evt.clientX - r.width - pad;
    if (y + r.height > window.innerHeight - 8) y = evt.clientY - r.height - pad;
    t.style.left = Math.max(8, x) + "px";
    t.style.top = Math.max(8, y) + "px";
  }
  function hideTip() { if (tipEl) tipEl.classList.remove("show"); }

  function tipRows(title, rows) {
    const body = rows
      .map((r) => `<div class="tooltip-row"><span class="k">${esc(r[0])}</span><span class="v">${esc(r[1])}</span></div>`)
      .join("");
    return `<div class="tooltip-title">${esc(title)}</div>${body}`;
  }

  /* ------------------------------------------------------------- helpers */
  function esc(s) {
    return String(s === null || s === undefined ? "" : s).replace(
      /[&<>"']/g,
      (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c])
    );
  }

  function el(name, attrs, parent) {
    const node = document.createElementNS(NS, name);
    for (const k in attrs) {
      if (attrs[k] !== null && attrs[k] !== undefined) node.setAttribute(k, attrs[k]);
    }
    if (parent) parent.appendChild(node);
    return node;
  }

  function fmt(n, digits) {
    if (n === null || n === undefined || Number.isNaN(n)) return "—";
    if (digits !== undefined) return Number(n).toFixed(digits);
    const a = Math.abs(n);
    if (a >= 1000) return Math.round(n).toLocaleString();
    if (a >= 10) return Number(n).toFixed(1);
    if (a >= 1) return Number(n).toFixed(2);
    return Number(n).toPrecision(2);
  }

  function cssVar(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  }

  /** Rounded-end bar: square at the baseline, 4px radius at the data end. */
  function barPath(x, y, w, h, r, horizontal) {
    const rad = Math.max(0, Math.min(r, horizontal ? w : h));
    if (rad <= 0.5) return `M${x},${y}h${w}v${h}h${-w}Z`;
    if (horizontal) {
      return `M${x},${y}h${w - rad}a${rad},${rad} 0 0 1 ${rad},${rad}v${h - 2 * rad}a${rad},${rad} 0 0 1 ${-rad},${rad}h${-(w - rad)}Z`;
    }
    return `M${x},${y + rad}a${rad},${rad} 0 0 1 ${rad},${-rad}h${w - 2 * rad}a${rad},${rad} 0 0 1 ${rad},${rad}v${h - rad}h${-w}Z`;
  }

  /* --------------------------------------------------- chart + table shell
   * Every chart is wrapped so it can flip to an accessible table. The table is
   * the WCAG-clean twin: it carries the same values without relying on colour,
   * which is what lets the light-mode aqua slot (sub-3:1 against the surface)
   * be used at all.
   */
  function shell(container, opts, draw, tableFn) {
    container.innerHTML = "";
    const wrap = document.createElement("div");

    if (tableFn) {
      const bar = document.createElement("div");
      bar.className = "chart-toolbar";
      bar.style.cssText = "justify-content:flex-end;margin-bottom:6px";
      const btn = document.createElement("button");
      btn.className = "btn btn-sm";
      btn.type = "button";
      btn.textContent = "Table view";
      btn.setAttribute("aria-pressed", "false");
      let showing = false;
      btn.addEventListener("click", () => {
        showing = !showing;
        btn.setAttribute("aria-pressed", String(showing));
        btn.textContent = showing ? "Chart view" : "Table view";
        plot.hidden = showing;
        tableHost.hidden = !showing;
        if (showing && !tableHost.dataset.built) {
          tableHost.innerHTML = tableFn();
          tableHost.dataset.built = "1";
        }
      });
      bar.appendChild(btn);
      wrap.appendChild(bar);
    }

    const plot = document.createElement("div");
    const tableHost = document.createElement("div");
    tableHost.hidden = true;
    tableHost.className = "table-wrap";
    wrap.appendChild(plot);
    wrap.appendChild(tableHost);
    container.appendChild(wrap);

    draw(plot);
    if (opts && opts.legend) plot.appendChild(legendEl(opts.legend));
    return wrap;
  }

  function legendEl(items) {
    const box = document.createElement("div");
    box.className = "legend";
    items.forEach((it) => {
      const span = document.createElement("span");
      span.className = "legend-item";
      span.innerHTML = `<span class="legend-swatch" style="background:${esc(it.color)}"></span>${esc(it.label)}`;
      box.appendChild(span);
    });
    return box;
  }

  function htmlTable(head, rows, caption) {
    const th = head.map((h) => `<th scope="col">${esc(h)}</th>`).join("");
    const body = rows
      .map(
        (r) =>
          "<tr>" +
          r
            .map((c, i) => (i === 0 ? `<th scope="row" style="font-weight:500">${esc(c)}</th>` : `<td class="num">${esc(c)}</td>`))
            .join("") +
          "</tr>"
      )
      .join("");
    return `<table>${caption ? `<caption>${esc(caption)}</caption>` : ""}<thead><tr>${th}</tr></thead><tbody>${body}</tbody></table>`;
  }

  /** Attach hover + keyboard focus to a mark, with an enlarged hit area. */
  function interactive(svg, hit, html) {
    hit.setAttribute("tabindex", "0");
    hit.setAttribute("role", "img");
    hit.setAttribute("aria-label", html.replace(/<[^>]+>/g, " ").replace(/\s+/g, " ").trim());
    hit.addEventListener("mouseenter", (e) => showTip(html, e));
    hit.addEventListener("mousemove", moveTip);
    hit.addEventListener("mouseleave", hideTip);
    hit.addEventListener("focus", () => {
      const r = hit.getBoundingClientRect();
      showTip(html, { clientX: r.left + r.width / 2, clientY: r.top });
    });
    hit.addEventListener("blur", hideTip);
  }

  /* ------------------------------------------------------- horizontal bar
   * Used for ordered magnitudes: the curation funnel and simple rankings.
   * `ordinal: true` steps the blue ramp with position (a legitimate ordered
   * ramp); otherwise every bar is slot 1, because a value-ramp on nominal
   * categories would double-encode length as hue.
   */
  function hbar(container, cfg) {
    const data = cfg.data || [];
    const labelW = cfg.labelWidth || 190;
    const rowH = cfg.rowHeight || 30;
    const gap = 8;
    const padR = cfg.valueLabels === false ? 16 : 68;
    const h = data.length * rowH + 14;
    const max = cfg.max || Math.max(...data.map((d) => d.value), 1);

    shell(
      container,
      { legend: cfg.legend },
      (host) => {
        const w = Math.max(host.clientWidth || 640, 320);
        const plotW = w - labelW - padR;
        const svg = el("svg", {
          class: "chart",
          viewBox: `0 0 ${w} ${h}`,
          height: h,
          role: "img",
          "aria-label": cfg.ariaLabel || cfg.title || "bar chart",
        }, host);

        // hairline grid at quarter steps, drawn under the marks
        for (let i = 0; i <= 4; i++) {
          const x = labelW + (plotW * i) / 4;
          el("line", { x1: x, y1: 0, x2: x, y2: h - 14, class: "gridline" }, svg);
        }

        data.forEach((d, i) => {
          const y = i * rowH + 4;
          const barH = rowH - gap;
          const bw = Math.max(2, (d.value / max) * plotW);
          const color = d.color || (cfg.ordinal
            ? cssVar(`--ramp-${Math.min(10, Math.round(1 + (i / Math.max(1, data.length - 1)) * 9))}`)
            : cssVar("--series-1"));

          el("text", {
            x: labelW - 10, y: y + barH / 2 + 4, "text-anchor": "end",
            class: d.emphasis ? "value-label-strong" : "axis-label",
          }, svg).textContent = d.label;

          el("path", {
            d: barPath(labelW, y, bw, barH, 4, true),
            fill: color,
            opacity: d.muted ? 0.45 : 1,
          }, svg);

          if (cfg.valueLabels !== false) {
            el("text", {
              x: labelW + bw + 8, y: y + barH / 2 + 4,
              class: d.emphasis ? "value-label-strong" : "value-label",
            }, svg).textContent = d.display || fmt(d.value);
          }

          // hit area spans the full row so the target clears 24px
          const hit = el("rect", {
            x: 0, y: y - gap / 2, width: w, height: rowH, class: "mark-hit",
          }, svg);
          interactive(svg, hit, tipRows(d.label, d.tip || [["value", fmt(d.value)]]));
        });
      },
      () =>
        htmlTable(
          [cfg.labelHeader || "Step", cfg.valueHeader || "Value"],
          data.map((d) => [d.label, d.display || fmt(d.value)]),
          cfg.title
        )
    );
  }

  /* ------------------------------------------------------- stacked bars
   * Part-to-whole across a handful of categories. A 2px surface gap separates
   * segments instead of a stroke.
   */
  function stacked(container, cfg) {
    const data = cfg.data || [];
    const series = cfg.series || [];
    const labelW = cfg.labelWidth || 120;
    const rowH = 36;
    const h = data.length * rowH + 16;

    shell(
      container,
      {
        legend: series.map((s) => ({ label: s.label, color: cssVar(s.colorVar) })),
      },
      (host) => {
        const w = Math.max(host.clientWidth || 640, 320);
        const plotW = w - labelW - 70;
        const svg = el("svg", {
          class: "chart", viewBox: `0 0 ${w} ${h}`, height: h,
          role: "img", "aria-label": cfg.ariaLabel || cfg.title || "stacked bar chart",
        }, host);
        const max = Math.max(...data.map((d) => series.reduce((a, s) => a + (d[s.key] || 0), 0)), 1);

        data.forEach((d, i) => {
          const y = i * rowH + 5;
          const barH = rowH - 14;
          let x = labelW;
          const total = series.reduce((a, s) => a + (d[s.key] || 0), 0);

          el("text", {
            x: labelW - 10, y: y + barH / 2 + 4, "text-anchor": "end", class: "axis-label",
          }, svg).textContent = d.label;

          series.forEach((s, si) => {
            const v = d[s.key] || 0;
            if (v <= 0) return;
            const segW = (v / max) * plotW;
            const drawW = Math.max(2, segW - (si < series.length - 1 ? 2 : 0)); // 2px surface gap
            el("path", {
              d: barPath(x, y, drawW, barH, si === series.length - 1 ? 4 : 0, true),
              fill: cssVar(s.colorVar),
            }, svg);
            const hit = el("rect", { x, y: y - 5, width: segW, height: rowH, class: "mark-hit" }, svg);
            interactive(svg, hit, tipRows(`${d.label} — ${s.label}`, [
              ["count", fmt(v)],
              ["share", total ? ((v / total) * 100).toFixed(1) + "%" : "—"],
            ]));
            x += segW;
          });

          el("text", {
            x: labelW + (total / max) * plotW + 8, y: y + barH / 2 + 4, class: "value-label",
          }, svg).textContent = fmt(total);
        });
      },
      () =>
        htmlTable(
          [cfg.labelHeader || "Group", ...series.map((s) => s.label), "Total"],
          data.map((d) => [
            d.label,
            ...series.map((s) => fmt(d[s.key] || 0)),
            fmt(series.reduce((a, s) => a + (d[s.key] || 0), 0)),
          ]),
          cfg.title
        )
    );
  }

  /* ---------------------------------------------------------- dumbbell
   * Two measurements of the same thing (random vs scaffold split; baseline vs
   * leave-one-out). A connector plus two markers shows the gap directly, which
   * a grouped bar pair does not.
   */
  function dumbbell(container, cfg) {
    const data = cfg.data || [];
    const labelW = cfg.labelWidth || 210;
    const rowH = 30;
    const h = data.length * rowH + 34;
    const lo = cfg.min !== undefined ? cfg.min : 0;
    const hi = cfg.max !== undefined ? cfg.max : 1;
    const scale = cfg.log
      ? (v) => (Math.log10(Math.max(v, lo)) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo))
      : (v) => (v - lo) / (hi - lo);

    shell(
      container,
      {
        legend: [
          { label: cfg.aLabel, color: cssVar(cfg.aColorVar || "--series-neutral") },
          { label: cfg.bLabel, color: cssVar(cfg.bColorVar || "--series-1") },
        ],
      },
      (host) => {
        const w = Math.max(host.clientWidth || 640, 320);
        const plotW = w - labelW - 56;
        const svg = el("svg", {
          class: "chart", viewBox: `0 0 ${w} ${h}`, height: h,
          role: "img", "aria-label": cfg.ariaLabel || cfg.title || "dumbbell chart",
        }, host);

        const ticks = cfg.ticks || [lo, lo + (hi - lo) / 2, hi];
        ticks.forEach((t) => {
          const x = labelW + scale(t) * plotW;
          el("line", { x1: x, y1: 0, x2: x, y2: h - 22, class: "gridline" }, svg);
          el("text", { x, y: h - 7, "text-anchor": "middle", class: "axis-label" }, svg).textContent =
            cfg.tickFmt ? cfg.tickFmt(t) : fmt(t, 2);
        });

        const cA = cssVar(cfg.aColorVar || "--series-neutral");
        const cB = cssVar(cfg.bColorVar || "--series-1");

        data.forEach((d, i) => {
          const y = i * rowH + 15;
          const xa = labelW + scale(d.a) * plotW;
          const xb = labelW + scale(d.b) * plotW;

          el("text", {
            x: labelW - 10, y: y + 4, "text-anchor": "end", class: "axis-label",
          }, svg).textContent = d.label;

          el("line", {
            x1: Math.min(xa, xb), y1: y, x2: Math.max(xa, xb), y2: y,
            stroke: cssVar("--baseline"), "stroke-width": 2, "stroke-linecap": "round",
          }, svg);

          // 2px surface ring keeps the markers readable where they overlap
          [[xa, cA], [xb, cB]].forEach(([x, c]) => {
            el("circle", { cx: x, cy: y, r: 5.5, fill: c, stroke: cssVar("--surface"), "stroke-width": 2 }, svg);
          });

          if (cfg.showDelta !== false) {
            const delta = d.b - d.a;
            el("text", {
              x: w - 8, y: y + 4, "text-anchor": "end",
              class: "value-label",
              fill: delta < 0 ? cssVar("--text-secondary") : cssVar("--good-text"),
            }, svg).textContent = (delta >= 0 ? "+" : "") + fmt(delta, cfg.deltaDigits || 3);
          }

          const hit = el("rect", { x: 0, y: y - rowH / 2, width: w, height: rowH, class: "mark-hit" }, svg);
          interactive(svg, hit, tipRows(d.label, [
            [cfg.aLabel, cfg.valueFmt ? cfg.valueFmt(d.a) : fmt(d.a, 3)],
            [cfg.bLabel, cfg.valueFmt ? cfg.valueFmt(d.b) : fmt(d.b, 3)],
            ["difference", (d.b - d.a >= 0 ? "+" : "") + fmt(d.b - d.a, 3)],
          ].concat(d.extra || [])));
        });
      },
      () =>
        htmlTable(
          [cfg.labelHeader || "Model", cfg.aLabel, cfg.bLabel, "Difference"],
          data.map((d) => [
            d.label,
            cfg.valueFmt ? cfg.valueFmt(d.a) : fmt(d.a, 3),
            cfg.valueFmt ? cfg.valueFmt(d.b) : fmt(d.b, 3),
            (d.b - d.a >= 0 ? "+" : "") + fmt(d.b - d.a, 3),
          ]),
          cfg.title
        )
    );
  }

  /* --------------------------------------------------------- histogram */
  function histogram(container, cfg) {
    const bins = cfg.bins || [];
    const h = cfg.height || 190;

    shell(
      container,
      null,
      (host) => {
        const w = Math.max(host.clientWidth || 640, 320);
        const padL = 52, padR = 12, padB = 30, padT = 8;
        const plotW = w - padL - padR;
        const plotH = h - padB - padT;
        const svg = el("svg", {
          class: "chart", viewBox: `0 0 ${w} ${h}`, height: h,
          role: "img", "aria-label": cfg.ariaLabel || cfg.title || "histogram",
        }, host);
        const max = Math.max(...bins.map((b) => b.count), 1);
        const bw = plotW / bins.length;

        [0, 0.5, 1].forEach((f) => {
          const y = padT + plotH * (1 - f);
          el("line", { x1: padL, y1: y, x2: w - padR, y2: y, class: "gridline" }, svg);
          el("text", { x: padL - 8, y: y + 4, "text-anchor": "end", class: "axis-label" }, svg).textContent =
            fmt(max * f);
        });

        bins.forEach((b, i) => {
          const bh = Math.max(1, (b.count / max) * plotH);
          const x = padL + i * bw;
          el("path", {
            d: barPath(x, padT + plotH - bh, Math.max(1, bw - 2), bh, 3, false),
            fill: cssVar("--series-1"),
          }, svg);
          const hit = el("rect", { x, y: padT, width: bw, height: plotH, class: "mark-hit" }, svg);
          interactive(svg, hit, tipRows(`${fmt(b.lo, 3)} – ${fmt(b.hi, 3)}`, [["compounds", fmt(b.count)]]));
        });

        el("line", { x1: padL, y1: padT + plotH, x2: w - padR, y2: padT + plotH, class: "baseline-rule" }, svg);
        [0, Math.floor(bins.length / 2), bins.length - 1].forEach((i) => {
          if (!bins[i]) return;
          el("text", {
            x: padL + i * bw + bw / 2, y: h - 9, "text-anchor": "middle", class: "axis-label",
          }, svg).textContent = fmt(bins[i].lo, 2);
        });
        if (cfg.xLabel) {
          el("text", { x: padL + plotW / 2, y: h - 9, "text-anchor": "middle", class: "axis-label" }, svg)
            .textContent = "";
        }
      },
      () =>
        htmlTable(
          ["Range", "Compounds"],
          bins.map((b) => [`${fmt(b.lo, 3)} – ${fmt(b.hi, 3)}`, fmt(b.count)]),
          cfg.title
        )
    );
  }

  /* ----------------------------------------------------------- funnel
   * A staged drop-off. Ordered stages, so the ordinal blue ramp is correct
   * here (and validated for the ordinal case).
   */
  function funnel(container, cfg) {
    const stages = cfg.stages || [];
    const first = stages.length ? stages[0].value : 1;
    hbar(container, {
      title: cfg.title,
      ariaLabel: cfg.ariaLabel,
      ordinal: true,
      labelWidth: cfg.labelWidth || 230,
      labelHeader: cfg.labelHeader || "Stage",
      valueHeader: "Compounds",
      data: stages.map((s, i) => ({
        label: s.label,
        value: s.value,
        display: fmt(s.value),
        emphasis: i === stages.length - 1,
        tip: [
          ["count", fmt(s.value)],
          ["of start", first ? ((s.value / first) * 100).toFixed(1) + "%" : "—"],
          ...(i > 0
            ? [["dropped here", fmt(stages[i - 1].value - s.value)]]
            : []),
          ...(s.note ? [["note", s.note]] : []),
        ],
      })),
    });
  }

  global.Charts = { hbar, stacked, dumbbell, histogram, funnel, fmt, esc, cssVar, htmlTable, hideTip };
})(window);
