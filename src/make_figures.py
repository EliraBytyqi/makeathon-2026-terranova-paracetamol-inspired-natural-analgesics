"""Presentation figures.

Three charts, each answering one question:

  fig1  Does the model actually generalise?       random vs scaffold PR-AUC
  fig2  How much reactive chemistry is in plants? alert prevalence
  fig3  Where does the library sit vs the models? applicability domain

Colors come from the reference categorical palette, assigned in fixed slot
order. Every chart is a magnitude comparison, so the form is a bar chart with
the category on the y-axis where labels are long.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
FIGS = RESULTS / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

# reference palette, fixed slot order
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
CRITICAL = "#d03b3b"

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK2,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def fig_generalisation():
    p = RESULTS / "model_metrics.csv"
    if not p.exists():
        return
    m = pd.read_csv(p)
    m = m[m["label"] != "label_10um"] if "label_10um" in set(m["label"]) else m
    piv = m.pivot_table(index=["target", "scope", "model"], columns="split",
                        values="pr_auc").reset_index()
    if "scaffold" not in piv.columns or "random" not in piv.columns:
        return
    piv = piv[piv["scope"] == "human"].sort_values("scaffold")
    if piv.empty:
        return

    labels = [f"{r.target} · {r.model.replace('_', ' ')}" for r in piv.itertuples()]
    y = np.arange(len(piv))
    h = 0.36

    fig, ax = plt.subplots(figsize=(8.4, 0.72 * len(piv) + 2.0))
    ax.barh(y + h / 2, piv["random"], height=h, color=ORANGE,
            label="Random split (flattering)", zorder=3)
    ax.barh(y - h / 2, piv["scaffold"], height=h, color=BLUE,
            label="Scaffold split (honest)", zorder=3)

    for yi, (rv, sv) in enumerate(zip(piv["random"], piv["scaffold"])):
        ax.text(rv + 0.012, yi + h / 2, f"{rv:.2f}", va="center",
                fontsize=9, color=INK2)
        ax.text(sv + 0.012, yi - h / 2, f"{sv:.2f}", va="center",
                fontsize=9, color=INK2, fontweight="bold")

    ax.set_yticks(y, labels, fontsize=9, color=INK2)
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("PR-AUC")
    ax.set_title("Generalisation gap: memorising series vs predicting new chemotypes",
                 fontsize=12, color=INK, pad=14, loc="left")
    ax.xaxis.grid(True, color=GRID, lw=1, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.text(0.01, 0.015,
             "Natural products share almost no scaffolds with ChEMBL, so the "
             "scaffold-split bar is what predicts screening performance.",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIGS / "fig1_generalisation.png", dpi=200)
    plt.close(fig)
    print("wrote fig1_generalisation.png")


def fig_alerts():
    p = PROC / "np_library.csv"
    if not p.exists():
        return
    df = pd.read_csv(p, low_memory=False)
    ex = df["reactive_alerts"].fillna("").str.split("|").explode()
    ex = ex[ex != ""]
    counts = ex.value_counts().head(10)[::-1]
    if counts.empty:
        return
    pct = counts / len(df) * 100

    fig, ax = plt.subplots(figsize=(8.4, 0.5 * len(counts) + 2.2))
    ax.barh(range(len(counts)), pct.values, color=AQUA, height=0.62, zorder=3)
    for i, v in enumerate(pct.values):
        ax.text(v + max(pct.values) * 0.012, i, f"{v:.1f}%", va="center",
                fontsize=9, color=INK2)
    ax.set_yticks(range(len(counts)),
                  [s.replace("_", " ") for s in counts.index], fontsize=9)
    ax.set_xlabel("% of the plant compound library")
    ax.set_title("Reactive-metabolite alerts across plant chemical space",
                 fontsize=12, color=INK, pad=14, loc="left")
    ax.xaxis.grid(True, color=GRID, lw=1, zorder=0)
    ax.set_axisbelow(True)
    fig.text(0.01, 0.015,
             "Plants make electrophiles deliberately. Flags are reported, not "
             "tuned away - paracetamol itself trips two severity-3 alerts.",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(FIGS / "fig2_alerts.png", dpi=200)
    plt.close(fig)
    print("wrote fig2_alerts.png")


def fig_domain():
    p = RESULTS / "screening_all.csv"
    if not p.exists():
        return
    df = pd.read_csv(p, low_memory=False)
    if "ad_best" not in df.columns:
        return
    vals = df["ad_best"].dropna()

    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.hist(vals, bins=60, color=BLUE, zorder=3)
    ax.axvline(0.25, color=CRITICAL, lw=2, ls="--", zorder=4)
    ax.text(0.26, ax.get_ylim()[1] * 0.92,
            "applicability-domain floor (0.25)", color=CRITICAL, fontsize=9)
    ax.set_xlabel("Max Tanimoto similarity to the model's training set")
    ax.set_ylabel("plant compounds")
    ax.set_title("Most of the plant library is genuine extrapolation",
                 fontsize=12, color=INK, pad=14, loc="left")
    ax.yaxis.grid(True, color=GRID, lw=1, zorder=0)
    ax.set_axisbelow(True)
    frac = float((vals < 0.25).mean())
    fig.text(0.01, 0.015,
             f"{frac:.0%} of compounds fall below the floor. Predictions there "
             "are flagged low-confidence rather than silently trusted.",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(FIGS / "fig3_applicability_domain.png", dpi=200)
    plt.close(fig)
    print("wrote fig3_applicability_domain.png")


def main():
    fig_generalisation()
    fig_alerts()
    fig_domain()
    print(f"figures in {FIGS}")


if __name__ == "__main__":
    main()
