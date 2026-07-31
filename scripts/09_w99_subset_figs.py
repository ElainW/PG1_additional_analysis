"""
W99 deep-dive, 23-model-subset version (same models/colors as Fig5_subset).

Left panel : heatmap of within-site normalized rank (0 = model thinks least fit,
             1 = most fit) for the 19 substitutions at W99. Rows = the 23 subset
             models grouped by model_type, plus the DMS experiment on top;
             columns ordered by the experiment (left = most fit).
Right panel: per-model statistic — Spearman correlation between the model's
             within-site ranking and the experiment's ranking at this site
             (the same 19 substitutions), as a horizontal bar colored by
             model_type. This quantifies how well each model recovers the local
             fitness ordering the assay actually measured.

Colors reuse 04_subset_figs.PALETTE (model_type tab10) and the RdYlBu_r heatmap.

Outputs per assay:
  results/site_W99_<assay>_subset_rank_matrix.csv
  figures/fig_site_W99_<assay>_subset.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import common as C
from importlib import import_module

sub = import_module("04_subset_figs")

SITE = 99
ASSAYS = ["HXK4_HUMAN_Gersing_2022_activity", "HXK4_HUMAN_Gersing_2023_abundance"]


def build(dms_id, site=SITE, results_dir=C.RESULTS):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
    present_cols = [c for c in sub.COLS if c in df.columns]
    rows = df[df["mutant"].str.match(rf"^[A-Z]{site}[A-Z]$")].copy()
    rows["alt"] = rows["mutant"].str[-1]
    rows = rows.set_index("alt")
    wt = rows["mutant"].iloc[0][0]

    cols = present_cols + ["DMS_score"]
    ranks = rows[cols].rank()
    norm = (ranks - 1) / (len(rows) - 1)
    dms_norm = norm["DMS_score"]

    # per-model Spearman vs experiment at this site
    spear = {c: spearmanr(rows[c], rows["DMS_score"]).statistic for c in present_cols}

    # row order: DMS, then by model_type, within type best-aligned first
    model_order = []
    for t in sub.MTYPES:
        members = [c for c in present_cols if sub.TYPE[sub.COL2DISP[c]] == t]
        members = sorted(members, key=lambda c: -spear[c])
        model_order += members
    col_order = dms_norm.sort_values(ascending=False).index.tolist()

    M = norm[["DMS_score"] + model_order].T[col_order]
    disp_index = ["DMS (experiment)"] + [sub.COL2DISP[c] for c in model_order]
    M.index = disp_index
    M.to_csv(os.path.join(results_dir, f"site_{wt}{site}_{dms_id}_subset_rank_matrix.csv"))

    stat = pd.Series({sub.COL2DISP[c]: spear[c] for c in model_order})
    row_type = {"DMS (experiment)": "DMS"}
    row_type.update({sub.COL2DISP[c]: sub.TYPE[sub.COL2DISP[c]] for c in model_order})
    return M, col_order, wt, stat, row_type


def figure(dms_id, M, col_order, wt, stat, row_type, site=SITE, outdir=C.FIGURES, region_note=None):
    palette = dict(sub.PALETTE)
    palette["DMS"] = (0.1, 0.1, 0.1)
    n = len(M)

    fig = plt.figure(figsize=(13.5, 10))
    # rows: [plots, spacer, colorbar];  cols: [heatmap, bar]
    gs = fig.add_gridspec(3, 2, width_ratios=[3.0, 1.15], height_ratios=[28, 1.5, 1.2],
                          wspace=0.06, hspace=0.05)
    axH = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1], sharey=axH)
    cax = fig.add_subplot(gs[2, 0])

    # ---- left: heatmap ----
    im = axH.imshow(M.astype(float).values, cmap="RdYlBu_r", vmin=0, vmax=1, aspect="auto")
    axH.set_xticks(np.arange(len(col_order)))
    axH.set_xticklabels([f"{wt}{site}{a}" for a in col_order], rotation=90, fontsize=8)
    axH.set_yticks(np.arange(n))
    axH.set_yticklabels(M.index, fontsize=8)
    for tick, r in zip(axH.get_yticklabels(), M.index):
        tick.set_color("black" if row_type[r] == "DMS" else palette[row_type[r]])
        if row_type[r] == "DMS":
            tick.set_fontweight("bold")
    axH.axhline(0.5, color="black", lw=2)          # separate DMS from models
    axH.set_xlabel("substitution (ordered by experiment, left = most fit)")
    title = f"{dms_id}\nWithin-site rank of the 19 substitutions at {wt}{site}"
    if region_note:
        title += f"\n[{region_note}]"
    axH.set_title(title, fontsize=11)

    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("within-site rank (0 = least fit, 1 = most fit)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # ---- right: Spearman vs experiment ----
    vals = [np.nan] + [stat.get(r, np.nan) for r in M.index[1:]]
    colors = [palette[row_type[r]] for r in M.index]
    ypos = np.arange(n)
    axB.barh(ypos, vals, color=colors, height=0.7)
    axB.axvline(0, color="gray", lw=0.8)
    axB.set_xlim(-1, 1)
    plt.setp(axB.get_yticklabels(), visible=False)
    for y, v in zip(ypos, vals):
        if not np.isnan(v):
            axB.text(v + (0.04 if v >= 0 else -0.04), y, f"{v:.2f}",
                     va="center", ha="left" if v >= 0 else "right", fontsize=6.5)
    axB.set_xlabel("Spearman vs experiment\n(within-site, 19 subs)", fontsize=9)
    axB.set_title("per-model rank corr.", fontsize=10, pad=8)
    axB.grid(axis="x", alpha=0.3)
    med = np.nanmedian(vals)
    axB.axvline(med, color="black", ls="--", lw=1)
    axB.text(med, n - 0.2, f"median {med:.2f}", fontsize=7, ha="center", va="top")

    # ---- legend: outside, to the right of the bar panel ----
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[t]) for t in sub.MTYPES]
    axB.legend(handles, sub.MTYPES, title="model_type", fontsize=7, title_fontsize=8,
               loc="upper left", bbox_to_anchor=(1.04, 1.0), frameon=False)

    out = os.path.join(outdir, f"fig_site_{wt}{site}_{dms_id}_subset.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


def main():
    for dms_id in ASSAYS:
        M, col_order, wt, stat, row_type = build(dms_id)
        figure(dms_id, M, col_order, wt, stat, row_type)
        print(f"  {dms_id}: median within-site Spearman vs DMS = {stat.median():.3f}; "
              f"best {stat.idxmax()} ({stat.max():.2f}), worst {stat.idxmin()} ({stat.min():.2f})")


if __name__ == "__main__":
    main()
