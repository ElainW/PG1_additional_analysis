"""
Deep dive on a single site (default W99 of HXK4_HUMAN), for both assays.

Shows how each model ranks the 19 substitutions at the site, as a heatmap of
within-site normalized rank (0 = model thinks least fit, 1 = most fit). Models
are grouped by model_type; the experiment (DMS) is shown as the top row, and
columns (the 19 alternative residues) are ordered by the DMS ranking so the
left side = experimentally most-tolerated substitutions.

Uses the 45-model family-best subset. Outputs per assay:
  results/site_<SITE>_<assay>_rank_matrix.csv
  figures/fig_site_<SITE>_<assay>_rankings.png
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import common as C
from importlib import import_module

sns.set_style("white")
fam = import_module("05_fig6_familybest")

SITE = 99
ASSAYS = ["HXK4_HUMAN_Gersing_2022_activity", "HXK4_HUMAN_Gersing_2023_abundance"]


def build(dms_id, chosen, type_of, mtypes):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
    present = [m for m in chosen if m in df.columns]
    sub = df[df["mutant"].str.match(rf"^[A-Z]{SITE}[A-Z]$")].copy()
    sub["alt"] = sub["mutant"].str[-1]
    sub = sub.set_index("alt")
    wt = sub["mutant"].iloc[0][0]

    # within-site normalized rank (0..1) per model + DMS
    cols = present + ["DMS_score"]
    ranks = sub[cols].rank()
    norm = (ranks - 1) / (len(sub) - 1)

    # order models by model_type (then by overall agreement w/ DMS at this site)
    dms_norm = norm["DMS_score"]
    model_order = []
    for t in mtypes:
        members = [m for m in present if type_of[m] == t]
        members = sorted(members, key=lambda m: -np.corrcoef(norm[m], dms_norm)[0, 1])
        model_order += members
    # columns ordered by DMS rank (most fit -> least)
    col_order = dms_norm.sort_values(ascending=False).index.tolist()

    M = norm[["DMS_score"] + model_order].T            # rows = DMS + models
    M = M[col_order]
    M.index = ["DMS (experiment)"] + model_order
    M.to_csv(C.rp(f"site_{wt}{SITE}_{dms_id}_rank_matrix.csv"))
    row_type = {"DMS (experiment)": "DMS"}
    row_type.update({m: type_of[m] for m in model_order})
    return M, col_order, wt, row_type, mtypes


def figure(dms_id, M, col_order, wt, row_type, mtypes):
    palette = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
    palette["DMS"] = (0.1, 0.1, 0.1)
    row_colors = [palette[row_type[r]] for r in M.index]

    fig, ax = plt.subplots(figsize=(11, 13))
    sns.heatmap(M.astype(float), cmap="RdYlBu_r", vmin=0, vmax=1, ax=ax,
                cbar_kws={"label": "within-site rank  (0 = least fit, 1 = most fit)", "shrink": 0.4},
                linewidths=0.3, linecolor="white")
    ax.set_xticks(np.arange(len(col_order)) + 0.5)
    ax.set_xticklabels([f"{wt}{SITE}{a}" for a in col_order], fontsize=8, rotation=90)
    ax.set_yticklabels(M.index, fontsize=7, rotation=0)
    # color the y tick labels by model_type and draw a separator under DMS
    for tick, r in zip(ax.get_yticklabels(), M.index):
        tick.set_color(palette[row_type[r]] if row_type[r] != "DMS" else "black")
        if row_type[r] == "DMS":
            tick.set_fontweight("bold")
    ax.axhline(1, color="black", lw=2)
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[t]) for t in mtypes]
    ax.legend(handles, mtypes, title="model_type", fontsize=7, loc="upper left",
              bbox_to_anchor=(1.02, 1.0))
    ax.set_title(f"{dms_id}\nHow each model ranks the 19 substitutions at {wt}{SITE} "
                 f"(columns ordered by experiment, left = most fit)", fontsize=11)
    ax.set_xlabel("substitution"); ax.set_ylabel("")
    plt.tight_layout()
    out = os.path.join(C.FIGURES, f"fig_site_{wt}{SITE}_{dms_id}_rankings.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


def main():
    chosen, type_of, _ = fam.build_family_best()
    mtypes = [t for t in C.MODEL_TYPE_ORDER if t in set(type_of.values())]
    for dms_id in ASSAYS:
        M, col_order, wt, row_type, mtypes = build(dms_id, chosen, type_of, mtypes)
        figure(dms_id, M, col_order, wt, row_type, mtypes)
        # quick text: experiment's preferred vs avoided, and which models invert it
        dms_row = M.loc["DMS (experiment)"]
        print(f"\n==== {dms_id}: site {wt}{SITE} ====")
        print("Experiment — most fit subs:", ", ".join([f"{wt}{SITE}{a}" for a in dms_row.sort_values(ascending=False).head(4).index]))
        print("Experiment — least fit subs:", ", ".join([f"{wt}{SITE}{a}" for a in dms_row.sort_values().head(4).index]))
        corr_dms = M.drop(index="DMS (experiment)").apply(lambda r: np.corrcoef(r.values, dms_row.values)[0, 1], axis=1)
        print("Models MOST opposite to experiment at this site:")
        print(corr_dms.sort_values().head(5).round(2).to_string())
        print("Models MOST aligned with experiment at this site:")
        print(corr_dms.sort_values(ascending=False).head(5).round(2).to_string())


if __name__ == "__main__":
    main()
