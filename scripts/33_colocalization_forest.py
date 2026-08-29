"""
Per-assay co-localization forest plot (companion to 32).

For every assay of every multi-assay protein we already computed, under the
permutation null of 32:
  exp_overlap_perm            expected # of this assay's high-disagreement sites
                              that also fall in the OTHER assays' flagged sets
  null_lo95 / null_hi95       2.5/97.5 percentiles of that permutation null
  obs_overlap_with_others     the actually observed overlap count

This makes one small-multiple panel per protein (mirroring the per-protein forest
plots elsewhere): each row is an assay; a dot marks the EXPECTED overlap with a
whisker for the 95% null band, and a separate diamond marks the OBSERVED overlap.
Observed points significant after BH-FDR across all assays (q_perm < 0.05) are
drawn filled/red, i.e. co-localizing more than chance; the rest are grey.

Input  -> results/high_disagreement_colocalization/colocalization_per_assay.csv
Output -> figures/fig_colocalization_forest.png
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import common as C

# "high" (default) or "low" -- picks the matching colocalization result dir from 32.
MODE = sys.argv[1] if len(sys.argv) > 1 else "high"
assert MODE in ("high", "low"), "usage: 33_colocalization_forest.py [high|low]"
LABEL = f"{MODE}-disagreement"

SIG = 0.05
EXP_COLOR = "#4C72B0"     # expected dot + null band
OBS_SIG = "#C44E52"       # observed, outside null band (co-localizes > chance)
OBS_NS = "#8C8C8C"        # observed, within null band


def short_assay(dms, uni):
    """Drop the protein prefix so rows read e.g. 'Gersing_2023_abundance'."""
    s = dms[len(uni) + 1:] if dms.startswith(uni + "_") else dms
    return s if len(s) <= 34 else s[:31] + "..."


def main():
    indir = C.rp(f"{MODE}_disagreement_colocalization")   # folder name -> its dir
    per = pd.read_csv(os.path.join(indir, "colocalization_per_assay.csv"))
    # order proteins by their co-localized fraction (same as the summary table)
    summ = pd.read_csv(os.path.join(indir, "colocalization_summary.csv"))
    order = summ.sort_values("frac_colocalized", ascending=False)["UniProt_ID"].tolist()
    order = [u for u in order if u in set(per.UniProt_ID)]

    ncol = 4
    nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 2.0 * nrow + 1.0))
    axes = np.atleast_1d(axes).ravel()

    for ax, uni in zip(axes, order):
        g = per[per.UniProt_ID == uni].reset_index(drop=True)
        y = np.arange(len(g))[::-1]            # first assay at the top

        # 95% null band as a horizontal line with end caps, expected mean as a dot
        # (drawn separately so a skewed discrete null whose mean sits just outside
        # the 2.5-97.5 band can't produce a negative whisker length)
        ax.hlines(y, g.null_lo95, g.null_hi95, color=EXP_COLOR, lw=1.5, zorder=2)
        ax.scatter(g.null_lo95, y, marker="|", s=55, color=EXP_COLOR, zorder=2)
        ax.scatter(g.null_hi95, y, marker="|", s=55, color=EXP_COLOR, zorder=2)
        ax.scatter(g.exp_overlap_perm, y, marker="o", s=28, color=EXP_COLOR,
                   zorder=3)
        # observed diamond, coloured by BH-FDR significance across all assays
        obs_col = np.where(g.q_perm < SIG, OBS_SIG, OBS_NS)
        ax.scatter(g.obs_overlap_with_others, y, marker="D", s=42,
                   color=obs_col, edgecolor="black", linewidth=0.5, zorder=3)

        ax.set_yticks(y)
        ax.set_yticklabels([short_assay(d, uni) for d in g.DMS_id], fontsize=7)
        ax.set_ylim(-0.6, len(g) - 0.4)
        cat = g.Neff_L_category.iloc[0]
        ax.set_title(f"{uni}  (Neff/L {cat})", fontsize=9)
        ax.tick_params(axis="x", labelsize=7)
        ax.grid(axis="x", ls=":", lw=0.5, alpha=0.6)
        left = min(0, g.null_lo95.min())
        right = max(g.obs_overlap_with_others.max(), g.null_hi95.max())
        pad = 0.10 * (right - left) + 0.5
        ax.set_xlim(left - pad, right + pad)

    for ax in axes[len(order):]:
        ax.set_visible(False)

    handles = [
        plt.Line2D([0], [0], marker="o", color=EXP_COLOR, ls="none", ms=6,
                   label="Expected overlap (dot) + 95% null band (whisker)"),
        plt.Line2D([0], [0], marker="D", color=OBS_SIG, mec="black", ls="none",
                   ms=7, label="Observed overlap — significant (BH q<0.05)"),
        plt.Line2D([0], [0], marker="D", color=OBS_NS, mec="black", ls="none",
                   ms=7, label="Observed overlap — n.s. (BH q≥0.05)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"{LABEL.capitalize()} site co-localization per assay: observed vs "
                 "permutation-null expectation", fontsize=12, y=0.995)
    fig.supxlabel(f"number of {LABEL} sites shared with the protein's other assays",
                  fontsize=9, y=0.035)
    fig.tight_layout(rect=[0, 0.05, 1, 0.98])

    suffix = "" if MODE == "high" else f"_{MODE}"
    out = os.path.join(C.FIGURES, f"fig_colocalization_forest{suffix}.png")
    os.makedirs(C.FIGURES, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"-> wrote {out}  ({len(order)} proteins, {len(per)} assays)")


if __name__ == "__main__":
    main()
