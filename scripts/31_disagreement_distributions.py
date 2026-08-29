"""
Per-protein disagreement distribution plots: one figure per UniProt_ID, with the
within-site disagreement distribution of each of its assays overlaid (so the 53
assays across 23 multi-assay proteins are shown as 23 plots).

Disagreement = 1 - mean pairwise within-site Spearman (45 family-best models),
read from results/generalization_functional_disagreement_sites.csv (script 15).
The 0.5 line marks the high-/low-disagreement boundary used by 27 / 30; the shaded
band is the high-disagreement region (> 0.5).

Output (one per protein) in figures/site_disagreement/distribution/:
  <common>_disagreement_dist.png     e.g. gck_disagreement_dist.png
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import common as C

# import 23's helpers (module name can't start with a digit)
_S23 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "23_structure_maps_all.py")
_spec = importlib.util.spec_from_file_location("structure_maps_all", _S23)
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)

OUTDIR = os.path.join(C.FIGURES, "site_disagreement", "distribution")
BINS = np.linspace(0.0, 1.0, 41)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    sites = pd.read_csv(S.SITES)
    names = dict(zip(pd.read_csv(S.CANDID)["UniProt_ID"], pd.read_csv(S.CANDID)["protein"]))
    ref = C.load_reference()
    cst = dict(zip(ref["DMS_id"], ref["coarse_selection_type"]))   # assay -> category (fig9 colours)

    for uni, g in sites.groupby("UniProt_ID"):
        cname = S.common_name(uni, names.get(uni, uni))
        assays = sorted(g["DMS_id"].unique())
        xmax = max(0.6, g["disagreement"].max() + 0.03)
        # one panel per assay, stacked, shared x
        fig, axes = plt.subplots(len(assays), 1, sharex=True, squeeze=False,
                                 figsize=(8.5, 1.9 * len(assays) + 1.0))
        axes = axes[:, 0]
        for ax, dms_id in zip(axes, assays):
            d = g.loc[g["DMS_id"] == dms_id, "disagreement"].dropna()
            cat = cst.get(dms_id, "")
            col = C.CATEGORY_COLORS.get(cat, "#555555")   # colour by assay category (matches fig9)
            short = dms_id.replace(uni + "_", "")
            # per-assay high/low selection cuts (same as 27/30)
            hi_thr = max(S.elbow_tail_threshold(d.values)[0], S.DISAGREE_FLOOR)
            lo_thr = min(S.elbow_low_tail_threshold(d.values)[0], S.DISAGREE_CEIL)
            # boundary line shared by 27/30
            ax.axvline(S.DISAGREE_FLOOR, color="#555", lw=0.9, ls="--", zorder=1)
            # histogram: opaque in the selected high/low tails, faded in the bulk
            _, edges, patches = ax.hist(d, bins=BINS, density=True, color=col,
                                        edgecolor=col, linewidth=0.4, zorder=2)
            half = (edges[1] - edges[0]) / 2
            for patch, left in zip(patches, edges[:-1]):
                c = left + half
                patch.set_alpha(0.85 if (c >= hi_thr or c <= lo_thr) else 0.15)
            ax.axvline(d.median(), color=col, lw=1.3, ls=":", alpha=0.9, zorder=3)
            ax.set_ylabel("density", fontsize=8)
            ax.text(0.99, 0.92, f"{short}\n{cat} · n={len(d)}, median={d.median():.2f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
                    color=col, fontweight="bold")
        axes[-1].set_xlim(0, xmax)
        axes[-1].set_xlabel("within-site disagreement  (1 - mean pairwise within-site Spearman)")
        axes[0].set_title(f"{names.get(uni, uni)} [{uni}] — per-site disagreement distribution\n"
                          f"{len(assays)} assay(s), {g['up_pos'].nunique()} sites  "
                          "(opaque = selected high/low-disagreement tails, 27/30; "
                          "dashed = 0.5 boundary; dotted = median)", fontsize=10)
        plt.tight_layout()
        out = os.path.join(OUTDIR, f"{cname}_disagreement_dist.png")
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"{uni:16s} -> {os.path.basename(out):40s} {len(assays)} assays, "
              f"{g['up_pos'].nunique()} sites")


if __name__ == "__main__":
    main()
