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
COLORS = plt.get_cmap("tab10").colors


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    sites = pd.read_csv(S.SITES)
    names = dict(zip(pd.read_csv(S.CANDID)["UniProt_ID"], pd.read_csv(S.CANDID)["protein"]))

    for uni, g in sites.groupby("UniProt_ID"):
        cname = S.common_name(uni, names.get(uni, uni))
        assays = sorted(g["DMS_id"].unique())
        fig, ax = plt.subplots(figsize=(8.5, 5))
        # high-disagreement region + boundary shared by 27/30
        ax.axvspan(S.DISAGREE_FLOOR, 1.0, color="#999", alpha=0.07, zorder=0)
        ax.axvline(S.DISAGREE_FLOOR, color="#555", lw=0.9, ls="--", zorder=1,
                   label="high/low boundary (0.5)")
        for i, dms_id in enumerate(assays):
            d = g.loc[g["DMS_id"] == dms_id, "disagreement"].dropna()
            col = COLORS[i % len(COLORS)]
            short = dms_id.replace(uni + "_", "")
            ax.hist(d, bins=BINS, density=True, histtype="step", lw=1.8,
                    color=col, zorder=3, label=f"{short}  (n={len(d)}, med={d.median():.2f})")
            ax.hist(d, bins=BINS, density=True, histtype="stepfilled", lw=0,
                    color=col, alpha=0.12, zorder=2)
            ax.axvline(d.median(), color=col, lw=1.1, ls=":", alpha=0.8, zorder=3)
        ax.set_xlim(0, max(0.6, g["disagreement"].max() + 0.03))
        ax.set_xlabel("within-site disagreement  (1 - mean pairwise within-site Spearman)")
        ax.set_ylabel("density")
        ax.set_title(f"{names.get(uni, uni)} [{uni}] — per-site disagreement distribution\n"
                     f"{len(assays)} assay(s), {g['up_pos'].nunique()} sites "
                     "(dotted = per-assay median)", fontsize=11)
        ax.legend(fontsize=8, title="assay")
        plt.tight_layout()
        out = os.path.join(OUTDIR, f"{cname}_disagreement_dist.png")
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"{uni:16s} -> {os.path.basename(out):40s} {len(assays)} assays, "
              f"{g['up_pos'].nunique()} sites")


if __name__ == "__main__":
    main()
