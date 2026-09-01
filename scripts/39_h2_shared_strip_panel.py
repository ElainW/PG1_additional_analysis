"""
Combined panel (ONE COLUMN PER UniProt: PTEN | HXK4/GCK | CP2C9; two assay rows)
of the per-site model Spearman at SHARED high-disagreement / low-conservation
sites, as box + strip plots.

Reads results/persite_analysis/h2_lowcons_shared_persite_spearman_by_site.csv (38).
Boxes are coloured by model_type class (standard palette). Each dot is one site,
drawn in a single-hue white->black greyscale encoding how far that site's DMS
signal sits from wild-type: |median(site substitution DMS percentiles) - WT DMS
percentile|, percentiles taken over the assay's entire DMS range. White = site
median at the WT percentile (close to WT), black = farthest. One shared colourbar.

Output: figures/h2_persite_spearman_density/fig_shared_lowcons_persite_spearman_panel.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import seaborn as sns
import common as C

sns.set_style("whitegrid")
rng = np.random.default_rng(0)

CSV = C.rp("h2_lowcons_shared_persite_spearman_by_site.csv")

# columns = UniProt; rows = that protein's assays (top -> bottom)
PROT_ASSAYS = {
    "PTEN_HUMAN":  ["PTEN_HUMAN_Matreyek_2021", "PTEN_HUMAN_Mighell_2018"],
    "HXK4_HUMAN":  ["HXK4_HUMAN_Gersing_2022_activity",
                    "HXK4_HUMAN_Gersing_2023_abundance"],
    "CP2C9_HUMAN": ["CP2C9_HUMAN_Amorosi_2021_abundance",
                    "CP2C9_HUMAN_Amorosi_2021_activity"],
}
COL_ORDER = ["PTEN_HUMAN", "HXK4_HUMAN", "CP2C9_HUMAN"]

SHORT = {
    "MSA": "MSA",
    "Single sequence": "SeqOnly",
    "Structure": "Struct",
    "Single sequence & Structure": "Seq+Struct",
    "Structure & MSA": "Struct+MSA",
    "Single sequence, Structure & Function annotations": "Seq+Struct+Fn",
}
CMAP = cm.get_cmap("RdBu")    # signed: red = below WT (deleterious), blue = above, white = at WT


def main():
    df = pd.read_csv(CSV)
    order = [t for t in C.MODEL_TYPE_ORDER if t in set(df.model_type)]
    palette = dict(zip(C.MODEL_TYPE_ORDER,
                       sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
    box_pal = {t: palette[t] for t in order}
    model_of = df.drop_duplicates("model_type").set_index("model_type")["model"].to_dict()
    # each chosen model's overall benchmark Spearman (mean over 217 assays)
    MT = pd.read_csv(C.rp("directionality_check.csv"))
    overall_rho = dict(zip(MT.model, MT.mean_spearman_vs_DMS))

    M = float(df.dms_pct_diff_from_wt.abs().max())
    norm = Normalize(vmin=-M, vmax=M)          # symmetric, 0 -> white

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), sharey=True)
    for c, uni in enumerate(COL_ORDER):
        for r, dms_id in enumerate(PROT_ASSAYS[uni]):
            ax = axes[r, c]
            d = df[df.DMS_id == dms_id]
            cat = d["category"].iloc[0]
            n_sites = d["up_pos"].nunique()
            wtp = d["wt_dms_pct"].iloc[0]

            sns.boxplot(data=d, x="model_type", y="persite_spearman", order=order,
                        hue="model_type", palette=box_pal, legend=False, ax=ax,
                        width=0.6, fliersize=0, linewidth=1,
                        boxprops=dict(alpha=0.35), showcaps=True)
            # dots: one per site per class, shaded by DMS-median distance to WT
            for xi, t in enumerate(order):
                dt = d[d.model_type == t]
                x = xi + rng.uniform(-0.16, 0.16, len(dt))
                ax.scatter(x, dt.persite_spearman,
                           c=CMAP(norm(dt.dms_pct_diff_from_wt.values)),
                           s=48, edgecolor="black", linewidth=0.6, zorder=3)

            ax.axhline(0, color="0.4", ls="--", lw=1)
            ax.set_title(f"{dms_id}\ncategory: {cat}   (n = {n_sites} sites; "
                         f"WT at pct {wtp:.0f})", fontsize=9)
            ax.set_xlabel("")
            ax.set_ylabel("per-site Spearman (model vs DMS)" if c == 0 else "")
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels([SHORT[t] for t in order], rotation=40,
                               ha="right", fontsize=8)
            ax.set_ylim(-1.05, 1.05)
            ax.yaxis.grid(False)     # drop horizontal gridlines
        axes[0, c].annotate(uni.replace("_HUMAN", ""), xy=(0.5, 1.22),
                            xycoords="axes fraction", ha="center",
                            fontsize=13, fontweight="bold")

    fig.suptitle("Per-site Spearman at SHARED high-disagreement, low-conservation "
                 "sites (cons pct < 20 or NaN; ≥5 subs)", fontsize=13, y=1.0)
    plt.tight_layout(rect=(0, 0.08, 0.9, 0.96))

    # shared colourbar on the far right (dot shade = DMS-percentile distance to WT)
    sm = cm.ScalarMappable(norm=norm, cmap=CMAP)
    sm.set_array([])
    cax = fig.add_axes([0.925, 0.32, 0.014, 0.42])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("median(site DMS percentile) − WT DMS percentile\n"
                   "red = below WT (deleterious), white = at WT, blue = above WT",
                   fontsize=9)

    # shared legend: class -> best-in-class model (box colours)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=palette[t], alpha=0.5, ec="none",
                             label=f"{SHORT[t]} = {t} "
                                   f"({model_of.get(t,'?')}, "
                                   f"overall ρ={overall_rho.get(model_of.get(t),float('nan')):.2f})")
               for t in order]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.0),
               title="box colour = model_type class "
                     "(best-in-class model, overall Spearman over 217 assays)")
    out = os.path.join(C.FIGURES, "h2_persite_spearman_density",
                       "fig_shared_lowcons_persite_spearman_panel.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("saved", out)


if __name__ == "__main__":
    main()
