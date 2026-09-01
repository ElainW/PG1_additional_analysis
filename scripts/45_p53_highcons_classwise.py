"""
P53 disagree-at-CONSERVED: do sequence/MSA models beat structure-based models at
P53's shared high-disagreement, HIGH-conservation sites?

Sites: per P53 assay, positions flagged high-disagreement in >=2 P53 assays
(colocalization_shared_sites.csv) with conservation_percentile > 80 and >=5
missense subs -- reused from script 40 (selected_high).
Pool: the 45 family-best models. For each model, its MEDIAN per-site Spearman over
those sites (one value per model). Then per assay, a box+strip over the 5 populated
model classes (each dot = one model), so we can compare classes directly. The 6th
class (Single sequence, Structure & Function annotations = ESM3 only) is dropped --
a single model can't form a box; its value is printed for reference.

Outputs:
  results/MSA_conservation/H2/h2_p53_highcons_classwise_persite.csv
  figures/h2_persite_spearman_density/fig_p53_highcons_classwise.png
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import common as C

sns.set_style("whitegrid")
rng = np.random.default_rng(0)

_spec = importlib.util.spec_from_file_location(
    "m40", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "40_h2_highcons_shared_strip_panel.py"))
M40 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M40)
M38 = M40.M38

H2DIR = M40.H2DIR if hasattr(M40, "H2DIR") else os.path.join(C.RESULTS, "MSA_conservation", "H2")
FB = pd.read_csv(os.path.join(C.RESULTS, "persite_analysis", "family_best_selection.csv"))
MODELS = list(FB.selected_model)
MT = pd.read_csv(C.rp("directionality_check.csv"))
M2T = dict(zip(MT.model, MT.model_type))

# sequence -> structure ordering; the Function class (ESM3) is a single model -> 1 dot
FUNC_CLASS = "Single sequence, Structure & Function annotations"
CLASS_ORDER = ["MSA", "Single sequence", "Single sequence & Structure",
               "Structure & MSA", "Structure", FUNC_CLASS]
SHORT = M40.SHORT
PALETTE = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
SEQ_MSA = {"MSA", "Single sequence"}
STRUCT_BASED = {"Structure", "Structure & MSA", "Single sequence & Structure"}
ASSAYS = M40.PROT_ASSAYS["P53_HUMAN"]


def main():
    hd = M40.shared_high_sites()
    rows = []
    for dms_id in ASSAYS:
        df, sites = M40.selected_high(dms_id, hd)
        if not sites:
            continue
        for m in [x for x in MODELS if x in df.columns]:
            rhos = [M38.persite_spearman(sub, m) for _, _, _, sub in sites]
            rhos = [r for r in rhos if pd.notna(r)]
            if rhos:
                rows.append({"DMS_id": dms_id, "n_sites": len(sites), "model": m,
                             "model_type": M2T[m],
                             "median_persite_spearman": round(float(np.median(rhos)), 4)})
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(H2DIR, "h2_p53_highcons_classwise_persite.csv"), index=False)

    present = [a for a in ASSAYS if a in set(tab.DMS_id)]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharey=True, squeeze=False)
    print("\nseq/MSA vs structure-based median per-site Spearman (per model), by assay:")
    for i, dms_id in enumerate(present):
        ax = axes[i // 2][i % 2]
        d = tab[(tab.DMS_id == dms_id) & (tab.model_type.isin(CLASS_ORDER))]
        sns.boxplot(data=d, x="model_type", y="median_persite_spearman", order=CLASS_ORDER,
                    hue="model_type", palette={c: PALETTE[c] for c in CLASS_ORDER},
                    legend=False, ax=ax, width=0.6, fliersize=0, linewidth=1,
                    boxprops=dict(alpha=0.35), showcaps=True)
        for xi, c in enumerate(CLASS_ORDER):
            dd = d[d.model_type == c]
            x = xi + rng.uniform(-0.15, 0.15, len(dd))
            ax.scatter(x, dd.median_persite_spearman, color=PALETTE[c], s=45,
                       edgecolor="black", linewidth=0.5, alpha=0.9, zorder=3)
        ax.axhline(0, color="0.4", ls="--", lw=1)
        n = tab[tab.DMS_id == dms_id].n_sites.iloc[0]
        ax.set_title(f"{dms_id}\n(high-cons shared sites: n={n})", fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("median per-site Spearman (per model)")
        ax.set_xticks(range(len(CLASS_ORDER)))
        ax.set_xticklabels([f"{SHORT[c]}\n(n={ (d.model_type==c).sum() })" for c in CLASS_ORDER],
                           fontsize=8)
        ax.set_ylim(-1.05, 1.05)
        ax.yaxis.grid(False)
        # seq/MSA vs structure-based comparison
        a = tab[(tab.DMS_id == dms_id) & (tab.model_type.isin(SEQ_MSA))].median_persite_spearman
        b = tab[(tab.DMS_id == dms_id) & (tab.model_type.isin(STRUCT_BASED))].median_persite_spearman
        p = mannwhitneyu(a, b, alternative="greater").pvalue if len(a) and len(b) else np.nan
        ax.text(0.02, 0.97, f"seq/MSA med {a.median():.2f}  vs  struct med {b.median():.2f}\n"
                            f"MWU seq>struct p={p:.3f}", transform=ax.transAxes,
                va="top", ha="left", fontsize=8, color="0.2")
        esm3 = tab[(tab.DMS_id == dms_id) & (tab.model_type == FUNC_CLASS)].median_persite_spearman
        print(f"  {dms_id}: seq/MSA med {a.median():.2f} (n={len(a)}) vs struct med "
              f"{b.median():.2f} (n={len(b)})  MWU p(seq>struct)={p:.3f}"
              + (f"  [ESM3={esm3.iloc[0]:.2f}]" if len(esm3) else ""))
    for j in range(len(present), 4):
        axes[j // 2][j % 2].axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=PALETTE[c], alpha=0.6, ec="none",
                             label=f"{SHORT[c]} = {c}") for c in CLASS_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.02), title="model_type class (box/dot colour); each dot = one model")
    fig.suptitle("P53 shared high-disagreement, HIGH-conservation sites (pct>80) — median "
                 "per-site Spearman by model class  (45 family-best; Seq+Struct+Fn = ESM3, "
                 "1 model)", fontsize=12, y=1.0)
    plt.tight_layout(rect=(0, 0.05, 1, 0.96))
    out = os.path.join(C.FIGURES, "h2_persite_spearman_density", "fig_p53_highcons_classwise.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("\nsaved", out)


if __name__ == "__main__":
    main()
