"""
H2: box+strip of per-site Spearman at each assay's assay-specific high-disagreement
sites, for 3 models chosen by that assay's WHOLE-ASSAY within-site ranking.

For each assay, rank the 45 family-best models by their median per-site Spearman
over ALL >=5-sub positions in the assay (whole-assay within-site accuracy). Pick:
  - Top #1
  - Rank 23 (middle)
  - Lowest (#N, N = models present in the assay)
Then plot each chosen model's per-site Spearman distribution over the assay-SPECIFIC
high-disagreement sites (script-27 high-disagreement minus shared positions, >=5
subs). A dotted mark per box shows the model's whole-assay within-site median (the
ranking basis) so you can see how the hard sites compare to the model's typical
within-site accuracy on that assay.

Outputs:
  results/MSA_conservation/H2/h2_assayspecific_rank_tiers.csv
  figures/h2_persite_spearman_density/fig_assayspecific_rank_tiers.pdf
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import common as C

sns.set_style("whitegrid")
rng = np.random.default_rng(0)

_spec = importlib.util.spec_from_file_location(
    "m41", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "41_assayspecific_category_best.py"))
M41 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M41)
M38 = M41.M38

H2DIR = M41.H2DIR
FB = pd.read_csv(os.path.join(C.RESULTS, "persite_analysis", "family_best_selection.csv"))
MODELS = list(FB.selected_model)                       # 45 family-best
MT = pd.read_csv(C.rp("directionality_check.csv"))
M2T = dict(zip(MT.model, MT.model_type))
OVERALL_RHO = dict(zip(MT.model, MT.mean_spearman_vs_DMS))

PALETTE_MT = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
SHORT = M41.SHORT


def whole_assay_medians(dms_id, present):
    """median per-site Spearman over ALL >=5-sub positions in the assay, per model."""
    df = M38._scores(dms_id)
    subs = [g for _, g in df.groupby("p") if len(g) >= M41.MIN_SUBS]
    med = {}
    for m in present:
        rhos = [M38.persite_spearman(s, m) for s in subs]
        rhos = [r for r in rhos if pd.notna(r)]
        if rhos:
            med[m] = float(np.median(rhos))
    return med


def main():
    shared = M41.shared_positions()
    order, rows, tier_pick = [], [], {}
    for uni in M41.COL_ORDER:
        hd = pd.read_csv(os.path.join(C.RESULTS, "high_disagreement_sites",
                                      f"{M41.HD_FILES[uni]}_high_disagreement.tsv"), sep="\t")
        for dms_id in sorted(hd.DMS_id.unique()):
            df, sites = M41.assay_specific_sites(uni, dms_id, shared)
            if not sites:
                continue
            cat = M38.REF.loc[dms_id, "coarse_selection_type"]
            present = [m for m in MODELS if m in df.columns]
            med = whole_assay_medians(dms_id, present)
            ranked = sorted(med, key=med.get, reverse=True)          # best -> worst
            n = len(ranked)
            picks = [("top", 1, ranked[0]),
                     ("mid", 23, ranked[22]),
                     ("low", n, ranked[-1])]
            tier_pick[dms_id] = picks
            order.append(dms_id)
            wtp = M41.wt_pct_of(df, dms_id)
            for tier, rank, m in picks:
                for up, wt, sub in sites:
                    med_pct = sub["dms_pct"].median()
                    rows.append({"UniProt_ID": uni, "DMS_id": dms_id, "category": cat,
                                 "n_sites": len(sites), "tier": tier, "assay_rank": rank,
                                 "model": m, "model_type": M2T[m],
                                 "whole_assay_median": round(med[m], 4),
                                 "overall_rho": round(OVERALL_RHO[m], 4),
                                 "up_pos": up, "residue": f"{wt}{up}",
                                 "persite_spearman": M38.persite_spearman(sub, m),
                                 "dms_pct_diff_from_wt": round(med_pct - wtp, 2)})
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(H2DIR, "h2_assayspecific_rank_tiers.csv"), index=False)

    tiers = ["top", "mid", "low"]
    ncols = 4
    nrows = int(np.ceil(len(order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.8 * nrows),
                             sharey=True, squeeze=False)
    used_mt = []
    for i, dms_id in enumerate(order):
        ax = axes[i // ncols][i % ncols]
        d = tab[tab.DMS_id == dms_id]
        tier_mt = {t: d[d.tier == t].model_type.iloc[0] for t in tiers}   # class per tier
        tier_pal = {t: PALETTE_MT[tier_mt[t]] for t in tiers}
        for mt in tier_mt.values():
            if mt not in used_mt:
                used_mt.append(mt)
        sns.boxplot(data=d, x="tier", y="persite_spearman", order=tiers,
                    hue="tier", palette=tier_pal, legend=False, ax=ax, width=0.55,
                    fliersize=0, linewidth=1, boxprops=dict(alpha=0.35), showcaps=True)
        for xi, t in enumerate(tiers):
            dd = d[d.tier == t]
            x = xi + rng.uniform(-0.15, 0.15, len(dd))
            ax.scatter(x, dd.persite_spearman, color=tier_pal[t], s=40,
                       edgecolor="black", linewidth=0.5, alpha=0.9, zorder=3)
            wam = dd.whole_assay_median.iloc[0]           # ranking basis (dotted mark)
            ax.hlines(wam, xi - 0.32, xi + 0.32, color="black", ls=":", lw=1.6, zorder=4)
        ax.axhline(0, color="0.4", ls="--", lw=1)
        cat = d.category.iloc[0]
        n = d.n_sites.iloc[0]
        ax.set_title(f"{dms_id}\ncategory: {cat}  (n={n} assay-specific sites)", fontsize=9)
        ax.set_ylim(-1.05, 1.05)
        ax.yaxis.grid(False)
        ax.set_ylabel("per-site Spearman (model vs DMS)" if i % ncols == 0 else "")
        ax.set_xlabel("")
        labs = []
        for tier, rank, m in tier_pick[dms_id]:
            head = {"top": "Top #1", "mid": "Rank 23", "low": f"Lowest #{rank}"}[tier]
            labs.append(f"{head}\n{M41.MODEL_ABBR.get(m, m)}\nassay med {tab.loc[(tab.DMS_id==dms_id)&(tab.tier==tier),'whole_assay_median'].iloc[0]:.2f}")
        ax.set_xticks(range(3))
        ax.set_xticklabels(labs, fontsize=7)
    for j in range(len(order), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=PALETTE_MT[t], alpha=0.6, ec="none",
                             label=f"{SHORT[t]} = {t}")
               for t in C.MODEL_TYPE_ORDER if t in used_mt]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.03), title="box/dot colour = model_type class")
    fig.suptitle("Assay-specific high-disagreement sites — per-site Spearman for the assay's "
                 "top / rank-23 / lowest within-site model  (45 family-best; dotted = model's "
                 "whole-assay within-site median)", fontsize=11, y=1.0)
    plt.tight_layout(rect=(0, 0.05, 1, 0.95))
    out = os.path.join(C.FIGURES, "h2_persite_spearman_density",
                       "fig_assayspecific_rank_tiers.pdf")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("saved", out)
    print("\nchosen models per assay (whole-assay within-site rank):")
    for dms_id in order:
        picks = tier_pick[dms_id]
        print(f"  {dms_id}: " + " | ".join(
            f"{t}#{r} {m} ({tab.loc[(tab.DMS_id==dms_id)&(tab.tier==t),'whole_assay_median'].iloc[0]:.2f})"
            for t, r, m in picks))


if __name__ == "__main__":
    main()
