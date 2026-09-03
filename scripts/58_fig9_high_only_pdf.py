"""fig9, high-disagreement row only, as a PDF (companion to `29`'s fig9 PNG).

`29_family_best_category_analysis.py` makes fig9 as a 2x2 PNG (rows = high/low
disagreement sites, cols = ranked by category Spearman / NDCG@10). This script
reuses its cached output (family_best_highdisagreement_within_site.csv) to render
just the top row (high-disagreement) as a 1x2 PDF -- same panel logic, same
colors/legend -- without touching or regenerating the original PNG.

Input:  results/persite_analysis/family_best_highdisagreement_within_site.csv
Output: figures/persite_analysis/fig9_highdisagreement_only.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import common as C

MIN_SITE_SUBS = 5  # must match 29's threshold used to build the cached CSV


def _panel(ax, ws, rank_col, tiers, cats, cat_label):
    width = 0.8 / len(cats)
    x = np.arange(len(tiers))
    for i, cat in enumerate(cats):
        g = ws[ws["coarse_selection_type"] == cat]
        sub = [g.loc[g[rank_col] <= thr, "within_site_spearman"] for _, thr in tiers]
        vals = [v.mean() for v in sub]
        errs = [v.std(ddof=1) / np.sqrt(v.count()) for v in sub]
        ax.bar(x + i * width, vals, width, yerr=errs, label=cat_label.get(cat, cat),
               color=C.CATEGORY_COLORS.get(cat), error_kw=dict(elinewidth=0.7, capsize=2, alpha=0.6))
    ax.set_xticks(x + width * (len(cats) - 1) / 2)
    ax.set_xticklabels([t for t, _ in tiers])
    ax.axhline(0, color="0.6", lw=0.8)


def main():
    ws_high = pd.read_csv(C.rp("family_best_highdisagreement_within_site.csv"))
    tiers = [("Top 1", 1), ("Top 3", 3), ("Top 5", 5), ("Top 10", 10), ("All", 10**9)]
    cats = sorted(ws_high["coarse_selection_type"].unique())
    n_models = ws_high["model"].nunique()

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.3), sharey=True)
    cat_label = {c: f"{c}  ({ws_high.loc[ws_high['coarse_selection_type'] == c, 'DMS_id'].nunique()} "
                    f"assays, {ws_high.loc[ws_high['coarse_selection_type'] == c, 'UniProt_ID'].nunique()} UniProt)"
                 for c in cats}
    for col, (rank_col, metric) in enumerate([("spearman_rank", "category Spearman"),
                                              ("ndcg_rank", "category NDCG@10")]):
        ax = axes[col]
        _panel(ax, ws_high, rank_col, tiers, cats, cat_label)
        ax.set_title(f"High-disagreement sites (>= {MIN_SITE_SUBS} subs)\nranked by {metric}", fontsize=11)
        ax.set_xlabel("models kept (by category rank)")
    axes[0].set_ylabel("mean within-site Spearman with DMS\n(error bars = SE across sites)")
    axes[1].legend(title="assay category", fontsize=8)
    fig.suptitle(f"Do better-ranked models win at high-disagreement sites?\n"
                 f"family-best subset: {n_models} models (one per family)", y=1.02, fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    figdir = os.path.join(C.FIGURES, "persite_analysis")
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, "fig9_highdisagreement_only.pdf")
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
