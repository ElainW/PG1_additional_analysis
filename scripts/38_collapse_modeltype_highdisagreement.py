"""Winners-collapse by model_type, generalized to ALL high-disagreement sites
(pooled across the multi-assay proteins) — companion to `16` (which is HXK4-only).

For every multi-assay assay that has labelled high-disagreement sites
(results/high_disagreement_sites/*.tsv), pool substitutions and compute the
'winners collapse': mean Spearman(model, DMS) vs experiment across fitness
subsets (all / top50 / top25 / top10), under a GLOBAL and a WITHIN-SITE fitness
split, averaged by model_type (45 family-best models).

Two site sets are overlaid so the collapse at hard sites is read against baseline:
  - HIGH  = the labelled high-disagreement positions (up_pos matches the score-file
            mutant position directly; verified incl. offset proteins)
  - ALL   = every position with >= MIN_SUBS substitutions

Aggregation pools across assays robustly: per-(assay, model) Spearman on each
subset, then averaged over all (assay, model) of a model_type (no cross-assay
DMS-scale mixing).

Outputs:
  results/persite_analysis/collapse_modeltype_highdisagreement.csv
  figures/persite_analysis/fig_collapse_modeltype_highdisagreement.png
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
import common as C
from importlib import import_module

fam = import_module("05_fig6_familybest")

SUBSETS = [("all", 0), ("top50", 50), ("top25", 75), ("top10", 90)]
SPLITS = [("global", "glob_pct"), ("within-site", "ws_pct")]
SITE_SETS = ["all", "high"]
MIN_SUBS = 5
MIN_N = 10          # min variants in an (assay, subset, site_set) cell to score
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def load_high_sites():
    """{DMS_id: set(up_pos)} of labelled high-disagreement positions."""
    frames = []
    for f in glob.glob(os.path.join(ROOT, "results/high_disagreement_sites/*_high_disagreement.tsv")):
        try:
            frames.append(pd.read_csv(f, sep="\t"))
        except pd.errors.EmptyDataError:
            continue
    hd = pd.concat(frames, ignore_index=True)
    return {k: set(g["up_pos"]) for k, g in hd.groupby("DMS_id")}


def main():
    chosen, type_of, _ = fam.build_family_best()
    mtypes = [t for t in C.MODEL_TYPE_ORDER if t in set(type_of.values())]
    high_by_assay = load_high_sites()

    rows = []
    for dms_id, high_pos in high_by_assay.items():
        sf = os.path.join(C.SCORES_DIR, dms_id + ".csv")
        if not os.path.exists(sf):
            continue
        df = pd.read_csv(sf)
        df = df[df["mutant"].str.match(r"^[A-Z]\d+[A-Z]$")].copy()
        df["pos"] = df["mutant"].str.extract(r"(\d+)").astype(int)
        df = df[df.groupby("pos")["pos"].transform("size") >= MIN_SUBS]
        if not len(df):
            continue
        present = [m for m in chosen if m in df.columns]

        g = df.groupby("pos")["DMS_score"]
        r = g.rank(method="average")
        n = g.transform("size")
        df["ws_pct"] = (r - 1) / (n - 1) * 100
        df["glob_pct"] = df["DMS_score"].rank(pct=True) * 100
        df["is_high"] = df["pos"].isin(high_pos)
        y = df["DMS_score"].values

        for site_set in SITE_SETS:
            base = np.ones(len(df), bool) if site_set == "all" else df["is_high"].values
            for split, pctcol in SPLITS:
                pct = df[pctcol].values
                for sublab, q in SUBSETS:
                    mask = base & (pct >= q if q > 0 else np.ones(len(df), bool))
                    if mask.sum() < MIN_N:
                        continue
                    for m in present:
                        s = df[m].values[mask]
                        if np.ptp(s) == 0:
                            continue
                        rho = spearmanr(y[mask], s).statistic
                        if np.isnan(rho):
                            continue
                        rows.append({"DMS_id": dms_id, "site_set": site_set,
                                     "split": split, "subset": sublab,
                                     "model": m, "model_type": type_of[m],
                                     "n_variants": int(mask.sum()), "spearman": rho})

    per = pd.DataFrame(rows)
    agg = (per.groupby(["site_set", "split", "subset", "model_type"])
              .agg(mean_spearman=("spearman", "mean"),
                   n_obs=("spearman", "size"),
                   n_assays=("DMS_id", "nunique"))
              .reset_index())
    dst = os.path.join(ROOT, "results/persite_analysis/collapse_modeltype_highdisagreement.csv")
    agg.to_csv(dst, index=False)
    print(f"wrote {dst}")
    print(f"assays pooled: {per['DMS_id'].nunique()}")

    # ---- figure: 2 panels (split); color=model_type; solid=high, dashed=all ----
    palette = dict(zip(C.MODEL_TYPE_ORDER, plt.cm.tab10(np.linspace(0, 1, len(C.MODEL_TYPE_ORDER)))))
    x = np.arange(len(SUBSETS))
    ls = {"high": "-", "all": "--"}
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for i, (split, _) in enumerate(SPLITS):
        ax = axes[i]
        for site_set in SITE_SETS:
            # within-site split: too few high-disagreement sites for a reliable
            # within-site quartile pool -> show all-sites only
            if split == "within-site" and site_set == "high":
                continue
            for t in mtypes:
                tt = (agg[(agg.split == split) & (agg.site_set == site_set) & (agg.model_type == t)]
                      .set_index("subset").reindex([s for s, _ in SUBSETS]))
                lw = 3.0 if t == "Structure" else 1.6
                ax.plot(x, tt["mean_spearman"].values, marker="o" if site_set == "high" else None,
                        color=palette[t], lw=lw, ls=ls[site_set],
                        alpha=1.0 if site_set == "high" else 0.55,
                        zorder=5 if t == "Structure" else 3)
        ax.axhline(0, color="gray", lw=0.6)
        ax.set_xticks(x); ax.set_xticklabels(["all", "top 50%", "top 25%", "top 10%"])
        ax.set_title(f"{split} fitness split", fontsize=12)
        ax.set_xlabel("fitness subset")
        if i == 0:
            ax.set_ylabel("mean Spearman vs experiment")

    type_handles = [Line2D([0], [0], color=palette[t], lw=3.0 if t == "Structure" else 1.8, label=t)
                    for t in mtypes]
    set_handles = [Line2D([0], [0], color="gray", lw=1.8, ls="-", marker="o", label="high-disagreement sites"),
                   Line2D([0], [0], color="gray", lw=1.8, ls="--", label="all sites (≥5 subs)")]
    leg1 = fig.legend(handles=type_handles, loc="upper center", ncol=6, fontsize=8,
                      bbox_to_anchor=(0.5, 1.02))
    fig.add_artist(leg1)
    fig.legend(handles=set_handles, loc="upper center", ncol=2, fontsize=9,
               bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("Winners-collapse by model_type — ALL high-disagreement sites pooled "
                 "(solid) vs all sites (dashed); Structure = bold", y=1.10, fontsize=13)
    plt.tight_layout()
    figdir = os.path.join(C.FIGURES, "persite_analysis")
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, "fig_collapse_modeltype_highdisagreement.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"saved {out}")

    # headline: high vs all collapse, per split, pooled over model types
    piv = (agg.groupby(["site_set", "split", "subset"])["mean_spearman"].mean()
              .reset_index()
              .pivot_table(index=["split", "site_set"], columns="subset", values="mean_spearman")
              [["all", "top50", "top25", "top10"]])
    print("\nmean Spearman (avg over model_types):")
    print(piv.round(3).to_string())


if __name__ == "__main__":
    main()
