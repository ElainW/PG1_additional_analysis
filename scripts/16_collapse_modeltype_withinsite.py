"""
'Winners collapse' by model_TYPE, pooled across ALL substitution assays.

For every assay we compute Spearman(model, DMS) over each fitness subset
(all / top50 / top25 / top10) under two fitness-split definitions:
  - GLOBAL      : percentile of a substitution's DMS score across the whole assay
  - WITHIN-SITE : percentile among the substitutions at its own position (so the
                  "winners" are the best substitution per position, spread over all
                  positions, rather than the globally most-fit variants).
Per-(assay, model) Spearmans are then averaged by model_type. Pooling gives each
assay equal weight, so each coarse_selection_type category contributes in
proportion to its number of assays ("weighted by #assays per category"). Uses the
45-model family-best set. Sites need >= MIN_SUBS substitutions.

Outputs:
  results/collapse_modeltype_withinsite.csv
  figures/persite_analysis/fig_collapse_modeltype_pooled.png
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import common as C
from importlib import import_module

fam = import_module("05_fig6_familybest")

SUBSETS = [("all", 0), ("top50", 50), ("top25", 75), ("top10", 90)]
SPLITS = [("global", "glob_pct"), ("within-site", "ws_pct")]
MIN_SUBS = 5
MIN_N = 10          # min variants in a subset to compute a Spearman


def main():
    chosen, type_of, _ = fam.build_family_best()
    mtypes = [t for t in C.MODEL_TYPE_ORDER if t in set(type_of.values())]

    ref = pd.read_csv(C.REF_PATH)
    cat_of = dict(zip(ref["DMS_id"], ref["coarse_selection_type"]))

    rows = []
    for sf in sorted(glob.glob(os.path.join(C.SCORES_DIR, "*.csv"))):
        dms_id = os.path.basename(sf)[:-4]
        df = pd.read_csv(sf)
        if "mutant" not in df.columns or "DMS_score" not in df.columns:
            continue
        df = df[df["mutant"].astype(str).str.match(r"^[A-Z]\d+[A-Z]$")].copy()
        if not len(df):
            continue
        df["pos"] = df["mutant"].str.extract(r"(\d+)").astype(int)
        df = df[df.groupby("pos")["pos"].transform("size") >= MIN_SUBS]
        present = [m for m in chosen if m in df.columns]
        if not len(df) or not present:
            continue

        g = df.groupby("pos")["DMS_score"]
        r = g.rank(method="average")
        n = g.transform("size")
        df["ws_pct"] = (r - 1) / (n - 1) * 100
        df["glob_pct"] = df["DMS_score"].rank(pct=True) * 100
        y = df["DMS_score"].values

        for split, pctcol in SPLITS:
            pct = df[pctcol].values
            for sublab, q in SUBSETS:
                mask = pct >= q if q > 0 else np.ones(len(df), bool)
                if mask.sum() < MIN_N:
                    continue
                for m in present:
                    s = df[m].values[mask]
                    if np.ptp(s) == 0:
                        continue
                    rho = spearmanr(y[mask], s).statistic
                    if np.isnan(rho):
                        continue
                    rows.append({"DMS_id": dms_id, "category": cat_of.get(dms_id),
                                 "split": split, "subset": sublab, "model": m,
                                 "model_type": type_of[m], "spearman": rho})

    per = pd.DataFrame(rows)
    agg = (per.groupby(["split", "subset", "model_type"])
              .agg(mean_spearman=("spearman", "mean"),
                   n_obs=("spearman", "size"),
                   n_assays=("DMS_id", "nunique"))
              .reset_index())
    agg.to_csv(C.rp("collapse_modeltype_withinsite.csv"), index=False)
    print(f"pooled over {per['DMS_id'].nunique()} assays")
    print("assays per category:",
          per.drop_duplicates("DMS_id")["category"].value_counts().to_dict())

    # ---- figure: 2 panels (split), model_type lines, Structure bold ----
    palette = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
    x = np.arange(len(SUBSETS))
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for i, (split, _) in enumerate(SPLITS):
        ax = axes[i]
        sub = agg[agg.split == split]
        for t in mtypes:
            tt = sub[sub.model_type == t].set_index("subset").reindex([s for s, _ in SUBSETS])
            lw = 3.0 if t == "Structure" else 1.6
            ax.plot(x, tt["mean_spearman"].values, marker="o", color=palette[t],
                    lw=lw, label=t if i == 0 else None,
                    zorder=5 if t == "Structure" else 3)
        ax.axhline(0, color="gray", lw=0.6)
        ax.set_xticks(x); ax.set_xticklabels(["all", "top 50%", "top 25%", "top 10%"])
        ax.set_title(f"{split} fitness split", fontsize=12)
        ax.set_xlabel("fitness subset")
        if i == 0:
            ax.set_ylabel("mean Spearman vs experiment")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=6, fontsize=8, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("Winners-collapse by model_type — pooled across all assays "
                 "(weighted by #assays per category); Structure = bold", y=1.045, fontsize=13)
    plt.tight_layout()
    figdir = os.path.join(C.FIGURES, "persite_analysis")
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, "fig_collapse_modeltype_pooled.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"saved {out}")

    piv = (agg.pivot_table(index=["split", "model_type"], columns="subset",
                           values="mean_spearman")[["all", "top50", "top25", "top10"]])
    print(piv.round(3).to_string())


if __name__ == "__main__":
    main()
