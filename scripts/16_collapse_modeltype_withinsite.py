"""
Two extensions of the 'winners collapse' analysis (HXK4, both assays, sites >=5 subs):

(A) broken down by model_TYPE  -- do the structure-only models (biggest dissenters)
    collapse on the high-fitness tail differently from other types?
(B) with the fitness split defined WITHIN-SITE (percentile of a substitution among
    the substitutions at its own position) instead of globally. Then the "winners"
    are the best substitution at each position, spread over all positions, rather
    than the globally most-fit variants (which concentrate at tolerant positions).

For each assay x split x fitness-subset x model we compute Spearman(model, DMS)
over the subset, then average across the models of each model_type. Uses the
45-model family-best set.

Outputs:
  results/collapse_modeltype_withinsite.csv
  figures/fig_collapse_modeltype_withinsite.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import common as C
from importlib import import_module

fam = import_module("05_fig6_familybest")

ASSAYS = {"activity": "HXK4_HUMAN_Gersing_2022_activity",
          "abundance": "HXK4_HUMAN_Gersing_2023_abundance"}
SUBSETS = [("all", 0), ("top50", 50), ("top25", 75), ("top10", 90)]
MIN_SUBS = 5


def main():
    chosen, type_of, _ = fam.build_family_best()
    mtypes = [t for t in C.MODEL_TYPE_ORDER if t in set(type_of.values())]

    rows = []
    for key, dms_id in ASSAYS.items():
        df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
        df = df[df["mutant"].str.match(r"^[A-Z]\d+[A-Z]$")].copy()
        df["pos"] = df["mutant"].str.extract(r"(\d+)").astype(int)
        df = df[df.groupby("pos")["pos"].transform("size") >= MIN_SUBS]
        present = [m for m in chosen if m in df.columns]

        # within-site fitness percentile (0..100, higher = fitter at its own site)
        g = df.groupby("pos")["DMS_score"]
        r = g.rank(method="average")
        n = g.transform("size")
        df["ws_pct"] = (r - 1) / (n - 1) * 100
        df["glob_pct"] = df["DMS_score"].rank(pct=True) * 100

        y = df["DMS_score"].values
        for split, pctcol in [("global", "glob_pct"), ("within-site", "ws_pct")]:
            for sublab, q in SUBSETS:
                mask = df[pctcol].values >= q if q > 0 else np.ones(len(df), bool)
                if mask.sum() < 20:
                    continue
                # per-model spearman on the subset, grouped to model_type
                by_type = {t: [] for t in mtypes}
                for m in present:
                    s = df[m].values[mask]
                    if np.ptp(s) == 0:
                        continue
                    rho = spearmanr(y[mask], s).statistic
                    by_type[type_of[m]].append(rho)
                for t in mtypes:
                    if by_type[t]:
                        rows.append({"assay": key, "split": split, "subset": sublab,
                                     "model_type": t, "n_models": len(by_type[t]),
                                     "n_variants": int(mask.sum()),
                                     "mean_spearman": float(np.nanmean(by_type[t]))})
    res = pd.DataFrame(rows)
    res.to_csv(C.rp("collapse_modeltype_withinsite.csv"), index=False)

    # ---- console: structure vs rest collapse ----
    print("mean Spearman vs experiment by model_type, subset, split:")
    piv = res.pivot_table(index=["assay", "split", "model_type"], columns="subset",
                          values="mean_spearman")[["all", "top50", "top25", "top10"]]
    print(piv.round(3).to_string())

    # ---- figure: 2 rows (split) x 2 cols (assay), model_type lines ----
    palette = dict(zip(C.MODEL_TYPE_ORDER, plt.cm.tab10(np.linspace(0, 1, len(C.MODEL_TYPE_ORDER)))))
    x = np.arange(len(SUBSETS))
    xt = [s for s, _ in SUBSETS]
    splits = ["global", "within-site"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=True)
    for i, split in enumerate(splits):
        for j, key in enumerate(ASSAYS):
            ax = axes[i, j]
            sub = res[(res.split == split) & (res.assay == key)]
            for t in mtypes:
                tt = sub[sub.model_type == t].set_index("subset").reindex([s for s, _ in SUBSETS])
                lw = 3.0 if t == "Structure" else 1.6
                ax.plot(x, tt["mean_spearman"].values, marker="o", color=palette[t],
                        lw=lw, label=t if (i == 0 and j == 0) else None,
                        zorder=5 if t == "Structure" else 3)
            ax.axhline(0, color="gray", lw=0.6)
            ax.set_xticks(x); ax.set_xticklabels(["all", "top 50%", "top 25%", "top 10%"])
            ax.set_title(f"{key}  —  {split} fitness split", fontsize=11)
            if j == 0:
                ax.set_ylabel("mean Spearman vs experiment")
    axes[1, 0].set_xlabel("fitness subset"); axes[1, 1].set_xlabel("fitness subset")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=6, fontsize=8, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("Winners-collapse by model_type — global vs within-site fitness split "
                 "(Structure = bold)", y=1.035, fontsize=13)
    plt.tight_layout()
    out = os.path.join(C.FIGURES, "fig_collapse_modeltype_withinsite.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")

    # headline contrasts
    def cell(assay, split, mt, sub):
        v = res[(res.assay == assay) & (res.split == split) & (res.model_type == mt) & (res.subset == sub)]
        return v["mean_spearman"].iloc[0] if len(v) else np.nan
    print("\nHeadline (abundance):")
    for mt in ["Structure", "MSA", "Structure & MSA"]:
        print(f"  {mt:16s} global all->top25: {cell('abundance','global',mt,'all'):.2f}->{cell('abundance','global',mt,'top25'):.2f}"
              f"   within-site all->top25: {cell('abundance','within-site',mt,'all'):.2f}->{cell('abundance','within-site',mt,'top25'):.2f}")


if __name__ == "__main__":
    main()
