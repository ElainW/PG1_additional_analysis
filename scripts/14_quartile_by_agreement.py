"""
Fitness-quartile (high-fitness-variant) analysis, restricted to the within-site set
(HXK4_HUMAN both assays, sites with >=5 subs, 23 subset models), split by whether
the SITE is high- vs low-agreement (agreement terciles from 13_agreement_vs_accuracy).

Question: models lose ranking power on the top-fitness tail (the "poor prioritizer
of high-fitness variants" effect). Is that collapse different at high-agreement vs
low-agreement sites?

For each assay, agreement group (low / high), and model we compute Spearman(model,
DMS) over: all variants in the group, then the top-50/25/10% by DMS_score within
the group. We aggregate the mean across the 23 subset models.

Outputs:
  results/hxk4_quartile_by_agreement.csv
  figures/GCK_pilot_analysis/fig_hxk4_quartile_by_agreement.png
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

sub = import_module("04_subset_figs")

ASSAYS = {"activity": "HXK4_HUMAN_Gersing_2022_activity",
          "abundance": "HXK4_HUMAN_Gersing_2023_abundance"}
QUANTILES = [("overall", 0), ("top50", 50), ("top25", 75), ("top10", 90)]


def model_spearmans(y, s):
    out = {}
    for lab, q in QUANTILES:
        m = y >= np.percentile(y, q) if q > 0 else np.ones(len(y), bool)
        out[lab] = (spearmanr(y[m], s[m]).statistic
                    if m.sum() >= 20 and np.ptp(s[m]) > 0 else np.nan)
    return out


def main():
    rows = []
    permodel = []
    for key, dms_id in ASSAYS.items():
        sites = pd.read_csv(C.rp(f"hxk4_agreement_vs_accuracy_sites_{key}.csv"))
        grp_of = dict(zip(sites["pos"], sites["agreement_group"]))
        df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
        df = df[df["mutant"].str.match(r"^[A-Z]\d+[A-Z]$")].copy()
        df["pos"] = df["mutant"].str.extract(r"(\d+)").astype(int)
        df["grp"] = df["pos"].map(grp_of)
        df = df[df["grp"].isin(["low", "high"])]
        models = [c for c in sub.COLS if c in df.columns]

        for grp in ["low", "high"]:
            g = df[df["grp"] == grp]
            y = g["DMS_score"].values
            per = {lab: [] for lab, _ in QUANTILES}
            for m in models:
                sp = model_spearmans(y, g[m].values)
                for lab, _ in QUANTILES:
                    per[lab].append(sp[lab])
                permodel.append({"assay": key, "agreement": grp, "model": sub.COL2DISP[m], **sp})
            row = {"assay": key, "agreement": grp, "n_sites": (sites["agreement_group"] == grp).sum(),
                   "n_variants": len(g)}
            for lab, _ in QUANTILES:
                row[lab] = float(np.nanmean(per[lab]))
            row["collapse_overall_to_top25"] = row["overall"] - row["top25"]
            rows.append(row)

    tab = pd.DataFrame(rows)
    tab.to_csv(C.rp("hxk4_quartile_by_agreement.csv"), index=False)
    pd.DataFrame(permodel).to_csv(
        C.rp("hxk4_quartile_by_agreement_permodel.csv"), index=False)

    print("Mean Spearman vs experiment across 23 models, by fitness quartile subset:")
    print(tab[["assay", "agreement", "n_sites", "n_variants",
               "overall", "top50", "top25", "top10", "collapse_overall_to_top25"]].round(3).to_string(index=False))

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), sharey=True)
    xlabs = [lab for lab, _ in QUANTILES]
    x = np.arange(len(xlabs))
    style = {"high": ("#1f77b4", "-", "high-agreement sites"),
             "low": ("#d62728", "--", "low-agreement sites")}
    for ax, key in zip(axes, ASSAYS):
        for grp in ["high", "low"]:
            r = tab[(tab.assay == key) & (tab.agreement == grp)].iloc[0]
            col, ls, lab = style[grp]
            ax.plot(x, [r[l] for l in xlabs], marker="o", color=col, ls=ls, label=lab)
            for xi, l in zip(x, xlabs):
                ax.annotate(f"{r[l]:.2f}", (xi, r[l]), fontsize=7, ha="center",
                            va="bottom", xytext=(0, 4), textcoords="offset points", color=col)
        ax.axhline(0, color="gray", lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(["all", "top 50%", "top 25%", "top 10%"])
        ax.set_xlabel("fitness subset (by DMS_score)")
        ax.set_title(key, fontsize=11)
        ax.legend(fontsize=9)
    axes[0].set_ylabel("mean Spearman vs experiment\n(23 subset models)")
    fig.suptitle("HXK4_HUMAN: models fail to rank high-fitness variants on the high-fitness tail, high- vs low-agreement sites",
                 fontsize=12)
    plt.tight_layout()
    figdir = os.path.join(C.FIGURES, "GCK_pilot_analysis")
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, "fig_hxk4_quartile_by_agreement.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
