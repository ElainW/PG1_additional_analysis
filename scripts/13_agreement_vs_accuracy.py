"""
HXK4_HUMAN, both assays: does within-site model *consensus* predict within-site
model *accuracy* vs the experiment? (23-model subset.)

For every site (>=5 single substitutions) we compute, among the 23 subset models:
  * among_model_agreement = mean pairwise within-site Spearman between models
    (the same quantity used for the disagreement analysis, here on the 23 subset).
  * per-model within-site Spearman vs the experiment (DMS), and their mean.

We then split sites into agreement terciles (low / mid / high) and compare each
model's mean within-site Spearman-vs-experiment on HIGH- vs LOW-agreement sites,
and test the site-level relationship between agreement and accuracy.

Outputs:
  results/hxk4_agreement_vs_accuracy_sites_<assay>.csv
  results/hxk4_agreement_vs_accuracy_bymodel.csv
  figures/GCK_pilot_analysis/fig_hxk4_agreement_vs_accuracy_scatter.png
  figures/GCK_pilot_analysis/fig_hxk4_agreement_vs_accuracy_bymodel.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import common as C
from importlib import import_module

sns.set_style("whitegrid")
sub = import_module("04_subset_figs")

ASSAYS = {"activity": "HXK4_HUMAN_Gersing_2022_activity",
          "abundance": "HXK4_HUMAN_Gersing_2023_abundance"}
MIN_SUBS = 5


def per_site(dms_id):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
    df = df[df["mutant"].str.match(r"^[A-Z]\d+[A-Z]$")].copy()
    df["pos"] = df["mutant"].str.extract(r"(\d+)").astype(int)
    df["wt"] = df["mutant"].str[0]
    models = [c for c in sub.COLS if c in df.columns]

    site_rows, permodel_rows = [], []
    for pos, g in df.groupby("pos"):
        if len(g) < MIN_SUBS:
            continue
        R = g[models].rank().values                 # subs x models
        cc = np.corrcoef(R.T)
        off = cc[~np.eye(cc.shape[0], dtype=bool)]
        agree = float(np.nanmean(off))
        dms_rank = g["DMS_score"].rank().values
        model_vs_dms = {}
        for j, m in enumerate(models):
            rj = R[:, j]
            model_vs_dms[m] = (float(np.corrcoef(rj, dms_rank)[0, 1])
                               if rj.std() > 0 and dms_rank.std() > 0 else np.nan)
        site_rows.append({"pos": pos, "wt": g["wt"].iloc[0], "n_subs": len(g),
                          "among_model_agreement": agree,
                          "mean_model_vs_dms": float(np.nanmean(list(model_vs_dms.values())))})
        for m, v in model_vs_dms.items():
            permodel_rows.append({"pos": pos, "model": sub.COL2DISP[m], "spearman_vs_dms": v,
                                  "among_model_agreement": agree})
    sites = pd.DataFrame(site_rows)
    # agreement terciles
    sites["agreement_group"] = pd.qcut(sites["among_model_agreement"], 3,
                                       labels=["low", "mid", "high"])
    pm = pd.DataFrame(permodel_rows).merge(sites[["pos", "agreement_group"]], on="pos")
    return sites, pm, models


def main():
    all_sites = {}
    bymodel_rows = []
    headline = []
    for key, dms_id in ASSAYS.items():
        sites, pm, models = per_site(dms_id)
        sites.to_csv(C.rp(f"hxk4_agreement_vs_accuracy_sites_{key}.csv"), index=False)
        all_sites[key] = sites

        rho = spearmanr(sites["among_model_agreement"], sites["mean_model_vs_dms"]).statistic
        lo = sites[sites.agreement_group == "low"]["mean_model_vs_dms"].mean()
        hi = sites[sites.agreement_group == "high"]["mean_model_vs_dms"].mean()
        headline.append((key, len(sites), rho, lo, hi))

        for disp in [sub.COL2DISP[m] for m in models]:
            sub_pm = pm[pm.model == disp]
            row = {"assay": key, "model": disp,
                   "model_type": sub.TYPE[disp]}
            for grp in ["low", "high"]:
                row[grp] = sub_pm[sub_pm.agreement_group == grp]["spearman_vs_dms"].mean()
            row["high_minus_low"] = row["high"] - row["low"]
            bymodel_rows.append(row)

    bym = pd.DataFrame(bymodel_rows)
    bym.to_csv(C.rp("hxk4_agreement_vs_accuracy_bymodel.csv"), index=False)

    print("=== headline: does among-model agreement predict accuracy vs experiment? ===")
    print(f"{'assay':10s} {'n_sites':>7s} {'rho(agree,acc)':>15s} {'acc@LOW':>9s} {'acc@HIGH':>9s}")
    for k, n, rho, lo, hi in headline:
        print(f"{k:10s} {n:7d} {rho:15.3f} {lo:9.3f} {hi:9.3f}")
    print("\nmean model-vs-experiment Spearman by agreement group (over 23 models):")
    print(bym.groupby("assay")[["low", "high", "high_minus_low"]].mean().round(3).to_string())

    # ---------- Figure 1: site scatter ----------
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True)
    for ax, key in zip(axes, ASSAYS):
        s = all_sites[key]
        colors = {"low": "#d62728", "mid": "#999999", "high": "#1f77b4"}
        for grp in ["low", "mid", "high"]:
            gg = s[s.agreement_group == grp]
            ax.scatter(gg["among_model_agreement"], gg["mean_model_vs_dms"],
                       s=16, c=colors[grp], alpha=0.6, label=f"{grp} agreement", edgecolor="none")
        rho = spearmanr(s["among_model_agreement"], s["mean_model_vs_dms"]).statistic
        ax.axhline(0, color="gray", lw=0.6)
        ax.set_xlabel("among-model agreement at site\n(mean pairwise within-site Spearman)")
        ax.set_title(f"{key}   (Spearman of the two = {rho:.2f})", fontsize=11)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("mean model-vs-experiment\nwithin-site Spearman")
    fig.suptitle("HXK4_HUMAN: do sites where models agree with each other also match the experiment?",
                 fontsize=12)
    plt.tight_layout()
    figdir = os.path.join(C.FIGURES, "GCK_pilot_analysis")
    os.makedirs(figdir, exist_ok=True)
    out1 = os.path.join(figdir, "fig_hxk4_agreement_vs_accuracy_scatter.png")
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out1}")

    # ---------- Figure 2: per-model low vs high ----------
    order = (bym[bym.assay == "abundance"].sort_values("high", ascending=False)["model"].tolist())
    palette = dict(sub.PALETTE)
    fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
    for ax, key in zip(axes, ASSAYS):
        b = bym[bym.assay == key].set_index("model").reindex(order)
        yy = np.arange(len(b))
        colors = [palette[t] for t in b["model_type"]]
        ax.barh(yy - 0.2, b["high"], height=0.4, color=colors, label="high-agreement sites")
        ax.barh(yy + 0.2, b["low"], height=0.4, color=colors, alpha=0.4, hatch="//",
                label="low-agreement sites")
        ax.set_yticks(yy); ax.set_yticklabels(b.index, fontsize=8)
        for t, mt in zip(ax.get_yticklabels(), b["model_type"]):
            t.set_color(palette[mt])
        ax.invert_yaxis()
        ax.axvline(0, color="gray", lw=0.6)
        ax.set_xlabel("mean within-site Spearman vs experiment")
        ax.set_title(key, fontsize=11)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("HXK4_HUMAN: model accuracy vs experiment on HIGH- vs LOW-agreement sites (23-model subset)",
                 fontsize=12)
    plt.tight_layout()
    out2 = os.path.join(figdir, "fig_hxk4_agreement_vs_accuracy_bymodel.png")
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out2}")


if __name__ == "__main__":
    main()
