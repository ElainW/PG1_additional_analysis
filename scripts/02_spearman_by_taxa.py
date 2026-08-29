"""
Analysis 1: Does each model's Spearman (score vs DMS_score) differ by taxon?
And what happens when models are grouped by model_type?

Inputs : results/spearman_long.csv
Outputs:
  results/spearman_model_by_taxon.csv        (model x taxon mean Spearman, wide)
  results/spearman_modeltype_by_taxon.csv    (model_type x taxon mean Spearman)
  results/kruskal_model_across_taxa.csv      (per-model Kruskal-Wallis across taxa)
  results/kruskal_modeltype_across_taxa.csv  (per-model_type Kruskal-Wallis)
  figures/fig1_model_by_taxon_heatmap.png
  figures/fig2_modeltype_by_taxon_box.png
  figures/fig3_modeltype_by_taxon_heatmap.png
  figures/fig4_taxon_overall_box.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import common as C

sns.set_style("whitegrid")


def main():
    df = pd.read_csv(C.rp("spearman_long.csv"))
    df = df.dropna(subset=["spearman"])
    taxa = [t for t in C.TAXA if t in df["taxon"].unique()]
    mtypes = [t for t in C.MODEL_TYPE_ORDER if t in df["model_type"].unique()]

    # ---- model x taxon mean Spearman (wide) ----
    mbt = (df.groupby(["model", "taxon"])["spearman"].mean()
             .unstack("taxon").reindex(columns=taxa))
    # order models by overall mean Spearman (desc)
    overall = df.groupby("model")["spearman"].mean().sort_values(ascending=False)
    mbt = mbt.reindex(overall.index)
    mtype_map = dict(zip(df["model"], df["model_type"]))
    mbt.insert(0, "model_type", [mtype_map[m] for m in mbt.index])
    mbt.insert(1, "overall_mean", overall.reindex(mbt.index).values)
    mbt.to_csv(C.rp("spearman_model_by_taxon.csv"))

    # ---- model_type x taxon mean Spearman (assay-model level mean) ----
    tbt = (df.groupby(["model_type", "taxon"])["spearman"].mean()
             .unstack("taxon").reindex(index=mtypes, columns=taxa))
    tbt["overall_mean"] = df.groupby("model_type")["spearman"].mean().reindex(mtypes).values
    tbt.to_csv(C.rp("spearman_modeltype_by_taxon.csv"))

    # ---- Kruskal-Wallis: does a model's per-assay Spearman differ across taxa? ----
    def kw_table(group_col, out_name):
        rows = []
        for key, g in df.groupby(group_col):
            samples = [g[g["taxon"] == t]["spearman"].values for t in taxa]
            samples = [s for s in samples if len(s) >= 2]
            if len(samples) >= 2:
                H, p = kruskal(*samples)
            else:
                H, p = np.nan, np.nan
            row = {group_col: key, "kruskal_H": H, "p_value": p}
            for t in taxa:
                row[f"mean_{t}"] = g[g["taxon"] == t]["spearman"].mean()
            row["range_across_taxa"] = (
                np.nanmax([row[f"mean_{t}"] for t in taxa]) -
                np.nanmin([row[f"mean_{t}"] for t in taxa]))
            rows.append(row)
        out = pd.DataFrame(rows).sort_values("p_value")
        # Benjamini-Hochberg FDR
        p = out["p_value"].values
        ok = ~np.isnan(p)
        q = np.full_like(p, np.nan, dtype=float)
        pv = p[ok]
        order = np.argsort(pv)
        m = len(pv)
        ranked = pv[order]
        qv = ranked * m / (np.arange(m) + 1)
        qv = np.minimum.accumulate(qv[::-1])[::-1]
        qtmp = np.empty(m)
        qtmp[order] = qv
        q[ok] = np.clip(qtmp, 0, 1)
        out["q_value_BH"] = q
        out.to_csv(C.rp(out_name), index=False)
        return out

    kw_model = kw_table("model", "kruskal_model_across_taxa.csv")
    kw_type = kw_table("model_type", "kruskal_modeltype_across_taxa.csv")

    # =================== FIGURES ===================
    os.makedirs(C.FIGURES, exist_ok=True)

    # Fig 1: model x taxon heatmap
    hm = mbt[taxa]
    fig, ax = plt.subplots(figsize=(7, 22))
    sns.heatmap(hm, cmap="viridis", annot=True, fmt=".2f", annot_kws={"size": 6},
                cbar_kws={"label": "Mean Spearman"}, ax=ax, vmin=0, vmax=0.7)
    ax.set_title("Mean Spearman (model vs DMS_score) by taxon\n(models ordered by overall mean)")
    ax.set_xlabel("Taxon"); ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(os.path.join(C.FIGURES, "fig1_model_by_taxon_heatmap.png"), dpi=150)
    plt.close()

    # Fig 2: model_type x taxon boxplots (distribution over models-within-type, per assay)
    fig, ax = plt.subplots(figsize=(13, 6))
    sns.boxplot(data=df[df["model_type"].isin(mtypes)], x="model_type", y="spearman",
                hue="taxon", order=mtypes, hue_order=taxa, showfliers=False, ax=ax)
    ax.set_title("Per-assay Spearman distribution by model_type and taxon")
    ax.set_xlabel(""); ax.set_ylabel("Spearman (model vs DMS_score)")
    ax.set_xticks(range(len(mtypes)))
    ax.set_xticklabels([t.replace(" & ", "\n& ").replace(", ", ",\n") for t in mtypes], fontsize=8)
    ax.legend(title="Taxon", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(C.FIGURES, "fig2_modeltype_by_taxon_box.png"), dpi=150)
    plt.close()

    # Fig 3: model_type x taxon heatmap
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.heatmap(tbt[taxa], cmap="viridis", annot=True, fmt=".3f",
                cbar_kws={"label": "Mean Spearman"}, ax=ax, vmin=0.2, vmax=0.55)
    ax.set_title("Mean Spearman by model_type and taxon")
    ax.set_xlabel("Taxon"); ax.set_ylabel("")
    ax.set_yticklabels([t.replace(" & ", " &\n") for t in tbt.index], rotation=0, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(C.FIGURES, "fig3_modeltype_by_taxon_heatmap.png"), dpi=150)
    plt.close()

    # Fig 4: overall per-assay Spearman by taxon (all models pooled)
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(data=df, x="taxon", y="spearman", order=taxa, showfliers=False, ax=ax,
                color="#69b3a2")
    sns.stripplot(data=df.groupby(["DMS_id", "taxon"])["spearman"].mean().reset_index(),
                  x="taxon", y="spearman", order=taxa, color="black", alpha=0.3, size=3, ax=ax)
    ax.set_title("Spearman by taxon (box=per model-assay; dots=per-assay mean over models)")
    ax.set_xlabel("Taxon"); ax.set_ylabel("Spearman")
    plt.tight_layout()
    plt.savefig(os.path.join(C.FIGURES, "fig4_taxon_overall_box.png"), dpi=150)
    plt.close()

    # ---- console summary ----
    print("=== Mean Spearman by model_type x taxon ===")
    print(tbt.round(3).to_string())
    print("\n=== Taxa overall (all models pooled) ===")
    print(df.groupby("taxon")["spearman"].agg(["mean", "median", "count"]).reindex(taxa).round(3).to_string())
    n_sig = (kw_model["q_value_BH"] < 0.05).sum()
    print(f"\nModels with significant taxon difference (Kruskal-Wallis, BH q<0.05): {n_sig}/{len(kw_model)}")
    print("\nTop 5 models most variable across taxa (by range of mean Spearman):")
    print(kw_model.sort_values("range_across_taxa", ascending=False)
          .head(5)[["model", "range_across_taxa", "p_value", "q_value_BH"]].to_string(index=False))
    print("\n=== model_type Kruskal-Wallis across taxa ===")
    print(kw_type[["model_type", "kruskal_H", "p_value", "q_value_BH", "range_across_taxa"]].round(4).to_string(index=False))
    print("\nSaved results + figures.")


if __name__ == "__main__":
    main()
