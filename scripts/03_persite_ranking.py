"""
Analysis 2: How does the ranking of mutations AT EACH SITE differ across models
and model types?

For every assay we keep single substitutions, group variants by residue position
(the "site"), and within each site convert every model's score to a normalized
within-site rank in [0,1] (0 = least fit at that site, 1 = most fit). This removes
the easy across-position signal and isolates "which substitution does the model
prefer at a given position".

Agreement metrics (per assay, then pooled across assays via Fisher-z):
  * model x model  : Pearson correlation of within-site normalized ranks
                     (pooled over all variants) = within-site ranking agreement.
  * model x DMS    : same, against the experiment's within-site ranking
                     => "site-controlled" recovery of the local fitness landscape.

Sites must have >=3 substitutions to contribute (need a non-trivial ranking).
Scores are used as-is (already higher=fitter; see 00_check_directionality.py),
so a positive value = the model ranks substitutions within a site the same way.

Inputs : O2/20250613/*.csv, config, reference
Outputs:
  results/persite_model_agreement_matrix.csv      (model x model)
  results/persite_modeltype_agreement_matrix.csv  (model_type x model_type)
  results/persite_dms_agreement_long.csv          (per assay: model within-site corr w/ DMS)
  results/persite_dms_agreement_by_modeltype_taxon.csv
  figures/fig5_persite_model_agreement_clustermap.png
  figures/fig6_persite_modeltype_agreement_heatmap.png
  figures/fig7_persite_dms_agreement_by_model.png
  figures/fig8_persite_dms_agreement_modeltype_taxon.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import common as C

sns.set_style("whitegrid")
MIN_SITE_VARIANTS = 3


def within_site_norm_ranks(df, cols, pos):
    """Return df of within-site normalized ranks in [0,1] for the given columns.
    Sites with < MIN_SITE_VARIANTS are dropped. NaNs preserved per column."""
    d = df[cols].copy()
    d["__pos__"] = pos
    sizes = d.groupby("__pos__")["__pos__"].transform("size")
    d = d[sizes >= MIN_SITE_VARIANTS]
    if d.empty:
        return d.drop(columns="__pos__")
    g = d.groupby("__pos__")
    out = {}
    for c in cols:
        r = g[c].rank(method="average")              # 1..k within site
        n = g[c].transform("count")                   # site size for non-NaN
        out[c] = (r - 1) / (n - 1).where(n > 1, np.nan)  # -> [0,1]
    return pd.DataFrame(out, index=d.index)


def main():
    models, mtype, _ = C.model_meta()
    tax = C.dms_to_taxon()
    files = C.list_score_files()
    cols_target = models + ["DMS_score"]
    idx = {c: i for i, c in enumerate(cols_target)}
    K = len(cols_target)

    # Fisher-z accumulators for the pooled correlation matrix
    z_sum = np.zeros((K, K))
    z_cnt = np.zeros((K, K))
    dms_records = []

    for path in files:
        dms_id = C.dms_id_from_path(path)
        df = pd.read_csv(path, usecols=lambda c: c in C.NON_MODEL_COLS or c in set(models))
        pos = df["mutant"].map(C.parse_position)
        keep = pos.notna()
        if keep.sum() == 0:
            print(f"{dms_id}: no single mutants, skipped")
            continue
        present = [m for m in models if m in df.columns]
        cols = present + ["DMS_score"]
        wr = within_site_norm_ranks(df[keep], cols, pos[keep])
        if wr.empty or len(wr) < 10:
            print(f"{dms_id}: too few multi-variant sites, skipped")
            continue

        cm = wr.corr(method="pearson", min_periods=10)  # within-site agreement
        # accumulate into global matrix (Fisher-z)
        for a in cols:
            ia = idx[a]
            for b in cols:
                ib = idx[b]
                v = cm.at[a, b]
                if pd.notna(v):
                    z_sum[ia, ib] += np.arctanh(np.clip(v, -0.999999, 0.999999))
                    z_cnt[ia, ib] += 1
        # per-model within-site agreement with DMS for this assay
        for m in present:
            v = cm.at[m, "DMS_score"]
            dms_records.append({
                "DMS_id": dms_id, "taxon": tax.get(dms_id, "Unknown"),
                "model": m, "model_type": mtype[m],
                "persite_corr_with_DMS": v,
                "n_sites_variants": len(wr),
            })
        print(f"{dms_id:45s} variants_in_multisites={len(wr):6d}")

    # ---- pooled model x model agreement matrix ----
    with np.errstate(invalid="ignore"):
        z_mean = np.where(z_cnt > 0, z_sum / np.maximum(z_cnt, 1), np.nan)
    corr = np.tanh(z_mean)
    cmat = pd.DataFrame(corr, index=cols_target, columns=cols_target)
    cmat.to_csv(C.rp("persite_model_agreement_matrix.csv"))

    # ---- model_type x model_type agreement ----
    mm = cmat.drop(index="DMS_score", columns="DMS_score")
    type_of = {m: mtype[m] for m in models}
    mtypes = [t for t in C.MODEL_TYPE_ORDER if t in set(type_of.values())]
    tt = pd.DataFrame(index=mtypes, columns=mtypes, dtype=float)
    for ta in mtypes:
        ma = [m for m in models if type_of[m] == ta]
        for tb in mtypes:
            mb = [m for m in models if type_of[m] == tb]
            block = mm.loc[ma, mb].values.astype(float)
            if ta == tb:  # exclude self-correlations (diagonal=1) for within-type
                iu = ~np.eye(len(ma), dtype=bool)
                vals = block[iu]
            else:
                vals = block.ravel()
            tt.at[ta, tb] = np.nanmean(vals)
    tt.to_csv(C.rp("persite_modeltype_agreement_matrix.csv"))

    # ---- per-model within-site DMS agreement (long + aggregates) ----
    dms = pd.DataFrame(dms_records)
    dms.to_csv(C.rp("persite_dms_agreement_long.csv"), index=False)
    taxa = [t for t in C.TAXA if t in dms["taxon"].unique()]
    by_mt_tax = (dms.groupby(["model_type", "taxon"])["persite_corr_with_DMS"].mean()
                   .unstack("taxon").reindex(index=mtypes, columns=taxa))
    by_mt_tax["overall"] = dms.groupby("model_type")["persite_corr_with_DMS"].mean().reindex(mtypes)
    by_mt_tax.to_csv(C.rp("persite_dms_agreement_by_modeltype_taxon.csv"))

    model_dms = (dms.groupby(["model", "model_type"])["persite_corr_with_DMS"]
                   .mean().reset_index().sort_values("persite_corr_with_DMS", ascending=False))
    model_dms.to_csv(C.rp("persite_dms_agreement_by_model.csv"), index=False)

    # =================== FIGURES ===================
    os.makedirs(C.FIGURES, exist_ok=True)

    # Fig 5: clustered model x model within-site agreement
    mm_f = mm.astype(float).fillna(0.0)
    type_codes = pd.Series({m: type_of[m] for m in mm_f.index})
    palette = dict(zip(mtypes, sns.color_palette("tab10", len(mtypes))))
    row_colors = type_codes.map(palette)
    g = sns.clustermap(mm_f, cmap="rocket", figsize=(20, 20),
                       row_colors=row_colors, col_colors=row_colors,
                       vmin=0, vmax=1, xticklabels=True, yticklabels=True,
                       cbar_kws={"label": "within-site rank agreement"})
    g.ax_heatmap.tick_params(labelsize=5)
    for label, color in palette.items():
        g.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    g.ax_col_dendrogram.legend(title="model_type", loc="center", ncol=3, fontsize=8)
    g.fig.suptitle("Within-site mutation-ranking agreement between models (clustered)", y=1.02)
    g.savefig(os.path.join(C.FIGURES, "fig5_persite_model_agreement_clustermap.png"), dpi=130, bbox_inches="tight")
    plt.close()

    # Fig 6: model_type x model_type agreement heatmap
    fig, ax = plt.subplots(figsize=(8, 6.5))
    sns.heatmap(tt.astype(float), cmap="rocket", annot=True, fmt=".2f", vmin=0, vmax=1,
                cbar_kws={"label": "mean within-site rank agreement"}, ax=ax)
    ax.set_title("Within-site mutation-ranking agreement between model TYPES\n(diagonal = within-type, self excluded)")
    lab = [t.replace(" & ", " &\n").replace(", ", ",\n") for t in mtypes]
    ax.set_xticklabels(lab, fontsize=7, rotation=30, ha="right")
    ax.set_yticklabels(lab, fontsize=7, rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(C.FIGURES, "fig6_persite_modeltype_agreement_heatmap.png"), dpi=150)
    plt.close()

    # Fig 7: per-model within-site agreement with DMS, colored by model_type
    md = model_dms.copy()
    fig, ax = plt.subplots(figsize=(9, 20))
    colors = md["model_type"].map(palette)
    ax.barh(md["model"], md["persite_corr_with_DMS"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("Mean within-site rank agreement with DMS")
    ax.set_title("How well each model ranks substitutions WITHIN a site (vs experiment)")
    ax.tick_params(axis="y", labelsize=5)
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[t]) for t in mtypes]
    ax.legend(handles, mtypes, title="model_type", fontsize=7, loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(C.FIGURES, "fig7_persite_dms_agreement_by_model.png"), dpi=150)
    plt.close()

    # Fig 8: model_type x taxon within-site DMS agreement
    fig, ax = plt.subplots(figsize=(7.5, 5))
    sns.heatmap(by_mt_tax[taxa].astype(float), cmap="mako", annot=True, fmt=".3f",
                cbar_kws={"label": "within-site agreement with DMS"}, ax=ax)
    ax.set_title("Within-site DMS agreement by model_type and taxon")
    ax.set_xlabel("Taxon"); ax.set_ylabel("")
    ax.set_yticklabels([t.replace(" & ", " &\n") for t in by_mt_tax.index], rotation=0, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(C.FIGURES, "fig8_persite_dms_agreement_modeltype_taxon.png"), dpi=150)
    plt.close()

    # ---- console summary ----
    print("\n=== model_type x model_type within-site ranking agreement ===")
    print(tt.astype(float).round(3).to_string())
    print("\n=== within-site DMS agreement by model_type x taxon ===")
    print(by_mt_tax.astype(float).round(3).to_string())
    print("\nTop 10 models by within-site DMS agreement:")
    print(model_dms.head(10).to_string(index=False))
    print("\nBottom 5 models by within-site DMS agreement:")
    print(model_dms.tail(5).to_string(index=False))
    print("\nSaved results + figures.")


if __name__ == "__main__":
    main()
