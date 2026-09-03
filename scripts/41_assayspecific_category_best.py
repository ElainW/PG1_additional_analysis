"""
Assay-SPECIFIC high-disagreement sites: does category-tuned model choice help?

Question (H2 follow-up): for high-disagreement sites flagged in only ONE of a
protein's assays (NOT shared across assays; complement of script 32's shared set),
is per-site model accuracy better with the model that's best for that assay's
category? I.e. is assay-specific disagreement partly because some models predict an
assay category better?

Sites: per assay, script-27 high-disagreement positions MINUS the protein's shared
positions (colocalization_shared_sites.csv, >=2 assays), keeping positions with
>=5 missense substitutions. Proteins: PTEN, HXK4/GCK, CP2C9, P53 (7 assays have
assay-specific sites).

Per-site Spearman = spearman(model_score, DMS_score) over a site's subs (oriented).
Dots coloured by signed median(site DMS pct) - WT pct (WT oriented by
raw_DMS_directionality; red = below WT/deleterious). On every panel a green dotted
line marks the assay's own-category top model's TYPICAL median per-site Spearman
(over all >=5-sub positions across that category's assays) -- a baseline to judge
whether these hard sites sit below the model's usual within-site accuracy.

Three figures (figures/h2_persite_spearman_density/):
  fig_assayspecific_category_best_spearman_panel.pdf       5 boxes = top model per category by Spearman
  fig_assayspecific_category_best_ndcg_panel.pdf  5 boxes = top model per category by NDCG@10
  fig_assayspecific_modelclass_best_panel.pdf     6 boxes = best model per model_type (by Spearman)
Tables (results/MSA_conservation/H2/):
  h2_assayspecific_persite_spearman_by_site.csv           (category-Spearman)
  h2_assayspecific_persite_spearman_by_site_ndcg.csv      (category-NDCG@10)
  h2_assayspecific_persite_spearman_by_modelclass.csv     (per model class)
  h2_assayspecific_category_best_models.csv               (all model picks)
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import seaborn as sns
import common as C

sns.set_style("whitegrid")
rng = np.random.default_rng(0)

_spec = importlib.util.spec_from_file_location(
    "m38", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "38_h2_lowcons_persite_spearman_density.py"))
M38 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M38)

H2DIR = os.path.join(C.RESULTS, "MSA_conservation", "H2")
MIN_SUBS = 5
SHARED = os.path.join(C.RESULTS, "high_disagreement_colocalization",
                      "colocalization_shared_sites.csv")
HD_FILES = {"PTEN_HUMAN": "pten", "HXK4_HUMAN": "gck", "CP2C9_HUMAN": "cyp2c9",
            "P53_HUMAN": "p53"}
COL_ORDER = ["PTEN_HUMAN", "HXK4_HUMAN", "CP2C9_HUMAN", "P53_HUMAN"]

SHORT = {
    "MSA": "MSA", "Single sequence": "SeqOnly", "Structure": "Struct",
    "Single sequence & Structure": "Seq+Struct", "Structure & MSA": "Struct+MSA",
    "Single sequence, Structure & Function annotations": "Seq+Struct+Fn",
}
CMAP = cm.get_cmap("RdBu")

CATEGORIES5 = ["Expression", "Activity", "Stability", "OrganismalFitness", "Binding"]
CAT_ABBR = {"Expression": "Expression", "Activity": "Activity", "Stability": "Stability",
            "OrganismalFitness": "OrgFitness", "Binding": "Binding"}
MODEL_ABBR = {"AIDO.Protein-RAG-16B": "AIDO-RAG-16B", "TranceptEVE_L": "TranceptEVE_L"}

LB = pd.read_csv(C.rp("family_best_category_leaderboard.csv"))
MT = pd.read_csv(C.rp("directionality_check.csv"))
LB["model_type"] = LB.model.map(dict(zip(MT.model, MT.model_type)))
OVERALL_RHO = dict(zip(MT.model, MT.mean_spearman_vs_DMS))   # overall Spearman (217 assays)


def assay_cols(dms_id):
    return set(pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"), nrows=1).columns)


def top_model(categories, present_in, metric="mean_spearman"):
    """Single top model over category/categories by `metric`, present everywhere."""
    if isinstance(categories, str):
        categories = [categories]
    cand = LB[LB.coarse_selection_type.isin(categories)].sort_values(metric, ascending=False)
    for _, r in cand.iterrows():
        if all(r.model in assay_cols(a) for a in present_in):
            return r.model_type, r.model, r[metric], r.coarse_selection_type
    return None


def best_per_class(category, present_in, metric="mean_spearman"):
    """Best model of each model_type for `category` by `metric`, present everywhere."""
    out = []
    for mt in C.MODEL_TYPE_ORDER:
        cand = LB[(LB.model_type == mt) & (LB.coarse_selection_type == category)]
        for _, r in cand.sort_values(metric, ascending=False).iterrows():
            if all(r.model in assay_cols(a) for a in present_in):
                out.append((mt, r.model))
                break
    return out


def shared_positions():
    sh = pd.read_csv(SHARED)
    return {u: set(g.up_pos.astype(int)) for u, g in sh.groupby("UniProt_ID")}


def assay_specific_sites(uni, dms_id, shared):
    hd = pd.read_csv(os.path.join(C.RESULTS, "high_disagreement_sites",
                                  f"{HD_FILES[uni]}_high_disagreement.tsv"), sep="\t")
    hd = hd[hd.DMS_id == dms_id]
    shp = shared.get(uni, set())
    off = M38._offset(dms_id)
    df = M38._scores(dms_id)
    out = []
    for _, r in hd.iterrows():
        up = int(r.up_pos)
        if up in shp:
            continue
        sub = df[df.p == up - off]
        if len(sub) >= MIN_SUBS:
            out.append((up, r.wt_aa, sub))
    return df, out


def wt_pct_of(df, dms_id):
    wt_raw = M38.WT.get(dms_id, np.nan)
    ddir = M38.REF.loc[dms_id, "raw_DMS_directionality"]
    wt_dms = ddir * wt_raw
    return 100.0 * (df["DMS_score"] < wt_dms).mean() if pd.notna(wt_dms) else np.nan


def typical_median(model, assays):
    """Median per-site Spearman of `model` over ALL >=5-sub positions across `assays`."""
    vals = []
    for a in assays:
        df = M38._scores(a)
        for _, sub in df.groupby("p"):
            if len(sub) >= MIN_SUBS:
                r = M38.persite_spearman(sub, model)
                if pd.notna(r):
                    vals.append(r)
    return float(np.median(vals)) if vals else np.nan


def build(models_for_assay, group_col, out_csv):
    """models_for_assay(assay_category) -> [(group_value, model), ...]; one row per
    (assay-specific site x group_value)."""
    shared = shared_positions()
    rows = []
    for uni in COL_ORDER:
        hd = pd.read_csv(os.path.join(C.RESULTS, "high_disagreement_sites",
                                      f"{HD_FILES[uni]}_high_disagreement.tsv"), sep="\t")
        for dms_id in sorted(hd.DMS_id.unique()):
            acat = M38.REF.loc[dms_id, "coarse_selection_type"]
            df, sites = assay_specific_sites(uni, dms_id, shared)
            print(f"{dms_id} [{acat}] assay-specific sites: {len(sites)}")
            if not sites:
                continue
            wtp = wt_pct_of(df, dms_id)
            for up, wt, sub in sites:
                med_pct = sub["dms_pct"].median()
                diff = med_pct - wtp
                for gval, mdl in models_for_assay(acat):
                    rows.append({"UniProt_ID": uni, "DMS_id": dms_id, "assay_category": acat,
                                 "up_pos": up, "residue": f"{wt}{up}", "n_subs": len(sub),
                                 group_col: gval, "model": mdl,
                                 "persite_spearman": M38.persite_spearman(sub, mdl),
                                 "wt_dms_pct": round(wtp, 2) if pd.notna(wtp) else np.nan,
                                 "dms_pct_diff_from_wt": round(diff, 2)})
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(H2DIR, out_csv), index=False)
    return tab


def make_layout(tab):
    layout = {}
    for c, uni in enumerate(COL_ORDER):
        for r, dms in enumerate(sorted(tab[tab.UniProt_ID == uni].DMS_id.unique())):
            layout[(r, c)] = dms
    return layout


def _panels(df, layout, suptitle, out, group_col, group_order, group_palette,
            xtick, legend_handles, legend_title, baseline, bottom_text=None,
            legend_ncol=3, rect_bottom=0.09):
    M = float(df.dms_pct_diff_from_wt.abs().max())
    norm = Normalize(vmin=-M, vmax=M)
    nrows = max(r for r, _ in layout) + 1
    ncols = max(c for _, c in layout) + 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 5.0 * nrows),
                             sharey=True, squeeze=False)
    used = set()
    for (r, c), dms in layout.items():
        used.add((r, c))
        ax = axes[r, c]
        d = df[df.DMS_id == dms]
        acat = d.assay_category.iloc[0]
        n = d.up_pos.nunique()
        wtp = d.wt_dms_pct.iloc[0]
        sns.boxplot(data=d, x=group_col, y="persite_spearman", order=group_order,
                    hue=group_col, palette=group_palette, legend=False, ax=ax, width=0.6,
                    fliersize=0, linewidth=1, boxprops=dict(alpha=0.35), showcaps=True)
        for xi, g in enumerate(group_order):
            dd = d[d[group_col] == g]
            x = xi + rng.uniform(-0.16, 0.16, len(dd))
            ax.scatter(x, dd.persite_spearman, c=CMAP(norm(dd.dms_pct_diff_from_wt.values)),
                       s=42, edgecolor="black", linewidth=0.5, zorder=3)
        ax.axhline(0, color="0.4", ls="--", lw=1)
        bl = baseline.get(acat, np.nan)
        if pd.notna(bl):
            ax.axhline(bl, color="darkgreen", ls=":", lw=1.7)
            ax.text(len(group_order) - 0.5, bl + 0.02,
                    f"own-cat model typ. median = {bl:.2f}", ha="right", va="bottom",
                    fontsize=6.5, color="darkgreen")
        ax.set_ylim(-1.05, 1.05)
        ax.yaxis.grid(False)
        ax.set_title(f"{dms}\ncategory: {acat}   (n={n} sites; WT pct {wtp:.0f})", fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("per-site Spearman (model vs DMS)" if c == 0 else "")
        ax.set_xticks(range(len(group_order)))
        labs, rot, fs = xtick(group_order)
        ax.set_xticklabels(labs, rotation=rot, ha=("right" if rot else "center"), fontsize=fs)
        for tick, g in zip(ax.get_xticklabels(), group_order):
            if g == acat:                       # bold the assay's own category (category figs)
                tick.set_fontweight("bold")
    for r in range(nrows):
        for c in range(ncols):
            if (r, c) not in used:
                axes[r, c].axis("off")
    fig.suptitle(suptitle, fontsize=12, y=1.0)
    plt.tight_layout(rect=(0, rect_bottom, 0.9, 0.96))
    sm = cm.ScalarMappable(norm=norm, cmap=CMAP); sm.set_array([])
    cax = fig.add_axes([0.925, 0.34, 0.014, 0.40])
    fig.colorbar(sm, cax=cax).set_label(
        "median(site DMS percentile) − WT DMS percentile\n"
        "red = below WT (deleterious), white = at WT, blue = above WT", fontsize=9)
    fig.legend(handles=legend_handles, loc="lower center", ncol=legend_ncol, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, rect_bottom - 0.075), title=legend_title)
    if bottom_text:
        fig.text(0.5, 0.005, bottom_text, ha="center", va="bottom", fontsize=7,
                 family="monospace", color="0.2")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("saved", out)



def main():
    figdir = os.path.join(C.FIGURES, "h2_persite_spearman_density")
    os.makedirs(figdir, exist_ok=True)
    palette_cat = C.CATEGORY_COLORS
    palette_mt = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))

    # assays / categories present
    all_assays = sorted({d for uni in COL_ORDER for d in pd.read_csv(os.path.join(
        C.RESULTS, "high_disagreement_sites",
        f"{HD_FILES[uni]}_high_disagreement.tsv"), sep="\t").DMS_id.unique()})
    cat_assays = {}
    for uni in COL_ORDER:
        for d in pd.read_csv(os.path.join(C.RESULTS, "high_disagreement_sites",
                                          f"{HD_FILES[uni]}_high_disagreement.tsv"),
                             sep="\t").DMS_id.unique():
            cat_assays.setdefault(M38.REF.loc[d, "coarse_selection_type"], []).append(d)

    pick_rows = []

    def xtick_cat(order):
        return ([f"{CAT_ABBR[g]}\n{MODEL_ABBR.get(ct[g][1], ct[g][1])}" for g in order], 0, 7.5)

    def cat_legend(ct, metric_label):
        return [plt.Rectangle((0, 0), 1, 1, fc=palette_cat[g], alpha=0.45, ec="none",
                              label=f"{g}: {MODEL_ABBR.get(ct[g][1], ct[g][1])} "
                                    f"({SHORT[ct[g][0]]}, ρ{OVERALL_RHO[ct[g][1]]:.2f})")
                for g in CATEGORIES5]

    # ---- Figure A: category top by SPEARMAN -----------------------------------
    ct = {c: top_model(c, all_assays, "mean_spearman") for c in CATEGORIES5}
    base = {c: typical_median(top_model(c, cat_assays[c], "mean_spearman")[1], cat_assays[c])
            for c in cat_assays}
    for c, tm in ct.items():
        pick_rows.append({"scheme": "category-top-spearman", "category": c, "model_type": tm[0],
                          "model": tm[1], "metric_value": round(tm[2], 4)})
    tabA = build(lambda a: [(ec, ct[ec][1]) for ec in CATEGORIES5], "eval_category",
                 "h2_assayspecific_persite_spearman_by_site.csv")
    _panels(tabA, make_layout(tabA),
            "Assay-SPECIFIC high-disagreement sites — top model per assay category (by Spearman)",
            os.path.join(figdir, "fig_assayspecific_category_best_spearman_panel.pdf"),
            "eval_category", CATEGORIES5, {cc: palette_cat[cc] for cc in CATEGORIES5},
            xtick_cat, cat_legend(ct, "Spearman"),
            "box = top model for this assay category (model class, overall ρ)  ·  bold x-tick "
            "= own category  ·  green line = own-cat model's typical within-site median", base)

    # ---- Figure A alt 1: category top by NDCG@10 ------------------------------
    ctn = {c: top_model(c, all_assays, "mean_ndcg10") for c in CATEGORIES5}
    basen = {c: typical_median(top_model(c, cat_assays[c], "mean_ndcg10")[1], cat_assays[c])
             for c in cat_assays}
    for c, tm in ctn.items():
        pick_rows.append({"scheme": "category-top-ndcg10", "category": c, "model_type": tm[0],
                          "model": tm[1], "metric_value": round(tm[2], 4)})

    def xtick_cat_n(order):
        return ([f"{CAT_ABBR[g]}\n{MODEL_ABBR.get(ctn[g][1], ctn[g][1])}" for g in order], 0, 7.5)

    tabAn = build(lambda a: [(ec, ctn[ec][1]) for ec in CATEGORIES5], "eval_category",
                  "h2_assayspecific_persite_spearman_by_site_ndcg.csv")
    legn = [plt.Rectangle((0, 0), 1, 1, fc=palette_cat[g], alpha=0.45, ec="none",
                          label=f"{g}: {MODEL_ABBR.get(ctn[g][1], ctn[g][1])} "
                                f"({SHORT[ctn[g][0]]}, ρ{OVERALL_RHO[ctn[g][1]]:.2f})")
            for g in CATEGORIES5]
    _panels(tabAn, make_layout(tabAn),
            "Assay-SPECIFIC high-disagreement sites — top model per assay category (by NDCG@10)",
            os.path.join(figdir, "fig_assayspecific_category_best_ndcg_panel.pdf"),
            "eval_category", CATEGORIES5, {cc: palette_cat[cc] for cc in CATEGORIES5},
            xtick_cat_n, legn,
            "box = top-NDCG@10 model for this assay category (model class, overall ρ)  ·  bold "
            "x-tick = own category  ·  green line = own-cat model's typical within-site median",
            basen)

    # ---- Figure A alt 2: best model per MODEL CLASS (by Spearman), 7 assays ----
    order_mt = list(C.MODEL_TYPE_ORDER)
    cat_class = {c: best_per_class(c, all_assays, "mean_spearman") for c in cat_assays}

    def models_mc(acat):
        return cat_class[acat]                      # [(model_type, model), ...]

    for c, picks in cat_class.items():
        for mt, mdl in picks:
            pick_rows.append({"scheme": "modelclass-best-spearman", "category": c,
                              "model_type": mt, "model": mdl, "metric_value": np.nan})
    tabC = build(models_mc, "model_type", "h2_assayspecific_persite_spearman_by_modelclass.csv")

    def xtick_mc(order):
        return ([SHORT[t] for t in order], 40, 8)

    leg_mc = [plt.Rectangle((0, 0), 1, 1, fc=palette_mt[t], alpha=0.5, ec="none",
                            label=f"{SHORT[t]} = {t}") for t in order_mt]
    tbl = ["Best model per model class per category  (overall ρ over 217 assays):"]
    for c in ["Expression", "Activity", "OrganismalFitness", "Binding"]:
        if c in cat_class:
            parts = [f"{SHORT[mt]} {MODEL_ABBR.get(mdl, mdl)}({OVERALL_RHO[mdl]:.2f})"
                     for mt, mdl in cat_class[c]]
            tbl.append(f"{c+':':18s} " + "  ".join(parts))
    _panels(tabC, make_layout(tabC),
            "Assay-SPECIFIC high-disagreement sites — best model in each MODEL CLASS "
            "(by Spearman for the assay's category; 7 assays)",
            os.path.join(figdir, "fig_assayspecific_modelclass_best_panel.pdf"),
            "model_type", order_mt, {t: palette_mt[t] for t in order_mt},
            xtick_mc, leg_mc, "box colour = model_type class", base,
            bottom_text="\n".join(tbl), legend_ncol=3, rect_bottom=0.20)

    pd.DataFrame(pick_rows).to_csv(
        os.path.join(H2DIR, "h2_assayspecific_category_best_models.csv"), index=False)
    print("\nwrote model picks table")


if __name__ == "__main__":
    main()
