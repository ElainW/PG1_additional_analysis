"""
Diagnostic: for every assay-specific high-disagreement site used in
fig_assayspecific_rank_tiers_bycategory.pdf / fig_assayspecific_modelclass_best_panel.pdf
(scripts 44 / 41), plot the DMS-score spread of that site's substitutions
relative to the assay's wild-type score and to 0 -- a visual check for
degenerate/floor-effect sites (all subs at ~identical score, as found for
CCDB_ECOLI's flagged sites -- see script 29's within-site table) vs sites with
real spread.

Part 1 (unchanged): per-assay spread plots -> figures/h2_persite_spearman_density/check_noise/

Part 2 (noise classification): a site is flagged NOISE if its substitutions'
DMS-score standard deviation is small relative to the assay's overall missense
DMS-score standard deviation (relative_std = site_std / assay_std < REL_STD_THRESHOLD).
Threshold (0.5) was picked by eye against CP2C9_HUMAN_Amorosi_2021_abundance,
where the manual call was: S24/P73/R105/C435 = real signal, A238/P221 = noise.
Their relative_std values are 0.664/0.573/0.551/0.818 (real) vs 0.341/0.197
(noise) -- a threshold of 0.5 reproduces that call exactly and is applied
uniformly to all 6 assays' assay-specific sites.
Outputs -> results/high_disagreement_sites/assay_specific_no_noise/
  assay_specific_site_noise_classification.csv   (all 52 sites, audit trail)
  <DMS_id>_assay_specific_no_noise.tsv            (per assay, KEPT sites only)

Part 3 (denoised figures): the same two panel figures as scripts 41 (modelclass
best) and 44 (rank tiers by category), rebuilt using ONLY the non-noise sites
from Part 2, saved to figures/h2_persite_spearman_density/assay_specific_denoised/.
Reuses the pure helper functions from 41/38 (best_per_class, typical_median,
wt_pct_of, persite_spearman, _panels) but re-derives the input tables locally
from the filtered site set instead of calling assay_specific_sites() again, so
nothing in the original tracked scripts/results is touched or overwritten.
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

rng = np.random.default_rng(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "m41", os.path.join(_HERE, "41_assayspecific_category_best.py"))
M41 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M41)
M38 = M41.M38

REF = C.load_reference().set_index("DMS_id")
WT = M38.WT

FIGDIR = os.path.join(C.FIGURES, "h2_persite_spearman_density")
CHECK_NOISE_DIR = os.path.join(FIGDIR, "check_noise")
DENOISED_FIGDIR = os.path.join(FIGDIR, "assay_specific_denoised")
NOISE_RESULTS_DIR = os.path.join(C.RESULTS, "high_disagreement_sites", "assay_specific_no_noise")

REL_STD_THRESHOLD = 0.5


def wt_oriented(dms_id):
    """WT DMS score in the same higher=fitter orientation as the DMS_score column."""
    wt_raw = WT.get(dms_id, np.nan)
    ddir = REF.loc[dms_id, "raw_DMS_directionality"] if dms_id in REF.index else np.nan
    if pd.isna(wt_raw) or pd.isna(ddir):
        return np.nan
    return float(ddir) * float(wt_raw)


def collect_assay_specific_sites():
    """(uni_of, order, sites_of, df_of) for every assay with >=1 assay-specific site."""
    shared = M41.shared_positions()
    order, sites_of, df_of, uni_of = [], {}, {}, {}
    for uni in M41.COL_ORDER:
        hd = pd.read_csv(os.path.join(C.RESULTS, "high_disagreement_sites",
                                      f"{M41.HD_FILES[uni]}_high_disagreement.tsv"), sep="\t")
        for dms_id in sorted(hd.DMS_id.unique()):
            df, sites = M41.assay_specific_sites(uni, dms_id, shared)
            if sites:
                order.append(dms_id)
                sites_of[dms_id] = sites
                df_of[dms_id] = df
                uni_of[dms_id] = uni
    return uni_of, order, sites_of, df_of


# --------------------------------------------------------------------------
# Part 1: per-site DMS-score spread plots (one PNG per assay)
# --------------------------------------------------------------------------
def plot_assay_spread(dms_id, sites, outdir):
    wt_y = wt_oriented(dms_id)
    n = len(sites)
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.0 * nrows), squeeze=False)

    for i, (up, wt_aa, sub) in enumerate(sites):
        ax = axes[i // ncols][i % ncols]
        y = sub["DMS_score"].values
        x = rng.uniform(-0.15, 0.15, len(y))
        ax.scatter(x, y, s=32, color="#4C72B0", edgecolor="black", linewidth=0.4, zorder=3)
        ax.axhline(0, color="0.4", ls="--", lw=1.2, zorder=1, label="0")
        if pd.notna(wt_y):
            ax.axhline(wt_y, color="#C44E52", ls=":", lw=1.6, zorder=2, label="WT")
        if len(sub) <= 20:
            for xi, yi, mut in zip(x, y, sub["mutant"]):
                ax.annotate(mut, (xi, yi), fontsize=5.5, xytext=(3, 0),
                           textcoords="offset points", va="center")
        spread = float(np.ptp(y))
        ax.set_title(f"{wt_aa}{up}  (n={len(sub)}, range={spread:.3g})", fontsize=8)
        ax.set_xlim(-0.6, 0.9)
        ax.set_xticks([])
        if i % ncols == 0:
            ax.set_ylabel("DMS score", fontsize=8)
        ax.tick_params(axis="y", labelsize=7)

    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    handles = [plt.Line2D([0], [0], color="0.4", ls="--", lw=1.2, label="0"),
              plt.Line2D([0], [0], color="#C44E52", ls=":", lw=1.6, label="WT (oriented)"),
              plt.Line2D([0], [0], marker="o", color="#4C72B0", ls="none", ms=6,
                         mec="black", label="substitution")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False,
              bbox_to_anchor=(0.5, -0.02 / nrows))
    fig.suptitle(f"{dms_id}: DMS-score spread at assay-specific high-disagreement sites",
                fontsize=10, y=1.0)
    plt.tight_layout(rect=(0, 0.04, 1, 0.95))
    out = os.path.join(outdir, f"{dms_id}_site_dms_spread.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"saved {out}  ({n} sites)")


# --------------------------------------------------------------------------
# Part 2: noise classification
# --------------------------------------------------------------------------
def classify_sites(uni_of, order, sites_of, df_of):
    rows = []
    for dms_id in order:
        uni = uni_of[dms_id]
        df = df_of[dms_id]
        assay_std = df["DMS_score"].std()
        for up, wt_aa, sub in sites_of[dms_id]:
            site_std = sub["DMS_score"].std()
            rel = site_std / assay_std if assay_std > 0 else np.nan
            rows.append({
                "UniProt_ID": uni, "DMS_id": dms_id, "up_pos": up, "wt_aa": wt_aa,
                "n_subs": len(sub),
                "dms_range": round(float(sub["DMS_score"].max() - sub["DMS_score"].min()), 4),
                "site_std": round(float(site_std), 4),
                "assay_std": round(float(assay_std), 4),
                "relative_std": round(float(rel), 4) if pd.notna(rel) else np.nan,
                "is_noise": bool(rel < REL_STD_THRESHOLD) if pd.notna(rel) else True,
            })
    return pd.DataFrame(rows)


def write_noise_outputs(classification):
    os.makedirs(NOISE_RESULTS_DIR, exist_ok=True)
    classification.to_csv(
        os.path.join(NOISE_RESULTS_DIR, "assay_specific_site_noise_classification.csv"),
        index=False)
    print(f"saved {NOISE_RESULTS_DIR}/assay_specific_site_noise_classification.csv "
          f"({len(classification)} sites, {classification.is_noise.sum()} flagged noise, "
          f"{(~classification.is_noise).sum()} kept)")
    for dms_id, g in classification.groupby("DMS_id"):
        kept = g[~g.is_noise]
        out = os.path.join(NOISE_RESULTS_DIR, f"{dms_id}_assay_specific_no_noise.tsv")
        kept.to_csv(out, sep="\t", index=False)
        print(f"  {dms_id}: kept {len(kept)}/{len(g)} -> {os.path.basename(out)}")


def denoised_sites(order, sites_of, classification):
    """Filter sites_of down to non-noise sites; drop assays left with zero sites."""
    noisy = set(zip(classification[classification.is_noise].DMS_id,
                    classification[classification.is_noise].up_pos))
    order2, sites_of2 = [], {}
    for dms_id in order:
        kept = [(up, wt, sub) for up, wt, sub in sites_of[dms_id] if (dms_id, up) not in noisy]
        if kept:
            order2.append(dms_id)
            sites_of2[dms_id] = kept
    return order2, sites_of2


# --------------------------------------------------------------------------
# Part 3a: denoised fig_assayspecific_modelclass_best_panel.pdf (mirrors script 41)
# --------------------------------------------------------------------------
def build_from_sites(order, sites_of, df_of, uni_of, models_for_assay, group_col):
    rows = []
    for dms_id in order:
        acat = M38.REF.loc[dms_id, "coarse_selection_type"]
        df = df_of[dms_id]
        wtp = M41.wt_pct_of(df, dms_id)
        for up, wt, sub in sites_of[dms_id]:
            med_pct = sub["dms_pct"].median()
            diff = med_pct - wtp
            for gval, mdl in models_for_assay(acat):
                rows.append({"UniProt_ID": uni_of[dms_id], "DMS_id": dms_id,
                            "assay_category": acat,
                            "up_pos": up, "residue": f"{wt}{up}", "n_subs": len(sub),
                            group_col: gval, "model": mdl,
                            "persite_spearman": M38.persite_spearman(sub, mdl),
                            "wt_dms_pct": round(wtp, 2) if pd.notna(wtp) else np.nan,
                            "dms_pct_diff_from_wt": round(diff, 2)})
    return pd.DataFrame(rows)


def _panels_grid(df, order, nrows, ncols, suptitle, out, group_col, group_order,
                 group_palette, xtick, legend_handles, legend_title, bottom_text=None,
                 legend_ncol=3, rect_bottom=0.09):
    """Like 41._panels, but (a) lays panels out sequentially in a fixed nrows x ncols
    grid (no per-protein column grouping / gaps) and (b) draws NO baseline dotted line."""
    M = float(df.dms_pct_diff_from_wt.abs().max())
    norm = Normalize(vmin=-M, vmax=M)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 5.0 * nrows),
                             sharey=True, squeeze=False)
    for i, dms in enumerate(order):
        r, c = i // ncols, i % ncols
        ax = axes[r][c]
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
            ax.scatter(x, dd.persite_spearman, c=cm.get_cmap("RdBu")(norm(dd.dms_pct_diff_from_wt.values)),
                       s=42, edgecolor="black", linewidth=0.5, zorder=3)
        ax.axhline(0, color="0.4", ls="--", lw=1)
        ax.set_ylim(-1.05, 1.05)
        ax.yaxis.grid(False)
        ax.set_title(f"{dms}\ncategory: {acat}   (n={n} sites; WT pct {wtp:.0f})", fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("per-site Spearman (model vs DMS)" if c == 0 else "")
        ax.set_xticks(range(len(group_order)))
        labs, rot, fs = xtick(group_order)
        ax.set_xticklabels(labs, rotation=rot, ha=("right" if rot else "center"), fontsize=fs)
        for tick, g in zip(ax.get_xticklabels(), group_order):
            if g == acat:
                tick.set_fontweight("bold")
    for j in range(len(order), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(suptitle, fontsize=12, y=1.0)
    plt.tight_layout(rect=(0, rect_bottom, 0.9, 0.96))
    sm = cm.ScalarMappable(norm=norm, cmap=cm.get_cmap("RdBu")); sm.set_array([])
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


def make_modelclass_figure(order, sites_of, df_of, uni_of):
    os.makedirs(DENOISED_FIGDIR, exist_ok=True)
    palette_mt = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
    order_mt = list(C.MODEL_TYPE_ORDER)

    cat_assays = {}
    for d in order:
        cat_assays.setdefault(M38.REF.loc[d, "coarse_selection_type"], []).append(d)

    cat_class = {c: M41.best_per_class(c, order, "mean_spearman") for c in cat_assays}

    def models_mc(acat):
        return cat_class[acat]

    tabC = build_from_sites(order, sites_of, df_of, uni_of, models_mc, "model_type")

    def xtick_mc(ord_):
        return ([M41.SHORT[t] for t in ord_], 40, 8)

    leg_mc = [plt.Rectangle((0, 0), 1, 1, fc=palette_mt[t], alpha=0.5, ec="none",
                            label=f"{M41.SHORT[t]} = {t}") for t in order_mt]
    tbl = ["Best model per model class per category  (overall ρ over 217 assays), DENOISED sites:"]
    for c in ["Expression", "Activity", "OrganismalFitness", "Binding"]:
        if c in cat_class:
            parts = [f"{M41.SHORT[mt]} {M41.MODEL_ABBR.get(mdl, mdl)}({M41.OVERALL_RHO[mdl]:.2f})"
                    for mt, mdl in cat_class[c]]
            tbl.append(f"{c+':':18s} " + "  ".join(parts))

    out = os.path.join(DENOISED_FIGDIR, "fig_assayspecific_modelclass_best_panel.pdf")
    _panels_grid(tabC, order, 2, 3,
               "Assay-SPECIFIC high-disagreement sites (NOISE-FILTERED) — best model in each "
               "MODEL CLASS (by Spearman for the assay's category)",
               out, "model_type", order_mt, {t: palette_mt[t] for t in order_mt},
               xtick_mc, leg_mc,
               "box colour = model_type class",
               bottom_text="\n".join(tbl), legend_ncol=3, rect_bottom=0.20)


# --------------------------------------------------------------------------
# Part 3a: denoised fig_assayspecific_category_best_spearman_panel.pdf
# (mirrors script 41's "Figure A" exactly -- same style as the original,
# including the green baseline dotted line -- built on the noise-filtered site
# set. Saved into assay_specific_denoised/ alongside its modelclass/rank-tiers
# siblings.)
# --------------------------------------------------------------------------
def make_spearman_figure(order, sites_of, df_of, uni_of):
    os.makedirs(DENOISED_FIGDIR, exist_ok=True)
    palette_cat = C.CATEGORY_COLORS
    CATEGORIES5 = M41.CATEGORIES5

    cat_assays = {}
    for d in order:
        cat_assays.setdefault(M38.REF.loc[d, "coarse_selection_type"], []).append(d)

    ct = {c: M41.top_model(c, order, "mean_spearman") for c in CATEGORIES5}
    base = {c: M41.typical_median(M41.top_model(c, cat_assays[c], "mean_spearman")[1], cat_assays[c])
            for c in cat_assays}

    def models_for_assay(acat):
        return [(ec, ct[ec][1]) for ec in CATEGORIES5]

    tabA = build_from_sites(order, sites_of, df_of, uni_of, models_for_assay, "eval_category")

    def xtick_cat(ord_):
        return ([f"{M41.CAT_ABBR[g]}\n{M41.MODEL_ABBR.get(ct[g][1], ct[g][1])}" for g in ord_], 0, 7.5)

    leg = [plt.Rectangle((0, 0), 1, 1, fc=palette_cat[g], alpha=0.45, ec="none",
                        label=f"{g}: {M41.MODEL_ABBR.get(ct[g][1], ct[g][1])} "
                              f"({M41.SHORT[ct[g][0]]}, ρ{M41.OVERALL_RHO[ct[g][1]]:.2f})")
           for g in CATEGORIES5]

    out = os.path.join(DENOISED_FIGDIR, "fig_assayspecific_category_best_spearman_panel.pdf")
    M41._panels(tabA, M41.make_layout(tabA),
               "Assay-SPECIFIC high-disagreement sites (NOISE-FILTERED) — top model per "
               "assay category (by Spearman)",
               out, "eval_category", CATEGORIES5, {cc: palette_cat[cc] for cc in CATEGORIES5},
               xtick_cat, leg,
               "box = top model for this assay category (model class, overall ρ)  ·  bold "
               "x-tick = own category  ·  green line = own-cat model's typical within-site median",
               base)


# --------------------------------------------------------------------------
# Part 3a-bis: denoised fig_assayspecific_category_best_ndcg_panel.pdf
# (mirrors script 41's "Figure A alt 1" exactly -- same style as the original,
# baseline dotted line and per-protein-column layout included -- just built on
# the noise-filtered site set. Saved alongside the ORIGINAL in
# figures/h2_persite_spearman_density/ under a _denoised suffix so it doesn't
# overwrite the tracked file.)
# --------------------------------------------------------------------------
def make_ndcg_figure(order, sites_of, df_of, uni_of):
    palette_cat = C.CATEGORY_COLORS
    CATEGORIES5 = M41.CATEGORIES5

    cat_assays = {}
    for d in order:
        cat_assays.setdefault(M38.REF.loc[d, "coarse_selection_type"], []).append(d)

    ctn = {c: M41.top_model(c, order, "mean_ndcg10") for c in CATEGORIES5}

    def models_for_assay(acat):
        return [(ec, ctn[ec][1]) for ec in CATEGORIES5]

    tabAn = build_from_sites(order, sites_of, df_of, uni_of, models_for_assay, "eval_category")

    def xtick_cat_n(ord_):
        return ([f"{M41.CAT_ABBR[g]}\n{M41.MODEL_ABBR.get(ctn[g][1], ctn[g][1])}" for g in ord_], 0, 7.5)

    legn = [plt.Rectangle((0, 0), 1, 1, fc=palette_cat[g], alpha=0.45, ec="none",
                          label=f"{g}: {M41.MODEL_ABBR.get(ctn[g][1], ctn[g][1])} "
                                f"({M41.SHORT[ctn[g][0]]}, ρ{M41.OVERALL_RHO[ctn[g][1]]:.2f})")
            for g in CATEGORIES5]

    out = os.path.join(FIGDIR, "fig_assayspecific_category_best_ndcg_panel_denoised.pdf")
    M41._panels(tabAn, M41.make_layout(tabAn),
               "Assay-SPECIFIC high-disagreement sites (NOISE-FILTERED) — top model per "
               "assay category (by NDCG@10)",
               out, "eval_category", CATEGORIES5, {cc: palette_cat[cc] for cc in CATEGORIES5},
               xtick_cat_n, legn,
               "box = top-NDCG@10 model for this assay category (model class, overall ρ)  ·  "
               "bold x-tick = own category",
               {})


# --------------------------------------------------------------------------
# Part 3b: denoised fig_assayspecific_rank_tiers_bycategory.pdf (mirrors script 44)
# --------------------------------------------------------------------------
def make_rank_tiers_figure(order, sites_of, df_of):
    os.makedirs(DENOISED_FIGDIR, exist_ok=True)
    PS_DIR = os.path.join(C.RESULTS, "persite_analysis")
    lb = pd.read_csv(os.path.join(PS_DIR, "persite_category_leaderboard.csv"))
    MT = pd.read_csv(C.rp("directionality_check.csv"))
    M2T = dict(zip(MT.model, MT.model_type))
    PALETTE_MT = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
    tiers = ["top", "mid", "low"]
    rng2 = np.random.default_rng(0)

    cat_assays = {}
    for d in order:
        cat_assays.setdefault(REF.loc[d, "coarse_selection_type"], []).append(d)
    cat_picks = {}
    for cat, assays in cat_assays.items():
        present = set.intersection(*[set(df_of[a].columns) for a in assays])
        sub = lb[(lb.coarse_selection_type == cat) & (lb.model.isin(present))] \
                .sort_values("persite_rank")
        ranked = sub.model.tolist()
        if len(ranked) < 23:
            print(f"  WARNING: only {len(ranked)} ranked models present for {cat}; "
                  f"'mid'=rank23 tier not meaningful, using last available")
            mid_idx = len(ranked) // 2
        else:
            mid_idx = 22
        cat_picks[cat] = [("top", ranked[0]), ("mid", ranked[mid_idx]), ("low", ranked[-1])]

    rows = []
    for dms_id in order:
        cat = REF.loc[dms_id, "coarse_selection_type"]
        wtp = M41.wt_pct_of(df_of[dms_id], dms_id)
        for tier, m in cat_picks[cat]:
            catmed = lb[(lb.coarse_selection_type == cat) & (lb.model == m)] \
                       .mean_persite_spearman.iloc[0]
            for up, wt, sub in sites_of[dms_id]:
                med_pct = sub["dms_pct"].median()
                rows.append({"DMS_id": dms_id, "category": cat, "n_sites": len(sites_of[dms_id]),
                            "tier": tier, "model": m, "model_type": M2T[m],
                            "cat_persite_mean": round(catmed, 4),
                            "persite_spearman": M38.persite_spearman(sub, m),
                            "dms_pct_diff_from_wt": round(med_pct - wtp, 2)})
    tab = pd.DataFrame(rows)

    ncols = 4
    nrows = int(np.ceil(len(order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.8 * nrows),
                             sharey=True, squeeze=False)
    used_mt = []
    for i, dms_id in enumerate(order):
        ax = axes[i // ncols][i % ncols]
        d = tab[tab.DMS_id == dms_id]
        tier_mt = {t: d[d.tier == t].model_type.iloc[0] for t in tiers}
        tier_pal = {t: PALETTE_MT[tier_mt[t]] for t in tiers}
        for mt in tier_mt.values():
            if mt not in used_mt:
                used_mt.append(mt)
        sns.boxplot(data=d, x="tier", y="persite_spearman", order=tiers, hue="tier",
                    palette=tier_pal, legend=False, ax=ax, width=0.55, fliersize=0,
                    linewidth=1, boxprops=dict(alpha=0.35), showcaps=True)
        for xi, t in enumerate(tiers):
            dd = d[d.tier == t]
            x = xi + rng2.uniform(-0.15, 0.15, len(dd))
            ax.scatter(x, dd.persite_spearman, color=tier_pal[t], s=40, edgecolor="black",
                       linewidth=0.5, alpha=0.9, zorder=3)
            ax.hlines(dd.cat_persite_mean.iloc[0], xi - 0.32, xi + 0.32, color="black",
                      ls=":", lw=1.6, zorder=4)
        ax.axhline(0, color="0.4", ls="--", lw=1)
        cat = d.category.iloc[0]
        ax.set_title(f"{dms_id}\ncategory: {cat}  (n={d.n_sites.iloc[0]} denoised assay-specific sites)",
                     fontsize=9)
        ax.set_ylim(-1.05, 1.05)
        ax.yaxis.grid(False)
        ax.set_ylabel("per-site Spearman (model vs DMS)" if i % ncols == 0 else "")
        ax.set_xlabel("")
        heads = {"top": "Top #1", "mid": "Rank 23", "low": "Lowest"}
        labs = [f"{heads[t]}\n{M41.MODEL_ABBR.get(m, m)}\ncat med "
                f"{lb[(lb.coarse_selection_type==cat)&(lb.model==m)].mean_persite_spearman.iloc[0]:.2f}"
                for t, m in cat_picks[cat]]
        ax.set_xticks(range(3))
        ax.set_xticklabels(labs, fontsize=7)
    for j in range(len(order), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=PALETTE_MT[t], alpha=0.6, ec="none",
                             label=f"{M41.SHORT[t]} = {t}")
               for t in C.MODEL_TYPE_ORDER if t in used_mt]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False,
              bbox_to_anchor=(0.5, -0.03), title="box/dot colour = model_type class")
    fig.suptitle("Assay-specific high-disagreement sites (NOISE-FILTERED) — per-site Spearman "
                "for the assay CATEGORY's top / rank-23 / lowest per-site model  (45 family-best; "
                "dotted = model's category mean per-site Spearman)", fontsize=11, y=1.0)
    plt.tight_layout(rect=(0, 0.05, 1, 0.95))
    out = os.path.join(DENOISED_FIGDIR, "fig_assayspecific_rank_tiers_bycategory.pdf")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("saved", out)


def main():
    os.makedirs(CHECK_NOISE_DIR, exist_ok=True)
    uni_of, order, sites_of, df_of = collect_assay_specific_sites()
    print(f"{len(order)} assays with assay-specific sites")

    for dms_id in order:
        plot_assay_spread(dms_id, sites_of[dms_id], CHECK_NOISE_DIR)

    classification = classify_sites(uni_of, order, sites_of, df_of)
    write_noise_outputs(classification)

    order2, sites_of2 = denoised_sites(order, sites_of, classification)
    print(f"\n{len(order2)}/{len(order)} assays retain >=1 site after noise filtering")
    for d in order2:
        print(f"  {d}: {len(sites_of2[d])} sites kept")

    make_modelclass_figure(order2, sites_of2, df_of, uni_of)
    make_rank_tiers_figure(order2, sites_of2, df_of)
    make_ndcg_figure(order2, sites_of2, df_of, uni_of)
    make_spearman_figure(order2, sites_of2, df_of, uni_of)


if __name__ == "__main__":
    main()
