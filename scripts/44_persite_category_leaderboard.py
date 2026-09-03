"""
Per-site Spearman leaderboard over ALL 217 assays + a category-ranked tier plot.

Part 1 (leaderboard). For every one of the 217 DMS assays and each of the 45
family-best models, compute the model's per-site Spearman at each position with
>=5 missense subs, then take the MEAN over positions -> one value per (assay,
model). Averaging those over the assays in each coarse_selection_type gives the
per-site-Spearman category leaderboard (mirrors family_best_category_leaderboard
but with WITHIN-site instead of whole-assay Spearman).
Outputs (results/persite_analysis/):
  persite_spearman_by_assay.csv        (217 assays x 45 models: mean_persite_spearman)
  persite_category_leaderboard.csv     (category x model: mean_persite_spearman, rank)

Part 2 (figure, side-by-side companion to fig_assayspecific_rank_tiers.pdf). Same
3-box/strip layout, but the top / rank-23 / lowest models are chosen from the
CATEGORY per-site leaderboard (not the individual assay's own within-site rank).
Each box = the chosen model's per-site Spearman at that assay's assay-specific
high-disagreement sites; dotted line = the model's category mean per-site Spearman.
Output: figures/h2_persite_spearman_density/fig_assayspecific_rank_tiers_bycategory.pdf
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import common as C

sns.set_style("whitegrid")
rng = np.random.default_rng(0)

_spec = importlib.util.spec_from_file_location(
    "m41", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "41_assayspecific_category_best.py"))
M41 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M41)
M38 = M41.M38

PS_DIR = os.path.join(C.RESULTS, "persite_analysis")
MIN_SUBS = M41.MIN_SUBS
FB = pd.read_csv(os.path.join(PS_DIR, "family_best_selection.csv"))
MODELS = list(FB.selected_model)
MT = pd.read_csv(C.rp("directionality_check.csv"))
M2T = dict(zip(MT.model, MT.model_type))
OVERALL_RHO = dict(zip(MT.model, MT.mean_spearman_vs_DMS))
REF = C.load_reference().set_index("DMS_id")
PALETTE_MT = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
SHORT = M41.SHORT
MISSENSE = r'^[A-Z]\d+[A-Z]$'


def model_persite_mean(d, m):
    """mean per-site Spearman of model m over >=5-sub positions in frame d."""
    dd = d[["p", "DMS_score", m]].dropna()
    sizes = dd.groupby("p")["p"].transform("size")
    dd = dd[sizes >= MIN_SUBS]
    if dd.empty:
        return np.nan, 0
    n = dd.groupby("p")["p"].transform("size")
    ra = dd.groupby("p")[m].rank() - (n + 1) / 2
    rd = dd.groupby("p")["DMS_score"].rank() - (n + 1) / 2
    agg = pd.DataFrame({"p": dd.p, "num": (ra * rd).values,
                        "a2": (ra * ra).values, "b2": (rd * rd).values}).groupby("p").sum()
    rho = agg["num"] / np.sqrt(agg["a2"] * agg["b2"])
    rho = rho[np.isfinite(rho)]
    return (float(rho.mean()) if len(rho) else np.nan), int(len(rho))


def build_leaderboard():
    rows = []
    for path in C.list_score_files():
        dms_id = C.dms_id_from_path(path)
        if dms_id not in REF.index:
            continue
        cat = REF.loc[dms_id, "coarse_selection_type"]
        cols = ["mutant", "DMS_score"] + MODELS
        df = pd.read_csv(path, usecols=lambda c: c in cols)
        df = df[df.mutant.str.match(MISSENSE)].copy()
        df["p"] = df.mutant.str.extract(r'^[A-Z](\d+)[A-Z]$')[0].astype(int)
        df = df[df.mutant.str[0] != df.mutant.str[-1]]
        for m in [x for x in MODELS if x in df.columns]:
            mean_ps, npos = model_persite_mean(df, m)
            if npos:
                rows.append({"DMS_id": dms_id, "coarse_selection_type": cat, "model": m,
                             "model_type": M2T[m], "mean_persite_spearman": round(mean_ps, 4),
                             "n_positions": npos})
    by_assay = pd.DataFrame(rows)
    by_assay.to_csv(os.path.join(PS_DIR, "persite_spearman_by_assay.csv"), index=False)

    lb = (by_assay.groupby(["coarse_selection_type", "model", "model_type"], as_index=False)
                  .agg(mean_persite_spearman=("mean_persite_spearman", "mean"),
                       n_assays=("DMS_id", "nunique")))
    lb["mean_persite_spearman"] = lb["mean_persite_spearman"].round(4)
    lb["persite_rank"] = lb.groupby("coarse_selection_type")["mean_persite_spearman"] \
                           .rank(ascending=False, method="first").astype(int)
    lb = lb.sort_values(["coarse_selection_type", "persite_rank"])
    lb.to_csv(os.path.join(PS_DIR, "persite_category_leaderboard.csv"), index=False)
    return by_assay, lb


def figure(lb):
    shared = M41.shared_positions()
    # assays with assay-specific sites, in protein order
    order, sites_of, df_of = [], {}, {}
    for uni in M41.COL_ORDER:
        hd = pd.read_csv(os.path.join(C.RESULTS, "high_disagreement_sites",
                                      f"{M41.HD_FILES[uni]}_high_disagreement.tsv"), sep="\t")
        for dms_id in sorted(hd.DMS_id.unique()):
            df, sites = M41.assay_specific_sites(uni, dms_id, shared)
            if sites:
                order.append(dms_id); sites_of[dms_id] = sites; df_of[dms_id] = df

    # category -> our assays of that category; tier picks from the category leaderboard
    cat_assays = {}
    for d in order:
        cat_assays.setdefault(REF.loc[d, "coarse_selection_type"], []).append(d)
    cat_picks = {}
    for cat, assays in cat_assays.items():
        present = set.intersection(*[set(df_of[a].columns) for a in assays])
        sub = lb[(lb.coarse_selection_type == cat) & (lb.model.isin(present))] \
                .sort_values("persite_rank")
        ranked = sub.model.tolist()
        cat_picks[cat] = [("top", ranked[0]), ("mid", ranked[22]), ("low", ranked[-1])]

    tiers = ["top", "mid", "low"]
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
            x = xi + rng.uniform(-0.15, 0.15, len(dd))
            ax.scatter(x, dd.persite_spearman, color=tier_pal[t], s=40, edgecolor="black",
                       linewidth=0.5, alpha=0.9, zorder=3)
            ax.hlines(dd.cat_persite_mean.iloc[0], xi - 0.32, xi + 0.32, color="black",
                      ls=":", lw=1.6, zorder=4)
        ax.axhline(0, color="0.4", ls="--", lw=1)
        cat = d.category.iloc[0]
        ax.set_title(f"{dms_id}\ncategory: {cat}  (n={d.n_sites.iloc[0]} assay-specific sites)",
                     fontsize=9)
        ax.set_ylim(-1.05, 1.05)
        ax.yaxis.grid(False)
        ax.set_ylabel("per-site Spearman (model vs DMS)" if i % ncols == 0 else "")
        ax.set_xlabel("")
        heads = {"top": "Top #1", "mid": "Rank 23", "low": "Lowest"}
        labs = [f"{heads[t]}\n{M41.MODEL_ABBR.get(m, m)}\ncat med {lb[(lb.coarse_selection_type==cat)&(lb.model==m)].mean_persite_spearman.iloc[0]:.2f}"
                for t, m in cat_picks[cat]]
        ax.set_xticks(range(3))
        ax.set_xticklabels(labs, fontsize=7)
    for j in range(len(order), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=PALETTE_MT[t], alpha=0.6, ec="none",
                             label=f"{SHORT[t]} = {t}")
               for t in C.MODEL_TYPE_ORDER if t in used_mt]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.03), title="box/dot colour = model_type class")
    fig.suptitle("Assay-specific high-disagreement sites — per-site Spearman for the assay "
                 "CATEGORY's top / rank-23 / lowest per-site model  (45 family-best; dotted = "
                 "model's category mean per-site Spearman)", fontsize=11, y=1.0)
    plt.tight_layout(rect=(0, 0.05, 1, 0.95))
    out = os.path.join(C.FIGURES, "h2_persite_spearman_density",
                       "fig_assayspecific_rank_tiers_bycategory.pdf")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("saved", out)
    print("\ncategory tier picks (per-site leaderboard):")
    for cat, picks in cat_picks.items():
        print(f"  {cat}: " + " | ".join(
            f"{t} {m} ({lb[(lb.coarse_selection_type==cat)&(lb.model==m)].mean_persite_spearman.iloc[0]:.2f})"
            for t, m in picks))


def main():
    print("computing per-site leaderboard over 217 assays x 45 models ...")
    _, lb = build_leaderboard()
    print("top per-site model per category:")
    print(lb[lb.persite_rank == 1][["coarse_selection_type", "model", "model_type",
                                    "mean_persite_spearman", "n_assays"]].to_string(index=False))
    figure(lb)


if __name__ == "__main__":
    main()
