"""
Per-category model performance for the family-best subset, and a within-site test
of whether category-level rank predicts accuracy at high-disagreement sites.

Motivation
----------
family_best_selection.csv picks one model per family by OVERALL mean Spearman
(a single number). Here we instead break performance down by assay CATEGORY
(coarse_selection_type) and by TWO metrics -- Spearman and NDCG@10 -- then test a
hypothesis:

    At a high-disagreement site, if we pick the model that tops the leaderboard
    for that assay's category, do better-ranked models actually achieve higher
    within-site correlation with the DMS scores?

High-disagreement sites are taken as-is from results/high_disagreement_sites/
(script 27's output). Only sites with >= MIN_SITE_SUBS single substitutions
contribute to the within-site test.

Metric conventions
------------------
Spearman: spearmanr(model_score, DMS_score) over all variants in an assay
(matches ProteinGym's performance_DMS_benchmarks.py; scores already higher=fitter,
see 00_check_directionality.py). NDCG@10 uses ProteinGym's own calc_ndcg with its
default (quantile=True, top=10 -> gains from the top 10% of variants), copied
verbatim below so the number matches their reported benchmark.

Outputs (results/persite_analysis/):
  family_best_category_leaderboard.csv               (model x category, both metrics + ranks)
  family_best_highdisagreement_within_site.csv       (per site x model within-site corr + rank)
  family_best_highdisagreement_rank_vs_accuracy.csv  (per category hypothesis summary)
Figure:
  figures/fig9_highdisagreement_rank_vs_within_site.png
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import common as C
from importlib import import_module

fam = import_module("05_fig6_familybest")

MIN_SITE_SUBS = 5   # a site must have >= this many single substitutions to be tested


# --- ProteinGym NDCG, copied verbatim from proteingym/performance_DMS_benchmarks.py ---
def minmax(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x))


def calc_ndcg(y_true, y_score, **kwargs):
    """ProteinGym NDCG. Default quantile=True, top=10 => top 10% gains."""
    kwargs.setdefault("quantile", True)
    kwargs.setdefault("top", 10)
    if kwargs["quantile"]:
        k = np.floor(y_true.shape[0] * (kwargs["top"] / 100)).astype(int)
    else:
        k = kwargs["top"]
    if k < 1:
        return np.nan                       # too few variants for a top-k gain
    if isinstance(y_true, pd.Series):
        y_true = y_true.values
    if isinstance(y_score, pd.Series):
        y_score = y_score.values
    if np.ptp(y_true) == 0:
        return np.nan                       # constant target -> gains undefined
    gains = minmax(y_true)
    ranks = np.argsort(np.argsort(-y_score)) + 1
    ranks_k = ranks[ranks <= k]
    gains_k = gains[ranks <= k]
    ranks_fil = ranks_k[gains_k != 0]
    gains_fil = gains_k[gains_k != 0]
    if len(ranks_fil) == 0:
        return 0.0
    dcg = np.sum([g / np.log2(r + 1) for r, g in zip(ranks_fil, gains_fil)])
    ideal_ranks = np.argsort(np.argsort(-gains)) + 1
    ideal_ranks_k = ideal_ranks[ideal_ranks <= k]
    ideal_gains_k = gains[ideal_ranks <= k]
    ideal_ranks_fil = ideal_ranks_k[ideal_gains_k != 0]
    ideal_gains_fil = ideal_gains_k[ideal_gains_k != 0]
    idcg = np.sum([g / np.log2(r + 1) for r, g in zip(ideal_ranks_fil, ideal_gains_fil)])
    return dcg / idcg


def safe_spearman(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return spearmanr(a, b).statistic


# ---------------------------------------------------------------------------
# Part A: per-category leaderboard for the family-best subset
# ---------------------------------------------------------------------------
def build_leaderboard(chosen, fam_of, cst):
    recs = []
    for path in C.list_score_files():
        dms = C.dms_id_from_path(path)
        cat = cst.get(dms)
        if cat is None:
            continue
        df = pd.read_csv(path, usecols=lambda c: c in C.NON_MODEL_COLS or c in set(chosen))
        present = [m for m in chosen if m in df.columns]
        y = df["DMS_score"]
        for m in present:
            s = df[m]
            mask = y.notna() & s.notna()
            n = int(mask.sum())
            if n < 3:
                continue
            recs.append({
                "DMS_id": dms, "coarse_selection_type": cat,
                "family": fam_of[m], "model": m,
                "spearman": safe_spearman(s[mask], y[mask]),
                "ndcg10": calc_ndcg(y[mask].values, s[mask].values),
                "n_variants": n,
            })
    per_assay = pd.DataFrame(recs)

    lead = (per_assay.groupby(["coarse_selection_type", "family", "model"])
            .agg(n_assays=("DMS_id", "nunique"),
                 mean_spearman=("spearman", "mean"),
                 mean_ndcg10=("ndcg10", "mean"))
            .reset_index())
    grp = lead.groupby("coarse_selection_type")
    lead["spearman_rank"] = grp["mean_spearman"].rank(ascending=False, method="min").astype(int)
    lead["ndcg_rank"] = grp["mean_ndcg10"].rank(ascending=False, method="min").astype(int)
    lead = lead.sort_values(["coarse_selection_type", "spearman_rank"]).reset_index(drop=True)
    lead.to_csv(C.rp("family_best_category_leaderboard.csv"), index=False)
    return lead


# ---------------------------------------------------------------------------
# Part B: within-site accuracy at high-disagreement sites vs category rank
# ---------------------------------------------------------------------------
def load_high_disagreement_sites():
    hd = []
    for f in sorted(glob.glob(os.path.join(C.rp("high_disagreement_sites"), "*_high_disagreement.tsv"))):
        try:
            d = pd.read_csv(f, sep="\t")
        except pd.errors.EmptyDataError:
            continue                        # some proteins have no labelled sites
        if len(d):
            hd.append(d[["UniProt_ID", "DMS_id", "up_pos", "residue"]])
    return pd.concat(hd, ignore_index=True).drop_duplicates()


def build_within_site(chosen, fam_of, cst, lead):
    hd = load_high_disagreement_sites()
    rows = []
    for dms, grp in hd.groupby("DMS_id"):
        cat = cst.get(dms)
        path = os.path.join(C.SCORES_DIR, dms + ".csv")
        if cat is None or not os.path.exists(path):
            continue
        uni_id = grp["UniProt_ID"].iloc[0]
        df = pd.read_csv(path, usecols=lambda c: c in C.NON_MODEL_COLS or c in set(chosen))
        df["__pos__"] = df["mutant"].map(C.parse_position)   # NaN for multi-mutants
        present = [m for m in chosen if m in df.columns]
        res_of = dict(zip(grp["up_pos"], grp["residue"]))
        for pos in sorted(grp["up_pos"].unique()):
            sub = df[df["__pos__"] == pos]                   # single subs at this site
            if len(sub) < MIN_SITE_SUBS:
                continue
            y = sub["DMS_score"]
            for m in present:
                s = sub[m]
                mask = y.notna() & s.notna()
                if int(mask.sum()) < MIN_SITE_SUBS:
                    continue
                rows.append({
                    "DMS_id": dms, "UniProt_ID": uni_id, "coarse_selection_type": cat,
                    "up_pos": int(pos), "residue": res_of.get(pos, ""),
                    "family": fam_of[m], "model": m,
                    "n_subs": int(mask.sum()),
                    "within_site_spearman": safe_spearman(s[mask], y[mask]),
                })
    ws = pd.DataFrame(rows)
    ws = ws.merge(
        lead[["coarse_selection_type", "model", "spearman_rank", "ndcg_rank",
              "mean_spearman", "mean_ndcg10"]],
        on=["coarse_selection_type", "model"], how="left")
    ws = ws.dropna(subset=["within_site_spearman"])
    ws.to_csv(C.rp("family_best_highdisagreement_within_site.csv"), index=False)
    return ws


def summarize(ws):
    """Per category: does a better category rank => higher within-site correlation?"""
    out = []
    for cat, g in ws.groupby("coarse_selection_type"):
        n_sites = g.groupby(["DMS_id", "up_pos"]).ngroups
        oracle = g.groupby(["DMS_id", "up_pos"])["within_site_spearman"].max().mean()
        out.append({
            "coarse_selection_type": cat,
            "n_assays": g["DMS_id"].nunique(),
            "n_uniprot": g["UniProt_ID"].nunique(),
            "n_sites": n_sites,
            "n_site_model_pairs": len(g),
            # negative rho => lower rank number (better model) has higher within-site corr
            "rho_spearmanrank_vs_withinsite": safe_spearman(g["spearman_rank"], g["within_site_spearman"]),
            "rho_ndcgrank_vs_withinsite": safe_spearman(g["ndcg_rank"], g["within_site_spearman"]),
            "mean_withinsite_top1_spearman": g.loc[g["spearman_rank"] == 1, "within_site_spearman"].mean(),
            "mean_withinsite_top3_spearman": g.loc[g["spearman_rank"] <= 3, "within_site_spearman"].mean(),
            "mean_withinsite_top1_ndcg": g.loc[g["ndcg_rank"] == 1, "within_site_spearman"].mean(),
            "mean_withinsite_top3_ndcg": g.loc[g["ndcg_rank"] <= 3, "within_site_spearman"].mean(),
            "mean_withinsite_all_models": g["within_site_spearman"].mean(),
            "mean_withinsite_oracle_persite": oracle,
        })
    summ = pd.DataFrame(out).sort_values("coarse_selection_type").reset_index(drop=True)
    summ.to_csv(C.rp("family_best_highdisagreement_rank_vs_accuracy.csv"), index=False)
    return summ


def make_figure(ws):
    """Mean within-site Spearman by category-rank tier, per category, for both metrics."""
    tiers = [("Top 1", 1), ("Top 3", 3), ("Top 5", 5), ("Top 10", 10), ("All", 10**9)]
    cats = sorted(ws["coarse_selection_type"].unique())
    # per-category assay / UniProt counts for the legend labels
    counts = {c: (ws.loc[ws["coarse_selection_type"] == c, "DMS_id"].nunique(),
                  ws.loc[ws["coarse_selection_type"] == c, "UniProt_ID"].nunique())
              for c in cats}
    cat_label = {c: f"{c}  ({counts[c][0]} assays, {counts[c][1]} UniProt)" for c in cats}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    for ax, rank_col, title in [(axes[0], "spearman_rank", "ranked by category Spearman"),
                                (axes[1], "ndcg_rank", "ranked by category NDCG@10")]:
        width = 0.8 / len(cats)
        x = np.arange(len(tiers))
        for i, cat in enumerate(cats):
            g = ws[ws["coarse_selection_type"] == cat]
            vals = [g.loc[g[rank_col] <= thr, "within_site_spearman"].mean() for _, thr in tiers]
            ax.bar(x + i * width, vals, width, label=cat_label[cat])
        ax.set_xticks(x + width * (len(cats) - 1) / 2)
        ax.set_xticklabels([t for t, _ in tiers])
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_title(f"High-disagreement sites (>= {MIN_SITE_SUBS} subs)\n{title}", fontsize=11)
        ax.set_xlabel("models kept (by category rank)")
    axes[0].set_ylabel("mean within-site Spearman with DMS")
    axes[1].legend(title="assay category", fontsize=8)
    fig.suptitle("Do better-ranked models win at high-disagreement sites?", y=1.02, fontsize=13)
    plt.tight_layout()
    os.makedirs(C.FIGURES, exist_ok=True)
    out = os.path.join(C.FIGURES, "fig9_highdisagreement_rank_vs_within_site.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    return out


def main():
    chosen, _type_of, sel = fam.build_family_best()
    fam_of = dict(zip(sel["selected_model"], sel["family"]))
    ref = C.load_reference()
    cst = dict(zip(ref["DMS_id"], ref["coarse_selection_type"]))

    print(f"Family-best subset: {len(chosen)} models")
    lead = build_leaderboard(chosen, fam_of, cst)
    print("\n=== per-category leaderboard (top 5 by Spearman) ===")
    for cat, g in lead.groupby("coarse_selection_type"):
        print(f"\n{cat}:")
        print(g.nsmallest(5, "spearman_rank")[
            ["model", "n_assays", "mean_spearman", "spearman_rank", "mean_ndcg10", "ndcg_rank"]
        ].to_string(index=False))

    ws = build_within_site(chosen, fam_of, cst, lead)
    summ = summarize(ws)
    print("\n=== hypothesis summary (high-disagreement sites, >= "
          f"{MIN_SITE_SUBS} subs) ===")
    print(summ.round(3).to_string(index=False))
    print("\nInterpretation: rho_*rank_vs_withinsite < 0 means better-ranked "
          "models (lower rank #) tend to have higher within-site Spearman.")

    out = make_figure(ws)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
