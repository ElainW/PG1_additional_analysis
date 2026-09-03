"""Condensed per-category leaderboard: the top Spearman- and top NDCG@10-ranked
model for each coarse_selection_type (paper Table T2, main-text cut).

Reads results/persite_analysis/family_best_category_leaderboard.csv (from `29`)
and, per category, keeps the rank-1 model by Spearman and the rank-1 model by
NDCG@10. `top_by` = 'both' when one model tops both metrics, else 'spearman' /
'ndcg'. Writes results/persite_analysis/family_best_category_leaderboard_top.csv.
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CAT_ORDER = ["Activity", "Binding", "Expression", "OrganismalFitness", "Stability"]
KEEP = ["coarse_selection_type", "family", "model", "model_type", "n_assays",
        "mean_spearman", "mean_ndcg10", "spearman_rank", "ndcg_rank"]


def main():
    src = os.path.join(ROOT, "results/persite_analysis/family_best_category_leaderboard.csv")
    d = pd.read_csv(src)
    # annotate model class (model_type) from the family-best selection, keyed by family
    fam = pd.read_csv(os.path.join(ROOT, "results/persite_analysis/family_best_selection.csv"))
    d = d.merge(fam[["family", "model_type"]].drop_duplicates("family"), on="family", how="left")
    rows = []
    for cat, g in d.groupby("coarse_selection_type"):
        top_sp = g.loc[g["spearman_rank"] == g["spearman_rank"].min()]
        top_nd = g.loc[g["ndcg_rank"] == g["ndcg_rank"].min()]
        sp_models, nd_models = set(top_sp["model"]), set(top_nd["model"])
        for _, r in top_sp.iterrows():
            rows.append({**r[KEEP].to_dict(),
                         "top_by": "both" if r["model"] in nd_models else "spearman"})
        for _, r in top_nd.iterrows():
            if r["model"] in sp_models:
                continue  # already emitted as 'both'
            rows.append({**r[KEEP].to_dict(), "top_by": "ndcg"})

    out = pd.DataFrame(rows)[["coarse_selection_type", "top_by", "model_type"] + [c for c in KEEP[1:] if c != "model_type"]]
    out["coarse_selection_type"] = pd.Categorical(out["coarse_selection_type"], CAT_ORDER, ordered=True)
    out = out.sort_values(["coarse_selection_type", "top_by"]).reset_index(drop=True)
    dst = os.path.join(ROOT, "results/persite_analysis/family_best_category_leaderboard_top.csv")
    out.to_csv(dst, index=False)
    print(f"wrote {dst} ({len(out)} rows)")
    same = out.groupby("coarse_selection_type", observed=True)["top_by"].apply(lambda s: "both" in set(s)).sum()
    print(f"categories where top-Spearman == top-NDCG model: {same}/{len(CAT_ORDER)}")


if __name__ == "__main__":
    main()
