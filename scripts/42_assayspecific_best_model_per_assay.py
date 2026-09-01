"""
H2 test: which single model has the best per-site Spearman on each assay's
assay-specific high-disagreement sites?

Pool: the 45 family-best models (family_best_selection.csv, script 05/29) -- the
same deduplicated set the rest of the H2 assay-specific analysis (41) uses; the
97-model directionality_check.csv only supplies model_type.
Sites: per assay, script-27 high-disagreement positions MINUS the protein's shared
positions (colocalization, >=2 assays), >=5 missense subs -- reused from script 41.
Metric: per model, the MEDIAN of its per-site Spearman over that assay's
assay-specific sites (median is robust to the small n; mean also written).

Outputs (results/MSA_conservation/H2/):
  h2_assayspecific_model_ranking_by_assay.csv  (all 45 models x 7 assays)
  h2_assayspecific_best_model_per_assay.csv    (winner per assay)
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import common as C

_spec = importlib.util.spec_from_file_location(
    "m41", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "41_assayspecific_category_best.py"))
M41 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M41)
M38 = M41.M38

H2DIR = M41.H2DIR

FB = pd.read_csv(os.path.join(C.RESULTS, "persite_analysis", "family_best_selection.csv"))
MODELS = list(FB.selected_model)                       # 45 family-best models
MT = pd.read_csv(C.rp("directionality_check.csv"))
M2T = dict(zip(MT.model, MT.model_type))
OVERALL_RHO = dict(zip(MT.model, MT.mean_spearman_vs_DMS))


def main():
    shared = M41.shared_positions()
    rows = []
    order = []                                          # assays with sites, in protein order
    for uni in M41.COL_ORDER:
        hd = pd.read_csv(os.path.join(C.RESULTS, "high_disagreement_sites",
                                      f"{M41.HD_FILES[uni]}_high_disagreement.tsv"), sep="\t")
        for dms_id in sorted(hd.DMS_id.unique()):
            df, sites = M41.assay_specific_sites(uni, dms_id, shared)
            if not sites:
                continue
            cat = M38.REF.loc[dms_id, "coarse_selection_type"]
            order.append(dms_id)
            present = [m for m in MODELS if m in df.columns]
            for m in present:
                rhos = [M38.persite_spearman(sub, m) for _, _, sub in sites]
                rhos = [r for r in rhos if pd.notna(r)]
                if not rhos:
                    continue
                rows.append({"UniProt_ID": uni, "DMS_id": dms_id, "category": cat,
                             "n_sites": len(sites), "model": m, "model_type": M2T[m],
                             "median_persite_spearman": round(float(np.median(rhos)), 4),
                             "mean_persite_spearman": round(float(np.mean(rhos)), 4),
                             "n_scored_sites": len(rhos),
                             "overall_rho": round(OVERALL_RHO[m], 4)})
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(H2DIR, "h2_assayspecific_model_ranking_by_assay.csv"), index=False)

    # winner per assay (by median, tie-break mean)
    best = (tab.sort_values(["median_persite_spearman", "mean_persite_spearman"],
                            ascending=False)
              .groupby("DMS_id", as_index=False).first())
    best = best.set_index("DMS_id").loc[order].reset_index()
    best.to_csv(os.path.join(H2DIR, "h2_assayspecific_best_model_per_assay.csv"), index=False)
    print("\nBest model per assay (by median per-site Spearman over assay-specific sites):")
    print(best[["DMS_id", "category", "n_sites", "model", "model_type",
                "median_persite_spearman", "overall_rho"]].to_string(index=False))


if __name__ == "__main__":
    main()
