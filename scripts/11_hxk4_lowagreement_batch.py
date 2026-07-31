"""
Batch the W99-style disagreement analysis over ALL low-model-agreement sites of
HXK4_HUMAN, excluding the N-/C-terminal tail sites (L2, D4, R5, L460-L463).

Site set = union of the top-15 most-disagreed sites across both assays, minus the
tails -> T49, R63, W99, H156, Y214, Q287, R345, Y349, R377, V401, M402.

T49 and R63 are labelled as lying in a predicted disordered region.

Produces, into figures/site_disagreement/HXK4_HUMAN/:
  * fig_hxk4_HXK4_HUMAN_<assay>_disagreement.png       (assay overview, both assays)
  * fig_site_<wt><pos>_HXK4_HUMAN_<assay>_subset.png    (per-site deep dive, all sites x both assays)
  * site_<wt><pos>_<assay>_subset_rank_matrix.csv       (underlying rank matrices)
"""
import os
from importlib import import_module
import common as C

mod06 = import_module("06_hxk4_site_disagreement")
mod09 = import_module("09_w99_subset_figs")
fam = import_module("05_fig6_familybest")

ASSAYS = ["HXK4_HUMAN_Gersing_2022_activity", "HXK4_HUMAN_Gersing_2023_abundance"]
SITES = [49, 63, 99, 156, 214, 287, 345, 349, 377, 401, 402]
DISORDERED = {49, 63}
OUTDIR = os.path.join(C.FIGURES, "site_disagreement", "HXK4_HUMAN")


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # ---- assay-level overview figures (family-best 45-model set) ----
    chosen, type_of, _ = fam.build_family_best()
    mtypes = [t for t in C.MODEL_TYPE_ORDER if t in set(type_of.values())]
    for dms_id in ASSAYS:
        res, members_by_type = mod06.analyze(dms_id, chosen, type_of, mtypes)
        mod06.make_figure(dms_id, res, mtypes, members_by_type,
                          outdir=OUTDIR, disordered=DISORDERED)

    # ---- per-site deep dives (23-model subset) ----
    for dms_id in ASSAYS:
        for site in SITES:
            M, col_order, wt, stat, row_type = mod09.build(dms_id, site, results_dir=OUTDIR)
            note = (f"{wt}{site} lies in a predicted disordered region"
                    if site in DISORDERED else None)
            mod09.figure(dms_id, M, col_order, wt, stat, row_type,
                         site=site, outdir=OUTDIR, region_note=note)
            print(f"  {dms_id} {wt}{site}: median model-vs-DMS Spearman = {stat.median():.2f}"
                  + ("  [disordered]" if site in DISORDERED else ""))

    n_png = len([f for f in os.listdir(OUTDIR) if f.endswith(".png")])
    print(f"\nDone. {n_png} PNGs in {OUTDIR}")


if __name__ == "__main__":
    main()
