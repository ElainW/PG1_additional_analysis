"""
3VEV vs 1V4S using each structure's OWN native ligands (no ligand transfer).

Both are "glucokinase + glucose + activator", closed/active:
  1V4S : glucose GLC + activator MRK  (res 14-461)
  3VEV : glucose GLC + activator 0H4  (res 4-458; loop 94-98 disordered)
The activator compounds differ (MRK vs 0H4) but occupy the same allosteric site.

For each structure we compute, per residue, the minimum heavy-atom distance to
its nearest native ligand (glucose or activator), then compare the two
independent measurements on shared residues and re-check the disagreement links.

Outputs:
  results/gck_native_3vev_vs_1v4s_perresidue.csv
  results/gck_native_comparison_summary.csv
  figures/fig_gck_native_3vev_vs_1v4s.png
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser
import common as C

warnings.simplefilter("ignore")
STRUCTS = {
    "1V4S": (os.path.join(C.ANALYSIS, "data", "1V4S.pdb"), "GLC", "MRK"),
    "3VEV": (os.path.join(C.ANALYSIS, "data", "3VEV.pdb"), "GLC", "0H4"),
}
CURATED = [49, 63, 99, 156, 214, 287, 345, 349, 377, 401, 402]


def residue_ligand_dist(path, glc_name, act_name):
    ch = PDBParser(QUIET=True).get_structure("s", path)[0]["A"]
    heavy, ca, glc, act = {}, {}, [], []
    for r in ch:
        if r.id[0] == " ":
            n = r.id[1]
            heavy[n] = np.array([a.coord for a in r if a.element != "H"])
            if "CA" in r:
                ca[n] = r["CA"].coord
        elif r.resname == glc_name:
            glc += [a.coord for a in r if a.element != "H"]
        elif r.resname == act_name:
            act += [a.coord for a in r if a.element != "H"]
    glc = np.array(glc); act = np.array(act)
    rows = []
    for n, atoms in heavy.items():
        dg = cdist(atoms, glc).min() if len(glc) else np.nan
        da = cdist(atoms, act).min() if len(act) else np.nan
        rows.append({"position": n, "dist_glucose": dg, "dist_activator": da,
                     "dist_nearest": np.nanmin([dg, da])})
    return pd.DataFrame(rows), ca


def main():
    (d1, ca1), (d3, ca3) = (residue_ligand_dist(*STRUCTS["1V4S"]),
                            residue_ligand_dist(*STRUCTS["3VEV"]))
    m = d1.merge(d3, on="position", suffixes=("_1V4S", "_3VEV"))
    m["abs_diff_nearest"] = (m["dist_nearest_1V4S"] - m["dist_nearest_3VEV"]).abs()
    m.to_csv(os.path.join(C.RESULTS, "gck_native_3vev_vs_1v4s_perresidue.csv"), index=False)

    rho, _ = spearmanr(m["dist_nearest_1V4S"], m["dist_nearest_3VEV"])
    pear = np.corrcoef(m["dist_nearest_1V4S"], m["dist_nearest_3VEV"])[0, 1]
    summ = [
        {"metric": "shared_residues", "value": len(m)},
        {"metric": "nearest_dist_Pearson", "value": round(pear, 3)},
        {"metric": "nearest_dist_Spearman", "value": round(rho, 3)},
        {"metric": "median_abs_diff_nearest_A", "value": round(m["abs_diff_nearest"].median(), 2)},
        {"metric": "glucose_dist_Pearson",
         "value": round(np.corrcoef(m["dist_glucose_1V4S"], m["dist_glucose_3VEV"])[0, 1], 3)},
        {"metric": "activator_dist_Pearson",
         "value": round(np.corrcoef(m["dist_activator_1V4S"], m["dist_activator_3VEV"])[0, 1], 3)},
    ]

    # re-check disagreement <-> distance and pocket enrichment, per structure
    for key, dms_id in [("activity", "HXK4_HUMAN_Gersing_2022_activity"),
                        ("abundance", "HXK4_HUMAN_Gersing_2023_abundance")]:
        dis = pd.read_csv(os.path.join(C.RESULTS, f"hxk4_{dms_id}_site_disagreement.csv"))
        mm = m.merge(dis[["position", "disagreement"]], on="position")
        for st in ["1V4S", "3VEV"]:
            r_, p_ = spearmanr(mm["disagreement"], mm[f"dist_nearest_{st}"])
            summ.append({"metric": f"{key}_spearman_disagree_vs_dist_{st}",
                         "value": f"{r_:.3f} (p={p_:.2f})"})
            base = (mm[f"dist_nearest_{st}"] <= 8).mean()
            top = mm.sort_values("disagreement", ascending=False).head(25)
            frac = (top[f"dist_nearest_{st}"] <= 8).mean()
            summ.append({"metric": f"{key}_{st}_top25_within8A_vs_base",
                         "value": f"{frac:.0%} vs {base:.0%}"})

    S = pd.DataFrame(summ)
    S.to_csv(os.path.join(C.RESULTS, "gck_native_comparison_summary.csv"), index=False)
    print("=== native 3VEV vs 1V4S (each with its own ligands) ===")
    print(S.to_string(index=False))

    hot = m[m["position"].isin(CURATED)]
    print("\n=== curated hotspots: dist to glucose / activator (Å), native each ===")
    print(hot[["position", "dist_glucose_1V4S", "dist_glucose_3VEV",
               "dist_activator_1V4S", "dist_activator_3VEV"]].round(1).to_string(index=False))

    # figure
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.6))
    ax[0].scatter(m["dist_nearest_1V4S"], m["dist_nearest_3VEV"], s=14, alpha=0.6, color="#4477aa")
    lim = [0, max(m[["dist_nearest_1V4S", "dist_nearest_3VEV"]].max()) + 2]
    ax[0].plot(lim, lim, "k--", lw=1, label="y = x")
    ax[0].set_xlabel("nearest-ligand dist — 1V4S native (Å)")
    ax[0].set_ylabel("nearest-ligand dist — 3VEV native (Å)")
    ax[0].set_title(f"Native per-residue distances\nPearson={pear:.3f}, Spearman={rho:.3f}, "
                    f"median |Δ|={m['abs_diff_nearest'].median():.2f} Å")
    ax[0].legend(fontsize=8)
    hb = hot.set_index("position")
    x = np.arange(len(hb))
    ax[1].bar(x - 0.2, hb["dist_activator_1V4S"], width=0.4, label="1V4S (MRK)", color="#4477aa")
    ax[1].bar(x + 0.2, hb["dist_activator_3VEV"], width=0.4, label="3VEV (0H4)", color="#ee6677")
    ax[1].set_xticks(x); ax[1].set_xticklabels(hb.index, fontsize=8)
    ax[1].set_xlabel("curated hotspot"); ax[1].set_ylabel("dist to ACTIVATOR (Å)")
    ax[1].set_title("Distance to (different) activator, by structure")
    ax[1].legend(fontsize=8)
    plt.tight_layout()
    out = os.path.join(C.FIGURES, "fig_gck_native_3vev_vs_1v4s.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
