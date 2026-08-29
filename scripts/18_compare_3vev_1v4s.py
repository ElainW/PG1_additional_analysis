"""
Does the structure choice (AF model vs 1V4S) change the GCK 3-D distance measurements?

AFmodel = data/structures/af2/HXK4_HUMAN.pdb : apo AlphaFold model of P35557 (res 1-465), chain A,
          NO ligands (its B-factor column is pLDDT, r=0.98 vs the AlphaFold DB).
1V4S    = data/structures/pdb/1V4S.pdb       : res 14-461, glucose (GLC) + activator (MRK) bound.
NOTE: for historical reasons the variables/outputs below tag the AF model "3VEV"
(e.g. dist_3VEV, fig_gck_3vev_vs_1v4s.png). It is NOT the 3VEV crystal -- the real
3VEV crystal (with native glucose + activator 0H4) is used separately in script 20.

Both are closed/active-like (global Calpha RMSD 1.70 A over 448 shared residues).
Because the AF model has no ligands, we superpose 1V4S onto it on shared Calpha and
transform the two ligands into the AF-model frame, then recompute per-residue minimum
heavy-atom distance to the nearest ligand exactly as script 17 did for 1V4S.

We report:
  * global and per-subdomain Calpha RMSD (is there a hidden domain rotation?)
  * per-residue distance-to-ligand: AF-model-frame vs 1V4S-native (corr, mean abs diff)
  * do the disagreement<->distance correlations change?
  * per-hotspot distances side by side
  * a LIGAND-FREE check: pairwise Calpha distances among curated hotspots in each
    structure (this needs no ligand transfer at all).

Outputs:
  results/gck_structure_comparison_summary.csv
  results/gck_3vev_vs_1v4s_perresidue.csv
  figures/fig_gck_3vev_vs_1v4s.png
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import three_to_index, index_to_one
import common as C

warnings.simplefilter("ignore")
VEV = os.path.join(C.ANALYSIS, "data", "structures", "af2", "HXK4_HUMAN.pdb")   # apo AF model (tagged "3VEV" in outputs)
V4S = os.path.join(C.ANALYSIS, "data", "structures", "pdb", "1V4S.pdb")
SMALL = range(67, 204)     # hexokinase small subdomain
LARGE = range(204, 444)    # hexokinase large subdomain
CURATED = [49, 63, 99, 156, 214, 287, 345, 349, 377, 401, 402]


def one(resname):
    try:
        return index_to_one(three_to_index(resname))
    except Exception:
        return "?"


def load(path):
    ch = PDBParser(QUIET=True).get_structure("s", path)[0]["A"]
    ca, heavy, glc, mrk = {}, {}, [], []
    for r in ch:
        if r.id[0] == " ":
            n = r.id[1]
            if "CA" in r:
                ca[n] = r["CA"].coord
            heavy[n] = np.array([a.coord for a in r if a.element != "H"])
        elif r.resname == "GLC":
            glc.append(np.array([a.coord for a in r if a.element != "H"]))
        elif r.resname == "MRK":
            mrk.append(np.array([a.coord for a in r if a.element != "H"]))
    glc = np.vstack(glc) if glc else np.empty((0, 3))
    mrk = np.vstack(mrk) if mrk else np.empty((0, 3))
    return ca, heavy, glc, mrk


def rmsd_on(fixed_ca, moving_ca, positions):
    pos = [n for n in positions if n in fixed_ca and n in moving_ca]
    from Bio.PDB.vectors import Vector
    # use Superimposer with Atom-like shims
    class A:
        def __init__(s, c): s.coord = np.array(c, float)
        def get_coord(s): return s.coord
        def transform(s, rot, tran): s.coord = np.dot(s.coord, rot) + tran
    f = [A(fixed_ca[n]) for n in pos]
    m = [A(moving_ca[n]) for n in pos]
    si = Superimposer(); si.set_atoms(f, m)
    return si.rms, si.rotran, len(pos)


def min_dist(atoms, ligand):
    if len(atoms) == 0 or len(ligand) == 0:
        return np.nan
    return float(cdist(atoms, ligand).min())


def main():
    ca_v, heavy_v, _, _ = load(VEV)          # 3VEV (no ligands)
    ca_s, heavy_s, glc_s, mrk_s = load(V4S)  # 1V4S (ligands)
    common = [n for n in ca_v if n in ca_s]

    # ---- RMSD: global + per subdomain ----
    g_rms, (g_rot, g_tran), g_n = rmsd_on(ca_v, ca_s, common)   # fixed=3VEV, moving=1V4S
    s_rms, _, s_n = rmsd_on(ca_v, ca_s, [n for n in common if n in SMALL])
    l_rms, _, l_n = rmsd_on(ca_v, ca_s, [n for n in common if n in LARGE])

    # ---- transfer 1V4S ligands into 3VEV frame (apply global rotran) ----
    glc_t = np.dot(glc_s, g_rot) + g_tran
    mrk_t = np.dot(mrk_s, g_rot) + g_tran

    # ---- per-residue distances in each structure ----
    rows = []
    for n in common:
        d_s = np.nanmin([min_dist(heavy_s[n], glc_s), min_dist(heavy_s[n], mrk_s)])
        d_v = np.nanmin([min_dist(heavy_v[n], glc_t), min_dist(heavy_v[n], mrk_t)])
        rows.append({"position": n, "dist_1V4S": d_s, "dist_3VEV": d_v})
    per = pd.DataFrame(rows)
    per["abs_diff"] = (per["dist_1V4S"] - per["dist_3VEV"]).abs()
    per.to_csv(os.path.join(C.RESULTS, "gck_3vev_vs_1v4s_perresidue.csv"), index=False)

    # ---- correlation of the two distance measurements ----
    rho_dist, _ = spearmanr(per["dist_1V4S"], per["dist_3VEV"])
    pear = np.corrcoef(per["dist_1V4S"], per["dist_3VEV"])[0, 1]

    # ---- does disagreement<->distance correlation change? ----
    summ = [{"metric": "global_CA_RMSD_A", "value": round(g_rms, 2), "n": g_n},
            {"metric": "small_subdomain_RMSD_A(67-203)", "value": round(s_rms, 2), "n": s_n},
            {"metric": "large_subdomain_RMSD_A(204-443)", "value": round(l_rms, 2), "n": l_n},
            {"metric": "dist_measurement_Spearman_1V4S_vs_3VEV", "value": round(rho_dist, 3), "n": len(per)},
            {"metric": "dist_measurement_Pearson_1V4S_vs_3VEV", "value": round(pear, 3), "n": len(per)},
            {"metric": "median_abs_dist_diff_A", "value": round(per["abs_diff"].median(), 2), "n": len(per)},
            {"metric": "max_abs_dist_diff_A", "value": round(per["abs_diff"].max(), 2), "n": len(per)}]

    for key, dms_id in [("activity", "HXK4_HUMAN_Gersing_2022_activity"),
                        ("abundance", "HXK4_HUMAN_Gersing_2023_abundance")]:
        dis = pd.read_csv(os.path.join(C.RESULTS, f"hxk4_{dms_id}_site_disagreement.csv"))
        mm = per.merge(dis[["position", "disagreement"]], on="position")
        r_s, p_s = spearmanr(mm["disagreement"], mm["dist_1V4S"])
        r_v, p_v = spearmanr(mm["disagreement"], mm["dist_3VEV"])
        summ.append({"metric": f"{key}_Spearman_disagree_vs_dist_1V4S", "value": round(r_s, 3), "n": len(mm)})
        summ.append({"metric": f"{key}_Spearman_disagree_vs_dist_3VEV", "value": round(r_v, 3), "n": len(mm)})

    # ---- ligand-free check: curated hotspot pairwise CA distances ----
    hs = [n for n in CURATED if n in ca_v and n in ca_s]
    Dv = cdist(np.array([ca_v[n] for n in hs]), np.array([ca_v[n] for n in hs]))
    Ds = cdist(np.array([ca_s[n] for n in hs]), np.array([ca_s[n] for n in hs]))
    iu = np.triu_indices(len(hs), k=1)
    rho_hs, _ = spearmanr(Ds[iu], Dv[iu])
    summ.append({"metric": "hotspot_pairwiseCA_Spearman_1V4S_vs_3VEV", "value": round(rho_hs, 3), "n": len(iu[0])})
    summ.append({"metric": "hotspot_pairwiseCA_median_abs_diff_A",
                 "value": round(float(np.median(np.abs(Ds[iu] - Dv[iu]))), 2), "n": len(iu[0])})

    S = pd.DataFrame(summ)
    S.to_csv(os.path.join(C.RESULTS, "gck_structure_comparison_summary.csv"), index=False)
    print("=== 3VEV vs 1V4S comparison ===")
    print(S.to_string(index=False))

    # ---- per-hotspot distances side by side ----
    hot = per[per["position"].isin(CURATED)].copy()
    print("\n=== curated hotspots: nearest-ligand distance (A) in each structure ===")
    print(hot[["position", "dist_1V4S", "dist_3VEV", "abs_diff"]].round(1).to_string(index=False))

    # =================== FIGURE ===================
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.6))
    ax[0].scatter(per["dist_1V4S"], per["dist_3VEV"], s=14, alpha=0.6, color="#4477aa")
    lim = [0, max(per[["dist_1V4S", "dist_3VEV"]].max()) + 2]
    ax[0].plot(lim, lim, "k--", lw=1, label="y = x")
    ax[0].set_xlabel("min dist to nearest ligand — 1V4S (Å)")
    ax[0].set_ylabel("min dist to nearest ligand — 3VEV frame (Å)")
    ax[0].set_title(f"Per-residue distance agreement\nPearson={pear:.3f}, Spearman={rho_dist:.3f}, "
                    f"median |Δ|={per['abs_diff'].median():.2f} Å")
    ax[0].legend(fontsize=8)

    ax[1].bar(np.arange(len(hot)) - 0.2, hot["dist_1V4S"], width=0.4, label="1V4S", color="#4477aa")
    ax[1].bar(np.arange(len(hot)) + 0.2, hot["dist_3VEV"], width=0.4, label="3VEV", color="#ee6677")
    ax[1].set_xticks(np.arange(len(hot)))
    ax[1].set_xticklabels([f"{int(p)}" for p in hot["position"]], fontsize=8)
    ax[1].set_xlabel("curated hotspot residue"); ax[1].set_ylabel("dist to nearest ligand (Å)")
    ax[1].set_title("Hotspot distances by structure")
    ax[1].legend(fontsize=8)
    fig.suptitle(f"GCK structure comparison — global Cα RMSD {g_rms:.2f} Å "
                 f"(small {s_rms:.2f} / large {l_rms:.2f})", fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(C.FIGURES, "fig_gck_3vev_vs_1v4s.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
