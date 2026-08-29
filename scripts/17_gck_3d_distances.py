"""
GCK (HXK4_HUMAN) 3-D distance analysis of per-site model disagreement.

Structure: PDB 1V4S chain A (human glucokinase, closed form) with glucose (GLC)
and the allosteric activator (MRK) bound. We use 1V4S (not 3VEV) because it is
the structure used in the published GCK literature, and because 3VEV's crystal
construct differs from the assay isoform at the N-terminus (residues 4-11),
adding avoidable wt/PDB mismatches. The two give near-identical distances anyway
(1V4S<->3VEV per-residue nearest-ligand Pearson 0.996; see 18/20). Numbering
matches UniProt P35557 / the assays (verified: T49,R63,W99,... all match).
Covers residues 14-461.

We test, in real 3-D space, two questions raised by the 1-D finding that
disagreement is elevated at binding-site residues:

  A. Are high-disagreement residues physically CLOSE to the ligand pockets?
     -> Spearman(disagreement, min heavy-atom distance to glucose / activator /
        nearest ligand) over all resolved residues; and a Mann-Whitney of
        ligand-distance for top-disagreement residues vs the rest.

  B. Are high-disagreement residues spatially CLUSTERED with each other?
     -> mean pairwise CA-CA distance of the top-N disagreement residues vs a
        permutation null drawn from resolved residues (structure-matched).

Distances are minimum heavy-atom (side-chain-aware) to each ligand's heavy atoms.
Disagreement = 1 - mean pairwise within-site Spearman (from 06 results).

Outputs:
  results/gck_3d_residue_distances.csv           (per residue: dists + disagreement)
  results/gck_3d_summary.csv                     (correlations + permutation tests)
  figures/fig_gck_3d_distance_vs_disagreement.png
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns  # registers the 'rocket' colormaps
from scipy.stats import spearmanr, mannwhitneyu
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one
import common as C

warnings.simplefilter("ignore")
PDB = os.path.join(C.ANALYSIS, "data", "structures", "pdb", "1V4S.pdb")
ACTIVATOR = "MRK"   # 1V4S's allosteric activator (3VEV uses 0H4)
ASSAYS = {
    "activity": "HXK4_HUMAN_Gersing_2022_activity",
    "abundance": "HXK4_HUMAN_Gersing_2023_abundance",
}
TOP_N = 25            # data-driven "hotspot" set size for enrichment/clustering
N_PERM = 20000


def one(resname):
    try:
        return index_to_one(three_to_index(resname))
    except Exception:
        return "?"


def load_structure():
    ch = PDBParser(QUIET=True).get_structure("g", PDB)[0]["A"]
    ca, heavy, resname = {}, {}, {}
    glc, mrk = [], []
    for r in ch:
        het = r.id[0]
        coords = np.array([a.coord for a in r if a.element != "H"])
        if het == " ":                                  # protein residue
            n = r.id[1]
            resname[n] = one(r.resname)
            if "CA" in r:
                ca[n] = r["CA"].coord
            heavy[n] = coords
        elif r.resname == "GLC":
            glc.append(coords)
        elif r.resname == ACTIVATOR:
            mrk.append(coords)
    glc = np.vstack(glc) if glc else np.empty((0, 3))
    mrk = np.vstack(mrk) if mrk else np.empty((0, 3))
    return ca, heavy, resname, glc, mrk


def min_dist(atoms, ligand):
    if len(atoms) == 0 or len(ligand) == 0:
        return np.nan
    return float(cdist(atoms, ligand).min())


def main():
    ca, heavy, resname, glc, mrk = load_structure()
    resolved = sorted(heavy.keys())
    rows = []
    for n in resolved:
        d_glc = min_dist(heavy[n], glc)
        d_mrk = min_dist(heavy[n], mrk)
        rows.append({
            "position": n, "pdb_aa": resname[n],
            "dist_glucose": d_glc, "dist_activator": d_mrk,
            "dist_nearest_ligand": np.nanmin([d_glc, d_mrk]),
        })
    dist = pd.DataFrame(rows)

    summary = []
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, (key, dms_id) in zip(axes, ASSAYS.items()):
        dis = pd.read_csv(os.path.join(C.RESULTS, f"hxk4_{dms_id}_site_disagreement.csv"))
        m = dist.merge(dis[["position", "wt_aa", "disagreement", "mean_pairwise_spearman"]],
                       on="position", how="inner")
        # wt sanity
        mism = (m["pdb_aa"] != m["wt_aa"]).sum()
        m.to_csv(os.path.join(C.RESULTS, f"gck_3d_residue_distances.csv"), index=False) if key == "activity" else None

        # ---- A. distance vs disagreement (all resolved+scored residues) ----
        for lig in ["dist_glucose", "dist_activator", "dist_nearest_ligand"]:
            rho, p = spearmanr(m["disagreement"], m[lig])
            summary.append({"assay": key, "test": f"spearman_disagree_vs_{lig}",
                            "value": round(rho, 3), "p_value": p, "n": len(m)})
        # ---- A'. top-N hotspots vs rest: ligand distance ----
        m_sorted = m.sort_values("disagreement", ascending=False)
        hot = m_sorted.head(TOP_N)
        rest = m_sorted.iloc[TOP_N:]
        U, p_mw = mannwhitneyu(hot["dist_nearest_ligand"], rest["dist_nearest_ligand"],
                               alternative="less")   # hot CLOSER than rest?
        summary.append({"assay": key, "test": "MW_top25_closer_to_ligand",
                        "value": round(hot["dist_nearest_ligand"].median(), 2),
                        "p_value": p_mw, "n": TOP_N,
                        "rest_median": round(rest["dist_nearest_ligand"].median(), 2)})

        # ---- B. spatial clustering of top-N hotspots (CA-CA) ----
        hot_pos = [pp for pp in hot["position"] if pp in ca]
        pool = [pp for pp in m["position"] if pp in ca]
        def mean_pair_ca(positions):
            P = np.array([ca[pp] for pp in positions])
            D = cdist(P, P)
            return D[np.triu_indices(len(P), k=1)].mean()
        obs = mean_pair_ca(hot_pos)
        rng = np.random.default_rng(0)
        null = np.array([mean_pair_ca(rng.choice(pool, size=len(hot_pos), replace=False))
                         for _ in range(N_PERM)])
        p_clust = (null <= obs).mean()               # smaller mean dist => clustered
        summary.append({"assay": key, "test": "perm_top25_CA_clustering",
                        "value": round(obs, 2), "p_value": p_clust, "n": len(hot_pos),
                        "rest_median": round(float(null.mean()), 2)})

        # ---- figure: disagreement vs nearest-ligand distance ----
        sc = ax.scatter(m["dist_nearest_ligand"], m["disagreement"], s=20,
                        c=m["disagreement"], cmap="rocket_r", edgecolor="none")
        for _, r in hot.head(12).iterrows():
            ax.annotate(f"{r['wt_aa']}{int(r['position'])}",
                        (r["dist_nearest_ligand"], r["disagreement"]),
                        fontsize=7, ha="center", va="bottom", xytext=(0, 3),
                        textcoords="offset points")
        rho_n, p_n = spearmanr(m["disagreement"], m["dist_nearest_ligand"])
        ax.set_title(f"{key}\nSpearman(disagree, dist-to-nearest-ligand) = {rho_n:.2f} (p={p_n:.1e})",
                     fontsize=11)
        ax.set_xlabel("min heavy-atom distance to nearest ligand (Å)\n[glucose or activator]")
        ax.set_ylabel("disagreement (1 - mean pairwise within-site Spearman)")
        ax.axvline(8, color="gray", ls="--", lw=1, label="8 Å (pocket shell)")
        ax.legend(fontsize=8)
        if mism:
            print(f"[{key}] WARNING: {mism} wt/pdb residue mismatches")

    plt.tight_layout()
    out = os.path.join(C.FIGURES, "fig_gck_3d_distance_vs_disagreement.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()

    S = pd.DataFrame(summary)
    S.to_csv(os.path.join(C.RESULTS, "gck_3d_summary.csv"), index=False)
    print("=== 3-D distance summary ===")
    print(S.to_string(index=False))
    print(f"\nSaved {out}")
    print("Saved results/gck_3d_residue_distances.csv, results/gck_3d_summary.csv")

    # curated hotspot distances to each pocket
    cur = [49, 63, 99, 156, 214, 287, 345, 349, 377, 401, 402]
    print("\n=== curated hotspots: distance to each ligand (Å) ===")
    cd = dist[dist["position"].isin(cur)][["position", "pdb_aa", "dist_glucose", "dist_activator"]]
    print(cd.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
