"""
2-D renders of GCK (1V4S_cleaned.pdb) colored per residue by
  (a) model disagreement = 1 - among_model_agreement
  (b) mean per-model Spearman with the DMS score (mean_model_vs_dms)
for HXK4_HUMAN_Gersing_2022_activity.

Residues are Calpha spheres on a thin backbone trace. The viewing direction is
chosen (by search) to keep the top-disagreement sites (incl. W99) unoccluded;
only high-disagreement sites (rule of 23/27) are coloured; they are drawn last, outlined, and labelled in BOTH panels
so the two images are directly comparable.

Outputs:
  figures/GCK_w99_deepdive/fig_gck_2d_disagreement.png   (panel a)
  figures/GCK_w99_deepdive/fig_gck_2d_model_vs_dms.png   (panel b)
  (both also copied to figures/GCK_3d/)
  results/gck_1V4S_bfactor_{disagreement,model_vs_dms}.pdb (for PyMOL/ChimeraX)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBIO
import common as C
import importlib.util

_sp = importlib.util.spec_from_file_location("structure_maps_all", os.path.join(os.path.dirname(os.path.abspath(__file__)), "23_structure_maps_all.py"))
S = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(S)
DMS_ID = "HXK4_HUMAN_Gersing_2022_activity"


def high_disagreement_sites():
    """Same selection as 23/27: disagreement > 0.5 and >= elbow threshold within the
    top 2.5-5% tail of the assay (default 5% mark if no clear drop)."""
    r = pd.read_csv(S.SITES)
    r = r[r["DMS_id"] == DMS_ID]
    thr, pct = S.elbow_tail_threshold(r["disagreement"].values)
    lab = r[(r["disagreement"] > S.DISAGREE_FLOOR) & (r["disagreement"] >= thr)]
    lab = lab.sort_values("disagreement", ascending=False).head(S.MAX_LABELS)
    return lab["up_pos"].astype(int).tolist(), thr, pct

PDB = os.path.join(C.ANALYSIS, "data", "structures", "pdb", "1V4S_cleaned.pdb")
CSV = os.path.join(C.ANALYSIS, "results", "hxk4_pilot_analysis", "hxk4_agreement_vs_accuracy_sites_activity.csv")
OUT = os.path.join(C.FIGURES, "GCK_w99_deepdive")
OUT3D = os.path.join(C.FIGURES, "GCK_3d")
MUST = [99]


def low_disagreement_sites():
    """Mirror of the above (23/30): disagreement < 0.5 and <= low-tail elbow threshold
    (95-97.5 percentile band, default lowest 5%) = high among-model agreement."""
    r = pd.read_csv(S.SITES)
    r = r[r["DMS_id"] == DMS_ID]
    thr, pct = S.elbow_low_tail_threshold(r["disagreement"].values)
    lab = r[(r["disagreement"] < S.DISAGREE_CEIL) & (r["disagreement"] <= thr)]
    lab = lab.sort_values("disagreement", ascending=True).head(S.MAX_LABELS)
    return lab["up_pos"].astype(int).tolist(), thr, pct


def main():
    d = pd.read_csv(CSV).set_index("pos")
    d["disagreement"] = 1 - d["among_model_agreement"]
    st = PDBParser(QUIET=True).get_structure("s", PDB)[0]
    ch = st["A"]
    pos, xyz = [], []
    for r in ch:
        if r.id[0] == " " and "CA" in r:
            pos.append(r.id[1]); xyz.append(r["CA"].coord)
    pos = np.array(pos); xyz = np.array(xyz)
    xyz = xyz - xyz.mean(0)
    glc = np.array([a.coord for r in ch if r.get_resname() == "GLC" for a in r]) - 0
    gc = np.array([r["CA"].coord for r in ch if r.id[0] == " " and "CA" in r]).mean(0)
    glc = glc - gc

    top, thr, pct = high_disagreement_sites()
    print(f"high-disagreement rule: >{S.DISAGREE_FLOOR} and >= elbow {thr:.3f} (top {pct:.1f}% tail): {len(top)} sites")
    hot = list(dict.fromkeys(top + MUST))
    absent = [f"{d.loc[p,'wt']}{p}" for p in hot if p not in pos]
    hot = [p for p in hot if p in pos]
    hot_idx = [int(np.where(pos == p)[0][0]) for p in hot]
    if absent:
        print("not resolved in structure (omitted):", absent)
    print("labelled hotspots:", [f"{d.loc[p,'wt']}{p}" for p in hot])

    # ---- pick a view that keeps hotspots unoccluded ----
    rng = np.random.default_rng(0)
    best, bestR = 1e9, None
    for _ in range(4000):
        v = rng.normal(size=3); v /= np.linalg.norm(v)            # viewing dir (toward viewer)
        a = np.cross(v, [0, 0, 1] if abs(v[2]) < .9 else [1, 0, 0]); a /= np.linalg.norm(a)
        b = np.cross(v, a)
        R = np.vstack([a, b, v])
        P = xyz @ R.T
        score = 0
        for i in hot_idx:
            near = (np.hypot(*(P[:, :2] - P[i, :2]).T) < 3.0) & (P[:, 2] > P[i, 2] + 0.5)
            score += near.sum() * (3 if pos[i] in MUST else 1)
        if score < best:
            best, bestR = score, R
    print(f"view occlusion score {best}")
    P = xyz @ bestR.T
    G = glc @ bestR.T

    panels = [
        ("disagreement", "disagreement",
         "Model disagreement  (1 − among-model agreement)", "RdBu_r", 0, 1,
         "fig_gck_2d_disagreement.png"),
        ("model_vs_dms", "mean_model_vs_dms",
         "Mean per-model Spearman with DMS", "RdBu_r", -1, 1, "fig_gck_2d_model_vs_dms.png"),
    ]
    os.makedirs(OUT, exist_ok=True); os.makedirs(OUT3D, exist_ok=True)
    for key, col, label, cmapname, lo, hi, fname in panels:
        vals = np.array([d[col].get(p, np.nan) for p in pos])
        cmap = plt.get_cmap(cmapname).copy(); cmap.set_bad("#cccccc")
        fig, ax = plt.subplots(figsize=(8, 7.5))
        ax.plot(P[:, 0], P[:, 1], color="#999", lw=0.5, alpha=0.5, zorder=0)
        order = np.argsort(P[:, 2])                      # back -> front
        order = [i for i in order if i not in hot_idx] + sorted(hot_idx, key=lambda i: P[i, 2])
        zr = (P[:, 2] - P[:, 2].min()) / np.ptp(P[:, 2])
        for rank, i in enumerate(order):
            is_hot = i in hot_idx
            if not is_hot:      # only high-disagreement sites are shown in colour
                ax.scatter(P[i, 0], P[i, 1], s=22, c="#d9d9d9", edgecolor="white", linewidth=0.3,
                           alpha=0.6, zorder=2 + rank)
                continue
            ax.scatter(P[i, 0], P[i, 1], s=(150 if is_hot else 38 + 25 * zr[i]),
                       c=[vals[i]], cmap=cmap, vmin=lo, vmax=hi,
                       edgecolor="black" if is_hot else "white", linewidth=1.3 if is_hot else 0.3,
                       alpha=1 if is_hot else 0.55 + 0.4 * zr[i], zorder=2 + rank)
        ax.scatter(G[:, 0], G[:, 1], marker="D", s=12, c="limegreen", edgecolor="k", linewidth=.3,
                   zorder=2 + len(order) + 1, label="glucose (GLC)")
        # labels, pushed radially outward from the structure centre
        for i in hot_idx:
            p = pos[i]
            dirn = P[i, :2] / (np.linalg.norm(P[i, :2]) + 1e-9)
            ax.annotate(f"{d.loc[p,'wt']}{p}", P[i, :2], xytext=P[i, :2] + dirn * 9, fontsize=11 if p in MUST else 9,
                        fontweight="bold" if p in MUST else "normal", ha="center", va="center",
                        arrowprops=dict(arrowstyle="-", color="black", lw=0.7), zorder=1000,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(lo, hi))
        cb = fig.colorbar(sm, ax=ax, orientation="horizontal", shrink=0.5, pad=0.03)
        cb.set_label(label)
        lp = np.array([P[i, :2] + P[i, :2] / (np.linalg.norm(P[i, :2]) + 1e-9) * 9 for i in hot_idx])
        allp = np.vstack([P[:, :2], lp])
        ax.set_xlim(allp[:, 0].min() - 6, allp[:, 0].max() + 6)
        ax.set_ylim(allp[:, 1].min() - 6, allp[:, 1].max() + 6)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"GCK (1V4S) — {label}\nHXK4_HUMAN_Gersing_2022_activity; coloured = high-disagreement sites (n={len(hot)}{'; not in 1V4S: ' + ', '.join(absent) if absent else ''})", fontsize=10)
        plt.tight_layout()
        out = os.path.join(OUT, fname)
        plt.savefig(out, dpi=200, bbox_inches="tight")
        out3d = os.path.join(OUT3D, fname)
        plt.savefig(out3d, dpi=200, bbox_inches="tight"); plt.close()
        print("Saved", out3d)
        print("Saved", out)

        # B-factor PDB (value x100; sites without data -> -999)
        s2 = PDBParser(QUIET=True).get_structure("s", PDB)
        for r in s2[0]["A"]:
            v = d[col].get(r.id[1], np.nan) if r.id[0] == " " else np.nan
            for a in r:
                a.set_bfactor(-999.0 if np.isnan(v) else float(v) * 100)
        io = PDBIO(); io.set_structure(s2)
        io.save(C.rp(f"gck_1V4S_bfactor_{key}.pdb"))


if __name__ == "__main__":
    main()
