"""
3-D view of GCK (HXK4_HUMAN) colored by per-site model disagreement.

Renders the 1V4S structure (the published GCK reference; has real glucose +
activator MRK) as a Calpha backbone trace with residues colored by disagreement,
both ligands drawn, and the top hotspots labelled, from two viewing angles per
assay. The backbone trace is broken at any residue-numbering gap.

Also writes B-factor-encoded PDB files (disagreement x 100 in the B-factor
column) so the structure can be opened in PyMOL / ChimeraX and colored by
B-factor for full interactive exploration:
  spectrum b, blue_white_red        # PyMOL example

Outputs:
  figures/fig_gck_3d_view.png
  results/gck_1V4S_bfactor_<assay>.pdb
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns  # registers the 'rocket' colormaps
from mpl_toolkits.mplot3d import Axes3D  # noqa
from Bio.PDB import PDBParser, PDBIO, Select
import common as C

warnings.simplefilter("ignore")
V4S = os.path.join(C.ANALYSIS, "data", "structures", "pdb", "1V4S.pdb")   # published GCK reference entry
ACTIVATOR = "MRK"
ASSAYS = {"activity": "HXK4_HUMAN_Gersing_2022_activity",
          "abundance": "HXK4_HUMAN_Gersing_2023_abundance"}
LABEL_TOP = 10


def load():
    ch = PDBParser(QUIET=True).get_structure("s", V4S)[0]["A"]
    ca, glc, mrk = {}, [], []
    for r in ch:
        if r.id[0] == " " and "CA" in r:
            ca[r.id[1]] = r["CA"].coord
        elif r.resname == "GLC":
            glc += [a.coord for a in r if a.element != "H"]
        elif r.resname == ACTIVATOR:
            mrk += [a.coord for a in r if a.element != "H"]
    return ca, np.array(glc), np.array(mrk)


def write_bfactor_pdb(dms_id, key, dis):
    """Copy 1V4S with B-factor = disagreement*100 (residues w/o value -> 0)."""
    struct = PDBParser(QUIET=True).get_structure("s", V4S)
    dmap = dict(zip(dis["position"], dis["disagreement"]))
    for r in struct[0]["A"]:
        if r.id[0] == " ":
            b = dmap.get(r.id[1], np.nan)
            val = 0.0 if np.isnan(b) else float(b) * 100
            for a in r:
                a.set_bfactor(val)
    io = PDBIO(); io.set_structure(struct)
    out = os.path.join(C.RESULTS, f"gck_1V4S_bfactor_{key}.pdb")
    io.save(out)
    return out


def broken_segments(positions, ca):
    """Yield contiguous Cα coordinate runs, split wherever residue numbering
    jumps (e.g. the disordered 94-98 loop), so the trace never draws across gaps."""
    seg = []
    prev = None
    for p in positions:
        if prev is not None and p != prev + 1 and seg:
            yield np.array(seg); seg = []
        seg.append(ca[p]); prev = p
    if seg:
        yield np.array(seg)


def main():
    ca, glc, mrk = load()
    fig = plt.figure(figsize=(15, 12))
    views = [(20, -60), (20, 120)]

    for i, (key, dms_id) in enumerate(ASSAYS.items()):
        dis = pd.read_csv(os.path.join(C.RESULTS, f"hxk4_{dms_id}_site_disagreement.csv"))
        dmap = dict(zip(dis["position"], dis["disagreement"]))
        wtmap = dict(zip(dis["position"], dis["wt_aa"]))
        pdb_path = write_bfactor_pdb(dms_id, key, dis)

        positions = sorted(ca.keys())
        segments = list(broken_segments(positions, ca))
        scored = [p for p in positions if p in dmap]
        pts = np.array([ca[p] for p in scored])
        cols = np.array([dmap[p] for p in scored])
        top = dis.sort_values("disagreement", ascending=False).head(LABEL_TOP)

        for j, (elev, azim) in enumerate(views):
            ax = fig.add_subplot(2, 2, i * 2 + j + 1, projection="3d")
            for seg in segments:
                ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color="0.8", lw=1, zorder=1)
            sc = ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=cols, cmap="rocket_r",
                            vmin=0.15, vmax=0.85, s=22, zorder=2, depthshade=True)
            if len(glc):
                ax.scatter(*glc.T, color="#2ca02c", s=45, marker="o", label="glucose", zorder=3)
            if len(mrk):
                ax.scatter(*mrk.T, color="#1f77b4", s=45, marker="^", label="activator", zorder=3)
            for _, r in top.iterrows():
                p = int(r["position"])
                if p in ca:
                    ax.text(*ca[p], f"{r['wt_aa']}{p}", fontsize=7, color="black", zorder=4)
            ax.set_title(f"{key}  (view {j+1})", fontsize=10)
            ax.view_init(elev=elev, azim=azim)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
            ax.set_box_aspect((1, 1, 1))
            if j == 0:
                ax.legend(fontsize=7, loc="upper left")
        print(f"[{key}] wrote {pdb_path}")

    cbar_ax = fig.add_axes([0.4, 0.06, 0.22, 0.015])
    fig.colorbar(sc, cax=cbar_ax, orientation="horizontal",
                 label="disagreement (1 - mean pairwise within-site Spearman)")
    fig.suptitle("GCK (1V4S) Cα trace colored by per-site model disagreement\n"
                 "green = glucose, blue ▲ = allosteric activator (MRK)", fontsize=13, y=0.98)
    plt.tight_layout(rect=[0, 0.09, 1, 0.95])
    out = os.path.join(C.FIGURES, "fig_gck_3d_view.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")
    print("Open the B-factor PDBs in PyMOL/ChimeraX and color by B-factor for interactive 3-D.")


if __name__ == "__main__":
    main()
