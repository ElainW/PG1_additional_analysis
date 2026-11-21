"""
Space-filling ("surface") 2-D renders of GCK (1V4S_cleaned.pdb), same two panels as
31_gck_2d_render.py but with the whole protein drawn as a molecular surface (every
heavy atom as a vdW sphere, z-buffered + Lambert-shaded).

Every residue is coloured by the panel's value. Two versions per panel:
  *_surface_<panel>_w99.png        only W99 outlined/labelled (view chosen to expose W99)
  *_surface_<panel>_highsites.png  all high-disagreement sites (rule of 23/27) outlined/labelled
  *_surface_<panel>_lowsites.png   high among-model-agreement sites = low-disagreement tail (rule of 23/30) Panels:
  (a) model disagreement = 1 - among_model_agreement          (PuOr_r, 0..1, white at 0.5)
  (b) mean per-model Spearman with DMS (mean_model_vs_dms)     (PuOr, -0.5..1, white at 0.25; low=orange, high=purple)
The viewing direction is searched to maximise the exposed surface area of the high-disagreement sites
(W99 weighted x3) so they are visible in the 2-D image.

Outputs (also copied to figures/GCK_3d/):
  figures/GCK_w99_deepdive/fig_gck_2d_surface_disagreement.png
  figures/GCK_w99_deepdive/fig_gck_2d_surface_model_vs_dms.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import binary_erosion
from Bio.PDB import PDBParser
import common as C
from importlib import import_module

R31 = import_module("31_gck_2d_render")
PDB, CSV, OUT, OUT3D = R31.PDB, R31.CSV, R31.OUT, R31.OUT3D
MUST = [99]
EDGE_PX = 5       # outline thickness in pixels (0.12 A/px)
MIN_PIX = 125    # sites with < MIN_PIX visible pixels (0.12 A/px) are treated as buried / not labelled
VDW = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}


def load():
    ch = PDBParser(QUIET=True).get_structure("s", PDB)[0]["A"]
    xyz, rad, resn = [], [], []
    for r in ch:
        if r.id[0] != " ":
            continue
        for a in r:
            xyz.append(a.coord); rad.append(VDW.get(a.element, 1.7)); resn.append(r.id[1])
    xyz = np.array(xyz); xyz -= xyz.mean(0)
    return xyz, np.array(rad), np.array(resn)


def zbuffer(P, rad, resn, res):
    """Height map H and per-pixel residue-number map ID (0 = background)."""
    x0, y0 = P[:, 0].min() - 3, P[:, 1].min() - 3
    W = int((P[:, 0].max() + 3 - x0) / res) + 1
    Hh = int((P[:, 1].max() + 3 - y0) / res) + 1
    H = np.full((Hh, W), -1e9); ID = np.zeros((Hh, W), dtype=np.int32)
    for (x, y, z), r, rn in zip(P, rad, resn):
        i0, i1 = int((y - r - y0) / res), int((y + r - y0) / res) + 1
        j0, j1 = int((x - r - x0) / res), int((x + r - x0) / res) + 1
        yy = y0 + (np.arange(i0, i1) + 0.5) * res - y
        xx = x0 + (np.arange(j0, j1) + 0.5) * res - x
        d2 = yy[:, None] ** 2 + xx[None, :] ** 2
        m = d2 < r * r
        h = z + np.sqrt(np.clip(r * r - d2, 0, None))
        win = H[i0:i1, j0:j1]; wid = ID[i0:i1, j0:j1]
        upd = m & (h > win)
        win[upd] = h[upd]; wid[upd] = rn
    return H, ID, (x0, y0)


def pick_view(xyz, rad, resn, hot, n=250, seed=0):
    rng = np.random.default_rng(seed)
    best, bestR = -1, None
    for _ in range(n):
        v = rng.normal(size=3); v /= np.linalg.norm(v)
        a = np.cross(v, [0, 0, 1] if abs(v[2]) < .9 else [1, 0, 0]); a /= np.linalg.norm(a)
        R = np.vstack([a, np.cross(v, a), v])
        _, ID, _ = zbuffer(xyz @ R.T, rad, resn, 0.5)
        score = sum((ID == p).sum() * (3 if p in MUST else 1) for p in hot)
        # require every highlighted site to be visible at all
        score *= 1 + sum((ID == p).sum() > 0 for p in hot)
        if score > best:
            best, bestR = score, R
    return bestR


def render(d, xyz, rad, resn, targets, absent, suffix, name="high-disagreement"):
    R = pick_view(xyz, rad, resn, targets)
    res = 0.12
    H, ID, (x0, y0) = zbuffer(xyz @ R.T, rad, resn, res)
    bg = H < -1e8
    Hf = np.where(bg, np.nanmin(np.where(bg, np.nan, H)) - 5, H)
    gy, gx = np.gradient(Hf, res)
    nrm = np.dstack([-gx, -gy, np.ones_like(gx)]); nrm /= np.linalg.norm(nrm, axis=2, keepdims=True)
    L = np.array([-0.4, 0.5, 0.75]); L /= np.linalg.norm(L)
    shade = 0.45 + 0.55 * np.clip(nrm @ L, 0, 1)
    vis = {p: int((ID == p).sum()) for p in targets}
    print("visible pixels per site:", vis)
    shown = [p for p in targets if vis[p] >= MIN_PIX]
    hidden = [f"{d.loc[p,'wt']}{p}" for p in targets if vis[p] < MIN_PIX]
    if hidden:
        print("buried in this view:", hidden)

    panels = [
        ("disagreement", "disagreement", "Model disagreement  (1 − among-model agreement)", 0, 1,
         "fig_gck_2d_surface_disagreement.png"),
        ("model_vs_dms", "mean_model_vs_dms", "Mean per-model Spearman with DMS", -0.5, 1.0,
         "fig_gck_2d_surface_model_vs_dms.png"),
    ]
    os.makedirs(OUT, exist_ok=True); os.makedirs(OUT3D, exist_ok=True)
    hot_mask = np.isin(ID, shown)
    edge = hot_mask & ~binary_erosion(hot_mask, iterations=EDGE_PX)
    ext = (x0, x0 + ID.shape[1] * res, y0, y0 + ID.shape[0] * res)
    for key, col, label, lo, hi, fname0 in panels:
        fname = fname0.replace('.png', suffix + '.png')
        # white at the midpoint: 0.5 (disagreement, high=orange) / 0.25 (Spearman, high=purple)
        cmap = plt.get_cmap("PuOr" if key == "model_vs_dms" else "PuOr_r")
        base = np.full(ID.shape + (3,), 0.80)            # residues without data stay grey
        for p in np.unique(ID[ID > 0]):
            v = d[col].get(int(p), np.nan)
            if not np.isnan(v):
                base[ID == p] = cmap((v - lo) / (hi - lo))[:3]
        img = base * shade[..., None]
        img[edge] = 0
        img[bg] = 1
        fig, ax = plt.subplots(figsize=(8, 7.5))
        ax.imshow(img, origin="lower", extent=ext, interpolation="bilinear")
        cx, cy = (ext[0] + ext[1]) / 2, (ext[2] + ext[3]) / 2
        pts = [(ext[0], ext[2]), (ext[1], ext[3])]
        for p in shown:
            yy, xx = np.nonzero(ID == p)
            px, py = x0 + (xx.mean() + 0.5) * res, y0 + (yy.mean() + 0.5) * res
            v = np.array([px - cx, py - cy]); v /= np.linalg.norm(v) + 1e-9
            pts.append((px + v[0] * 12, py + v[1] * 12))
            ax.annotate(f"{d.loc[p,'wt']}{p}", (px, py), xytext=(px + v[0] * 12, py + v[1] * 12),
                        fontsize=11 if p in MUST else 9, fontweight="bold" if p in MUST else "normal",
                        ha="center", va="center", arrowprops=dict(arrowstyle="-", color="black", lw=0.7),
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
        pts = np.array(pts)     # limits hug the structure + labels (no dead space above the colorbar)
        ax.set_xlim(pts[:, 0].min() - 6, pts[:, 0].max() + 6); ax.set_ylim(pts[:, 1].min() - 6, pts[:, 1].max() + 6)
        ax.set_aspect("equal"); ax.axis("off")
        multi = len(targets) > 1
        what = (f"outlined = {name} sites (n={len(shown)}"
                + ("; not in 1V4S: " + ", ".join(absent) if absent else "")
                + ("; buried here: " + ", ".join(hidden) if hidden else "") + ")") if multi else "W99 outlined"
        ax.set_title(f"GCK (1V4S) surface — {label}\nHXK4_HUMAN_Gersing_2022_activity; all residues coloured; {what}", fontsize=10)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(lo, hi))
        plt.tight_layout()
        ax.apply_aspect(); pos = ax.get_position()
        cax = fig.add_axes([pos.x0 + 0.25 * pos.width, pos.y0 - 0.035, 0.5 * pos.width, 0.022])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(label)
        if key == "model_vs_dms":
            cb.set_ticks(np.arange(-0.5, 1.01, 0.25))
        for o in (OUT, OUT3D):
            out = os.path.join(o, fname)
            plt.savefig(out, dpi=200, bbox_inches="tight")
            print("Saved", out)
        plt.close()


def main():
    d = pd.read_csv(CSV).set_index("pos")
    d["disagreement"] = 1 - d["among_model_agreement"]
    xyz, rad, resn = load()
    present = set(resn.tolist())
    top, thr, pct = R31.high_disagreement_sites()
    hot_all = list(dict.fromkeys(top + MUST))
    absent = [f"{d.loc[p,'wt']}{p}" for p in hot_all if p not in present]
    hot = [p for p in hot_all if p in present]
    # W99-only version (default names) and high-disagreement-sites version (_highsites)
    render(d, xyz, rad, resn, MUST, absent, "_w99")
    render(d, xyz, rad, resn, hot, absent, "_highsites")
    # high among-model agreement = LOW-disagreement tail (mirror rule of 23/30)
    low, lthr, lpct = R31.low_disagreement_sites()
    labsent = [f"{d.loc[p,'wt']}{p}" for p in low if p not in present]
    lows = [p for p in low if p in present]
    print(f"low-disagreement rule: <{R31.S.DISAGREE_CEIL} and <= elbow {lthr:.3f} (lowest {lpct:.1f}%): "
          f"{len(low)} sites: {[f'{d.loc[p,chr(119)+chr(116)]}{p}' for p in low]}")
    render(d, xyz, rad, resn, lows, labsent, "_lowsites", name="high-agreement (low-disagreement)")


if __name__ == "__main__":
    main()
