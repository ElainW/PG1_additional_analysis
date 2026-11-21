"""
Build two PyMOL session files for GCK (1V4S_cleaned.pdb), one per panel:
  (a) model disagreement   (1 - among_model_agreement)
  (b) mean per-model Spearman with DMS (mean_model_vs_dms)
for HXK4_HUMAN_Gersing_2022_activity.

Needs the B-factor PDBs from 31_gck_2d_render.py (value x100 in the B column,
-999 = no data). Only the high-disagreement sites (>0.5 and >= elbow in the top 2.5-5% tail, as in 23/27) + W99 are coloured/shown as CA spheres
and labelled in BOTH sessions (W99 enlarged), so the two views are comparable.

Run with a PyMOL-enabled python (pymol-open-source / conda pymol):
  python 32_gck_pymol_sessions.py
Outputs:
  results/gck_3d/gck_1V4S_disagreement.pse
  results/gck_3d/gck_1V4S_model_vs_dms.pse
  figures/GCK_animation/gck_1V4S_render.pml   (same commands, loadable in the PyMOL GUI)
"""
import os
import pandas as pd
import common as C
from importlib import import_module

CSV = os.path.join(C.ANALYSIS, "results", "hxk4_pilot_analysis", "hxk4_agreement_vs_accuracy_sites_activity.csv")
OUTDIR = os.path.join(C.RESULTS, "gck_3d")
PMLDIR = os.path.join(C.FIGURES, "GCK_animation")
PANELS = {
    "disagreement": dict(pal="blue_white_red", lo=0, hi=100),
    "model_vs_dms": dict(pal="blue_white_red", lo=-100, hi=100),
}
MUST = [99]


def commands(key, sel, pdb, cfg):
    return f"""reinitialize
load {pdb}, {key}
hide everything, {key}
show cartoon, {key} and polymer
set cartoon_transparency, 0.25
color grey85, {key} and polymer
show sticks, {key} and organic
color green, {key} and resn GLC
color magenta, {key} and resn MRK
select hot, {key} and polymer and resi {sel}
spectrum b, {cfg['pal']}, hot, minimum={cfg['lo']}, maximum={cfg['hi']}
show spheres, hot and name CA
set sphere_scale, 1.3, hot and name CA
select w99, {key} and resi 99 and name CA
set sphere_scale, 2.0, w99
label hot and name CA, resn+resi
set label_size, 18
set label_color, black
bg_color white
orient {key}
"""


def main():
    top, thr, pct = import_module("31_gck_2d_render").high_disagreement_sites()
    sel = "+".join(map(str, dict.fromkeys(top + MUST)))
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(PMLDIR, exist_ok=True)

    pml = ""
    for key, cfg in PANELS.items():
        pdb = os.path.join(OUTDIR, f"gck_1V4S_bfactor_{key}.pdb")
        if not os.path.exists(pdb):
            raise SystemExit(f"missing {pdb}; run 31_gck_2d_render.py first")
        cmds = commands(key, sel, pdb, cfg)
        pml += f"# ---- {key} ----\n{cmds}\nsave {os.path.join(OUTDIR, f'gck_1V4S_{key}.pse')}\n\n"
    open(os.path.join(PMLDIR, "gck_1V4S_render.pml"), "w").write(pml)

    try:
        from pymol import cmd
    except ImportError:
        print("PyMOL python module not found: wrote .pml only. Run it with:\n"
              f"  pymol -cq {os.path.join(PMLDIR, 'gck_1V4S_render.pml')}")
        return
    for key, cfg in PANELS.items():
        cmd.do(commands(key, sel, os.path.join(OUTDIR, f"gck_1V4S_bfactor_{key}.pdb"), cfg))
        out = os.path.join(OUTDIR, f"gck_1V4S_{key}.pse")
        cmd.save(out)
        print("Saved", out)


if __name__ == "__main__":
    main()
