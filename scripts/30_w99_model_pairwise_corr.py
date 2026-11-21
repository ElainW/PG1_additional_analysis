"""
W99 (HXK4_HUMAN, Gersing 2022 activity): model-vs-model agreement.

For the 23 subset models (same as Fig5_subset / 09_w99_subset_figs), compute the
pairwise Spearman correlation between their rankings of the 19 substitutions at
W99. Shows a 23x23 heatmap (rows/cols grouped by model_type, diagonal masked =
self excluded) and prints the mean off-diagonal Spearman.

Outputs:
  results/site_W99_<assay>_model_pairwise_spearman.csv
  figures/GCK_w99_deepdive/fig_site_W99_<assay>_model_pairwise_corr.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import common as C
from importlib import import_module

sub = import_module("04_subset_figs")

SITE = 99
DMS_ID = "HXK4_HUMAN_Gersing_2022_activity"


def main(dms_id=DMS_ID, site=SITE):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
    present = [c for c in sub.COLS if c in df.columns]
    rows = df[df["mutant"].str.match(rf"^[A-Z]{site}[A-Z]$")].copy()
    wt = rows["mutant"].iloc[0][0]

    # order by model_type (keeps subset order within type)
    order = []
    for t in sub.MTYPES:
        order += [c for c in present if sub.TYPE[sub.COL2DISP[c]] == t]
    names = [sub.COL2DISP[c] for c in order]
    types = [sub.TYPE[n] for n in names]

    X = rows[order].copy()
    X.columns = names
    corr = X.corr(method="spearman")                # pairwise Spearman over the 19 subs
    corr.to_csv(C.rp(f"site_{wt}{site}_{dms_id}_model_pairwise_spearman.csv"))

    n = len(names)
    off = corr.values[~np.eye(n, dtype=bool)]
    print(f"{dms_id} {wt}{site}: {n} models, {len(rows)} substitutions")
    print(f"Mean pairwise Spearman (off-diagonal, self excluded): {np.nanmean(off):.3f}")
    print(f"  median {np.nanmedian(off):.3f}; min {np.nanmin(off):.3f}; max {np.nanmax(off):.3f}")

    # within vs between model_type
    same = np.array([[types[i] == types[j] for j in range(n)] for i in range(n)]) & ~np.eye(n, dtype=bool)
    diff = np.array([[types[i] != types[j] for j in range(n)] for i in range(n)])
    print(f"  within-model_type mean: {corr.values[same].mean():.3f}; between-type mean: {corr.values[diff].mean():.3f}")

    # ---- figure ----
    palette = dict(sub.PALETTE)
    M = np.ma.masked_array(corr.values, mask=np.eye(n, dtype=bool))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e6e6e6")
    fig, ax = plt.subplots(figsize=(11.5, 9.5))
    im = ax.imshow(M, cmap=cmap, vmin=-1, vmax=1)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(names, rotation=90, fontsize=8)
    ax.set_yticklabels(names, fontsize=8)
    for lab_ax in (ax.get_xticklabels(), ax.get_yticklabels()):
        for tick, t in zip(lab_ax, types):
            tick.set_color(palette[t])
    # separators between model_type groups
    bounds = [i for i in range(1, n) if types[i] != types[i - 1]]
    for b in bounds:
        ax.axhline(b - 0.5, color="black", lw=0.8)
        ax.axvline(b - 0.5, color="black", lw=0.8)
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.5, pad=0.22)
    cb.set_label("Spearman (19 substitutions)")
    ax.set_title(f"{dms_id}\nModel-vs-model Spearman of rankings of 19 substitutions at {wt}{site}\n"
                 f"mean off-diagonal = {np.nanmean(off):.3f}", fontsize=11)
    plt.tight_layout()
    outdir = os.path.join(C.FIGURES, "GCK_w99_deepdive")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"fig_site_{wt}{site}_{dms_id}_model_pairwise_corr.pdf")
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
