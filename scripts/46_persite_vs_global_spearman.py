"""
Per-site Spearman (weighted by #mutations/site) vs overall Spearman.

For every (assay, model) among the 23-model subset we compute two accuracies:

  * overall_spearman  : Spearman(DMS_score, model_score) over ALL variants in the
                        assay -- the standard ProteinGym per-assay benchmark number
                        (matches 01_compute_spearman.py).
  * persite_spearman  : within-site accuracy. Restrict to single substitutions,
                        group by residue position, and for each site with >= 3
                        mutations compute Spearman(DMS_score, model_score) among
                        the mutations AT THAT SITE. Aggregate to one number per
                        (assay, model) as a weighted mean of the site Spearmans,
                        weighting each site by its number of mutations. Sites with
                        no variance (undefined Spearman) are dropped from both the
                        weighted sum and the weight total.

overall mixes the easy across-position signal with the hard within-site signal;
persite isolates the within-site landscape ("which substitution at this position").
The scatter asks: does whole-assay accuracy predict within-site accuracy?

Scores are used as-is (already higher = fitter; see 00_check_directionality.py).

Outputs
  results/persite_vs_global/persite_vs_global_long.csv     (per assay x model)
  results/persite_vs_global/persite_vs_global_summary.csv  (per-assay correlations)
  figures/persite_vs_global/aggregate_persite_vs_global.png/.pdf
  figures/persite_vs_global/<DMS_id>.png                   (one per target assay)
"""
import os
from importlib import import_module

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import common as C

sub = import_module("04_subset_figs")   # 23-model SUBSET, palette, type map

MIN_SITE = 3                            # min mutations/site for a site Spearman
OUTFIG = os.path.join(C.FIGURES, "persite_vs_global")
OUTRES = os.path.join(C.RESULTS, "persite_vs_global")
os.makedirs(OUTFIG, exist_ok=True)
os.makedirs(OUTRES, exist_ok=True)

# Target proteins -> their assays (one scatter per assay)
TARGET_ASSAYS = [
    "PTEN_HUMAN_Matreyek_2021",
    "PTEN_HUMAN_Mighell_2018",
    "P53_HUMAN_Giacomelli_2018_Null_Etoposide",
    "P53_HUMAN_Giacomelli_2018_Null_Nutlin",
    "P53_HUMAN_Giacomelli_2018_WT_Nutlin",
    "P53_HUMAN_Kotler_2018",
    "HXK4_HUMAN_Gersing_2022_activity",
    "HXK4_HUMAN_Gersing_2023_abundance",
    "CP2C9_HUMAN_Amorosi_2021_abundance",
    "CP2C9_HUMAN_Amorosi_2021_activity",
]

# subset columns actually present as scored model outputs
SUBSET_COLS = sub.COLS
COL2DISP = sub.COL2DISP
COL2TYPE = {s[1]: s[2] for s in sub.SUBSET}
PALETTE = sub.PALETTE          # model_type -> colour
MTYPES = sub.MTYPES


def spearman_all(x, y):
    m = x.notna() & y.notna()
    if m.sum() < 3:
        return np.nan
    xr, yr = x[m], y[m]
    if xr.nunique() < 2 or yr.nunique() < 2:
        return np.nan
    return float(spearmanr(xr, yr).correlation)


def persite_weighted_spearman(dms_col, model_col, pos):
    """Weighted mean of per-site Spearman(dms, model), weight = #mutations/site.
    Returns (weighted_rho, n_sites_used, n_muts_used)."""
    d = pd.DataFrame({"dms": dms_col.values, "mdl": model_col.values, "pos": pos.values})
    d = d.dropna(subset=["pos"])
    num = 0.0
    wsum = 0.0
    n_sites = 0
    n_muts = 0
    for _, g in d.groupby("pos"):
        gg = g.dropna(subset=["dms", "mdl"])
        n = len(gg)
        if n < MIN_SITE:
            continue
        if gg["dms"].nunique() < 2 or gg["mdl"].nunique() < 2:
            continue
        rho = spearmanr(gg["dms"], gg["mdl"]).correlation
        if not np.isfinite(rho):
            continue
        num += n * rho
        wsum += n
        n_sites += 1
        n_muts += n
    if wsum == 0:
        return np.nan, 0, 0
    return num / wsum, n_sites, n_muts


def compute_long():
    tax = C.dms_to_taxon()
    files = C.list_score_files()
    rows = []
    for path in files:
        dms_id = C.dms_id_from_path(path)
        df = pd.read_csv(path)
        present = [c for c in SUBSET_COLS if c in df.columns]
        if not present:
            continue
        pos = df["mutant"].map(C.parse_position)
        dms = df["DMS_score"]
        for c in present:
            overall = spearman_all(dms, df[c])
            ps, n_sites, n_muts = persite_weighted_spearman(dms, df[c], pos)
            rows.append({
                "DMS_id": dms_id,
                "taxon": tax.get(dms_id, "Unknown"),
                "model": COL2DISP[c],
                "model_col": c,
                "model_type": COL2TYPE[c],
                "overall_spearman": overall,
                "persite_spearman": ps,
                "n_sites_used": n_sites,
                "n_muts_in_sites": n_muts,
                "n_variants": len(df),
            })
        print(f"{dms_id:48s} models={len(present)}")
    return pd.DataFrame(rows)


def corr_stats(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    r = pearsonr(x[m], y[m])[0]
    rho = spearmanr(x[m], y[m]).correlation
    return float(r), float(rho), int(m.sum())


def scatter(ax, d, title, annotate_models=False):
    for mt in MTYPES:
        sd = d[d["model_type"] == mt]
        if sd.empty:
            continue
        ax.scatter(sd["overall_spearman"], sd["persite_spearman"],
                   s=(38 if annotate_models else 10),
                   c=[PALETTE[mt]], label=mt, alpha=0.75,
                   edgecolors="none")
    if annotate_models:
        for _, r in d.iterrows():
            if np.isfinite(r["overall_spearman"]) and np.isfinite(r["persite_spearman"]):
                ax.annotate(r["model"], (r["overall_spearman"], r["persite_spearman"]),
                            fontsize=5.5, alpha=0.7,
                            xytext=(2, 2), textcoords="offset points")
    lo = np.nanmin([d["overall_spearman"].min(), d["persite_spearman"].min(), 0])
    hi = np.nanmax([d["overall_spearman"].max(), d["persite_spearman"].max(), 0.1])
    pad = 0.05 * (hi - lo)
    lims = [lo - pad, hi + pad]
    ax.plot(lims, lims, ls="--", c="0.6", lw=1, zorder=0)
    ax.axhline(0, c="0.85", lw=0.8, zorder=0)
    ax.axvline(0, c="0.85", lw=0.8, zorder=0)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_aspect("equal")
    r, rho, n = corr_stats(d["overall_spearman"].values, d["persite_spearman"].values)
    ax.set_xlabel("Overall Spearman (whole-assay)")
    ax.set_ylabel("Per-site Spearman (weighted by #mut/site)")
    ax.set_title(f"{title}\nPearson r={r:.2f}, Spearman ρ={rho:.2f}  (n={n})", fontsize=10)
    return r, rho, n


def main():
    long = compute_long()
    long.to_csv(os.path.join(OUTRES, "persite_vs_global_long.csv"), index=False)
    print(f"\nwrote long table: {len(long)} rows, "
          f"{long['DMS_id'].nunique()} assays x {long['model'].nunique()} models")

    # ---- aggregate scatter across all assays ----
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    r, rho, n = scatter(ax, long, "Per-site vs overall Spearman — all 23 models × 217 assays")
    ax.legend(fontsize=7, loc="lower right", framealpha=0.9, title="model type")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTFIG, "aggregate_persite_vs_global.png"), dpi=200)
    fig.savefig(os.path.join(OUTFIG, "aggregate_persite_vs_global.pdf"))
    plt.close(fig)
    print(f"aggregate: Pearson r={r:.3f}, Spearman rho={rho:.3f}, n={n}")

    # ---- per-assay summary + one scatter per target assay ----
    summ = []
    for dms_id, d in long.groupby("DMS_id"):
        r, rho, n = corr_stats(d["overall_spearman"].values, d["persite_spearman"].values)
        summ.append({"DMS_id": dms_id, "pearson_r": r, "spearman_rho": rho, "n_models": n,
                     "mean_overall": d["overall_spearman"].mean(),
                     "mean_persite": d["persite_spearman"].mean()})
    pd.DataFrame(summ).sort_values("DMS_id").to_csv(
        os.path.join(OUTRES, "persite_vs_global_summary.csv"), index=False)

    for dms_id in TARGET_ASSAYS:
        d = long[long["DMS_id"] == dms_id]
        if d.empty:
            print(f"WARN: no rows for {dms_id}")
            continue
        fig, ax = plt.subplots(figsize=(6.8, 6.8))
        scatter(ax, d, dms_id, annotate_models=True)
        handles = [mlines.Line2D([], [], marker="o", ls="", color=PALETTE[mt], label=mt)
                   for mt in MTYPES if (d["model_type"] == mt).any()]
        ax.legend(handles=handles, fontsize=6.5, loc="lower right",
                  framealpha=0.9, title="model type")
        fig.tight_layout()
        fig.savefig(os.path.join(OUTFIG, f"{dms_id}.png"), dpi=200)
        plt.close(fig)
        print(f"assay fig: {dms_id}")


if __name__ == "__main__":
    main()
