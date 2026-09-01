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
import sys
from importlib import import_module

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from scipy.optimize import curve_fit
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


def _r2(y, yhat):
    y = np.asarray(y, float); yhat = np.asarray(yhat, float)
    ss_res = np.nansum((y - yhat) ** 2)
    ss_tot = np.nansum((y - np.nanmean(y)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else np.nan


def fit_simple_curves(x, y):
    """Fit three simple monotone/concave forms to (x, y) and rank by R^2.

    logarithmic : y = a + b*ln(x)          (x>0 subset only)
    exp_neg     : y = c - a*exp(-b*x)       saturating, defined everywhere
    quad_left   : y = a + b*x + c*x^2, c<0  fit on the left half of the x-range

    Returns list of dicts sorted by R^2 desc: {name, func, params, r2, xdomain}.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    xmin, xmax = x.min(), x.max()
    fits = []

    # --- logarithmic (x>0) ---
    mpos = x > 0
    if mpos.sum() > 5:
        f = lambda xx, a, b: a + b * np.log(xx)
        try:
            p, _ = curve_fit(f, x[mpos], y[mpos], p0=[0.3, 0.1], maxfev=10000)
            r2 = _r2(y[mpos], f(x[mpos], *p))
            # draw only where ln(x) is on-scale (it dives to -inf as x->0+)
            fits.append({"name": "logarithmic", "func": f, "params": p, "r2": r2,
                         "xdomain": (0.05, xmax),
                         "label": f"y = {p[0]:.3f} + {p[1]:.3f}·ln(x)"})
        except Exception as e:
            print("log fit failed:", e)

    # --- saturating exponential: c - a*exp(-b*x) ---
    f = lambda xx, a, b, c: c - a * np.exp(-b * xx)
    try:
        p, _ = curve_fit(f, x, y, p0=[0.4, 2.0, 0.5], maxfev=20000)
        r2 = _r2(y, f(x, *p))
        fits.append({"name": "exp_neg", "func": f, "params": p, "r2": r2,
                     "xdomain": (xmin, xmax),
                     "label": f"y = {p[2]:.3f} − {p[0]:.3f}·e^(−{p[1]:.2f}x)"})
    except Exception as e:
        print("exp fit failed:", e)

    # --- concave quadratic on the left half ---
    xcut = np.median(x)
    ml = x <= xcut
    if ml.sum() > 5:
        p = np.polyfit(x[ml], y[ml], 2)          # [c, b, a] for c*x^2+b*x+a
        f = lambda xx, c, b, a: c * xx ** 2 + b * xx + a
        r2 = _r2(y[ml], np.polyval(p, x[ml]))
        fits.append({"name": "quad_left", "func": f, "params": p, "r2": r2,
                     "xdomain": (xmin, xcut),
                     "label": f"y = {p[0]:.3f}x² + {p[1]:.3f}x + {p[2]:.3f}  (x≤{xcut:.2f})"})

    fits.sort(key=lambda d: (np.nan_to_num(d["r2"], nan=-1)), reverse=True)
    return fits


def curve_xy(fit, n=200):
    lo, hi = fit["xdomain"]
    xg = np.linspace(lo, hi, n)
    return xg, fit["func"](xg, *fit["params"])


def place_labels_leader(ax, xs, ys, texts, fontsize=7):
    """Force-directed label placement with thin leader lines (adjustText-lite).
    Repels labels from one another and from every data point, on a short leash
    to their own anchor; positions solved in pixel space, drawn in data space."""
    fig = ax.figure
    fig.canvas.draw()
    trans = ax.transData
    inv = trans.inverted()
    P = np.array([trans.transform((x, y)) for x, y in zip(xs, ys)], float)
    dpi = fig.dpi
    ch_w = fontsize * dpi / 72 * 0.60
    h = fontsize * dpi / 72 * 1.25
    W = np.array([max(len(t), 1) * ch_w for t in texts], float)
    # seed each label offset radially outward from the point cloud centroid, so
    # dense clusters explode into the surrounding empty space instead of piling up
    ctr = P.mean(axis=0)
    L = P.copy()
    d = L - ctr
    norm = np.hypot(d[:, 0], d[:, 1])[:, None]
    unit = np.where(norm > 1e-6, d / np.maximum(norm, 1e-6), np.array([0.0, 1.0]))
    L = P + unit * (h * 2.2)
    rng = np.random.default_rng(0)
    L += rng.uniform(-1.0, 1.0, L.shape)
    bb = ax.get_window_extent()
    for _ in range(800):
        for i in range(len(L)):
            for j in range(i + 1, len(L)):
                dx = L[i, 0] - L[j, 0]; dy = L[i, 1] - L[j, 1]
                ox = (W[i] + W[j]) / 2 + 2 - abs(dx)
                oy = h * 1.1 - abs(dy)
                if ox > 0 and oy > 0:
                    if ox < oy:
                        s = ox / 2 * (1 if dx >= 0 else -1)
                        L[i, 0] += s; L[j, 0] -= s
                    else:
                        s = oy / 2 * (1 if dy >= 0 else -1)
                        L[i, 1] += s; L[j, 1] -= s
        for i in range(len(L)):
            for k in range(len(P)):
                dx = L[i, 0] - P[k, 0]; dy = L[i, 1] - P[k, 1]
                ox = W[i] / 2 + 3 - abs(dx); oy = h / 2 + 3 - abs(dy)
                if ox > 0 and oy > 0:
                    if ox < oy:
                        L[i, 0] += ox * 0.5 * (1 if dx >= 0 else -1)
                    else:
                        L[i, 1] += oy * 0.5 * (1 if dy >= 0 else -1)
        for i in range(len(L)):
            vx = L[i, 0] - P[i, 0]; vy = L[i, 1] - P[i, 1]
            dist = np.hypot(vx, vy) + 1e-9
            # weak leash toward a modest standoff; lets crowded labels migrate to
            # open space, only reeling back leaders that grow very long
            pull = (dist - h * 2.5) / dist * 0.04
            L[i, 0] -= vx * pull; L[i, 1] -= vy * pull
        L[:, 0] = np.clip(L[:, 0], bb.x0 + W / 2, bb.x1 - W / 2)
        L[:, 1] = np.clip(L[:, 1], bb.y0 + h / 2, bb.y1 - h / 2)
    for i, t in enumerate(texts):
        lx, ly = inv.transform(L[i])
        ax.annotate(t, xy=(xs[i], ys[i]), xytext=(lx, ly), textcoords="data",
                    fontsize=fontsize, ha="center", va="center", color="0.12",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="0.55",
                                    shrinkA=1, shrinkB=3))


def corr_stats(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    r = pearsonr(x[m], y[m])[0]
    rho = spearmanr(x[m], y[m]).correlation
    return float(r), float(rho), int(m.sum())


def scatter(ax, d, title, annotate_models=False):
    ax.grid(False)
    for mt in MTYPES:
        sd = d[d["model_type"] == mt]
        if sd.empty:
            continue
        ax.scatter(sd["overall_spearman"], sd["persite_spearman"],
                   s=(38 if annotate_models else 10),
                   c=[PALETTE[mt]], label=mt, alpha=0.75,
                   edgecolors="none")
    lo = np.nanmin([d["overall_spearman"].min(), d["persite_spearman"].min(), 0])
    hi = np.nanmax([d["overall_spearman"].max(), d["persite_spearman"].max(), 0.1])
    pad = 0.05 * (hi - lo)
    lims = [lo - pad, hi + pad]
    ax.plot(lims, lims, ls="--", c="0.6", lw=1, zorder=0)
    ax.axhline(0, c="0.9", lw=0.8, zorder=0)
    ax.axvline(0, c="0.9", lw=0.8, zorder=0)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_aspect("equal")
    if annotate_models:
        fin = d[np.isfinite(d["overall_spearman"]) & np.isfinite(d["persite_spearman"])]
        place_labels_leader(ax, fin["overall_spearman"].values,
                            fin["persite_spearman"].values, list(fin["model"]))
    r, rho, n = corr_stats(d["overall_spearman"].values, d["persite_spearman"].values)
    ax.set_xlabel("Overall Spearman (whole-assay)")
    ax.set_ylabel("Per-site Spearman (weighted by #mut/site)")
    ax.set_title(f"{title}\nPearson r={r:.2f}, Spearman ρ={rho:.2f}  (n={n})", fontsize=10)
    return r, rho, n


def main():
    long_path = os.path.join(OUTRES, "persite_vs_global_long.csv")
    # Recompute is expensive (~8 min); reuse the cached long table for plot-only
    # reruns. Pass --recompute to force recomputation.
    if "--recompute" not in sys.argv and os.path.exists(long_path):
        long = pd.read_csv(long_path)
        print(f"loaded cached long table: {len(long)} rows "
              f"(pass --recompute to regenerate)")
    else:
        long = compute_long()
        long.to_csv(long_path, index=False)
        print(f"\nwrote long table: {len(long)} rows, "
              f"{long['DMS_id'].nunique()} assays x {long['model'].nunique()} models")

    # ---- aggregate scatter across all assays ----
    fig, ax = plt.subplots(figsize=(8.6, 7.0))
    r, rho, n = scatter(ax, long, "Per-site vs overall Spearman — all 23 models × 217 assays")
    # simple parametric trend: fit log / saturating-exp / concave-quad, keep best R^2
    m = np.isfinite(long["overall_spearman"].values) & np.isfinite(long["persite_spearman"].values)
    fits = fit_simple_curves(long["overall_spearman"].values[m], long["persite_spearman"].values[m])
    print("\ncurve fits (by R^2):")
    for fdt in fits:
        print(f"  {fdt['name']:12s} R2={fdt['r2']:.4f}  {fdt['label']}")
    pd.DataFrame([{"model": f["name"], "r2": f["r2"], "equation": f["label"]}
                  for f in fits]).to_csv(
        os.path.join(OUTRES, "aggregate_curve_fits.csv"), index=False)
    styles = {"exp_neg": dict(c="k", ls="-", lw=2.4),
              "logarithmic": dict(c="#08519c", ls="--", lw=2.0),
              "quad_left": dict(c="#b30000", ls="-.", lw=2.0)}
    nice = {"exp_neg": "saturating exp", "logarithmic": "logarithmic",
            "quad_left": "concave quad (left half)"}
    for fdt in sorted(fits, key=lambda d: d["name"]):
        xg, yg = curve_xy(fdt)
        ax.plot(xg, yg, zorder=5, **styles[fdt["name"]],
                label=f"{nice[fdt['name']]}  (R²={fdt['r2']:.2f})")
    ax.legend(fontsize=7.5, loc="center left", bbox_to_anchor=(1.02, 0.5),
              framealpha=0.95, title="model type / fit", borderaxespad=0)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTFIG, "aggregate_persite_vs_global.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(OUTFIG, "aggregate_persite_vs_global.pdf"), bbox_inches="tight")
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
        fig, ax = plt.subplots(figsize=(8.4, 6.8))
        scatter(ax, d, dms_id, annotate_models=True)
        handles = [mlines.Line2D([], [], marker="o", ls="", color=PALETTE[mt], label=mt)
                   for mt in MTYPES if (d["model_type"] == mt).any()]
        ax.legend(handles=handles, fontsize=7, loc="center left",
                  bbox_to_anchor=(1.02, 0.5), framealpha=0.95,
                  title="model type", borderaxespad=0)
        fig.tight_layout()
        fig.savefig(os.path.join(OUTFIG, f"{dms_id}.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"assay fig: {dms_id}")


if __name__ == "__main__":
    main()
