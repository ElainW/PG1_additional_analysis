"""
Merged per-assay conservation + DMS panels (consolidates the plotting of scripts 35
and 36). One 2x2 four-panel figure per assay, across the four categories:
    shared high-disagreement (= "low-agreement"), shared low-disagreement
    (= "high-agreement"), highly conserved (conservation_percentile >= 80), all.

Panels (per assay):
  (1) conservation percentile [0-100] density, one value per residue. For the
      low-agreement category, residues with a real percentile (ranked KDE) are
      separated from those with NaN conservation (gappy/lower-case column) shown as
      a marker. [from script 35]
  (2) raw DMS-score density over ALL single-substitution variants at each category's
      residues; wild-type DMS score dashed, x-axis labelled with the fitness
      direction (raw_DMS_directionality). [from script 35]
  (3) residue-level DMS-score percentile [0-100] density: each residue = mean
      DMS_score over its variants, ranked across residues; wild-type marked at the
      percentile of its DMS score among the residue means. [from script 36]
  (4) per-residue DMS-score percentile SPREAD [0-100] density = (max - min) of the
      substitution-level DMS percentiles at each residue. [from script 36]

Shared residues (per protein): shared HIGH from high_disagreement_colocalization/
colocalization_shared_sites.csv (flagged high-disagreement in >=2 assays); shared
LOW computed here (flagged low-disagreement in >=2 of the protein's assays).

Scope: PTEN_HUMAN, CP2C9_HUMAN, HXK4_HUMAN, P53_HUMAN.

Note: the H1 tables (H1_conservation_by_category*, H1_highcons_middisagreement_range*)
are still produced by script 35; this script only draws the merged figures and a
residue percentile/spread summary.

Outputs:
  figures/conservation_dms_panels/<DMS_id>.png   the 2x2 four-panel figure.
  results/conservation_dms_panels/dms_percentile_summary.csv
      per (assay, category): n residues, mean/median/min/max residue DMS percentile,
      and mean/median per-residue percentile spread.
"""
import os
import glob
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import common as C

PROTEINS = ["PTEN_HUMAN", "CP2C9_HUMAN", "HXK4_HUMAN", "P53_HUMAN"]
GENE = {"PTEN_HUMAN": "pten", "CP2C9_HUMAN": "cyp2c9",
        "HXK4_HUMAN": "gck", "P53_HUMAN": "p53"}  # tsv basename
HI_PCT, LO_PCT = 80, 20

OUT = C.rp("MSA_conservation")
PPDIR = os.path.join(OUT, "per_protein_conservation")  # 28's per-protein tails
FIGDIR = os.path.join(C.FIGURES, "conservation_dms_panels")
RESDIR = os.path.join(C.RESULTS, "conservation_dms_panels")
COLOC = os.path.join(C.rp("high_disagreement_colocalization"),
                     "colocalization_shared_sites.csv")

_S15 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "15_functional_disagreement_generalization.py")
_spec = importlib.util.spec_from_file_location("func_disagree_gen", _S15)
S15 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S15)

REF = C.load_reference().set_index("DMS_id")
CONS = pd.read_csv(os.path.join(OUT, "msa_conservation_by_site.csv"))
WT_DMS = C.load_wt_dms()

DMS_XMAX = {"PTEN_HUMAN_Matreyek_2021": 2.25}         # cap raw-DMS panel x upper limit


def clabel(pct):
    if pd.isna(pct) or pct <= LO_PCT:
        return "low"
    if pct >= HI_PCT:
        return "high"
    return "mid"


def _offset(dms_id):
    return (S15.features_for(REF.loc[dms_id].UniProt_ID, REF.loc[dms_id].target_seq) or [0])[0]


def assay_positions(dms_id):
    """{up_pos: wt aa} for a DMS (from single-substitution mutants)."""
    off = _offset(dms_id)
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"), usecols=["mutant"])
    df = df[df.mutant.str.match(r'^[A-Z]\d+[A-Z]$')]
    pos = df.mutant.str.extract(r'^([A-Z])(\d+)[A-Z]$')
    return dict(zip(pos[1].astype(int) + off, pos[0]))


def assay_dms_sub(dms_id):
    """DataFrame(up_pos, DMS_score) over single-substitution variants (per mutation)."""
    off = _offset(dms_id)
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"),
                     usecols=["mutant", "DMS_score"])
    df = df[df.mutant.str.match(r'^[A-Z]\d+[A-Z]$')].copy()
    df["up_pos"] = df.mutant.str.extract(r'^[A-Z](\d+)[A-Z]$')[0].astype(int) + off
    return df[["up_pos", "DMS_score"]].dropna(subset=["DMS_score"])


def assay_dms_res(dms_id):
    """DataFrame(up_pos, dms_score, dms_percentile, pct_spread) per assayed residue.

    dms_score      = mean DMS_score over the residue's variants.
    dms_percentile = percentile rank of that residue mean across residues, [0,100].
    pct_spread     = (max - min) of the substitution-level DMS percentiles at the
                     residue (each substitution ranked across all the assay's
                     substitutions); 0 for a single-variant residue."""
    df = assay_dms_sub(dms_id).copy()
    df["sub_pct"] = df["DMS_score"].rank(pct=True) * 100.0
    g = df.groupby("up_pos")
    res = g["DMS_score"].mean().rename("dms_score").reset_index()
    res["pct_spread"] = (g["sub_pct"].max() - g["sub_pct"].min()).values
    res["dms_percentile"] = res["dms_score"].rank(pct=True) * 100.0
    return res[["up_pos", "dms_score", "dms_percentile", "pct_spread"]]


def tail_positions(gene, tag):
    """{DMS_id: set(up_pos)} for a disagreement tail."""
    out = {}
    for f in glob.glob(os.path.join(PPDIR, f"{gene}_{tag}_disagreement_conservation.tsv")):
        t = pd.read_csv(f, sep="\t")
        for d, g in t.groupby("DMS_id"):
            out.setdefault(d, set()).update(g.up_pos.astype(int))
    return out


def build_protein(uni, gene):
    assays = sorted(REF.index[(REF.UniProt_ID == uni)].intersection(CONS.DMS_id.unique()))
    hi_tail = tail_positions(gene, "high")
    lo_tail = tail_positions(gene, "low")
    coloc = pd.read_csv(COLOC)
    shared_high = set(coloc.loc[coloc.UniProt_ID == uni, "up_pos"].astype(int))
    low_counts = {}
    for d in assays:
        for up in lo_tail.get(d, set()):
            low_counts[up] = low_counts.get(up, 0) + 1
    shared_low = {up for up, n in low_counts.items() if n >= 2}

    per_assay = {}
    for d in assays:
        wt = assay_positions(d)
        ck = CONS[CONS.DMS_id == d].set_index("up_pos")
        rows = []
        for up, aa in wt.items():
            pct = ck.conservation_percentile.get(up, np.nan)
            rows.append({"up_pos": up, "residue": f"{aa}{up}",
                         "relfrac_cons": 0.0 if pd.isna(pct) else float(pct) / 100.0,
                         "cons_pct": pct, "cons_category": clabel(pct)})
        per_assay[d] = pd.DataFrame(rows)
    return assays, per_assay, shared_high, shared_low


CATS = ["shared high-disagreement", "shared low-disagreement",
        "highly conserved", "all residues"]
COLORS = {"shared high-disagreement": "#c1121f", "shared low-disagreement": "#457b9d",
          "highly conserved": "#2a9d8f", "all residues": "#888888"}
# legend display names; high-disagreement reads as low-AGREEMENT and vice versa.
LABELS = {"shared high-disagreement": "shared low-agreement",
          "shared low-disagreement": "shared high-agreement",
          "highly conserved": "highly conserved",
          "all residues": "all residues"}


def _cat_positions(df, shared_high, shared_low):
    return {
        "shared high-disagreement": set(df.loc[df.up_pos.isin(shared_high), "up_pos"]),
        "shared low-disagreement":  set(df.loc[df.up_pos.isin(shared_low), "up_pos"]),
        "highly conserved":         set(df.loc[df.cons_category == "high", "up_pos"]),
        "all residues":             set(df.up_pos),
    }


def _kde(ax, name, vals, clip=None, label=None):
    """Smooth density; fall back to a rug if too few points for a KDE."""
    disp = label if label is not None else name
    vals = np.asarray(vals, float)
    vals = vals[~np.isnan(vals)]
    if len(vals) >= 2 and np.ptp(vals) > 0:
        sns.kdeplot(x=vals, ax=ax, color=COLORS[name], fill=True, alpha=0.25,
                    linewidth=1.8, bw_adjust=0.9, clip=clip, cut=0,
                    label=f"{disp} (n={len(vals)})")
    elif len(vals):
        ax.plot(vals, np.zeros_like(vals), marker="|", linestyle="none",
                color=COLORS[name], markersize=26, markeredgewidth=2.6,
                label=f"{disp} (n={len(vals)})")


def _style_legend(ax, title):
    ax.grid(False)
    leg = ax.legend(fontsize=7.5, title=title, title_fontsize=7.5, frameon=True,
                    facecolor="white", edgecolor="0.3", framealpha=1.0, loc="best")
    leg.get_frame().set_boxstyle("square")
    for h in getattr(leg, "legend_handles", getattr(leg, "legendHandles", [])):
        if hasattr(h, "get_marker") and h.get_marker() == "|":
            h.set_markersize(9)


def plot_assay(d, df, dms_sub, dms_res, shared_high, shared_low):
    catpos = _cat_positions(df, shared_high, shared_low)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(13.5, 9))

    # (1) conservation percentile density (one value per residue)
    for name in CATS:
        sdf = df[df.up_pos.isin(catpos[name])]
        if name == "shared high-disagreement":     # = "shared low-agreement"
            real = sdf.loc[sdf.cons_pct.notna(), "cons_pct"].values
            n_nan = int(sdf.cons_pct.isna().sum())
            _kde(ax1, name, real, clip=(0, 100), label=f"{LABELS[name]} · ranked")
            if n_nan:
                ax1.plot([0], [0], marker="x", linestyle="none", color=COLORS[name],
                         markersize=13, markeredgewidth=2.4,
                         label=f"{LABELS[name]} · NaN cons (n={n_nan})")
        else:
            _kde(ax1, name, sdf["relfrac_cons"].values * 100.0, clip=(0, 100),
                 label=LABELS[name])
    ax1.set_xlabel("conservation percentile  (0 if NaN)")
    ax1.set_ylabel("density")
    ax1.set_xlim(0, 100)
    _style_legend(ax1, "categories (number of positions)")

    # (2) raw DMS-score density (all substitutions at each category's residues)
    for name in CATS:
        sub = dms_sub[dms_sub.up_pos.isin(catpos[name])]["DMS_score"].values
        disp = "all mutations" if name == "all residues" else LABELS[name]
        _kde(ax2, name, sub, label=disp)
    wt = WT_DMS.get(d)
    if wt is not None:
        ax2.axvline(wt, ls="--", color="black", lw=1.3, label=f"wild-type = {wt:g}")
    direction = REF.loc[d, "raw_DMS_directionality"] if "raw_DMS_directionality" in REF.columns else None
    fit = ("higher is higher fitness" if direction == 1
           else "lower is higher fitness" if direction == -1 else None)
    ax2.set_xlabel(f"DMS score ({fit})" if fit else "DMS score")
    ax2.set_ylabel("density")
    if d in DMS_XMAX:
        ax2.set_xlim(right=DMS_XMAX[d])
    _style_legend(ax2, "categories (number of mutations)")

    # (3) residue-level DMS-score percentile density (residue mean, ranked)
    for name in CATS:
        sub = dms_res[dms_res.up_pos.isin(catpos[name])]["dms_percentile"].values
        _kde(ax3, name, sub, clip=(0, 100), label=LABELS[name])
    if wt is not None:
        wt_pct = percentileofscore(dms_res["dms_score"].values, wt, kind="mean")
        ax3.axvline(wt_pct, ls="--", color="black", lw=1.3,
                    label=f"wild-type (pct={wt_pct:.0f})")
    ax3.set_xlabel("DMS-score percentile  (residue mean, per-assay rank, 0–100"
                   + (f"; {fit})" if fit else ")"))
    ax3.set_ylabel("density")
    ax3.set_xlim(0, 100)
    _style_legend(ax3, "categories (number of residues)")

    # (4) per-residue DMS-score percentile spread density (max - min over substitutions)
    for name in CATS:
        sub = dms_res[dms_res.up_pos.isin(catpos[name])]["pct_spread"].values
        _kde(ax4, name, sub, clip=(0, 100), label=LABELS[name])
    ax4.set_xlabel("DMS-score percentile spread  (max − min over substitutions, 0–100)")
    ax4.set_ylabel("density")
    ax4.set_xlim(0, 100)
    _style_legend(ax4, "categories (number of residues)")

    fig.suptitle(d, y=1.0, fontsize=12)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    p = os.path.join(FIGDIR, f"{d}.png")
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return p, catpos


def summarize(d, dms_res, catpos):
    rows = []
    for name in CATS:
        sdf = dms_res[dms_res.up_pos.isin(catpos[name])]
        vals, spr = sdf["dms_percentile"].values, sdf["pct_spread"].values
        rows.append({"DMS_id": d, "category": name, "n_residues": len(vals),
                     "pct_mean": round(float(np.mean(vals)), 2) if len(vals) else np.nan,
                     "pct_median": round(float(np.median(vals)), 2) if len(vals) else np.nan,
                     "pct_min": round(float(np.min(vals)), 2) if len(vals) else np.nan,
                     "pct_max": round(float(np.max(vals)), 2) if len(vals) else np.nan,
                     "spread_mean": round(float(np.mean(spr)), 2) if len(spr) else np.nan,
                     "spread_median": round(float(np.median(spr)), 2) if len(spr) else np.nan})
    return rows


def main():
    summary = []
    for uni in PROTEINS:
        gene = GENE[uni]
        assays, per_assay, sh_hi, sh_lo = build_protein(uni, gene)
        print(f"{uni}: {len(assays)} assays; shared_high={len(sh_hi)} shared_low={len(sh_lo)}")
        for d in assays:
            dms_res = assay_dms_res(d)
            p, catpos = plot_assay(d, per_assay[d], assay_dms_sub(d), dms_res, sh_hi, sh_lo)
            summary.extend(summarize(d, dms_res, catpos))
            print("  wrote", p)
    os.makedirs(RESDIR, exist_ok=True)
    sdf = pd.DataFrame(summary)
    sp = os.path.join(RESDIR, "dms_percentile_summary.csv")
    sdf.to_csv(sp, index=False)
    print("\nwrote", sp, f"({len(sdf)} rows)")


if __name__ == "__main__":
    main()
