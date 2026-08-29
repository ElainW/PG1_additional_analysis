"""
H1 follow-up (P53 for now): per-assay distribution of conservation across four
residue categories, plus the underlying table.

"Relative fraction of conservation" for a residue = `conservation_percentile` / 100
(its conservation rank among the assay's match columns, in [0,1]), with NaN (gappy /
insertion column, no match-state value) mapped to 0. This is one density x-axis; the
second panel uses the raw per-substitution DMS score.

Conservation category (per assay, from `conservation_percentile`, script 28):
    high = percentile >= 80 ,  low = percentile <= 20 or NaN ,  mid = between.
Disagreement label (per assay): high = in <gene>_high_disagreement, low = in
<gene>_low_disagreement, else mid.

Shared residues (across a protein's assays):
    shared HIGH  = high_disagreement_colocalization/colocalization_shared_sites.csv
                   (positions flagged high-disagreement in >=2 assays).
    shared LOW   = positions flagged low-disagreement in >=2 of the protein's assays
                   (computed here; no precomputed file).

Scope: P53_HUMAN, PTEN_HUMAN, HXK4_HUMAN (GCK).

Outputs:
  results/MSA_conservation/H1_conservation_by_category.csv
      one row per (assay, shared high/low-disagreement residue): up_pos, residue,
      disagreement, relative_fraction_of_conservation, cons_category.
  results/MSA_conservation/H1_highcons_middisagreement_range.csv
      per assay: n and min/max relative-fraction for highly-conserved (cons high)
      residues with MID disagreement.
  figures/conservation_hist/<DMS_id>.png
      per-assay density of relative-fraction-of-conservation for the 4 categories
      (shared high-disagreement, shared low-disagreement, highly conserved, all),
      legend annotated with each category's residue count.
"""
import os
import glob
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import common as C

PROTEINS = ["P53_HUMAN", "PTEN_HUMAN", "HXK4_HUMAN"]
GENE = {"P53_HUMAN": "p53", "PTEN_HUMAN": "pten", "HXK4_HUMAN": "gck"}  # tsv/fig basename
HI_PCT, LO_PCT = 80, 20

OUT = C.rp("MSA_conservation")
FIGDIR = os.path.join(C.FIGURES, "conservation_hist")
COLOC = os.path.join(C.rp("high_disagreement_colocalization"),
                     "colocalization_shared_sites.csv")

_S15 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "15_functional_disagreement_generalization.py")
_spec = importlib.util.spec_from_file_location("func_disagree_gen", _S15)
S15 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S15)

REF = C.load_reference().set_index("DMS_id")
CONS = pd.read_csv(os.path.join(OUT, "msa_conservation_by_site.csv"))


def clabel(pct):
    if pd.isna(pct) or pct <= LO_PCT:
        return "low"
    if pct >= HI_PCT:
        return "high"
    return "mid"


# wild-type DMS score to mark on the DMS-score panel (dashed line), per assay family
WT_DMS = {
    "P53_HUMAN_Giacomelli_2018_Null_Etoposide": 0.0,
    "P53_HUMAN_Giacomelli_2018_Null_Nutlin": 0.0,
    "P53_HUMAN_Giacomelli_2018_WT_Nutlin": 0.0,
    "P53_HUMAN_Kotler_2018": -2.4,
    "PTEN_HUMAN_Matreyek_2021": 1.0,
    "PTEN_HUMAN_Mighell_2018": 0.0,
    "HXK4_HUMAN_Gersing_2022_activity": 1.0,
    "HXK4_HUMAN_Gersing_2023_abundance": 1.0,
}

# per-assay display overrides
DMS_XMAX = {"PTEN_HUMAN_Matreyek_2021": 2.25}         # cap DMS-panel x upper limit
# legend location per assay as (conservation-panel, DMS-panel); else matplotlib "best"
LEGEND_LOC = {"PTEN_HUMAN_Matreyek_2021": ("upper center", "upper right")}


def assay_positions(dms_id):
    """All assayed up_pos for a DMS (from single-substitution mutants)."""
    off = (S15.features_for(REF.loc[dms_id].UniProt_ID, REF.loc[dms_id].target_seq) or [0])[0]
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"),
                     usecols=["mutant"])
    df = df[df.mutant.str.match(r'^[A-Z]\d+[A-Z]$')]
    pos = df.mutant.str.extract(r'^([A-Z])(\d+)[A-Z]$')
    wt = dict(zip(pos[1].astype(int) + off, pos[0]))          # up_pos -> wt aa
    return wt


def assay_dms(dms_id):
    """DataFrame(up_pos, DMS_score) over single-substitution variants."""
    off = (S15.features_for(REF.loc[dms_id].UniProt_ID, REF.loc[dms_id].target_seq) or [0])[0]
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"),
                     usecols=["mutant", "DMS_score"])
    df = df[df.mutant.str.match(r'^[A-Z]\d+[A-Z]$')].copy()
    df["up_pos"] = df.mutant.str.extract(r'^[A-Z](\d+)[A-Z]$')[0].astype(int) + off
    return df[["up_pos", "DMS_score"]]


def tail_positions(gene, tag):
    """{DMS_id: set(up_pos)} for a disagreement tail."""
    out = {}
    for f in glob.glob(os.path.join(OUT, f"{gene}_{tag}_disagreement_conservation.tsv")):
        t = pd.read_csv(f, sep="\t")
        for d, g in t.groupby("DMS_id"):
            out.setdefault(d, set()).update(g.up_pos.astype(int))
    return out


def build_protein(uni, gene):
    assays = sorted(REF.index[(REF.UniProt_ID == uni)].intersection(CONS.DMS_id.unique()))
    hi_tail = tail_positions(gene, "high")
    lo_tail = tail_positions(gene, "low")
    # shared HIGH from colocalization; shared LOW computed from low tails
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
            dis = ("high" if up in hi_tail.get(d, set())
                   else "low" if up in lo_tail.get(d, set()) else "mid")
            rows.append({"DMS_id": d, "up_pos": up, "residue": f"{aa}{up}",
                         "relfrac_cons": 0.0 if pd.isna(pct) else float(pct) / 100.0,
                         "cons_pct": pct, "cons_category": clabel(pct),
                         "disagreement": dis})
        per_assay[d] = pd.DataFrame(rows)
    return assays, per_assay, shared_high, shared_low


def make_table(assays, per_assay, shared_high, shared_low):
    rows = []
    for d in assays:
        df = per_assay[d].set_index("up_pos")
        for up in sorted(shared_high):
            if up in df.index:
                r = df.loc[up]
                rows.append({"DMS_id": d, "up_pos": up, "residue": r.residue,
                             "disagreement": "high",
                             "relative_fraction_of_conservation": round(r.relfrac_cons, 4),
                             "cons_category": r.cons_category})
        for up in sorted(shared_low):
            if up in df.index:
                r = df.loc[up]
                rows.append({"DMS_id": d, "up_pos": up, "residue": r.residue,
                             "disagreement": "low",
                             "relative_fraction_of_conservation": round(r.relfrac_cons, 4),
                             "cons_category": r.cons_category})
    return pd.DataFrame(rows)


def highcons_middis_range(assays, per_assay):
    rows = []
    for d in assays:
        df = per_assay[d]
        sub = df[(df.cons_category == "high") & (df.disagreement == "mid")]["relfrac_cons"]
        rows.append({"DMS_id": d, "n": len(sub),
                     "relfrac_min": round(sub.min(), 4) if len(sub) else np.nan,
                     "relfrac_max": round(sub.max(), 4) if len(sub) else np.nan})
    return pd.DataFrame(rows)


CATS = ["shared high-disagreement", "shared low-disagreement",
        "highly conserved", "all residues"]
COLORS = {"shared high-disagreement": "#c1121f", "shared low-disagreement": "#457b9d",
          "highly conserved": "#2a9d8f", "all residues": "#888888"}
# legend display names (internal category keys unchanged); high-disagreement reads
# as low-AGREEMENT and vice versa. "all residues" -> "all mutations" on the DMS panel.
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
        # too few points for a KDE: draw a tall, visible marker line
        ax.plot(vals, np.zeros_like(vals), marker="|", linestyle="none",
                color=COLORS[name], markersize=26, markeredgewidth=2.6,
                label=f"{disp} (n={len(vals)})")


def plot_assay(d, df, dms, shared_high, shared_low):
    catpos = _cat_positions(df, shared_high, shared_low)
    fig, (axc, axd) = plt.subplots(1, 2, figsize=(12, 4.4))

    # left: conservation percentile (0-100) density (one value per residue).
    # For the low-agreement category, separate residues with an actual percentile
    # from those with NaN conservation (gappy column, no match-state value).
    for name in CATS:
        sdf = df[df.up_pos.isin(catpos[name])]
        if name == "shared high-disagreement":     # = "shared low-agreement"
            real = sdf.loc[sdf.cons_pct.notna(), "cons_pct"].values
            n_nan = int(sdf.cons_pct.isna().sum())
            _kde(axc, name, real, clip=(0, 100),
                 label=f"{LABELS[name]} · ranked")
            if n_nan:
                axc.plot([0], [0], marker="x", linestyle="none", color=COLORS[name],
                         markersize=13, markeredgewidth=2.4,
                         label=f"{LABELS[name]} · NaN cons (n={n_nan})")
        else:
            sub = sdf["relfrac_cons"].values * 100.0
            _kde(axc, name, sub, clip=(0, 100), label=LABELS[name])
    axc.set_xlabel("conservation percentile  (0 if NaN)")
    axc.set_ylabel("density")
    axc.set_xlim(0, 100)

    # right: DMS-score density (all substitutions at each category's residues)
    for name in CATS:
        sub = dms[dms.up_pos.isin(catpos[name])]["DMS_score"].values
        disp = "all mutations" if name == "all residues" else LABELS[name]
        _kde(axd, name, sub, label=disp)
    wt = WT_DMS.get(d)
    if wt is not None:
        axd.axvline(wt, ls="--", color="black", lw=1.3, label=f"wild-type = {wt:g}")
    direction = REF.loc[d, "raw_DMS_directionality"] if "raw_DMS_directionality" in REF.columns else None
    fit = ("higher is more fit" if direction == 1
           else "lower is more fit" if direction == -1 else None)
    axd.set_xlabel(f"DMS score ({fit})" if fit else "DMS score")
    axd.set_ylabel("density")
    if d in DMS_XMAX:
        axd.set_xlim(right=DMS_XMAX[d])

    locs = LEGEND_LOC.get(d, ("best", "best"))
    for ax, title, loc in ((axc, "categories (number of positions)", locs[0]),
                           (axd, "categories (number of mutations)", locs[1])):
        ax.grid(False)
        leg = ax.legend(fontsize=8, title=title, title_fontsize=8, frameon=True,
                        facecolor="white", edgecolor="0.3", framealpha=1.0, loc=loc)
        leg.get_frame().set_boxstyle("square")
        # shrink the tall rug-marker handle in the legend only (keeps plot tick tall)
        for h in getattr(leg, "legend_handles", getattr(leg, "legendHandles", [])):
            if hasattr(h, "get_marker") and h.get_marker() == "|":
                h.set_markersize(9)

    fig.suptitle(d, y=1.02, fontsize=11)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    p = os.path.join(FIGDIR, f"{d}.png")
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return p


def main():
    tables, ranges = [], []
    for uni in PROTEINS:
        gene = GENE[uni]
        assays, per_assay, sh_hi, sh_lo = build_protein(uni, gene)
        print(f"{uni}: {len(assays)} assays; shared_high={len(sh_hi)} shared_low={len(sh_lo)}")
        tables.append(make_table(assays, per_assay, sh_hi, sh_lo))
        ranges.append(highcons_middis_range(assays, per_assay))
        for d in assays:
            p = plot_assay(d, per_assay[d], assay_dms(d), sh_hi, sh_lo)
            print("  wrote", p)
    tab = pd.concat(tables, ignore_index=True)
    tab.to_csv(os.path.join(OUT, "H1_conservation_by_category.csv"), index=False)
    rng = pd.concat(ranges, ignore_index=True)
    rng.to_csv(os.path.join(OUT, "H1_highcons_middisagreement_range.csv"), index=False)
    print("\nwrote H1_conservation_by_category.csv (%d rows)" % len(tab))
    print("wrote H1_highcons_middisagreement_range.csv")
    print("\nhighly-conserved + mid-disagreement relfrac range per assay:")
    print(rng.to_string(index=False))


if __name__ == "__main__":
    main()
