"""
H2 "disagree-at-CONSERVED" counterpart of 38/39, for P53 and HXK4/GCK.

Same construction as 38/39 but the site filter is flipped to HIGH conservation:
per assay, take the protein's SHARED high-disagreement positions (flagged in >=2
of its assays, colocalization_shared_sites.csv), keep those with
conservation_percentile >= 80 (NaN/gappy columns are excluded -- they are the
low-conservation regime), require >= 5 missense substitutions.

Per-site model Spearman = spearman(model_score, DMS_score) over the site's subs
(scores as-is, already oriented). One box per model_type class (best-in-class
model, standard palette); each dot is one site, greyscaled by
|median(site DMS percentile) - WT DMS percentile| (white = close to WT, black =
far). Columns = UniProt (P53 | HXK4); rows = that protein's assays.

Outputs:
  results/persite_analysis/h2_highcons_shared_persite_spearman_by_site.csv
  figures/h2_persite_spearman_density/fig_shared_highcons_persite_spearman_panel.png
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import seaborn as sns
import common as C

sns.set_style("whitegrid")
rng = np.random.default_rng(0)

# reuse the helpers/paths from script 38 (offset, scores, models, WT, CONS, REF)
_spec = importlib.util.spec_from_file_location(
    "m38", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "38_h2_lowcons_persite_spearman_density.py"))
M38 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M38)

HIGH_PCT = 80
MIN_SUBS = 5

PROT_ASSAYS = {
    "P53_HUMAN": ["P53_HUMAN_Giacomelli_2018_Null_Etoposide",
                  "P53_HUMAN_Giacomelli_2018_Null_Nutlin",
                  "P53_HUMAN_Giacomelli_2018_WT_Nutlin",
                  "P53_HUMAN_Kotler_2018"],
    "HXK4_HUMAN": ["HXK4_HUMAN_Gersing_2022_activity",
                   "HXK4_HUMAN_Gersing_2023_abundance"],
}
COL_ORDER = ["P53_HUMAN", "HXK4_HUMAN"]

SHORT = {
    "MSA": "MSA", "Single sequence": "SeqOnly", "Structure": "Struct",
    "Single sequence & Structure": "Seq+Struct", "Structure & MSA": "Struct+MSA",
    "Single sequence, Structure & Function annotations": "Seq+Struct+Fn",
}
CMAP = cm.get_cmap("RdBu")    # signed: red = below WT (deleterious), blue = above, white = at WT


def shared_high_sites():
    """(DMS_id, up_pos, wt_aa) rows for P53/HXK4 shared positions."""
    sh = pd.read_csv(M38.SHARED)
    sh = sh[sh.UniProt_ID.isin(COL_ORDER)]
    rows = []
    for _, r in sh.iterrows():
        for dms in str(r.assays_flagging).split(";"):
            rows.append({"DMS_id": dms, "up_pos": int(r.up_pos),
                         "wt_aa": str(r.residue)[0]})
    return pd.DataFrame(rows)


def selected_high(dms_id, hd):
    off = M38._offset(dms_id)
    df = M38._scores(dms_id)
    ck = M38.CONS[M38.CONS.DMS_id == dms_id].set_index("up_pos")["conservation_percentile"]
    out = []
    for _, r in hd[hd.DMS_id == dms_id].iterrows():
        up = int(r.up_pos)
        cp = ck.get(up, np.nan)
        if not (pd.notna(cp) and cp >= HIGH_PCT):          # HIGH conservation only
            continue
        sub = df[df.p == up - off]
        if len(sub) < MIN_SUBS:
            continue
        out.append((up, r.wt_aa, cp, sub))
    return df, out


def build_table():
    models = M38.best_models_per_type()
    hd = shared_high_sites()
    rows = []
    for uni in COL_ORDER:
        for dms_id in PROT_ASSAYS[uni]:
            cat = M38.REF.loc[dms_id, "coarse_selection_type"]
            df, sites = selected_high(dms_id, hd)
            print(f"{dms_id}  [{cat}]  high-cons shared sites: {len(sites)}")
            if not sites:
                continue
            # WT scores are RAW paper values; orient by raw_DMS_directionality to
            # match the oriented DMS_score before ranking (flips the -1 P53 assays).
            wt_raw = M38.WT.get(dms_id, np.nan)
            ddir = M38.REF.loc[dms_id, "raw_DMS_directionality"]
            wt_dms = ddir * wt_raw
            wt_pct = 100.0 * (df["DMS_score"] < wt_dms).mean() if pd.notna(wt_dms) else np.nan
            for up, wt, cp, sub in sites:
                med_pct = sub["dms_pct"].median()
                diff = med_pct - wt_pct
                for mtype, m in models:
                    rows.append({"UniProt_ID": uni, "DMS_id": dms_id, "category": cat,
                                 "up_pos": up, "residue": f"{wt}{up}",
                                 "cons_pct": round(cp, 1), "n_subs": len(sub),
                                 "model_type": mtype, "model": m,
                                 "persite_spearman": M38.persite_spearman(sub, m),
                                 "site_median_dms_pct": round(med_pct, 2),
                                 "wt_dms_raw": wt_raw,
                                 "raw_DMS_directionality": int(ddir),
                                 "wt_dms_oriented": round(wt_dms, 4) if pd.notna(wt_dms) else np.nan,
                                 "wt_dms_pct": round(wt_pct, 2) if pd.notna(wt_pct) else np.nan,
                                 "dms_pct_diff_from_wt": round(diff, 2),
                                 "abs_dms_pct_diff_from_wt": round(abs(diff), 2)})
    tab = pd.DataFrame(rows)
    tab.to_csv(C.rp("h2_highcons_shared_persite_spearman_by_site.csv"), index=False)
    return tab, models


def main():
    df, models = build_table()
    order = [t for t in C.MODEL_TYPE_ORDER if t in set(df.model_type)]
    palette = dict(zip(C.MODEL_TYPE_ORDER,
                       sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))
    box_pal = {t: palette[t] for t in order}
    model_of = df.drop_duplicates("model_type").set_index("model_type")["model"].to_dict()
    MT = pd.read_csv(C.rp("directionality_check.csv"))
    overall_rho = dict(zip(MT.model, MT.mean_spearman_vs_DMS))
    M = float(df.dms_pct_diff_from_wt.abs().max())
    norm = Normalize(vmin=-M, vmax=M)          # symmetric, 0 -> white

    # P53 fills the first two columns (both rows); HXK4 fills the third column
    p53, hxk4 = PROT_ASSAYS["P53_HUMAN"], PROT_ASSAYS["HXK4_HUMAN"]
    placement = {
        (0, 0): p53[0], (0, 1): p53[1], (0, 2): hxk4[0],
        (1, 0): p53[2], (1, 1): p53[3], (1, 2): hxk4[1],
    }
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), sharey=True)
    for (r, c), dms_id in placement.items():
        ax = axes[r, c]
        d = df[df.DMS_id == dms_id]
        cat = d["category"].iloc[0]
        n_sites = d["up_pos"].nunique()
        wtp = d["wt_dms_pct"].iloc[0]
        sns.boxplot(data=d, x="model_type", y="persite_spearman", order=order,
                    hue="model_type", palette=box_pal, legend=False, ax=ax,
                    width=0.6, fliersize=0, linewidth=1,
                    boxprops=dict(alpha=0.35), showcaps=True)
        for xi, t in enumerate(order):
            dt = d[d.model_type == t]
            x = xi + rng.uniform(-0.16, 0.16, len(dt))
            ax.scatter(x, dt.persite_spearman,
                       c=CMAP(norm(dt.dms_pct_diff_from_wt.values)),
                       s=48, edgecolor="black", linewidth=0.6, zorder=3)
        ax.axhline(0, color="0.4", ls="--", lw=1)
        ax.set_title(f"{dms_id}\ncategory: {cat}   (n = {n_sites} sites; "
                     f"WT at pct {wtp:.0f})", fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("per-site Spearman (model vs DMS)")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([SHORT[t] for t in order], rotation=40,
                           ha="right", fontsize=8)
        ax.set_ylim(-1.05, 1.05)
        ax.yaxis.grid(False)

    fig.suptitle("Per-site Spearman at SHARED high-disagreement, HIGH-conservation "
                 "sites (cons pct ≥ 80; ≥5 subs)", fontsize=13, y=1.0)
    # DMS_score is pre-oriented higher=fitter for every assay (raw_DMS_directionality
    # already applied upstream), so the signed dot metric handles assay direction.
    # Caveat: WT is the fit anchor only where its percentile is high (PTEN/HXK4/CP2C9,
    # 69-88); in these P53 selection assays WT sits low (pct 0-48), so "below WT"
    # (red) does NOT read as deleterious for P53 -- see per-panel "WT at pct".
    fig.text(0.5, 0.01,
             "DMS_score is oriented higher = fitter; the raw paper WT score is oriented by "
             "raw_DMS_directionality to match before ranking.\nCaveat: WT is the clean fit "
             "anchor only where 'WT at pct' is high (Kotler 80); the Giacomelli selection "
             "screens put WT mid (pct 35–48), so read red/blue there relative to WT, not as "
             "absolute deleteriousness.",
             ha="center", va="bottom", fontsize=8, style="italic", color="0.25")
    plt.tight_layout(rect=(0, 0.17, 0.9, 0.96))

    sm = cm.ScalarMappable(norm=norm, cmap=CMAP)
    sm.set_array([])
    cax = fig.add_axes([0.925, 0.32, 0.014, 0.42])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("median(site DMS percentile) − WT DMS percentile\n"
                   "red = below WT (deleterious), white = at WT, blue = above WT",
                   fontsize=9)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=palette[t], alpha=0.5, ec="none",
                             label=f"{SHORT[t]} = {t} ({model_of.get(t,'?')}, "
                                   f"overall ρ={overall_rho.get(model_of.get(t),float('nan')):.2f})")
               for t in order]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.08),
               title="box colour = model_type class "
                     "(best-in-class model, overall Spearman over 217 assays)")
    out = os.path.join(C.FIGURES, "h2_persite_spearman_density",
                       "fig_shared_highcons_persite_spearman_panel.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("saved", out)


if __name__ == "__main__":
    main()
