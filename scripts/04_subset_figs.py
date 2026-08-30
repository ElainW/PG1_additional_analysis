"""
Subset versions of fig5 and fig6, driven by the user-supplied 23-model dict.

Fig5_subset : clustered model x model within-site ranking-agreement heatmap,
              restricted to the 23 listed models. Reuses the global agreement
              matrix from 03_persite_ranking.py (results/persite_model_agreement_matrix.csv).

Fig6 (per selection type): model_type x model_type within-site agreement, built
              from ONLY the 23 listed models, computed SEPARATELY for each coarse
              selection type (OrganismalFitness, Stability, Activity, Expression,
              Binding). Each panel's title shows the number of assays.

Model_type and display names come from the user-supplied dict (not config). The
three names that differ from the CSV columns are mapped:
    ESM-IF1_142M -> ESM-IF1 ; ESM3_1.4B_open -> ESM3 ; Progen3_3B -> Progen3_3b
Scores are used as-is (already higher=fitter; see 00_check_directionality.py).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import common as C
from importlib import import_module

sns.set_style("whitegrid")
persite = import_module("03_persite_ranking")  # reuse within_site_norm_ranks / MIN_SITE_VARIANTS

FIGDIR = os.path.join(C.FIGURES, "persite_analysis")

# (display_name, csv_column, model_type)  -- order preserved for plots
SUBSET = [
    ("Site_Independent",           "Site_Independent",          "MSA"),
    ("EVmutation",                 "EVmutation",                "MSA"),
    ("EVE_ensemble",               "EVE_ensemble",              "MSA"),
    ("MSA_Transformer_ensemble",   "MSA_Transformer_ensemble",  "MSA"),
    ("ESM2_650M",                  "ESM2_650M",                 "Single sequence"),
    ("GEMME",                      "GEMME",                     "MSA"),
    ("VespaG",                     "VespaG",                    "Single sequence"),
    ("Tranception_L_no_retrieval", "Tranception_L_no_retrieval","Single sequence"),
    ("Tranception_L",              "Tranception_L",             "MSA"),
    ("TranceptEVE_L",              "TranceptEVE_L",             "MSA"),
    ("ESM-IF1_142M",               "ESM-IF1",                   "Structure"),
    ("ProteinMPNN",                "ProteinMPNN",               "Structure"),
    ("SaProt_650M_AF2",            "SaProt_650M_AF2",           "Single sequence & Structure"),
    ("PoET",                       "PoET",                      "MSA"),
    ("ESCOTT",                     "ESCOTT",                    "Structure & MSA"),
    ("VenusREM",                   "VenusREM",                  "Structure & MSA"),
    ("ProSST-2048",                "ProSST-2048",               "Single sequence & Structure"),
    ("RSALOR",                     "RSALOR",                    "Structure & MSA"),
    ("S3F",                        "S3F",                       "Single sequence & Structure"),
    ("S3F_MSA",                    "S3F_MSA",                   "Structure & MSA"),
    ("ESM3_1.4B_open",             "ESM3",                      "Single sequence, Structure & Function annotations"),
    ("Progen3_3B",                 "Progen3_3b",                "Single sequence"),
    ("AIDO.Protein-RAG-16B",       "AIDO.Protein-RAG-16B",      "Structure & MSA"),
]

DISP = [s[0] for s in SUBSET]
COLS = [s[1] for s in SUBSET]
TYPE = {s[0]: s[2] for s in SUBSET}
COL2DISP = {s[1]: s[0] for s in SUBSET}
MTYPES = [t for t in C.MODEL_TYPE_ORDER if t in {s[2] for s in SUBSET}]
PALETTE = dict(zip(C.MODEL_TYPE_ORDER, sns.color_palette("tab10", len(C.MODEL_TYPE_ORDER))))


def fig5_subset():
    mat = pd.read_csv(C.rp("persite_model_agreement_matrix.csv"), index_col=0)
    sub = mat.loc[COLS, COLS].copy()
    sub.index = [COL2DISP[c] for c in sub.index]
    sub.columns = [COL2DISP[c] for c in sub.columns]
    sub = sub.astype(float).fillna(0.0)

    row_colors = pd.Series({d: PALETTE[TYPE[d]] for d in sub.index})
    g = sns.clustermap(sub, cmap="rocket", figsize=(11, 11),
                       row_colors=row_colors, col_colors=row_colors,
                       vmin=0, vmax=1, annot=True, fmt=".2f", annot_kws={"size": 5},
                       xticklabels=True, yticklabels=True,
                       cbar_kws={"label": "within-site rank agreement"})
    g.ax_heatmap.tick_params(labelsize=7)
    for t in MTYPES:
        g.ax_col_dendrogram.bar(0, 0, color=PALETTE[t], label=t, linewidth=0)
    g.ax_col_dendrogram.legend(title="model_type", loc="center", ncol=2, fontsize=7)
    g.fig.suptitle("Within-site mutation-ranking agreement (23-model subset, clustered)", y=1.02)
    os.makedirs(FIGDIR, exist_ok=True)
    out = os.path.join(FIGDIR, "Fig5_subset.png")
    g.savefig(out, dpi=150, bbox_inches="tight")
    g.savefig(os.path.join(FIGDIR, "Fig5_subset.pdf"), bbox_inches="tight")
    plt.close()
    print(f"Saved {out} (+ .pdf)")


def modeltype_matrix(model_corr):
    """Aggregate a model x model (display-named) matrix to model_type x model_type,
    excluding self-correlations on the within-type diagonal."""
    tt = pd.DataFrame(index=MTYPES, columns=MTYPES, dtype=float)
    members = {t: [d for d in DISP if TYPE[d] == t] for t in MTYPES}
    for ta in MTYPES:
        for tb in MTYPES:
            block = model_corr.loc[members[ta], members[tb]].values.astype(float)
            if ta == tb:
                mask = ~np.eye(len(members[ta]), dtype=bool)
                vals = block[mask]
            else:
                vals = block.ravel()
            tt.at[ta, tb] = np.nanmean(vals) if np.isfinite(vals).any() else np.nan
    return tt


def fig6_by_selection_type():
    ref = C.load_reference()
    sel_of = dict(zip(ref["DMS_id"], ref["coarse_selection_type"]))
    sel_types = (ref["coarse_selection_type"].value_counts().index.tolist())

    K = len(COLS) + 1
    cols_t = COLS + ["DMS_score"]
    idx = {c: i for i, c in enumerate(cols_t)}
    z_sum = {s: np.zeros((K, K)) for s in sel_types}
    z_cnt = {s: np.zeros((K, K)) for s in sel_types}
    assay_n = {s: 0 for s in sel_types}

    for path in C.list_score_files():
        dms_id = C.dms_id_from_path(path)
        sel = sel_of.get(dms_id)
        if sel not in z_sum:
            continue
        df = pd.read_csv(path, usecols=lambda c: c in C.NON_MODEL_COLS or c in set(COLS))
        pos = df["mutant"].map(C.parse_position)
        keep = pos.notna()
        if keep.sum() == 0:
            continue
        wr = persite.within_site_norm_ranks(df[keep], cols_t, pos[keep])
        if wr.empty or len(wr) < 10:
            continue
        cm = wr.corr(method="pearson", min_periods=10)
        contributed = False
        for a in cols_t:
            for b in cols_t:
                v = cm.at[a, b]
                if pd.notna(v):
                    z_sum[sel][idx[a], idx[b]] += np.arctanh(np.clip(v, -0.999999, 0.999999))
                    z_cnt[sel][idx[a], idx[b]] += 1
                    contributed = True
        if contributed:
            assay_n[sel] += 1

    # one heatmap per selection type, ordered by assay count (desc)
    order = sorted(sel_types, key=lambda s: -assay_n[s])
    ncol = 3
    nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 5.2 * nrow))
    axes = np.atleast_1d(axes).ravel()
    lab = [t.replace(" & ", " &\n").replace(", ", ",\n") for t in MTYPES]
    saved_tables = {}
    for k, sel in enumerate(order):
        with np.errstate(invalid="ignore"):
            zmean = np.where(z_cnt[sel] > 0, z_sum[sel] / np.maximum(z_cnt[sel], 1), np.nan)
        corr = pd.DataFrame(np.tanh(zmean), index=cols_t, columns=cols_t)
        corr_disp = corr.loc[COLS, COLS]
        corr_disp.index = [COL2DISP[c] for c in COLS]
        corr_disp.columns = [COL2DISP[c] for c in COLS]
        tt = modeltype_matrix(corr_disp)
        saved_tables[sel] = tt
        ax = axes[k]
        sns.heatmap(tt.astype(float), cmap="rocket", annot=True, fmt=".2f", vmin=0, vmax=1,
                    cbar_kws={"label": "within-site agreement"}, ax=ax, annot_kws={"size": 7})
        ax.set_title(f"{sel}  (n = {assay_n[sel]} assays)", fontsize=11)
        ax.set_xticklabels(lab, fontsize=6, rotation=30, ha="right")
        ax.set_yticklabels(lab, fontsize=6, rotation=0)
    for j in range(len(order), len(axes)):
        axes[j].axis("off")
    fig.suptitle("Within-site mutation-ranking agreement between model TYPES, by selection type\n"
                 "(23-model subset; diagonal = within-type, self excluded)", fontsize=13, y=1.0)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(FIGDIR, exist_ok=True)
    out = os.path.join(FIGDIR, "Fig6_by_selection_type.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")

    # also persist the per-selection-type tables
    big = pd.concat({s: saved_tables[s] for s in order}, names=["selection_type", "model_type"])
    big.to_csv(C.rp("persite_modeltype_agreement_by_selectiontype.csv"))
    print("assay counts per selection type:", {s: assay_n[s] for s in order})


if __name__ == "__main__":
    fig5_subset()
    fig6_by_selection_type()
