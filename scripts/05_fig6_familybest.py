"""
Fig6 (per coarse_selection_type) rebuilt from a "family-best" model subset.

Instead of the 23 hand-picked models, we collapse every set of models that are
the same base method differing only by parameter size or single/ensemble variant
into one "family", and keep the family member with the highest OVERALL mean
Spearman (the Fig1 row ordering, results/spearman_model_by_taxon.csv:overall_mean).
Families never cross model_type boundaries (e.g. Tranception with vs without
retrieval, the *_MSA structure models, and xTrimoPGLM MLM vs CLM stay separate),
so the model_type aggregation stays clean.

Outputs:
  results/family_best_selection.csv
  figures/Fig6_by_selection_type_familybest.png
  results/persite_modeltype_agreement_by_selectiontype_familybest.csv
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
persite = import_module("03_persite_ranking")

# family_name -> member model columns (multi-member only; everything else is its own family)
FAMILIES = {
    "DeepSequence":             ["DeepSequence_single", "DeepSequence_ensemble"],
    "EVE":                      ["EVE_single", "EVE_ensemble"],
    "MSA_Transformer":          ["MSA_Transformer_single", "MSA_Transformer_ensemble"],
    "ESM1v":                    ["ESM1v_single", "ESM1v_ensemble"],
    "ESM2":                     ["ESM2_8M", "ESM2_35M", "ESM2_150M", "ESM2_650M", "ESM2_3B", "ESM2_15B"],
    "RITA":                     ["RITA_s", "RITA_m", "RITA_l", "RITA_xl"],
    "Progen2":                  ["Progen2_small", "Progen2_medium", "Progen2_base", "Progen2_large", "Progen2_xlarge"],
    "VESPA":                    ["VESPA", "VESPAl"],
    "Tranception_no_retrieval": ["Tranception_S_no_retrieval", "Tranception_M_no_retrieval", "Tranception_L_no_retrieval"],
    "Tranception_retrieval":    ["Tranception_S", "Tranception_M", "Tranception_L"],
    "TranceptEVE":              ["TranceptEVE_S", "TranceptEVE_M", "TranceptEVE_L"],
    "CARP":                     ["CARP_38M", "CARP_76M", "CARP_600K", "CARP_640M"],
    "ProtSSN":                  ["ProtSSN_k10_h512", "ProtSSN_k10_h768", "ProtSSN_k10_h1280",
                                 "ProtSSN_k20_h512", "ProtSSN_k20_h768", "ProtSSN_k20_h1280",
                                 "ProtSSN_k30_h512", "ProtSSN_k30_h768", "ProtSSN_k30_h1280", "ProtSSN_ensemble"],
    "SaProt":                   ["SaProt_650M_AF2", "SaProt_35M_AF2"],
    "ProSST":                   ["ProSST-20", "ProSST-128", "ProSST-512", "ProSST-1024", "ProSST-2048", "ProSST-4096"],
    "ESMC":                     ["ESMC-300M", "ESMC-600M"],
    "xTrimoPGLM-MLM":           ["xTrimoPGLM-1B-MLM", "xTrimoPGLM-3B-MLM", "xTrimoPGLM-10B-MLM"],
    "xTrimoPGLM-CLM":           ["xTrimoPGLM-1B-CLM", "xTrimoPGLM-3B-CLM", "xTrimoPGLM-7B-CLM", "xTrimoPGLM-100B-int4"],
    "Progen3":                  ["Progen3_112m", "Progen3_219m", "Progen3_339m", "Progen3_762m", "Progen3_1b", "Progen3_3b"],
}


def build_family_best():
    models, mtype, _ = C.model_meta()
    overall = (pd.read_csv(C.rp("spearman_model_by_taxon.csv"), index_col=0)
                 ["overall_mean"])
    fam_of = {}
    for fam, members in FAMILIES.items():
        for m in members:
            fam_of[m] = fam
    for m in models:                       # singletons are their own family
        fam_of.setdefault(m, m)

    # sanity: families must not span model types
    rows = []
    chosen = []
    for fam in dict.fromkeys(fam_of.values()):
        members = [m for m in models if fam_of[m] == fam]
        types = {mtype[m] for m in members}
        assert len(types) == 1, f"family {fam} spans model_types {types}"
        best = max(members, key=lambda m: overall.get(m, np.nan))
        chosen.append(best)
        rows.append({"family": fam, "model_type": mtype[best], "n_members": len(members),
                     "selected_model": best, "selected_overall_spearman": round(float(overall[best]), 3),
                     "members": ";".join(members)})
    sel = pd.DataFrame(rows).sort_values(["model_type", "selected_overall_spearman"],
                                         ascending=[True, False])
    sel.to_csv(C.rp("family_best_selection.csv"), index=False)
    return chosen, {m: mtype[m] for m in chosen}, sel


def main():
    chosen, type_of, sel = build_family_best()
    MTYPES = [t for t in C.MODEL_TYPE_ORDER if t in set(type_of.values())]
    members_by_type = {t: [m for m in chosen if type_of[m] == t] for t in MTYPES}
    print(f"Family-best subset: {len(chosen)} models from {sel['family'].nunique()} families")
    print(sel[["family", "selected_model", "model_type", "selected_overall_spearman", "n_members"]]
          .to_string(index=False))
    print("\nmodels per type:", {t: len(v) for t, v in members_by_type.items()})

    ref = C.load_reference()
    sel_of = dict(zip(ref["DMS_id"], ref["coarse_selection_type"]))
    sel_types = ref["coarse_selection_type"].value_counts().index.tolist()

    cols_t = chosen + ["DMS_score"]
    idx = {c: i for i, c in enumerate(cols_t)}
    K = len(cols_t)
    z_sum = {s: np.zeros((K, K)) for s in sel_types}
    z_cnt = {s: np.zeros((K, K)) for s in sel_types}
    assay_n = {s: 0 for s in sel_types}

    for path in C.list_score_files():
        dms_id = C.dms_id_from_path(path)
        s = sel_of.get(dms_id)
        if s not in z_sum:
            continue
        df = pd.read_csv(path, usecols=lambda c: c in C.NON_MODEL_COLS or c in set(chosen))
        pos = df["mutant"].map(C.parse_position)
        keep = pos.notna()
        if keep.sum() == 0:
            continue
        present = [m for m in chosen if m in df.columns] + ["DMS_score"]
        wr = persite.within_site_norm_ranks(df[keep], present, pos[keep])
        if wr.empty or len(wr) < 10:
            continue
        cm = wr.corr(method="pearson", min_periods=10)
        contributed = False
        for a in present:
            for b in present:
                v = cm.at[a, b]
                if pd.notna(v):
                    z_sum[s][idx[a], idx[b]] += np.arctanh(np.clip(v, -0.999999, 0.999999))
                    z_cnt[s][idx[a], idx[b]] += 1
                    contributed = True
        if contributed:
            assay_n[s] += 1

    def modeltype_matrix(corr_disp):
        tt = pd.DataFrame(index=MTYPES, columns=MTYPES, dtype=float)
        for ta in MTYPES:
            for tb in MTYPES:
                block = corr_disp.loc[members_by_type[ta], members_by_type[tb]].values.astype(float)
                if ta == tb:
                    mask = ~np.eye(len(members_by_type[ta]), dtype=bool)
                    vals = block[mask]
                else:
                    vals = block.ravel()
                tt.at[ta, tb] = np.nanmean(vals) if np.isfinite(vals).any() else np.nan
        return tt

    order = sorted(sel_types, key=lambda s: -assay_n[s])
    ncol = 3
    nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 5.2 * nrow))
    axes = np.atleast_1d(axes).ravel()
    lab = [t.replace(" & ", " &\n").replace(", ", ",\n") for t in MTYPES]
    tables = {}
    for k, s in enumerate(order):
        with np.errstate(invalid="ignore"):
            zmean = np.where(z_cnt[s] > 0, z_sum[s] / np.maximum(z_cnt[s], 1), np.nan)
        corr = pd.DataFrame(np.tanh(zmean), index=cols_t, columns=cols_t).loc[chosen, chosen]
        tt = modeltype_matrix(corr)
        tables[s] = tt
        ax = axes[k]
        sns.heatmap(tt.astype(float), cmap="rocket", annot=True, fmt=".2f", vmin=0, vmax=1,
                    cbar_kws={"label": "within-site agreement"}, ax=ax, annot_kws={"size": 7})
        ax.set_title(f"{s}  (n = {assay_n[s]} assays)", fontsize=11)
        ax.set_xticklabels(lab, fontsize=6, rotation=30, ha="right")
        ax.set_yticklabels(lab, fontsize=6, rotation=0)
    for j in range(len(order), len(axes)):
        axes[j].axis("off")
    fig.suptitle("Within-site mutation-ranking agreement between model TYPES, by selection type\n"
                 f"(family-best subset: {len(chosen)} models, one per family; diagonal = within-type, self excluded)",
                 fontsize=13, y=1.0)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(C.FIGURES, "Fig6_by_selection_type_familybest.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")
    big = pd.concat({s: tables[s] for s in order}, names=["selection_type", "model_type"])
    big.to_csv(C.rp("persite_modeltype_agreement_by_selectiontype_familybest.csv"))
    print("assay counts per selection type:", {s: assay_n[s] for s in order})


if __name__ == "__main__":
    main()
