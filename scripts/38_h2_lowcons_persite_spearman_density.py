"""
H2 plots: per-site model Spearman at high-disagreement / low-conservation sites,
one density figure PER ASSAY, for PTEN, GCK (HXK4), and CP2C9.

Site selection (per assay), mirroring the H2 low-conservation regime (README):
  1. Start from the protein's SHARED high-disagreement sites -- positions flagged
     in >=2 of its assays (script 32,
     high_disagreement_colocalization/colocalization_shared_sites.csv). For these
     2-assay proteins a shared position is flagged by BOTH assays, so both assays
     of a protein draw from the same position list.
  2. Keep only LOW-conservation sites: conservation_percentile < 20 (percentile
     among the assay's MSA match columns, script 28) OR NaN (gappy/insertion
     column, i.e. not a match column at all).
  3. Require >= 5 missense substitutions at the site (alt != wt), as in the
     HXK4 W99 deep-dive (07).

Per-site model Spearman = spearman(model_score, DMS_score) over the site's
substitutions. Scores are used AS-IS (already oriented higher = fitter for every
model; directionality was folded in upstream -- see common.py / 00), matching
ProteinGym's performance_DMS_benchmarks.py and script 34's H2.

Curves: the single HIGHEST-performing model in each model_type class
(by overall mean_spearman_vs_DMS in directionality_check.csv). MSA uses PoET
(Protriever has the highest MSA mean but only over 200 assays and is absent from
the HXK4 files, so PoET -- best over the full 217 and present everywhere -- keeps
the MSA class comparable across assays).

Colours: the standard model_type palette (tab10 in common.MODEL_TYPE_ORDER),
matching scripts 04/07.

Outputs:
  figures/h2_persite_spearman_density/fig_<DMS_id>_lowcons_persite_spearman.png
  results/MSA_conservation/H2/h2_lowcons_shared_persite_spearman_by_site.csv
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import common as C

sns.set_style("white")

LOW_PCT = 20          # low-conservation cutoff (percentile < 20, or NaN)
MIN_SUBS = 5          # >= 5 missense substitutions per site

PROTEINS = ["PTEN_HUMAN", "HXK4_HUMAN", "CP2C9_HUMAN"]   # GCK == HXK4
SHARED = os.path.join(C.RESULTS, "high_disagreement_colocalization",
                      "colocalization_shared_sites.csv")

# ---- shared helpers (offset mapping) reused from script 15 ------------------
_S15 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "15_functional_disagreement_generalization.py")
_spec = importlib.util.spec_from_file_location("func_disagree_gen", _S15)
S15 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S15)

REF = C.load_reference().set_index("DMS_id")
CONS = pd.read_csv(os.path.join(C.RESULTS, "MSA_conservation",
                                "msa_conservation_by_site.csv"))
MT = pd.read_csv(C.rp("directionality_check.csv"))
WT = pd.read_csv(C.WT_DMS_PATH).set_index("DMS_id")["wt_score"]


def best_models_per_type():
    """Highest overall-Spearman model per model_type; force MSA->PoET so the
    MSA class is measured on the full 217 assays and is present everywhere."""
    idx = MT.groupby("model_type")["mean_spearman_vs_DMS"].idxmax()
    best = MT.loc[idx].set_index("model_type")["model"].to_dict()
    best["MSA"] = "PoET"
    # order by the canonical model_type order
    return [(t, best[t]) for t in C.MODEL_TYPE_ORDER if t in best]


def _offset(dms_id):
    r = REF.loc[dms_id]
    f = S15.features_for(r.UniProt_ID, r.target_seq)
    return f[0] if f else 0


def _scores(dms_id):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
    df = df[df.mutant.str.match(r'^[A-Z]\d+[A-Z]$')].copy()
    df["p"] = df.mutant.str.extract(r'^[A-Z](\d+)[A-Z]$')[0].astype(int)
    df["wt"] = df.mutant.str[0]
    df["alt"] = df.mutant.str[-1]
    df = df[df.alt != df.wt].copy()            # missense only
    # DMS-score percentile over the entire assay's DMS range (0..100)
    df["dms_pct"] = df["DMS_score"].rank(pct=True) * 100
    return df


def selected_sites(dms_id, hd):
    """Low-conservation (pct<20 or NaN) high-disagreement sites with >=5 subs."""
    off = _offset(dms_id)
    df = _scores(dms_id)
    ck = CONS[CONS.DMS_id == dms_id].set_index("up_pos")["conservation_percentile"]
    g = hd[hd.DMS_id == dms_id]
    out = []
    for _, r in g.iterrows():
        up = int(r.up_pos)
        cp = ck.get(up, np.nan)
        if not (pd.isna(cp) or cp < LOW_PCT):
            continue
        sub = df[df.p == up - off]
        if len(sub) < MIN_SUBS:
            continue
        out.append((up, r.wt_aa, cp, sub))
    return df, out


def persite_spearman(sub, model):
    if sub[model].nunique() < 2 or sub["DMS_score"].nunique() < 2:
        return np.nan
    return spearmanr(sub[model], sub["DMS_score"])[0]


def shared_sites_table():
    """Explode colocalization_shared_sites.csv into (DMS_id, up_pos, wt_aa) rows
    for our 3 proteins -- one row per (assay, shared position)."""
    sh = pd.read_csv(SHARED)
    sh = sh[sh.UniProt_ID.isin(PROTEINS)].copy()
    rows = []
    for _, r in sh.iterrows():
        for dms in str(r.assays_flagging).split(";"):
            rows.append({"UniProt_ID": r.UniProt_ID, "DMS_id": dms,
                         "up_pos": int(r.up_pos), "wt_aa": str(r.residue)[0]})
    return pd.DataFrame(rows)


def main():
    """Build the per-site table (one row per selected site x best-in-class model),
    including each site's DMS-median distance to WT. The 2x3 panel is drawn by
    script 39 from this CSV."""
    models = best_models_per_type()
    all_rows = []

    hd = shared_sites_table()
    for uni in PROTEINS:
        prot_hd = hd[hd.UniProt_ID == uni]
        for dms_id in sorted(prot_hd.DMS_id.unique()):
            cat = REF.loc[dms_id, "coarse_selection_type"]
            df, sites = selected_sites(dms_id, hd)
            print(f"\n{dms_id}  [{cat}]  selected sites: {len(sites)}")
            if not sites:
                continue
            # WT scores in DMS_wt_score.csv are the RAW paper values; DMS_score is
            # already oriented higher=fitter, so orient the WT by raw_DMS_directionality
            # before ranking it (a no-op for +1 assays; flips the -1 ones).
            wt_raw = WT.get(dms_id, np.nan)
            ddir = REF.loc[dms_id, "raw_DMS_directionality"]
            wt_dms = ddir * wt_raw
            # WT DMS-score percentile within the assay's whole DMS range (0..100)
            wt_pct = 100.0 * (df["DMS_score"] < wt_dms).mean() if pd.notna(wt_dms) else np.nan

            for up, wt, cp, sub in sites:
                med_pct = sub["dms_pct"].median()        # site median DMS percentile
                diff = med_pct - wt_pct                   # signed (site median - WT)
                for mtype, m in models:
                    rho = persite_spearman(sub, m)
                    all_rows.append({"DMS_id": dms_id, "category": cat,
                                     "up_pos": up, "residue": f"{wt}{up}",
                                     "cons_pct": round(cp, 1) if pd.notna(cp) else np.nan,
                                     "n_subs": len(sub), "model_type": mtype,
                                     "model": m, "persite_spearman": rho,
                                     "site_median_dms_pct": round(med_pct, 2),
                                     "wt_dms_raw": wt_raw,
                                     "raw_DMS_directionality": int(ddir),
                                     "wt_dms_oriented": round(wt_dms, 4) if pd.notna(wt_dms) else np.nan,
                                     "wt_dms_pct": round(wt_pct, 2) if pd.notna(wt_pct) else np.nan,
                                     "dms_pct_diff_from_wt": round(diff, 2),
                                     "abs_dms_pct_diff_from_wt": round(abs(diff), 2)})

    tab = pd.DataFrame(all_rows)
    out_csv = os.path.join(C.RESULTS, "MSA_conservation", "H2",
                           "h2_lowcons_shared_persite_spearman_by_site.csv")
    tab.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}  ({len(tab)} rows)")
    # quick class summary per assay
    print("\nmean per-site Spearman by assay x model_type:")
    print(tab.groupby(["DMS_id", "model_type"])["persite_spearman"]
             .mean().round(3).to_string())


if __name__ == "__main__":
    main()
