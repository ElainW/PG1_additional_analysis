"""
Mechanism of the per-assay disagreement <-> conservation sign: H1 (experimental
DMS-score spread) and H2 (per-site model-class accuracy), with a STRICT
conservation split.

High/low conservation are defined on `conservation_percentile` (percentile of a
column's weighted conservation among that assay's match columns, script 28):
    high  = conservation_percentile >= HIGH_PCT
    low   = conservation_percentile <= LOW_PCT  OR  NaN (gappy/insertion column)
    mid   = everything in between (excluded from the strict contrast)
Run for two threshold pairs: (80, 20) and (90, 10).

H1 (P53 Giacomelli, the disagree-at-CONSERVED case). Per site, DMS "relative
fraction" = max - min of its substitutions' DMS-percentiles within the assay
(fraction of the assay's DMS range the site's mutations span; ProteinGym assays
are pre-selected for a decent overall range). Compare relfrac at high- vs
low-disagreement selected sites, split by conservation, against the assay-wide
baseline over all high- / low-conservation match sites.

H2 (per-site model-class accuracy). For shared high-disagreement sites (flagged in
>=2 assays, `high_disagreement_colocalization/colocalization_shared_sites.csv`),
per-site Spearman of each model's scores vs DMS across a position's substitutions,
grouped into sequence+MSA vs structure-only classes
(`results/qc_directionality/directionality_check.csv:model_type`). Contrast the
disagree-at-conserved proteins (P53) with the disagree-at-variable ones
(PTEN, HXK4).

Outputs (results/MSA_conservation/), tagged by threshold, e.g. _p80_20 / _p90_10:
  H1_P53_dms_relfrac_by_site_<tag>.csv
  H1_P53_relfrac_summary_<tag>.csv
  H2_persite_spearman_by_class_<tag>.csv
"""
import os
import importlib.util
import pandas as pd
import numpy as np
import common as C
from scipy.stats import spearmanr

OUT = C.rp("MSA_conservation")
COLOC = os.path.join(C.rp("high_disagreement_colocalization"),
                     "colocalization_shared_sites.csv")

_S15 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "15_functional_disagreement_generalization.py")
_spec = importlib.util.spec_from_file_location("func_disagree_gen", _S15)
S15 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S15)

REF = C.load_reference().set_index("DMS_id")
CONS = pd.read_csv(os.path.join(OUT, "msa_conservation_by_site.csv"))
MT = pd.read_csv(C.rp("directionality_check.csv"))
M2T = dict(zip(MT.model, MT.model_type))
SEQMSA = {"Single sequence", "MSA"}
STRUCT = {"Structure"}

P53_GIAC = ["P53_HUMAN_Giacomelli_2018_Null_Etoposide",
            "P53_HUMAN_Giacomelli_2018_Null_Nutlin",
            "P53_HUMAN_Giacomelli_2018_WT_Nutlin"]


def clabel(pct, hi, lo):
    if pd.isna(pct) or pct <= lo:
        return "low"
    if pct >= hi:
        return "high"
    return "mid"


def _offset(dms_id):
    r = REF.loc[dms_id]
    f = S15.features_for(r.UniProt_ID, r.target_seq)
    return (f[0] if f else 0)


def _scores(dms_id):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
    df = df[df.mutant.str.match(r'^[A-Z]\d+[A-Z]$')].copy()
    df["p"] = df.mutant.str.extract(r'^[A-Z](\d+)[A-Z]$')[0].astype(int)
    return df


def dms_relfrac(dms_id):
    off = _offset(dms_id)
    df = _scores(dms_id)
    df["pct"] = df.DMS_score.rank(pct=True) * 100
    g = df.groupby("p")["pct"]
    o = pd.DataFrame({"relfrac": (g.max() - g.min())}).reset_index()
    o["up_pos"] = o["p"] + off
    return o[["up_pos", "relfrac"]]


def persite_class_spearman(dms_id, ups):
    off = _offset(dms_id)
    df = _scores(dms_id)
    models = [c for c in df.columns if c in M2T]
    sm = [m for m in models if M2T[m] in SEQMSA]
    st = [m for m in models if M2T[m] in STRUCT]
    ck = CONS[CONS.DMS_id == dms_id].set_index("up_pos")["conservation_percentile"]
    rows = []
    for up in ups:
        sub = df[df.p == up - off]
        if len(sub) < 8:
            continue

        def cm(mods):
            a = [spearmanr(sub[m], sub.DMS_score)[0] for m in mods if sub[m].nunique() > 1]
            return np.nanmean(a) if a else np.nan
        rows.append({"DMS_id": dms_id, "up_pos": up,
                     "cons_pct": round(ck.get(up, np.nan), 1) if up in ck.index else np.nan,
                     "seqmsa": round(cm(sm), 3), "struct": round(cm(st), 3)})
    return pd.DataFrame(rows)


def run_H1(hi, lo, tag):
    rows, summ = [], []
    for dms in P53_GIAC:
        rf = dms_relfrac(dms)
        cc = CONS[CONS.DMS_id == dms][["up_pos", "conservation_percentile"]].merge(rf, on="up_pos", how="left")
        cc["cl"] = cc.conservation_percentile.apply(lambda p: clabel(p, hi, lo))
        short = dms.replace("P53_HUMAN_Giacomelli_2018_", "")
        summ.append({"assay": short, "threshold": tag,
                     "all_high_n": int((cc.cl == "high").sum()),
                     "all_high_relfrac": round(cc.loc[cc.cl == "high", "relfrac"].mean(), 1),
                     "all_low_n": int((cc.cl == "low").sum()),
                     "all_low_relfrac": round(cc.loc[cc.cl == "low", "relfrac"].mean(), 1)})
        for dtail in ("high", "low"):
            ts = pd.read_csv(os.path.join(OUT, f"p53_{dtail}_disagreement_conservation.tsv"), sep="\t")
            ts = ts[ts.DMS_id == dms].merge(rf, on="up_pos", how="left")
            for _, r in ts.iterrows():
                rows.append({"assay": short, "up_pos": int(r.up_pos), "residue": r.residue,
                             "dis_tail": dtail, "disagreement": round(r.disagreement, 3),
                             "cons_pct": round(r.conservation_percentile, 1) if pd.notna(r.conservation_percentile) else np.nan,
                             "cons_label": clabel(r.conservation_percentile, hi, lo),
                             "dms_relfrac": round(r.relfrac, 1) if pd.notna(r.relfrac) else np.nan})
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, f"H1_P53_dms_relfrac_by_site_{tag}.csv"), index=False)
    pd.DataFrame(summ).to_csv(os.path.join(OUT, f"H1_P53_relfrac_summary_{tag}.csv"), index=False)
    sel = tab.groupby(["dis_tail", "cons_label"])["dms_relfrac"].agg(["mean", "count"]).round(1)
    return tab, pd.DataFrame(summ), sel


def run_H2(hi, lo, tag):
    coloc = pd.read_csv(COLOC)
    rows = []
    for uni, focus in [("P53_HUMAN", "conserved-regime"),
                       ("PTEN_HUMAN", "variable-regime"),
                       ("HXK4_HUMAN", "variable-regime")]:
        sub = coloc[coloc.UniProt_ID == uni]
        for _, s in sub.iterrows():
            up = int(s.up_pos)
            for dms in s.assays_flagging.split(";"):
                if dms not in REF.index:
                    continue
                r = persite_class_spearman(dms, [up])
                if len(r):
                    rr = r.iloc[0].to_dict()
                    rr.update({"UniProt_ID": uni, "regime": focus, "residue": s.residue,
                               "cons_label": clabel(rr["cons_pct"], hi, lo)})
                    rows.append(rr)
    h2 = pd.DataFrame(rows)
    h2.to_csv(os.path.join(OUT, f"H2_persite_spearman_by_class_{tag}.csv"), index=False)
    return h2


def main():
    for hi, lo in [(80, 20), (90, 10)]:
        tag = f"p{hi}_{lo}"
        print(f"\n{'='*78}\nTHRESHOLDS  high>={hi}  low<={lo} or NaN   ({tag})\n{'='*78}")
        tab, summ, sel = run_H1(hi, lo, tag)
        print("\n[H1] P53 selected sites: mean DMS relfrac by disagreement tail x conservation label")
        print(sel.to_string())
        print("\n[H1] baseline mean DMS relfrac over ALL match sites, by conservation class")
        print(summ.to_string(index=False))
        h2 = run_H2(hi, lo, tag)
        print("\n[H2] per-site Spearman by class at shared high-disagreement sites, by regime x conservation")
        agg = h2.groupby(["regime", "cons_label"])[["seqmsa", "struct"]].agg(["mean", "count"]).round(3)
        print(agg.to_string())
    print("\nwrote H1_*/H2_* CSVs (tagged _p80_20 and _p90_10) to", OUT)


if __name__ == "__main__":
    main()
