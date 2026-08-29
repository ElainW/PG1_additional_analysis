"""
Per-column MSA conservation for each DMS assay, joined onto the per-site
disagreement tables by up_pos. Conservation is theta-REWEIGHTED using the
ProteinGym sequence weights (data/msa/DMS_msa_weights/*.npy).

Motivation
----------
Within-site model disagreement is per-residue; evolutionary conservation is
per-MSA-column. This script computes column conservation from each assay's own
ProteinGym MSA (data/msa/DMS_msa_files/, referenced by DMS_substitutions.csv:
MSA_filename) and maps every column onto the SAME up_pos coordinate the site
tables use, so conservation can be correlated with disagreement / the
binding / disorder / disease annotations already on the structure maps.

Coordinate alignment (the crux)
-------------------------------
  * The MSA "focus" sequence (record 1, header >UNIPROT/START-END) equals the
    assay target_seq. a2m match columns are UPPERCASE; lowercase are insertion
    (unmodeled) columns.
  * target_seq index k (1-based, counting non-gap focus residues) is exactly the
    standardized `mutant` position in the score files -- verified match=1.000
    across all 53 assays, INCLUDING the three with a non-null raw_mut_offset
    (HSP82 Cote-Hammarlof=290, Q53Z42 McShan=25, RPC1 Li=18): ProteinGym has
    already folded raw_mut_offset into the standardized mutant column, so we do
    NOT re-apply it here.
  * up_pos = k + offset, where offset is script 15's features_for() alignment of
    target_seq into the UniProt sequence (up_pos = dms_pos + offset). Reusing
    that exact offset guarantees the conservation lands on the same up_pos the
    disagreement sites carry -- even for the isoform / off-by-one proteins
    (e.g. KRAS 4A target_seq vs 4B UniProt).
  * Conservation is reported only on match (uppercase) columns; lowercase focus
    positions are insertion columns and get no value.

Sequence weighting (theta reweighting)
--------------------------------------
The .npy weight vectors are indexed against the EVE-FILTERED sequence set, NOT
the raw a2m records. We reproduce ProteinGym's MSA_processing.preprocess_msa
verbatim (proteingym/utils/msa_utils.py) so weight[i] pairs with retained seq[i]:
  1. name-keyed dict (record 1 = focus); '.'->'-' and upper-case everything.
  2. keep only columns where the focus/WT is not a gap.
  3. drop sequences whose gap fraction over those columns > 0.5.
  4. drop sequences with an indeterminate residue (not in ACDEFGHIKLMNPQRSTVWY,
     not '-') at any kept column.
This reproduces every weight-vector length exactly (all 29 MSAs, diff = 0).

Per match column, over the retained sequences with weights w_s:
  occupancy      = sum(w_s over non-gap) / sum(w_s)
  freq_i         = sum(w_s : residue_s == i) / sum(w_s over the 20 std AAs)
  entropy_bits   = -sum_i freq_i log2 freq_i
  conservation   = 1 - entropy_bits / log2(20)     in [0, 1], 1 = fully conserved
  top_aa/_freq   = most frequent residue (weighted) and its frequency
(The effective depth n_eff = sum(w_s) is printed to the run log per assay; it
equals the reference MSA_N_eff and is not written as a column.)

Outputs (results/MSA_conservation/)
-----------------------------------
  msa_conservation_by_site.csv
      one row per (DMS_id, up_pos) match column, for every assay.
  generalization_functional_disagreement_sites_conservation.csv
      the sites table + conservation columns (left join on DMS_id, up_pos).
  conservation_disagreement_spearman_by_assay.csv
      per-assay Spearman(disagreement, conservation) with p-value.
  <name>_{high,low}_disagreement_conservation.tsv
      each per-protein high-/low-disagreement site table (27 / low-tail selection)
      + conservation columns (does not clobber the originals).
  conservation_high_vs_low_by_assay.csv
      per-assay mean conservation (+percentile) at the high- vs low-disagreement
      tails and their difference (delta_high_minus_low).
"""
import os
import glob
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import common as C

MSADIR = os.path.join(C.ANALYSIS, "data", "msa", "DMS_msa_files")
WDIR = os.path.join(C.ANALYSIS, "data", "msa", "DMS_msa_weights")
SITES = C.rp("generalization_functional_disagreement_sites.csv")
HIDIR = C.rp("high_disagreement_sites")   # 27's high-disagreement inputs
LODIR = C.rp("low_disagreement_sites")    # low-disagreement tail inputs
OUTDIR = C.rp("MSA_conservation")          # 28's per-assay outputs

AA = "ACDEFGHIKLMNPQRSTVWY"
LOG2_20 = np.log2(20)
SEQ_FRAC_GAPS = 0.5          # ProteinGym MSA_processing default

# residue -> code: 0 = gap, 1..20 = AA, 21 = indeterminate/other
_LUT = np.full(256, 21, dtype=np.uint8)
_LUT[ord("-")] = 0
for _i, _a in enumerate(AA):
    _LUT[ord(_a)] = _i + 1

# reuse script 15's offset logic verbatim (module name can't start with a digit)
_S15 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "15_functional_disagreement_generalization.py")
_spec = importlib.util.spec_from_file_location("func_disagree_gen", _S15)
S15 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S15)


def read_a2m_ordered(path):
    """EVE gen_alignment semantics: name-keyed dict, duplicate headers concatenate.
    Returns (focus_name, focus_raw, [(name, raw_seq), ...]) in first-seen order."""
    d, order, focus_name, cur = {}, [], None, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                cur = line
                if focus_name is None:
                    focus_name = cur
                if cur not in d:
                    d[cur] = []
                    order.append(cur)
            else:
                d[cur].append(line.strip())
    seqs = [(nm, "".join(d[nm])) for nm in order]
    return focus_name, "".join(d[focus_name]), seqs


def build_retained_matrix(focus_raw, seqs):
    """Reproduce preprocess_msa filtering; return the retained sequences' MATCH
    columns as an (n_ret, n_match) code matrix, plus per-match-column metadata.

    Filtering columns = positions where the '.'->'-' upper-cased focus is not '-'
    (EVE focus columns, threshold_focus_cols_frac_gaps=1.0 keeps them all).
    Match columns (reported) = the subset whose ORIGINAL focus char is uppercase.
    """
    focus_up = focus_raw.replace(".", "-").upper()
    ncol = len(focus_up)
    keep = np.array([c != "-" for c in focus_up])          # EVE focus columns
    keep_idx = np.where(keep)[0]
    n_keep = keep_idx.size
    # among kept columns, which are a2m match state (focus uppercase)?
    match_in_keep = np.array([focus_raw[c].isupper() for c in keep_idx])
    # up_pos for each kept column: k-th kept residue -> k (1-based)
    k_keep = np.arange(1, n_keep + 1)
    match_aa = np.array([focus_raw[c].upper() for c in keep_idx])[match_in_keep]

    rows = []
    for _nm, seq in seqs:
        su = seq.replace(".", "-").upper()
        b = np.frombuffer(su.encode("ascii", "replace"), dtype=np.uint8)
        if b.size != ncol:                                 # aligned rows only
            b = (b[:ncol] if b.size > ncol
                 else np.concatenate([b, np.full(ncol - b.size, 45, np.uint8)]))
        sub = b[keep_idx]
        if (sub == 45).mean() > SEQ_FRAC_GAPS:             # too-gappy fragment
            continue
        codes = _LUT[sub]
        if (codes == 21).any():                            # indeterminate AA
            continue
        rows.append(codes[match_in_keep])
    R = np.vstack(rows) if rows else np.empty((0, match_in_keep.sum()), np.uint8)
    return R, k_keep[match_in_keep], match_aa


def _col_stats(col, w):
    """Weighted conservation for one match column of codes."""
    c = np.bincount(col, weights=w, minlength=22)
    wsum = w.sum()
    aa = c[1:21]
    denom = aa.sum()
    occ = 1.0 - c[0] / wsum if wsum else np.nan
    if denom <= 0:
        return occ, np.nan, np.nan, "", np.nan
    p = aa / denom
    nz = p[p > 0]
    ent = float(-(nz * np.log2(nz)).sum())
    cons = 1.0 - ent / LOG2_20
    top = int(aa.argmax())
    return occ, ent, cons, AA[top], float(aa[top] / denom)


def conservation_for_assay(dms_id, uni, target_seq, msa_path, weight_path):
    offset = 0
    feats = S15.features_for(uni, target_seq)
    if feats is not None:
        offset = feats[0]
    focus_name, focus_raw, seqs = read_a2m_ordered(msa_path)
    if focus_raw.upper() != target_seq:
        print(f"  !! {dms_id}: focus != target_seq "
              f"({len(focus_raw)} vs {len(target_seq)}), skipped")
        return None
    R, k_match, match_aa = build_retained_matrix(focus_raw, seqs)
    n_ret = R.shape[0]

    w = np.load(weight_path)
    if w.shape[0] != n_ret:
        print(f"  !! {dms_id}: weight len {w.shape[0]} != retained {n_ret}, "
              f"falling back to unweighted")
        w = np.ones(n_ret)
    n_eff = float(w.sum())

    rows = []
    for j in range(R.shape[1]):
        occ, ent, cons, top, topf = _col_stats(R[:, j], w)
        rows.append({
            "DMS_id": dms_id, "UniProt_ID": uni,
            "up_pos": int(k_match[j] + offset),
            "wt_target": match_aa[j],
            "occupancy": round(occ, 4),
            "entropy_bits": round(ent, 4) if ent == ent else np.nan,
            "conservation": round(cons, 4) if cons == cons else np.nan,
            "top_aa": top,
            "top_aa_freq": round(topf, 4) if topf == topf else np.nan,
        })
    # n_eff / n_ret retained for the run log only (not emitted as columns)
    return pd.DataFrame(rows), n_ret, n_eff


def main():
    ref = C.load_reference().set_index("DMS_id")
    sites = pd.read_csv(SITES)
    dms_ids = sorted(sites["DMS_id"].unique())

    tables = []
    for d in dms_ids:
        if d not in ref.index:
            print(f"  !! {d}: not in reference, skipped");  continue
        r = ref.loc[d]
        msa_path = os.path.join(MSADIR, r["MSA_filename"])
        wpath = os.path.join(WDIR, str(r.get("weight_file_name", "")))
        if not os.path.exists(msa_path):
            print(f"  !! {d}: MSA {r['MSA_filename']} missing, skipped");  continue
        if not os.path.exists(wpath):
            print(f"  !! {d}: weights {r.get('weight_file_name')} missing -> unweighted")
        res = conservation_for_assay(d, r["UniProt_ID"], r["target_seq"], msa_path,
                                     wpath if os.path.exists(wpath) else None)
        if res is not None:
            t, n_ret, n_eff = res
            tables.append(t)
            print(f"{d:45s} {r['MSA_filename']:44s} "
                  f"Neff={n_eff:9.1f}/{n_ret:7d}  "
                  f"{len(t):4d} cols  mean_cons={t['conservation'].mean():.3f}")

    os.makedirs(OUTDIR, exist_ok=True)
    cons = pd.concat(tables, ignore_index=True)
    # per-assay percentile of each column's conservation (0-100, across all
    # match columns of that assay; higher = more conserved than the rest)
    cons["conservation_percentile"] = (
        cons.groupby("DMS_id")["conservation"].rank(pct=True) * 100).round(2)
    out1 = os.path.join(OUTDIR, "msa_conservation_by_site.csv")
    cons.to_csv(out1, index=False)
    print(f"\nwrote {out1}  ({len(cons)} rows)")

    keep = ["DMS_id", "up_pos", "occupancy", "entropy_bits", "conservation",
            "conservation_percentile", "top_aa", "top_aa_freq", "wt_target"]
    merged = sites.merge(cons[keep], on=["DMS_id", "up_pos"], how="left")
    for c in ("disagreement", "z_disagree"):
        if c in merged.columns:
            merged[c] = merged[c].round(4)
    out2 = os.path.join(OUTDIR,
                        "generalization_functional_disagreement_sites_conservation.csv")
    merged.to_csv(out2, index=False)
    cov = merged["conservation"].notna().mean()
    print(f"wrote {out2}  ({len(merged)} rows, conservation present for {cov:.1%})")

    # per-assay Spearman(disagreement, conservation)
    srows = []
    for d, g in merged.dropna(subset=["conservation"]).groupby("DMS_id"):
        rho, p = spearmanr(g["disagreement"], g["conservation"])
        srows.append({"DMS_id": d, "UniProt_ID": g["UniProt_ID"].iloc[0],
                      "n_sites": len(g),
                      "spearman_rho": round(float(rho), 4),
                      "p_value": float(f"{p:.3g}")})
    spear = pd.DataFrame(srows).sort_values("DMS_id")
    out3 = os.path.join(OUTDIR, "conservation_disagreement_spearman_by_assay.csv")
    spear.to_csv(out3, index=False)
    print(f"wrote {out3}  ({len(spear)} assays)")

    # augment per-protein high- and low-disagreement site tables with conservation
    tails = {}
    for indir, pattern in ((HIDIR, "*_high_disagreement.tsv"),
                           (LODIR, "*_low_disagreement.tsv")):
        tail = pattern.split("_")[1]   # 'high' or 'low'
        frames, n_tsv = [], 0
        for tsv in sorted(glob.glob(os.path.join(indir, pattern))):
            try:
                h = pd.read_csv(tsv, sep="\t")
            except pd.errors.EmptyDataError:
                continue
            if h.empty or "up_pos" not in h.columns:
                continue
            hm = h.merge(cons[keep].drop(columns=["wt_target"]),
                         on=["DMS_id", "up_pos"], how="left")
            base = os.path.basename(tsv).replace(".tsv", "_conservation.tsv")
            hm.to_csv(os.path.join(OUTDIR, base), sep="\t", index=False)
            frames.append(hm)
            n_tsv += 1
        tails[tail] = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        print(f"wrote {n_tsv} augmented {tail}-disagreement *_conservation.tsv files")

    # per-assay comparison of conservation at the two disagreement tails
    def _assay_means(df):
        g = df.dropna(subset=["conservation"]).groupby("DMS_id")
        return g.agg(n=("conservation", "size"),
                     mean_cons=("conservation", "mean"),
                     mean_pct=("conservation_percentile", "mean"))
    mh, ml = _assay_means(tails["high"]), _assay_means(tails["low"])
    cmp = mh.join(ml, lsuffix="_high", rsuffix="_low", how="outer")
    uni = dict(zip(cons["DMS_id"], cons["UniProt_ID"]))
    cmp.insert(0, "UniProt_ID", [uni.get(d, "") for d in cmp.index])
    cmp["delta_high_minus_low"] = (cmp["mean_cons_high"] - cmp["mean_cons_low"]).round(4)
    cmp["delta_pct_high_minus_low"] = (cmp["mean_pct_high"] - cmp["mean_pct_low"]).round(2)
    for c in ("mean_cons_high", "mean_cons_low", "mean_pct_high", "mean_pct_low"):
        cmp[c] = cmp[c].round(4)
    cmp = cmp.reset_index().sort_values("DMS_id")
    out4 = os.path.join(OUTDIR, "conservation_high_vs_low_by_assay.csv")
    cmp.to_csv(out4, index=False)
    print(f"wrote {out4}  ({len(cmp)} assays)")


if __name__ == "__main__":
    main()
