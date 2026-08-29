"""
Per-protein tables of LOW-disagreement sites -- the mirror image of
27_export_high_disagreement_sites.py.

Where 27 keeps the extreme HIGH tail (disagreement > 0.5 AND in the top 2.5-5%
elbow), this keeps the extreme LOW tail: within each assay a site is recorded iff

  * disagreement < DISAGREE_CEIL (0.5), AND
  * disagreement <= the low-tail elbow threshold: an elbow for the RISE in
    disagreement within the 95-97.5 percentile band (100 = lowest disagreement);
    with no pronounced rise, default to the 95-100 percentile band (lowest 5%).

Selection reuses 23's helpers verbatim (elbow_low_tail_threshold, DISAGREE_CEIL),
so the two exports are exact reverses of one another. Each site is recorded with
the same functional context (binding/disordered/disease) as 27.

Outputs (one per protein) in results/low_disagreement_sites/:
  <common>_low_disagreement.tsv     e.g. gck_low_disagreement.tsv
"""
import os
import importlib.util
import pandas as pd
import common as C

# import 23's helpers (module name can't start with a digit)
_S23 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "23_structure_maps_all.py")
_spec = importlib.util.spec_from_file_location("structure_maps_all", _S23)
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)

OUTDIR = C.rp("low_disagreement_sites")


def disease_set(uni, seq, g):
    """Same dz_set 23/27 use for disease annotation."""
    cv = S.clinvar_positions(uni, seq)
    if cv is not None:
        return cv | S.uniprot_disease_supplement(uni)
    return set(g.loc[g["disease"], "up_pos"].astype(int))


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    sites = pd.read_csv(S.SITES)
    names = dict(zip(pd.read_csv(S.CANDID)["UniProt_ID"], pd.read_csv(S.CANDID)["protein"]))

    for uni, g in sites.groupby("UniProt_ID"):
        cname = S.common_name(uni, names.get(uni, uni))
        seq = S.uniprot_seq(uni)
        dz_set = disease_set(uni, seq, g)
        rows = []
        for dms_id in sorted(g["DMS_id"].unique()):
            r = g[g["DMS_id"] == dms_id].sort_values("up_pos")
            thr, pct = S.elbow_low_tail_threshold(r["disagreement"].values)
            lab = r[(r["disagreement"] < S.DISAGREE_CEIL) & (r["disagreement"] <= thr)]
            lab = lab.sort_values("disagreement", ascending=True).head(S.MAX_LABELS)
            for rank, (_, row) in enumerate(lab.iterrows(), start=1):
                pos = int(row["up_pos"])
                wt = seq[pos - 1] if seq and 1 <= pos <= len(seq) else ""
                rows.append({
                    "UniProt_ID": uni,
                    "DMS_id": dms_id,
                    "up_pos": pos,
                    "wt_aa": wt,
                    "residue": f"{wt}{pos}",
                    "disagreement": round(float(row["disagreement"]), 4),
                    "rank_in_assay": rank,                 # 1 = lowest disagreement
                    "elbow_threshold": round(float(thr), 4),
                    "effective_ceiling": round(min(float(thr), S.DISAGREE_CEIL), 4),
                    "tail_percentile": round(float(pct), 2),
                    "near_binding": bool(row["near_binding"]),
                    "disordered": bool(row["disordered"]),
                    "disease": pos in dz_set,
                })
        df = pd.DataFrame(rows)
        out = os.path.join(OUTDIR, f"{cname}_low_disagreement.tsv")
        df.to_csv(out, sep="\t", index=False)
        n_assays = g["DMS_id"].nunique()
        print(f"{uni:16s} -> {os.path.basename(out):40s} {len(df):3d} low-disagreement sites "
              f"across {n_assays} assays")


if __name__ == "__main__":
    main()
