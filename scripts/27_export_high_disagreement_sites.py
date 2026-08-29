"""
Write down the exact sites that 23_structure_maps_all.py labels as high-disagreement
(the black dots on each structure map), one table per protein.

A site is labelled by 23 iff, within its assay, disagreement > DISAGREE_FLOOR (0.5)
AND disagreement >= the elbow threshold in the top 2.5-5% tail (capped at MAX_LABELS).
This reproduces that selection verbatim by importing 23's helpers, and records each
labelled site with its functional context (binding/disordered/disease) and the
per-assay threshold that admitted it.

Outputs (one per protein, named like the figures) in
figures/site_disagreement/high_disagreement_sites/:
  <common>_high_disagreement.tsv     e.g. gck_high_disagreement.tsv
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

OUTDIR = os.path.join(C.FIGURES, "site_disagreement", "high_disagreement_sites")


def disease_set(uni, seq, g):
    """Same dz_set 23 uses for its rug / bold labels."""
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
            thr, pct = S.elbow_tail_threshold(r["disagreement"].values)
            lab = r[(r["disagreement"] > S.DISAGREE_FLOOR) & (r["disagreement"] >= thr)]
            lab = lab.sort_values("disagreement", ascending=False).head(S.MAX_LABELS)
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
                    "rank_in_assay": rank,
                    "elbow_threshold": round(float(thr), 4),
                    "effective_floor": round(max(float(thr), S.DISAGREE_FLOOR), 4),
                    "tail_percentile": round(float(pct), 2),
                    "near_binding": bool(row["near_binding"]),
                    "disordered": bool(row["disordered"]),
                    "disease": pos in dz_set,     # bold in the figure
                })
        df = pd.DataFrame(rows)
        out = os.path.join(OUTDIR, f"{cname}_high_disagreement.tsv")
        df.to_csv(out, sep="\t", index=False)
        n_assays = g["DMS_id"].nunique()
        print(f"{uni:16s} -> {os.path.basename(out):40s} {len(df):3d} labelled sites "
              f"across {n_assays} assays")


if __name__ == "__main__":
    main()
