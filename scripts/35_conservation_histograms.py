"""
H1 follow-up (P53/PTEN/GCK): per-assay conservation across four residue categories.

Builds the underlying H1 tables. The per-assay density FIGURES that used to live here
are now drawn by 37_conservation_dms_panels.py (merged conservation + DMS panels);
this script only writes the tables.

"Relative fraction of conservation" for a residue = `conservation_percentile` / 100
(its conservation rank among the assay's match columns, in [0,1]), with NaN (gappy /
insertion column, no match-state value) mapped to 0.

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
  results/MSA_conservation/H1/H1_conservation_by_category.csv
      one row per (assay, shared high/low-disagreement residue): up_pos, residue,
      disagreement, relative_fraction_of_conservation, cons_category.
  results/MSA_conservation/H1/H1_highcons_middisagreement_range.csv
      per assay: n and min/max relative-fraction for highly-conserved (cons high)
      residues with MID disagreement.
"""
import os
import glob
import importlib.util
import numpy as np
import pandas as pd
import common as C

PROTEINS = ["P53_HUMAN", "PTEN_HUMAN", "HXK4_HUMAN"]
GENE = {"P53_HUMAN": "p53", "PTEN_HUMAN": "pten", "HXK4_HUMAN": "gck"}  # tsv basename
HI_PCT, LO_PCT = 80, 20

OUT = C.rp("MSA_conservation")
H1DIR = os.path.join(OUT, "H1")     # H1 tables
PPDIR = os.path.join(OUT, "per_protein_conservation")  # 28's per-protein tails
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


def assay_positions(dms_id):
    """All assayed up_pos for a DMS (from single-substitution mutants)."""
    off = (S15.features_for(REF.loc[dms_id].UniProt_ID, REF.loc[dms_id].target_seq) or [0])[0]
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"),
                     usecols=["mutant"])
    df = df[df.mutant.str.match(r'^[A-Z]\d+[A-Z]$')]
    pos = df.mutant.str.extract(r'^([A-Z])(\d+)[A-Z]$')
    wt = dict(zip(pos[1].astype(int) + off, pos[0]))          # up_pos -> wt aa
    return wt


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


def main():
    tables, ranges = [], []
    for uni in PROTEINS:
        gene = GENE[uni]
        assays, per_assay, sh_hi, sh_lo = build_protein(uni, gene)
        print(f"{uni}: {len(assays)} assays; shared_high={len(sh_hi)} shared_low={len(sh_lo)}")
        tables.append(make_table(assays, per_assay, sh_hi, sh_lo))
        ranges.append(highcons_middis_range(assays, per_assay))
    os.makedirs(H1DIR, exist_ok=True)
    tab = pd.concat(tables, ignore_index=True)
    tab.to_csv(os.path.join(H1DIR, "H1_conservation_by_category.csv"), index=False)
    rng = pd.concat(ranges, ignore_index=True)
    rng.to_csv(os.path.join(H1DIR, "H1_highcons_middisagreement_range.csv"), index=False)
    print("\nwrote H1/H1_conservation_by_category.csv (%d rows)" % len(tab))
    print("wrote H1/H1_highcons_middisagreement_range.csv")
    print("\nhighly-conserved + mid-disagreement relfrac range per assay:")
    print(rng.to_string(index=False))


if __name__ == "__main__":
    main()
