"""Build the enriched selected-assay table (paper Table T1).

Combines the selected-assay list with per-protein annotations and the
per-assay wild-type DMS score:

  base : results/assay_selection/assay_position_overlap_selected_assays.csv
  +conformational_switch, +taxon   <- results/generalization_candidate_proteins.csv (by UniProt_ID)
  +selected_pdb                    <- data/structures/pdb/README.md "Selected structures"
                                      (curated experimental structures; only the 6 proteins
                                       with a curated pick — blank otherwise)
  +DMS_wt_score                    <- data/DMS_wt_score.csv (by DMS_id; blank where absent)

Writes results/assay_selection/selected_assays_annotated.csv.
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def p(*a):
    return os.path.join(ROOT, *a)


# Curated experimental structures, read from data/structures/pdb/README.md
# "Selected structures" (one per protein, or per conformational state for switches).
# Only the 6 proteins with a curated pick appear there; all others are blank.
SELECTED_PDB = {
    "HXK4_HUMAN": "1V4T (super-open); 1V4S (closed); open state unassigned",
    "RASK_HUMAN": "KRAS set: 6VJJ (ON/RAF1), 5O2S (DARPin K27), 1NVW (SOS)",
    "CP2C9_HUMAN": "1R9O",
    "PTEN_HUMAN": "1D5R",
    "SRC_HUMAN": "2SRC (inactive); 1Y57 (active)",
    "VKOR1_HUMAN": "6WVH",
}


def main():
    base = pd.read_csv(p("results/assay_selection/assay_position_overlap_selected_assays.csv"))
    cand = pd.read_csv(p("results/generalization_candidate_proteins.csv"))
    wt = pd.read_csv(p("data/DMS_wt_score.csv"))

    # by UniProt_ID: taxon
    ann = cand[["UniProt_ID", "taxon"]].drop_duplicates("UniProt_ID")
    out = base.merge(ann, on="UniProt_ID", how="left")

    # by UniProt_ID: curated selected PDB(s)
    out["selected_pdb"] = out["UniProt_ID"].map(SELECTED_PDB).fillna("")

    # by DMS_id: wild-type DMS score
    out = out.merge(
        wt[["DMS_id", "wt_score"]].rename(columns={"wt_score": "DMS_wt_score"}),
        on="DMS_id", how="left",
    )

    cols = ["UniProt_ID", "DMS_index", "DMS_id", "coarse_selection_type", "taxon",
            "selected_pdb", "MSA_Neff_L", "MSA_perc_cov",
            "seq_len", "DMS_wt_score"]
    out = out[cols]
    dst = p("results/assay_selection/selected_assays_annotated.csv")
    out.to_csv(dst, index=False)

    n_pdb = (out["selected_pdb"] != "").sum()
    n_wt = out["DMS_wt_score"].notna().sum()
    print(f"wrote {dst}")
    print(f"  {len(out)} assays, {out['UniProt_ID'].nunique()} proteins")
    print(f"  selected_pdb populated: {n_pdb}/{len(out)} rows "
          f"({out.loc[out['selected_pdb'] != '', 'UniProt_ID'].nunique()} proteins)")
    print(f"  DMS_wt_score populated: {n_wt}/{len(out)} rows")
    print(f"  taxon missing: {out['taxon'].isna().sum()}")


if __name__ == "__main__":
    main()
