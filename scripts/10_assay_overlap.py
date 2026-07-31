"""
Which proteins have >=2 DMS assays that overlap in mutated residue positions?

For each assay we take the set of mutated residue positions from the `mutant`
column (splitting multi-mutants on ':'), grouped by protein (UniProt_ID). For
every protein with >=2 assays we test all assay pairs for a shared position; a
protein "overlaps" if at least one pair intersects. We also flag whether the
assays share the same target_seq (same coordinate frame) so overlaps are
meaningful.

Outputs:
  results/assay_position_overlap_proteins.csv  (one row per overlapping protein)
  results/assay_position_overlap_pairs.csv     (one row per overlapping assay pair)
"""
import os
import re
import itertools
import numpy as np
import pandas as pd
import common as C

POS_RE = re.compile(r"[A-Z](\d+)[A-Z]")


def positions_for(dms_id):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"), usecols=["mutant"])
    pos = set()
    for mut in df["mutant"].astype(str):
        for tok in mut.split(":"):
            m = POS_RE.fullmatch(tok)
            if m:
                pos.add(int(m.group(1)))
    return pos


def main():
    ref = C.load_reference()
    idx_of = dict(zip(ref["DMS_id"], ref["DMS_index"]))
    seq_of = dict(zip(ref["DMS_id"], ref["target_seq"]))
    meta = ref.set_index("DMS_id")[["coarse_selection_type", "MSA_Neff_L", "MSA_perc_cov", "seq_len"]]

    multi = ref.groupby("UniProt_ID")["DMS_id"].apply(list)
    multi = multi[multi.apply(len) >= 2]

    pos_cache = {}
    prot_rows, pair_rows, assay_rows = [], [], []
    for uni, dms_ids in multi.items():
        for d in dms_ids:
            if d not in pos_cache:
                pos_cache[d] = positions_for(d)
        same_seq_pairs, dropped_pairs = [], []
        for a, b in itertools.combinations(dms_ids, 2):
            inter = pos_cache[a] & pos_cache[b]
            if inter:
                same_seq = seq_of[a] == seq_of[b]
                (same_seq_pairs if same_seq else dropped_pairs).append((a, b, len(inter)))
                pair_rows.append({
                    "UniProt_ID": uni,
                    "DMS_id_A": a, "DMS_index_A": idx_of[a],
                    "DMS_id_B": b, "DMS_index_B": idx_of[b],
                    "n_shared_positions": len(inter),
                    "same_target_seq": same_seq,
                    "kept": same_seq,
                    "n_pos_A": len(pos_cache[a]), "n_pos_B": len(pos_cache[b]),
                    "shared_positions": ";".join(map(str, sorted(inter))),
                })
        # selected = assays connected by an overlap on the SAME target_seq
        if same_seq_pairs:
            involved = sorted(set(itertools.chain.from_iterable((a, b) for a, b, _ in same_seq_pairs)),
                              key=lambda d: idx_of[d])
            # DOUBLE CHECK: every selected assay must share one identical target_seq
            uniq_seqs = {seq_of[d] for d in involved}
            assert len(uniq_seqs) == 1, f"{uni}: selected assays span {len(uniq_seqs)} target_seqs!"
            dropped = sorted({d for d in dms_ids if d not in involved}, key=lambda d: idx_of[d])
            prot_rows.append({
                "UniProt_ID": uni,
                "n_assays_total": len(dms_ids),
                "n_assays_selected": len(involved),
                "n_overlapping_pairs_kept": len(same_seq_pairs),
                "selected_target_seq_unique": len(uniq_seqs) == 1,
                "DMS_indices_selected": "; ".join(str(idx_of[d]) for d in involved),
                "DMS_ids_selected": "; ".join(involved),
                "DMS_ids_dropped_diff_seq": "; ".join(dropped) if dropped else "",
            })
            for d in involved:
                m = meta.loc[d]
                assay_rows.append({
                    "UniProt_ID": uni, "DMS_index": idx_of[d], "DMS_id": d,
                    "coarse_selection_type": m["coarse_selection_type"],
                    "MSA_Neff_L": m["MSA_Neff_L"], "MSA_perc_cov": m["MSA_perc_cov"],
                    "seq_len": m["seq_len"],
                })

    prot = pd.DataFrame(prot_rows).sort_values("UniProt_ID")
    pairs = pd.DataFrame(pair_rows).sort_values(["UniProt_ID", "DMS_index_A", "DMS_index_B"])
    assays = pd.DataFrame(assay_rows).sort_values(["UniProt_ID", "DMS_index"])
    prot.to_csv(os.path.join(C.RESULTS, "assay_position_overlap_proteins.csv"), index=False)
    pairs.to_csv(os.path.join(C.RESULTS, "assay_position_overlap_pairs.csv"), index=False)
    assays.to_csv(os.path.join(C.RESULTS, "assay_position_overlap_selected_assays.csv"), index=False)

    print(f"Proteins with >=2 assays: {len(multi)}")
    print(f"Proteins with >=1 same-target_seq overlapping pair (selected): {len(prot)}")
    print(f"All proteins' selected assays share one identical target_seq: {bool(prot['selected_target_seq_unique'].all())}")
    dropped_any = prot[prot['DMS_ids_dropped_diff_seq'] != '']
    print(f"Proteins where an assay was dropped for different target_seq: "
          f"{dropped_any['UniProt_ID'].tolist()} -> dropped {dropped_any['DMS_ids_dropped_diff_seq'].tolist()}")
    print()
    with pd.option_context("display.max_colwidth", 120, "display.width", 240):
        print(assays.to_string(index=False))
    print(f"\nSaved proteins, pairs, and selected_assays CSVs to results/")


if __name__ == "__main__":
    main()
