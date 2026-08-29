"""
Compute Spearman correlation between each model's score and DMS_score, for every
assay. Uses all variants in each assay (matching ProteinGym's
performance_DMS_benchmarks.py). Scores are used as-is (already higher=fitter; see
00_check_directionality.py).

Output: results/spearman_long.csv with columns
    DMS_id, taxon, model, model_type, spearman, n_mutants
"""
import os
import numpy as np
import pandas as pd
import common as C


def main():
    models, mtype, _ = C.model_meta()
    tax = C.dms_to_taxon()
    files = C.list_score_files()

    records = []
    for path in files:
        dms_id = C.dms_id_from_path(path)
        df = pd.read_csv(path, usecols=lambda c: c in C.NON_MODEL_COLS or c in set(models))
        present = [m for m in models if m in df.columns]
        sp = C.spearman_with_target(df, "DMS_score", present)
        n = len(df)
        for m in present:
            records.append({
                "DMS_id": dms_id,
                "taxon": tax.get(dms_id, "Unknown"),
                "model": m,
                "model_type": mtype[m],
                "spearman": sp[m],
                "n_mutants": n,
            })
        print(f"{dms_id:45s} n={n:6d} models={len(present)}")

    out = pd.DataFrame(records)
    os.makedirs(C.RESULTS, exist_ok=True)
    out_path = C.rp("spearman_long.csv")
    out.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}  ({len(out)} rows = {out['DMS_id'].nunique()} assays x {out['model'].nunique()} models)")


if __name__ == "__main__":
    main()
