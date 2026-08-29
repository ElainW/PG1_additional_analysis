"""
Verify how the score columns in O2/20250613/*.csv relate to config directionality.

For each model we compute the mean Spearman against DMS_score across all assays.
If the files were stored as RAW model output, directionality=-1 models would
correlate NEGATIVELY with DMS_score (higher score = lower fitness). If the files
are already oriented "higher = fitter", every reasonable model -- including the
directionality=-1 ones -- correlates POSITIVELY.

Output: results/directionality_check.csv and a printed verdict.
"""
import os
import numpy as np
import pandas as pd
import common as C


def main():
    models, mtype, direction = C.model_meta()
    files = C.list_score_files()
    acc = {m: [] for m in models}
    for path in files:
        df = pd.read_csv(path, usecols=lambda c: c in C.NON_MODEL_COLS or c in set(models))
        sp = C.spearman_with_target(df, "DMS_score", [m for m in models if m in df.columns])
        for m, r in sp.items():
            if not np.isnan(r):
                acc[m].append(r)

    rows = []
    for m in models:
        rows.append({
            "model": m,
            "model_type": mtype[m],
            "config_directionality": direction[m],
            "mean_spearman_vs_DMS": np.mean(acc[m]) if acc[m] else np.nan,
            "n_assays": len(acc[m]),
        })
    res = pd.DataFrame(rows)
    os.makedirs(C.RESULTS, exist_ok=True)
    res.to_csv(C.rp("directionality_check.csv"), index=False)

    neg = res[res["config_directionality"] == -1]
    print("directionality=-1 models (raw config), mean Spearman vs DMS_score:")
    print(neg[["model", "mean_spearman_vs_DMS"]].to_string(index=False))
    frac_pos = (neg["mean_spearman_vs_DMS"] > 0).mean()
    print(f"\nFraction of directionality=-1 models with POSITIVE mean Spearman: {frac_pos:.0%}")
    if frac_pos == 1.0:
        print("VERDICT: columns are ALREADY oriented higher=fitter -> do NOT re-apply directionality.")
    else:
        print("VERDICT: some directionality=-1 models are negative -> columns may be RAW; re-apply directionality.")
    print(f"\nSaved {C.rp('directionality_check.csv')}")


if __name__ == "__main__":
    main()
