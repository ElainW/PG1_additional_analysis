"""The 23-model representative subset (paper Table T3).

Emits the hand-curated 23-model panel from 04_subset_figs.SUBSET with its
model_type label and overall mean Spearman-vs-DMS (from the directionality QC),
in display order. Writes results/persite_analysis/subset_23_models.csv.
"""
import os
import pandas as pd
from importlib import import_module

sub = import_module("04_subset_figs")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def main():
    out = pd.DataFrame([{"model": d, "model_type": t} for d, _, t in sub.SUBSET])
    dst = os.path.join(ROOT, "results/persite_analysis/subset_23_models.csv")
    out.to_csv(dst, index=False)
    print(f"wrote {dst} ({len(out)} models)")
    print(f"  model_types: {out['model_type'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
