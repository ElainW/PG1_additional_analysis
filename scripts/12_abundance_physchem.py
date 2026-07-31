"""
HXK4_HUMAN_Gersing_2023_abundance: is model performance explained by residue
hydrophobicity and net-charge differences?

The abundance assay is a stability/folding readout, so simple physicochemical
changes should carry much of the signal. We test two things for the 23-model
subset:

1. Features-only baseline. Per single mutant compute
     d_hydro   = KyteDoolittle(mut) - KyteDoolittle(wt)
     d_charge  = charge(mut) - charge(wt)          (D,E=-1; K,R=+1; H=+0.1)
   plus their magnitudes |d_hydro|, |d_charge| (stability cares about the size
   of the change). Controls X = [d_hydro, |d_hydro|, d_charge, |d_charge|].
   Baseline performance = Spearman(DMS, OLS fit of DMS ~ X).

2. Do models just recapitulate these features? For each model:
     raw       = Spearman(DMS, model)
     partial   = partial Spearman(DMS, model | X)   (residualize ranks on X)
     frac_expl = 1 - partial/raw  (share of performance attributable to X)
     model_R2  = R^2 of OLS(model_score ~ X)        (how physicochemical the model is)
   A large raw->partial drop means the model's performance is largely explained
   by hydrophobicity + charge.

Outputs:
  results/abundance_physchem_model_table.csv
  figures/fig_abundance_physchem.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata
import common as C
from importlib import import_module

sub = import_module("04_subset_figs")

DMS_ID = "HXK4_HUMAN_Gersing_2023_abundance"

KD = {  # Kyte-Doolittle hydropathy
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHG = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}


def charge(aa):
    return CHG.get(aa, 0.0)


def partial_spearman(y, x, Z):
    """Partial Spearman corr of y,x controlling for columns of Z (linear)."""
    ry, rx = rankdata(y), rankdata(x)
    Zc = np.column_stack([np.ones(len(Z))] + [Z[c].values for c in Z.columns])
    resy = ry - Zc @ np.linalg.lstsq(Zc, ry, rcond=None)[0]
    resx = rx - Zc @ np.linalg.lstsq(Zc, rx, rcond=None)[0]
    if resy.std() == 0 or resx.std() == 0:
        return np.nan
    return float(np.corrcoef(resy, resx)[0, 1])


def r2(y, Z):
    Zc = np.column_stack([np.ones(len(Z))] + [Z[c].values for c in Z.columns])
    beta = np.linalg.lstsq(Zc, y, rcond=None)[0]
    resid = y - Zc @ beta
    ss_tot = ((y - y.mean()) ** 2).sum()
    return 1 - (resid ** 2).sum() / ss_tot if ss_tot > 0 else np.nan


def main():
    df = pd.read_csv(os.path.join(C.SCORES_DIR, DMS_ID + ".csv"))
    df = df[df["mutant"].str.match(r"^[A-Z]\d+[A-Z]$")].copy()
    df["wt"] = df["mutant"].str[0]
    df["alt"] = df["mutant"].str[-1]
    df["d_hydro"] = df["alt"].map(KD) - df["wt"].map(KD)
    df["d_charge"] = df["alt"].map(charge) - df["wt"].map(charge)
    df["abs_d_hydro"] = df["d_hydro"].abs()
    df["abs_d_charge"] = df["d_charge"].abs()
    X = df[["d_hydro", "abs_d_hydro", "d_charge", "abs_d_charge"]]
    y = df["DMS_score"].values

    # feature characterization
    print(f"assay: {DMS_ID}   n single mutants: {len(df)}")
    print("Spearman(DMS, feature):")
    for f in ["d_hydro", "abs_d_hydro", "d_charge", "abs_d_charge"]:
        print(f"  {f:14s} {spearmanr(y, df[f]).statistic:+.3f}")

    # features-only baseline
    Zc = np.column_stack([np.ones(len(X))] + [X[c].values for c in X.columns])
    fit = Zc @ np.linalg.lstsq(Zc, y, rcond=None)[0]
    baseline_perf = spearmanr(y, fit).statistic
    baseline_r2 = r2(y, X)
    print(f"\nfeatures-only baseline: Spearman(DMS, fit) = {baseline_perf:.3f}  (R^2 = {baseline_r2:.3f})")

    present = [c for c in sub.COLS if c in df.columns]
    rows = []
    for c in present:
        raw = spearmanr(y, df[c]).statistic
        part = partial_spearman(y, df[c].values, X)
        rows.append({
            "model": sub.COL2DISP[c],
            "model_type": sub.TYPE[sub.COL2DISP[c]],
            "spearman_raw": raw,
            "spearman_partial": part,
            "frac_explained_by_physchem": 1 - part / raw if raw not in (0, np.nan) else np.nan,
            "model_R2_on_physchem": r2(df[c].values, X),
            "model_corr_d_hydro": spearmanr(df[c], df["d_hydro"]).statistic,
            "model_corr_abs_d_charge": spearmanr(df[c], df["abs_d_charge"]).statistic,
        })
    tab = pd.DataFrame(rows).sort_values("spearman_raw", ascending=False).reset_index(drop=True)
    tab.to_csv(os.path.join(C.RESULTS, "abundance_physchem_model_table.csv"), index=False)

    print("\nper-model raw vs partial Spearman (controlling hydrophobicity + charge):")
    print(tab[["model", "spearman_raw", "spearman_partial",
               "frac_explained_by_physchem", "model_R2_on_physchem"]].round(3).to_string(index=False))
    print(f"\nmean raw = {tab['spearman_raw'].mean():.3f}, mean partial = {tab['spearman_partial'].mean():.3f}, "
          f"mean fraction explained = {tab['frac_explained_by_physchem'].mean():.1%}")

    # ---------------- figure ----------------
    palette = dict(sub.PALETTE)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 8), gridspec_kw={"width_ratios": [1.5, 1]})

    yy = np.arange(len(tab))
    colors = tab["model_type"].map(palette)
    axL.barh(yy - 0.2, tab["spearman_raw"], height=0.4, color=colors, label="raw")
    axL.barh(yy + 0.2, tab["spearman_partial"], height=0.4, color=colors, alpha=0.45,
             hatch="//", label="partial | hydro+charge")
    axL.axvline(baseline_perf, color="black", ls="--", lw=1.2)
    axL.text(baseline_perf, -1.1, f"features-only\nbaseline = {baseline_perf:.2f}",
             fontsize=8, ha="center", va="bottom")
    axL.set_yticks(yy)
    axL.set_yticklabels(tab["model"], fontsize=8)
    for t, mt in zip(axL.get_yticklabels(), tab["model_type"]):
        t.set_color(palette[mt])
    axL.invert_yaxis()
    axL.set_xlabel("Spearman vs abundance DMS")
    axL.set_title(f"{DMS_ID}\nmodel performance: raw (solid) vs partial after\n"
                  "controlling for hydrophobicity + charge (hatched)", fontsize=11)
    axL.axvline(0, color="gray", lw=0.6)
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[t]) for t in sub.MTYPES]
    axL.legend([*handles], [*sub.MTYPES], fontsize=7, title="model_type", loc="lower right")

    # right: how physicochemical is each model (R2) vs its performance
    axR.scatter(tab["model_R2_on_physchem"], tab["spearman_raw"],
                c=[palette[t] for t in tab["model_type"]], s=45, edgecolor="none")
    for _, r in tab.iterrows():
        axR.annotate(r["model"], (r["model_R2_on_physchem"], r["spearman_raw"]),
                     fontsize=5.5, xytext=(3, 2), textcoords="offset points")
    axR.axhline(baseline_perf, color="black", ls="--", lw=1, label=f"features baseline ({baseline_perf:.2f})")
    axR.set_xlabel("model score variance explained by\nhydrophobicity + charge  (R²)")
    axR.set_ylabel("model Spearman vs abundance DMS")
    axR.set_title("more physicochemical models\nvs their performance", fontsize=11)
    axR.legend(fontsize=7)
    plt.tight_layout()
    out = os.path.join(C.FIGURES, "fig_abundance_physchem.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
