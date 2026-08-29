"""
HXK4_HUMAN (glucokinase / GCK): per-site model disagreement.

Two assays: HXK4_HUMAN_Gersing_2022_activity (OrganismalFitness) and
HXK4_HUMAN_Gersing_2023_abundance (Expression).

For each residue position (site) with >=5 single substitutions we ask: how much
do models disagree on *which* substitutions are more vs less fit at that site?

Per site:
  * within-site rank each model's score (higher = fitter; data already oriented).
  * mean_pairwise_spearman = mean Spearman over all model pairs of those
    within-site rankings. LOW value = models disagree the most.
  * consensus ranking = mean within-site rank across models; we then correlate
    each model_TYPE's mean ranking with the consensus to see who dissents.
  * dms_consensus_spearman = how the model consensus ranks vs the experiment.
  * dms_within_site_std = spread of DMS_score at the site (context: flat sites
    carry little real signal).

Models: the 45-model family-best subset (one model per family) from
05_fig6_familybest.py, so redundant parameter variants don't inflate agreement.

Outputs (per assay):
  results/hxk4_<assay>_site_disagreement.csv
  figures/fig_hxk4_<assay>_disagreement.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import common as C
from importlib import import_module

sns.set_style("whitegrid")
fam = import_module("05_fig6_familybest")

ASSAYS = ["HXK4_HUMAN_Gersing_2022_activity", "HXK4_HUMAN_Gersing_2023_abundance"]
MIN_SUBS = 5


def analyze(dms_id, chosen, type_of, mtypes):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
    present = [m for m in chosen if m in df.columns]
    members_by_type = {t: [m for m in present if type_of[m] == t] for t in mtypes}
    members_by_type = {t: v for t, v in members_by_type.items() if len(v) >= 1}

    df["pos"] = df["mutant"].map(C.parse_position)
    df["wt"] = df["mutant"].str[0]
    df = df[df["pos"].notna()].copy()

    rows = []
    for pos, g in df.groupby("pos"):
        if len(g) < MIN_SUBS:
            continue
        R = g[present].rank()                       # within-site ranks (variants x models)
        Rv = R.values
        # mean pairwise Spearman across models = mean off-diagonal of corr of ranks
        cc = np.corrcoef(Rv.T)
        m = cc.shape[0]
        off = cc[~np.eye(m, dtype=bool)]
        mean_pair = float(np.nanmean(off))
        consensus = R.mean(axis=1)                  # crowd ranking
        # model_type vs consensus
        type_agree = {}
        for t, members in members_by_type.items():
            tr = R[members].mean(axis=1)
            if tr.std() > 0 and consensus.std() > 0:
                type_agree[t] = float(np.corrcoef(tr.values, consensus.values)[0, 1])
            else:
                type_agree[t] = np.nan
        dms_rank = g["DMS_score"].rank()
        dms_cons = (float(np.corrcoef(consensus.values, dms_rank.values)[0, 1])
                    if dms_rank.std() > 0 else np.nan)
        row = {
            "position": int(pos), "wt_aa": g["wt"].iloc[0], "n_subs": len(g),
            "mean_pairwise_spearman": mean_pair, "disagreement": 1 - mean_pair,
            "dms_consensus_spearman": dms_cons,
            "dms_within_site_std": float(g["DMS_score"].std()),
        }
        for t in mtypes:
            row[f"agree_{t}"] = type_agree.get(t, np.nan)
        rows.append(row)

    res = pd.DataFrame(rows).sort_values("disagreement", ascending=False).reset_index(drop=True)
    res.to_csv(C.rp(f"hxk4_{dms_id}_site_disagreement.csv"), index=False)
    return res, members_by_type


def make_figure(dms_id, res, mtypes, members_by_type, outdir=C.FIGURES, disordered=frozenset()):
    fig = plt.figure(figsize=(15, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.3], wspace=0.32)

    # Panel A: position vs consensus (mean pairwise spearman), annotate top disagreers
    axA = fig.add_subplot(gs[0, 0])
    axA.scatter(res["position"], res["mean_pairwise_spearman"], s=18,
                c=res["disagreement"], cmap="rocket_r", edgecolor="none")
    axA.axhline(res["mean_pairwise_spearman"].median(), color="gray", ls="--", lw=1,
                label=f"median = {res['mean_pairwise_spearman'].median():.2f}")
    top = res.head(12)
    for _, r in top.iterrows():
        lab = f"{r['wt_aa']}{int(r['position'])}"
        if int(r["position"]) in disordered:
            lab += " (disordered)"
        axA.annotate(lab, (r["position"], r["mean_pairwise_spearman"]),
                     fontsize=7, ha="center", va="top", xytext=(0, -4), textcoords="offset points",
                     color="tab:purple" if int(r["position"]) in disordered else "black")
    axA.set_xlabel("residue position"); axA.set_ylabel("mean pairwise within-site Spearman\n(low = models disagree)")
    axA.set_title(f"{dms_id}\nper-site model agreement (n={len(res)} sites, >= {MIN_SUBS} subs)")
    axA.legend(fontsize=8)

    # Panel B: top-15 disagreement sites, model_type vs consensus agreement
    axB = fig.add_subplot(gs[0, 1])
    tcols = [f"agree_{t}" for t in mtypes if any(len(members_by_type.get(t, [])) >= 1 for _ in [0])]
    tcols = [c for c in tcols if c in res.columns]
    H = res.head(15).set_index(res.head(15).apply(lambda r: f"{r['wt_aa']}{int(r['position'])}", axis=1))
    sns.heatmap(H[tcols].astype(float), cmap="vlag", center=0.5, vmin=-0.2, vmax=1,
                annot=True, fmt=".2f", annot_kws={"size": 6}, ax=axB,
                cbar_kws={"label": "agreement with model consensus"})
    axB.set_title("Top-15 disagreed sites:\nwhich model_type dissents?", fontsize=10)
    axB.set_xticklabels([t.replace("agree_", "").replace(" & ", " &\n").replace(", ", ",\n") for t in tcols],
                        fontsize=6, rotation=30, ha="right")
    axB.set_ylabel("site (wt+pos)"); axB.set_xlabel("")
    plt.setp(axB.get_yticklabels(), fontsize=7, rotation=0)

    out = os.path.join(outdir, f"fig_hxk4_{dms_id}_disagreement.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


def main():
    chosen, type_of, sel = fam.build_family_best()
    mtypes = [t for t in C.MODEL_TYPE_ORDER if t in set(type_of.values())]
    for dms_id in ASSAYS:
        res, members_by_type = analyze(dms_id, chosen, type_of, mtypes)
        make_figure(dms_id, res, mtypes, members_by_type)
        print(f"\n==== {dms_id} ====")
        print(f"sites analyzed: {len(res)}   median pairwise Spearman: {res['mean_pairwise_spearman'].median():.3f}")
        show = ["position", "wt_aa", "n_subs", "mean_pairwise_spearman",
                "dms_consensus_spearman", "dms_within_site_std"]
        print("Top 15 sites where models disagree MOST:")
        print(res.head(15)[show].to_string(index=False))
        print("Bottom 5 (models agree most):")
        print(res.tail(5)[show].to_string(index=False))
        print(f"Saved results/hxk4_{dms_id}_site_disagreement.csv")


if __name__ == "__main__":
    main()
