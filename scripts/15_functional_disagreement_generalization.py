"""
Generalization test: is within-site model disagreement elevated near FUNCTIONAL
sites (binding/active-site residues +-2, and disease natural-variant positions)
and in INTERNAL DISORDERED regions (loops >10 aa from both termini), across ALL
UniProt IDs that have >1 assay?

Per protein we:
  * map DMS positions to UniProt coordinates (offset from target_seq alignment),
  * pull binding/active/site residues, disease variant positions, and disordered
    regions from that protein's own UniProt entry (analysis/data/uniprot/*.json),
  * compute per-site disagreement (1 - mean pairwise within-site Spearman over the
    45 family-best models; sites with >=5 substitutions) for each of its assays,
  * z-score disagreement within each (protein, assay) so assays are comparable,
  * test functional vs rest and disordered-internal vs rest (Mann-Whitney).

Outputs:
  results/generalization_functional_disagreement_sites.csv
  results/generalization_functional_disagreement_per_protein.csv
  figures/fig_generalization_functional_disagreement.png
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, zscore
import common as C
from importlib import import_module

fam = import_module("05_fig6_familybest")
persite = import_module("03_persite_ranking")

MIN_SUBS = 5
BIND_WINDOW = 2
TERM_PAD = 10  # ignore residues within 10 aa of either terminus for disorder
DIS_KW = ["mody", "hyperinsulin", "diabetes", "cancer", "carcinoma", "syndrome",
          "deafness", "disease", "phhi", "hhf", "encephalopath", "arrhythmi",
          "lqt", "myopathy", "anemia", "epilep", "lethal", "pathogenic",
          "susceptibility", "hypertension"]
UNIPROT_DIR = os.path.join(C.ANALYSIS, "data", "uniprot", "raw_json")


def features_for(uni, target_seq):
    """Return (offset, binding_pos, disease_pos, disordered_internal_pos, seq_len).
    offset maps DMS position -> UniProt position: up_pos = dms_pos + offset."""
    path = os.path.join(UNIPROT_DIR, uni + ".json")
    if not os.path.exists(path):
        return None
    e = json.load(open(path))
    up = e.get("sequence", {}).get("value", "")
    if not up:
        return None
    if up == target_seq or len(up) == len(target_seq):
        offset = 0
    elif target_seq in up:
        offset = up.find(target_seq)
    elif up in target_seq:
        offset = -target_seq.find(up)
    else:
        k = target_seq[:25]
        i = up.find(k)
        offset = i if i >= 0 else 0
    L = len(up)
    binding, disease, disordered = set(), set(), set()
    for f in e.get("features", []):
        t = f["type"]
        s = f["location"]["start"]["value"]; en = f["location"]["end"]["value"]
        if s is None or en is None:
            continue
        if t in ("Binding site", "Active site", "Site", "Metal binding"):
            binding.update(range(s, en + 1))
        elif t == "Natural variant":
            desc = (f.get("description", "") or "").lower()
            if any(kw in desc for kw in DIS_KW):
                disease.update(range(s, en + 1))
        elif t == "Region" and "disorder" in (f.get("description", "") or "").lower():
            for p in range(s, en + 1):
                if TERM_PAD < p <= L - TERM_PAD:
                    disordered.add(p)
    binding_prox = set()
    for b in binding:
        binding_prox.update(range(b - BIND_WINDOW, b + BIND_WINDOW + 1))
    return offset, binding_prox, disease, disordered, L


def site_disagreement(dms_id, chosen):
    df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
    df = df[df["mutant"].str.match(r"^[A-Z]\d+[A-Z]$")].copy()
    if df.empty:
        return pd.DataFrame()
    df["pos"] = df["mutant"].str.extract(r"(\d+)").astype(int)
    present = [m for m in chosen if m in df.columns]
    rows = []
    for pos, g in df.groupby("pos"):
        if len(g) < MIN_SUBS:
            continue
        R = g[present].rank().values
        cc = np.corrcoef(R.T)
        off = cc[~np.eye(cc.shape[0], dtype=bool)]
        rows.append({"pos": int(pos), "disagreement": 1 - float(np.nanmean(off))})
    return pd.DataFrame(rows)


def main():
    chosen, _, _ = fam.build_family_best()
    ref = C.load_reference()
    # Drop coordinate-frame-incompatible assays (e.g. P53 Kotler on a different
    # target_seq) before enumerating per-protein assay groups, so they never enter
    # the master site table or the shared-coordinate 'seqmap' offset. See
    # common.EXCLUDE_DMS.
    ref = C.drop_excluded_dms(ref)
    seqmap = ref.groupby("UniProt_ID")["target_seq"].first()
    multi = ref.groupby("UniProt_ID")["DMS_id"].apply(list)
    multi = multi[multi.apply(len) >= 2]

    site_rows = []
    for uni, dms_ids in multi.items():
        feats = features_for(uni, seqmap[uni])
        if feats is None:
            print(f"{uni}: no UniProt sequence, skipped")
            continue
        offset, binding_prox, disease, disordered, L = feats
        for dms_id in dms_ids:
            sd = site_disagreement(dms_id, chosen)
            if sd.empty:
                continue
            sd["UniProt_ID"] = uni
            sd["DMS_id"] = dms_id
            sd["up_pos"] = sd["pos"] + offset
            sd["near_binding"] = sd["up_pos"].isin(binding_prox)
            sd["disease"] = sd["up_pos"].isin(disease)
            sd["functional"] = sd["near_binding"] | sd["disease"]
            sd["disordered"] = sd["up_pos"].isin(disordered)
            sd["z_disagree"] = zscore(sd["disagreement"], nan_policy="omit")
            site_rows.append(sd)
        print(f"{uni:16s} binding+-2:{len(binding_prox):4d} disease:{len(disease):4d} disorder_int:{len(disordered):4d}")

    sites = pd.concat(site_rows, ignore_index=True)
    sites.to_csv(C.rp("generalization_functional_disagreement_sites.csv"), index=False)

    def per_protein(flag):
        rows = []
        for uni, g in sites.groupby("UniProt_ID"):
            a = g[g[flag]]["z_disagree"].dropna()
            b = g[~g[flag]]["z_disagree"].dropna()
            if len(a) >= 5 and len(b) >= 5:
                p = mannwhitneyu(a, b, alternative="greater").pvalue  # H1: functional MORE disagreement
                rows.append({"UniProt_ID": uni, "flag": flag, "n_flag": len(a), "n_rest": len(b),
                             "mean_z_flag": a.mean(), "mean_z_rest": b.mean(),
                             "delta_z": a.mean() - b.mean(), "p_greater": p})
        return pd.DataFrame(rows)

    func_tab = per_protein("functional")
    dis_tab = per_protein("disordered")
    per_prot = pd.concat([func_tab, dis_tab], ignore_index=True)
    per_prot.to_csv(C.rp("generalization_functional_disagreement_per_protein.csv"), index=False)

    def summarize(tab, label):
        n = len(tab)
        pos = (tab["delta_z"] > 0).sum()
        sig = (tab["p_greater"] < 0.05).sum()
        print(f"\n{label}: testable in {n} proteins")
        print(f"  higher disagreement in group (delta_z>0): {pos}/{n}")
        print(f"  significantly higher (Mann-Whitney p<0.05): {sig}/{n}")
        # pooled test
        a = sites[sites[tab['flag'].iloc[0]]]["z_disagree"].dropna() if n else []
        b = sites[~sites[tab['flag'].iloc[0]]]["z_disagree"].dropna() if n else []
        if n:
            pp = mannwhitneyu(a, b, alternative="greater").pvalue
            print(f"  POOLED across all proteins: mean z {a.mean():+.3f} (group) vs {b.mean():+.3f} (rest), p={pp:.2e}")
        print(tab.sort_values("delta_z", ascending=False)
              [["UniProt_ID", "n_flag", "n_rest", "delta_z", "p_greater"]].round(3).to_string(index=False))
        return tab

    print("\n" + "=" * 70)
    summarize(func_tab, "FUNCTIONAL sites (binding+-2 OR disease variant)")
    summarize(dis_tab, "INTERNAL DISORDERED regions (>10 aa from termini)")

    # ---------------- figure: per-protein delta_z ----------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, tab, title in [(axes[0], func_tab, "functional sites\n(binding ±2 or disease)"),
                           (axes[1], dis_tab, "internal disordered regions")]:
        t = tab.sort_values("delta_z")
        colors = ["#d62728" if p < 0.05 else "#9ecae1" for p in t["p_greater"]]
        ax.barh(t["UniProt_ID"], t["delta_z"], color=colors)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel("Δ mean z-disagreement (group − rest)")
        ax.set_title(f"{title}\n(red = Mann-Whitney p<0.05, higher in group)", fontsize=10)
        ax.tick_params(axis="y", labelsize=8)
    fig.suptitle("Is model disagreement elevated at functional / disordered sites? (multi-assay proteins)",
                 fontsize=12)
    plt.tight_layout()
    out = os.path.join(C.FIGURES, "fig_generalization_functional_disagreement.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
