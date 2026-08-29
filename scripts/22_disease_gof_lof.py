"""
Systematic GoF vs LoF test (W99C/R377H generalized) across disease/characterized
variants in HXK4_HUMAN, KCNE1_HUMAN, PTEN_HUMAN, RASK_HUMAN, SRC_HUMAN, CP2C9_HUMAN.

Claim under test: models "miss gain-of-function" (predict activating variants as
deleterious) and "under-call kinetic loss-of-function" (predict LoF variants as
tolerated). Direction is read from each UniProt Natural-variant description
(catalytic/pathway/current language), with a RAS-specific rule (KRAS disease
variants are activating even when the text says "impairs GTPase-activating
protein" -- that mechanism keeps RAS-GTP on).

For every classified single-substitution variant that is measured in an assay we
compute its WITHIN-SITE percentile (0 = worst, 1 = best substitution at that
position) for the experiment and for each of the 23 subset models. Model
percentile low = model predicts deleterious.

Outputs:
  results/disease_gof_lof_variants.csv      (per protein x assay x variant)
  results/disease_gof_lof_summary.csv       (GoF vs LoF aggregates)
  figures/fig_disease_gof_lof.png
"""
import os
import re
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import common as C
from importlib import import_module

sub = import_module("04_subset_figs")

PROTEINS = ["HXK4_HUMAN", "KCNE1_HUMAN", "PTEN_HUMAN", "RASK_HUMAN", "SRC_HUMAN", "CP2C9_HUMAN"]
UNIPROT_DIR = os.path.join(C.ANALYSIS, "data", "uniprot", "raw_json")


# Drop variants annotated benign / likely benign / uncertain significance UNLESS
# they also carry an explicit functional-effect characterization (a change in
# stability / activity / catalytic efficiency / binding, or an activating mechanism).
BENIGN_TERMS = ("benign", "uncertain significance")
FUNC_CHANGE_TERMS = ("increas", "decreas", "reduc", "loss of", "abolish", "impair",
                     "diminish", "constitutively activ", "hyperactiv", "destabil",
                     "thermostab", "gain-of-function", "gain of function", "activates",
                     "activation of the", "nucleotide exchange", "transform", "oncogen")


def classify(desc, protein):
    """Direction by CLINICAL disease phenotype (protein-aware), which is the axis
    that matters for flagging pathogenicity. Cancer = GoF for the RAS oncogene but
    LoF for the PTEN tumor-suppressor; GCK is the only bidirectional case
    (HHF3 hyperinsulinism = GoF; MODY2/PNDM1 diabetes = LoF)."""
    d = (desc or "").lower()
    ev = desc or ""                          # full description (not truncated)
    if any(t in d for t in BENIGN_TERMS) and not any(t in d for t in FUNC_CHANGE_TERMS):
        return "benign_no_function", ev      # benign/uncertain with no functional evidence
    if protein == "HXK4_HUMAN":
        if any(k in d for k in ["hhf", "hyperinsulin", "phhi"]):
            return "GoF", ev
        if any(k in d for k in ["mody", "pndm", "diabet"]):
            return "LoF", ev
        return "unclassified", ev
    if protein == "RASK_HUMAN":       # oncogene / RASopathy -> activating
        if any(k in d for k in ["ns3", "cfc", "costello", "carcinoma", "cancer", "astrocytoma",
                                "aml", "oncogen", "gasc", "activat", "transform", "leukemia",
                                "glm", "melanoma", "nucleotide exchange", "gtpase-activating"]):
            return "GoF", ev
        return "unclassified", ev
    if protein == "SRC_HUMAN":        # proto-oncogene -> activating
        if any(k in d for k in ["thc", "thrombocyt", "activat", "oncogen", "cancer", "carcinoma"]):
            return "GoF", ev
        return "unclassified", ev
    if protein == "PTEN_HUMAN":       # tumor suppressor -> LoF
        if "retains" in d and not re.search(r"(loss|reduc|decreas|sever)\w*[^.;]*phosphatase", d):
            return "neutral", ev
        if any(k in d for k in ["cws", "cowden", "glioma", "glioblastoma", "cancer", "carcinoma",
                                "tumor", "endometrial", "ldd", "loss of phosphatase",
                                "loss of protein phosphatase", "reduced phosphatase",
                                "severely reduced"]):
            return "LoF", ev
        return "unclassified", ev
    if protein == "KCNE1_HUMAN":      # channel LoF (LQT/JLNS)
        if any(k in d for k in ["lqt", "jlns", "deafness", "long qt", "jervell",
                                "reduces i(ks)", "impair"]):
            return "LoF", ev
        return "unclassified", ev
    if protein == "CP2C9_HUMAN":      # pharmacogenetic poor-metabolizer -> LoF
        if any(k in d for k in ["decreas", "reduced", "poor metaboli", "loss", "impair", "abolish"]):
            return "LoF", ev
        return "unclassified", ev
    return "unclassified", ev


def uniprot_offset_and_variants(protein, target_seq):
    e = json.load(open(os.path.join(UNIPROT_DIR, protein + ".json")))
    up = e.get("sequence", {}).get("value", "")
    if up == target_seq or len(up) == len(target_seq):
        offset = 0
    elif target_seq in up:
        offset = up.find(target_seq)
    elif up in target_seq:
        offset = -target_seq.find(up)
    else:
        offset = 0
    out = []
    for f in e.get("features", []):
        if f["type"] != "Natural variant":
            continue
        alt = f.get("alternativeSequence", {})
        orig = alt.get("originalSequence", "")
        alts = alt.get("alternativeSequences", []) or []
        pos = f["location"]["start"]["value"]
        if len(orig) != 1 or not alts or f["location"]["end"]["value"] != pos:
            continue                       # single substitutions only
        label, ev = classify(f.get("description", ""), protein)
        for a in alts:
            if len(a) == 1 and a != orig:
                out.append({"up_pos": pos, "orig": orig, "alt": a,
                            "label": label, "evidence": ev})
    return offset, pd.DataFrame(out)


def within_site_pct(series):
    return (series.rank() - 1) / (len(series) - 1) if len(series) > 1 else series * np.nan


def main():
    ref = C.load_reference()
    rows = []
    for protein in PROTEINS:
        g = ref[ref["UniProt_ID"] == protein]
        target = g["target_seq"].iloc[0]
        offset, variants = uniprot_offset_and_variants(protein, target)
        variants = variants[variants["label"].isin(["GoF", "LoF"])]
        if variants.empty:
            print(f"{protein}: no classified GoF/LoF variants"); continue
        for _, arow in g.iterrows():
            dms_id = arow["DMS_id"]; sel = arow["coarse_selection_type"]
            df = pd.read_csv(os.path.join(C.SCORES_DIR, dms_id + ".csv"))
            df = df[df["mutant"].str.match(r"^[A-Z]\d+[A-Z]$")].copy()
            df["pos"] = df["mutant"].str.extract(r"(\d+)").astype(int)
            models = [c for c in sub.COLS if c in df.columns]
            for _, v in variants.iterrows():
                dms_pos = v["up_pos"] - offset
                mut = f"{v['orig']}{dms_pos}{v['alt']}"
                site = df[df["pos"] == dms_pos]
                if mut not in set(site["mutant"]) or len(site) < 5:
                    continue
                dms_p = within_site_pct(site["DMS_score"])
                idx = site["mutant"] == mut
                mps = [float(within_site_pct(site[m])[idx].iloc[0]) for m in models]
                rows.append({
                    "protein": protein, "variant": mut, "label": v["label"],
                    "assay": dms_id.replace(protein + "_", ""), "selection": sel,
                    "dms_pct": float(dms_p[idx].iloc[0]),
                    "dms_bin": int(site.loc[idx, "DMS_score_bin"].iloc[0]) if "DMS_score_bin" in site else -1,
                    "model_median_pct": float(np.median(mps)),
                    "model_frac_bottom_half": float(np.mean(np.array(mps) < 0.5)),
                    "evidence": v["evidence"],
                })
    res = pd.DataFrame(rows)
    res.to_csv(C.rp("disease_gof_lof_variants.csv"), index=False)

    # unique variants (dedup across assays for the headline counts)
    uniq = res.drop_duplicates(["protein", "variant", "label"])
    print(f"classified measured variants: {len(uniq)} unique "
          f"({(uniq.label=='GoF').sum()} GoF, {(uniq.label=='LoF').sum()} LoF) "
          f"across {res['protein'].nunique()} proteins")

    summ = (res.groupby("label")
            .agg(n=("variant", "size"),
                 model_median_pct=("model_median_pct", "median"),
                 model_frac_bottom_half=("model_frac_bottom_half", "mean"),
                 dms_median_pct=("dms_pct", "median")).reset_index())
    summ.to_csv(C.rp("disease_gof_lof_summary.csv"), index=False)
    print("\n=== GoF vs LoF: how models rank them within-site (low = predicted deleterious) ===")
    print(summ.round(3).to_string(index=False))
    gof = res[res.label == "GoF"]["model_median_pct"]; lof = res[res.label == "LoF"]["model_median_pct"]
    if len(gof) >= 3 and len(lof) >= 3:
        p = mannwhitneyu(gof, lof, alternative="less").pvalue
        print(f"\nMann-Whitney model_median_pct(GoF) < (LoF): p={p:.4f} "
              f"(GoF median {gof.median():.2f} vs LoF {lof.median():.2f})")

    # focus: functional (activity) assays only, where GoF should read as fit
    func = res[res.selection.isin(["Activity", "OrganismalFitness"])]
    print("\n=== on functional (Activity/OrganismalFitness) assays only ===")
    if not func.empty:
        print(func.groupby("label").agg(n=("variant", "size"),
              model_median_pct=("model_median_pct", "median"),
              dms_median_pct=("dms_pct", "median")).round(3).to_string())

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    order = ["GoF", "LoF"]
    cols = {"GoF": "#d1495b", "LoF": "#2e6f95"}
    # A: model vs DMS percentile, colored by label
    ax = axes[0]
    for lab in order:
        d = res[res.label == lab]
        ax.scatter(d["dms_pct"], d["model_median_pct"], c=cols[lab], s=45, alpha=0.8,
                   label=f"{lab} (n={len(d)})", edgecolor="none")
    ax.plot([0, 1], [0, 1], color="gray", ls=":", lw=1)
    ax.axhline(0.5, color="gray", lw=0.5); ax.axvline(0.5, color="gray", lw=0.5)
    ax.set_xlabel("experiment within-site percentile"); ax.set_ylabel("model median within-site percentile")
    ax.set_title("Disease variants: model vs experiment\n(below diagonal = models rank it worse than the assay)")
    ax.legend(fontsize=9)
    # B: distribution of model percentile by label
    ax = axes[1]
    data = [res[res.label == lab]["model_median_pct"].values for lab in order]
    parts = ax.violinplot(data, showmedians=True, showextrema=False)
    for pc, lab in zip(parts["bodies"], order):
        pc.set_facecolor(cols[lab]); pc.set_alpha(0.6)
    for i, lab in enumerate(order, 1):
        vals = res[res.label == lab]["model_median_pct"]
        ax.scatter(np.full(len(vals), i) + np.random.uniform(-.05, .05, len(vals)), vals,
                   c=cols[lab], s=18, alpha=0.7, edgecolor="none", zorder=3)
    ax.axhline(0.5, color="gray", ls="--", lw=1)
    ax.set_xticks([1, 2]); ax.set_xticklabels([f"GoF (n={len(data[0])})", f"LoF (n={len(data[1])})"])
    ax.set_ylabel("model median within-site percentile")
    ax.set_title("Do models rank GoF lower than LoF?\n(low = predicted deleterious)")
    fig.suptitle("Models vs GoF/LoF disease variants — HXK4, KCNE1, PTEN, RASK, SRC, CP2C9", fontsize=13)
    plt.tight_layout()
    out = os.path.join(C.FIGURES, "fig_disease_gof_lof.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")
    print("\nper-variant table (dedup, sorted by label, model_median_pct):")
    with pd.option_context("display.max_colwidth", 55, "display.width", 200):
        print(uniq.sort_values(["label", "model_median_pct"])
              [["protein", "variant", "label", "assay", "dms_pct", "model_median_pct", "evidence"]]
              .to_string(index=False))


if __name__ == "__main__":
    main()
