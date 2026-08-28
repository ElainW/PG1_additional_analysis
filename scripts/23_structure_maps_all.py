"""
Per-protein "structure map" of per-site model disagreement for every protein in
results/generalization_functional_disagreement_sites.csv (the multi-assay set).

Generalizes 08_gck_structure_map.py to all proteins: for each protein, one panel
per assay showing the disagreement profile (1 - mean pairwise within-site Spearman,
45 family-best models) along the UniProt sequence, with binding/active-site
residues (blue lines), disordered regions (shaded), disease-variant positions
(orange rug; ClinVar pathogenic/LP for GCK via data/clinvar/, UniProt natural-
variant flag elsewhere), and the top-12 disagreement sites labelled.

Output (one per protein, named by common/gene name):
  figures/site_disagreement/fig_<common>_structure_map.png
  e.g. fig_gck_structure_map.png (HXK4_HUMAN), fig_kras_..., fig_pten_..., etc.
"""
import os
import re
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import common as C

SITES = os.path.join(C.RESULTS, "generalization_functional_disagreement_sites.csv")
CANDID = os.path.join(C.RESULTS, "generalization_candidate_proteins.csv")
OUTDIR = os.path.join(C.FIGURES, "site_disagreement")
UNIPROT_DIR = os.path.join(C.ANALYSIS, "data", "uniprot")
CLINVAR_DIR = os.path.join(C.ANALYSIS, "data", "clinvar")

# clean filename stems where the bare gene symbol is ambiguous/absent
NAME_OVERRIDE = {"SPIKE_SARS2": "spike", "P53_HUMAN": "p53",
                 "A0A2Z5U3Z0_9INFA": "flu_ha", "F7YBW8_MESOW": "f7ybw8"}

# proteins for which disease positions come from a curated ClinVar export
# (pathogenic + likely pathogenic, missense, SNV, expert-panel/practice-guideline;
# see README for the download recipe) instead of the UniProt natural-variant flag.
CLINVAR_TSV = {"HXK4_HUMAN": "GCK.tsv"}
# The GCK ClinVar query is for diabetes (MODY/PNDM, loss-of-function); it cannot
# contain the hyperinsulinism gain-of-function variants (e.g. W99C/HHF3, whose only
# ClinVar entry is the non-missense Trp99Ter). Supplement with those GoF variants
# from UniProt natural-variant annotations, matched by keyword (as in script 08).
CLINVAR_SUPPLEMENT_KW = {"HXK4_HUMAN": ("hhf", "hyperinsulin", "phhi")}
_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
        "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
        "Tyr": "Y", "Val": "V"}
_PCHANGE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")


def clinvar_positions(uni):
    """Residue positions of ClinVar pathogenic/LP missense variants, parsed from
    the canonical p.(3-letter) change in the Name column. Returns None if the
    protein has no ClinVar export configured."""
    fn = CLINVAR_TSV.get(uni)
    if not fn:
        return None
    df = pd.read_csv(os.path.join(CLINVAR_DIR, fn), sep="\t", index_col=False)
    pos = set()
    for nm in df["Name"]:
        m = _PCHANGE.search(str(nm))
        if m:
            pos.add(int(m.group(2)))
    return pos


def uniprot_disease_supplement(uni):
    """Positions of UniProt natural-variant disease annotations matching the
    supplement keywords (the GoF/hyperinsulinism variants a diabetes-only ClinVar
    query misses). Empty if none configured."""
    kws = CLINVAR_SUPPLEMENT_KW.get(uni)
    if not kws:
        return set()
    e = json.load(open(os.path.join(UNIPROT_DIR, uni + ".json")))
    pos = set()
    for f in e.get("features", []):
        if f["type"] == "Natural variant":
            d = (f.get("description", "") or "").lower()
            if any(k in d for k in kws):
                pos.add(f["location"]["start"]["value"])
    return pos


def common_name(uni, protein):
    if uni in NAME_OVERRIDE:
        return NAME_OVERRIDE[uni]
    gene = str(protein).split(" (")[0].strip()
    g = re.sub(r"[^A-Za-z0-9]+", "_", gene).strip("_").lower()
    return g if 2 <= len(g) <= 12 else uni.split("_")[0].lower()


def uniprot_seq(uni):
    p = os.path.join(UNIPROT_DIR, uni + ".json")
    if os.path.exists(p):
        return json.load(open(p)).get("sequence", {}).get("value", "")
    return ""


def contiguous(positions):
    """Yield (start, end) for runs of consecutive ints."""
    positions = sorted(positions)
    if not positions:
        return
    s = prev = positions[0]
    for p in positions[1:]:
        if p == prev + 1:
            prev = p
        else:
            yield s, prev
            s = prev = p
    yield s, prev


def main():
    sites = pd.read_csv(SITES)
    names = dict(zip(pd.read_csv(CANDID)["UniProt_ID"], pd.read_csv(CANDID)["protein"]))
    os.makedirs(OUTDIR, exist_ok=True)

    for uni, g in sites.groupby("UniProt_ID"):
        cname = common_name(uni, names.get(uni, uni))
        seq = uniprot_seq(uni)
        cv = clinvar_positions(uni)                     # None unless a ClinVar export exists
        if cv is not None:
            supp = uniprot_disease_supplement(uni)      # GoF variants ClinVar diabetes query misses
            dz_set = cv | supp
            dz_src = "ClinVar P/LP + UniProt GoF" if supp else "ClinVar P/LP"
        else:
            dz_set = set(g.loc[g["disease"], "up_pos"])
            dz_src = "UniProt"
        assays = sorted(g["DMS_id"].unique())
        fig, axes = plt.subplots(len(assays), 1, figsize=(16, 2.9 * len(assays) + 1.4),
                                 sharex=True, squeeze=False)
        axes = axes[:, 0]
        for ax, dms_id in zip(axes, assays):
            r = g[g["DMS_id"] == dms_id].sort_values("up_pos")
            ax.plot(r["up_pos"], r["disagreement"], color="#444", lw=0.8, zorder=3)
            ax.fill_between(r["up_pos"], r["disagreement"], color="#bbb", alpha=0.4, zorder=2)
            ymax = max(0.85, r["disagreement"].max() + 0.08)
            # disordered regions -> shaded spans
            for s0, e0 in contiguous(r.loc[r["disordered"], "up_pos"].tolist()):
                ax.axvspan(s0 - 0.5, e0 + 0.5, color="#8172b3", alpha=0.12, zorder=0)
            # binding / active-site residues -> blue lines
            for p in r.loc[r["near_binding"], "up_pos"]:
                ax.axvline(p, color="#1f77b4", lw=0.6, alpha=0.4, zorder=1)
            # disease positions -> orange rug at bottom (ClinVar where available)
            dz = r.loc[r["up_pos"].isin(dz_set), "up_pos"]
            if len(dz):
                ax.plot(dz, np.full(len(dz), 0.012 * ymax), "|", color="#e08214",
                        markersize=5, mew=0.8, zorder=4)
            # top-12 disagreement sites
            top = r.sort_values("disagreement", ascending=False).head(12)
            ax.scatter(top["up_pos"], top["disagreement"], color="black", s=22, zorder=5)
            for _, row in top.iterrows():
                pos = int(row["up_pos"])
                wt = seq[pos - 1] if seq and 1 <= pos <= len(seq) else ""
                tags = [t for t, on in [("bind", row["near_binding"]),
                                        ("dz", pos in dz_set),
                                        ("dis", row["disordered"])] if on]
                lab = f"{wt}{pos}" + (f"\n({'/'.join(tags)})" if tags else "")
                ax.annotate(lab, (pos, row["disagreement"]), fontsize=6, ha="center",
                            va="bottom", xytext=(0, 3), textcoords="offset points")
            ax.set_ylim(0, ymax)
            ax.set_ylabel(f"disagreement\n[{dms_id.replace(uni + '_', '')}]", fontsize=8)
        axes[-1].set_xlabel(f"{uni} residue position (UniProt numbering)")
        legend = [
            Patch(color="#1f77b4", label="binding / active-site residue (±2)"),
            Patch(color="#8172b3", alpha=0.3, label="disordered region"),
            plt.Line2D([0], [0], marker="|", color="#e08214", lw=0, label=f"disease variant ({dz_src})"),
            plt.Line2D([0], [0], marker="o", color="black", lw=0, label="top-12 disagreed sites"),
        ]
        axes[0].legend(handles=legend, fontsize=7.5, ncol=4, loc="upper center",
                       bbox_to_anchor=(0.5, 1.30))
        axes[0].set_title(f"{names.get(uni, uni)} [{uni}] — per-site model disagreement vs functional map",
                          fontsize=11, pad=26)
        plt.tight_layout()
        out = os.path.join(OUTDIR, f"fig_{cname}_structure_map.png")
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"{uni:18s} -> {os.path.basename(out)}  ({len(assays)} assays, {g['up_pos'].nunique()} sites)")


if __name__ == "__main__":
    main()
