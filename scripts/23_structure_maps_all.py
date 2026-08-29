"""
Per-protein "structure map" of per-site model disagreement for every protein in
results/generalization_functional_disagreement_sites.csv (the multi-assay set).

Generalizes 08_gck_structure_map.py to all proteins: for each protein, one panel
per assay showing the disagreement profile (1 - mean pairwise within-site Spearman,
45 family-best models) along the UniProt sequence. As in 08 (which does this only
for GCK's ATP/glucose sites and small/large subdomains), the functional map is
resolved per feature rather than lumped together:
  * binding residues coloured *separately per ligand* (one colour + legend entry
    per distinct UniProt ligand, e.g. ATP vs substrate vs heme/FAD/GTP/...),
    plus catalytic "Active site" residues;
  * subdomains / domains shaded *separately* (one light span colour per
    UniProt Region/Domain, e.g. GCK small vs large subdomain);
  * disordered regions shaded and disease-variant positions marked exactly as
    before (purple spans; orange rug — ClinVar pathogenic/LP + a UniProt
    gain-of-function supplement for the 10 human proteins with a data/clinvar/
    export, UniProt natural-variant flag for the rest).

Sites are *labelled* (not a fixed top-12) when they are both strongly disagreed
(disagreement > 0.5) and in the extreme tail: an elbow within the top 2.5-5% of
each assay's disagreement distribution (the largest sudden drop in that band sets
the exact percentile). Labels for disease-variant sites are drawn in bold.

Output (one per protein, named by common/gene name):
  figures/site_disagreement/fig_<common>_structure_map.png
  e.g. fig_gck_structure_map.png (HXK4_HUMAN), fig_kras_..., fig_pten_..., etc.
"""
import os
import re
import json
from collections import defaultdict
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
UNIPROT_DIR = os.path.join(C.ANALYSIS, "data", "uniprot", "raw_json")
CLINVAR_DIR = os.path.join(C.ANALYSIS, "data", "clinvar")

# clean filename stems where the bare gene symbol is ambiguous/absent
NAME_OVERRIDE = {"SPIKE_SARS2": "spike", "P53_HUMAN": "p53",
                 "A0A2Z5U3Z0_9INFA": "flu_ha", "F7YBW8_MESOW": "f7ybw8"}

# --- labelling rule -------------------------------------------------------
# A site is labelled iff its disagreement is BOTH > DISAGREE_FLOOR and in the
# extreme tail of its assay: an elbow located within the top [TAIL_LO, TAIL_HI]
# fraction of the disagreement distribution (largest sudden drop in that band).
DISAGREE_FLOOR = 0.50
TAIL_LO, TAIL_HI = 0.025, 0.05          # top 2.5%..5%
MAX_LABELS = 30                          # legibility guard for very large proteins

# --- functional-map styling ----------------------------------------------
# distinct line colours for binding residues, assigned per distinct ligand name.
LIGAND_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b",
                 "#17becf", "#bcbd22", "#e377c2"]
ACTIVE_COLOR = "#000000"                 # catalytic "Active site" residues
# distinct light fills for subdomains / domains.
SUBDOM_COLORS = ["#ffb74d", "#81c784", "#64b5f6", "#ba68b8", "#a1887f", "#4db6ac"]
DISORDER_COLOR = "#8172b3"               # unchanged
DISEASE_COLOR = "#e08214"                # unchanged
# UniProt calls some substrates generically "substrate"; prettify a few for the
# legend without hard-coding per protein (falls back to the raw ligand name).
LIGAND_ALIAS = {"substrate": "substrate"}
SUBDOM_MAX_FRAC = 0.85                   # skip whole-protein-spanning domains
SUBDOM_MIN_REGION = 15                   # skip tiny motifs annotated as "Region"

# Human proteins with a curated ClinVar export in data/clinvar/ (pathogenic +
# likely pathogenic missense SNVs; per-protein review-status filter documented in
# data/clinvar/README.md). Used instead of the UniProt natural-variant flag.
CLINVAR_TSV = {
    "HXK4_HUMAN": "GCK.tsv",   "PTEN_HUMAN": "PTEN.tsv",  "P53_HUMAN": "TP53.tsv",
    "SRC_HUMAN": "SRC.tsv",    "CP2C9_HUMAN": "CYP2C9.tsv", "KCNE1_HUMAN": "KCNE1.tsv",
    "CAR11_HUMAN": "CARD11.tsv", "VKOR1_HUMAN": "VKORC1.tsv", "RASK_HUMAN": "KRAS.tsv",
    "NPC1_HUMAN": "NPC1.tsv",
}
# ClinVar pathogenic queries are usually dominated by loss-of-function conditions
# and can miss gain-of-function disease variants (e.g. GCK hyperinsulinism/HHF3,
# CARD11 BENTA, KRAS activating). Supplement with UniProt natural-variant positions
# whose disease description carries GoF language.
GOF_KW = ("gain-of-function", "gain of function", "constitutively activ", "activating",
          "increased activ", "hyperactiv", "increases activ", "increase in activ",
          "hhf", "hyperinsulin", "phhi", "benta")
_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
        "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
        "Tyr": "Y", "Val": "V"}
_PCHANGE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")


def clinvar_positions(uni, seq):
    """Residue positions of ClinVar pathogenic/LP missense variants, parsed from
    the canonical p.(3-letter) change in the Name column and kept only when the
    wild-type matches our assay-isoform sequence (drops splice-isoform mismatches,
    e.g. KRAS D153V). Returns None if the protein has no ClinVar export."""
    fn = CLINVAR_TSV.get(uni)
    if not fn:
        return None
    df = pd.read_csv(os.path.join(CLINVAR_DIR, fn), sep="\t", index_col=False)
    pos = set()
    for nm in df["Name"]:
        m = _PCHANGE.search(str(nm))
        if not m:
            continue
        wt, p = _AA3.get(m.group(1), "?"), int(m.group(2))
        if seq and 1 <= p <= len(seq) and seq[p - 1] == wt:
            pos.add(p)
    return pos


def uniprot_disease_supplement(uni):
    """Positions of UniProt natural-variant disease annotations whose description
    carries gain-of-function language (variants a LoF-dominated ClinVar query misses)."""
    path = os.path.join(UNIPROT_DIR, uni + ".json")
    if not os.path.exists(path):
        return set()
    e = json.load(open(path))
    pos = set()
    for f in e.get("features", []):
        if f["type"] == "Natural variant":
            d = (f.get("description", "") or "").lower()
            if any(k in d for k in GOF_KW):
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


def structural_features(uni):
    """Parse the UniProt JSON into per-feature functional annotations, in UniProt
    numbering (== up_pos in the sites CSV). Generalizes 08_gck_structure_map's
    ATP/glucose/subdomain split to arbitrary ligands and regions.

    Returns (ligand_pos, active_pos, subdomains):
      ligand_pos : {ligand_name -> set(positions)}   one entry per distinct ligand
      active_pos : set(positions)                     catalytic Active-site residues
      subdomains : [(name, start, end)]               Region/Domain spans to shade
    """
    path = os.path.join(UNIPROT_DIR, uni + ".json")
    if not os.path.exists(path):
        return {}, set(), []
    e = json.load(open(path))
    L = len(e.get("sequence", {}).get("value", "")) or 10**9
    ligand_pos, active_pos, subdomains = defaultdict(set), set(), []
    for f in e.get("features", []):
        t = f["type"]
        s = f["location"]["start"]["value"]
        en = f["location"]["end"]["value"]
        if s is None or en is None:
            continue
        if t == "Binding site":
            lig = f.get("ligand", {}).get("name", "") or "ligand"
            ligand_pos[LIGAND_ALIAS.get(lig, lig)].update(range(s, en + 1))
        elif t == "Active site":
            active_pos.update(range(s, en + 1))
        elif t in ("Region", "Domain"):
            desc = (f.get("description", "") or "").strip()
            if "disorder" in desc.lower():          # disorder handled as its own track
                continue
            if (en - s + 1) >= SUBDOM_MAX_FRAC * L:  # whole-protein domain -> uninformative
                continue
            if t == "Region" and (en - s + 1) < SUBDOM_MIN_REGION:
                continue                             # drop tiny sequence motifs
            subdomains.append((desc or t, s, en))
    return dict(ligand_pos), active_pos, subdomains


def elbow_tail_threshold(vals, lo=TAIL_LO, hi=TAIL_HI):
    """Disagreement threshold for the extreme tail: keep the top sites down to an
    elbow located within the top [lo, hi] fraction (2.5-5%). Within that rank band
    we cut just above the largest consecutive drop; with no pronounced drop this
    defaults to the wider (hi) end. Returns (threshold, percentile_kept)."""
    v = np.sort(np.asarray(vals, dtype=float))[::-1]     # descending
    n = len(v)
    if n < 3:
        return (float(v.min()) if n else np.inf), (100.0 if n else 0.0)
    i_lo = max(1, int(np.floor(lo * n)))
    i_hi = min(max(i_lo + 1, int(np.ceil(hi * n))), n - 1)
    best_i, best_drop = i_hi, -1.0                       # default: wider (5%) end
    for i in range(i_lo, i_hi + 1):
        drop = v[i - 1] - v[i]
        if drop > best_drop:
            best_drop, best_i = drop, i
    return float(v[best_i - 1]), 100.0 * best_i / n      # keep ranks 0..best_i-1


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
        cv = clinvar_positions(uni, seq)                # None unless a ClinVar export exists
        if cv is not None:
            supp = uniprot_disease_supplement(uni)      # GoF variants ClinVar diabetes query misses
            dz_set = cv | supp
            dz_src = "ClinVar P/LP + UniProt GoF" if supp else "ClinVar P/LP"
        else:
            dz_set = set(g.loc[g["disease"], "up_pos"])
            dz_src = "UniProt"

        # per-ligand binding residues, catalytic sites, subdomains (UniProt coords)
        ligand_pos, active_pos, subdomains = structural_features(uni)
        lig_color = {lig: LIGAND_COLORS[i % len(LIGAND_COLORS)]
                     for i, lig in enumerate(sorted(ligand_pos))}
        sub_color = {i: SUBDOM_COLORS[i % len(SUBDOM_COLORS)]
                     for i in range(len(subdomains))}

        assays = sorted(g["DMS_id"].unique())
        fig, axes = plt.subplots(len(assays), 1, figsize=(16, 2.9 * len(assays) + 1.4),
                                 sharex=True, squeeze=False)
        axes = axes[:, 0]
        for ax, dms_id in zip(axes, assays):
            r = g[g["DMS_id"] == dms_id].sort_values("up_pos")
            ax.plot(r["up_pos"], r["disagreement"], color="#444", lw=0.8, zorder=3)
            ax.fill_between(r["up_pos"], r["disagreement"], color="#bbb", alpha=0.4, zorder=2)
            ymax = max(0.85, r["disagreement"].max() + 0.08)
            # subdomains / domains -> separate light shaded spans
            for i, (_, s0, e0) in enumerate(subdomains):
                ax.axvspan(s0 - 0.5, e0 + 0.5, color=sub_color[i], alpha=0.10, zorder=0)
            # disordered regions -> shaded spans (unchanged)
            for s0, e0 in contiguous(r.loc[r["disordered"], "up_pos"].tolist()):
                ax.axvspan(s0 - 0.5, e0 + 0.5, color=DISORDER_COLOR, alpha=0.12, zorder=0)
            # binding residues -> vertical lines coloured separately per ligand
            for lig, positions in ligand_pos.items():
                for p in positions:
                    ax.axvline(p, color=lig_color[lig], lw=0.6, alpha=0.45, zorder=1)
            # catalytic Active-site residues -> dashed dark lines
            for p in active_pos:
                ax.axvline(p, color=ACTIVE_COLOR, lw=0.7, alpha=0.6,
                           ls=(0, (2, 2)), zorder=1)
            # disease positions -> orange rug at bottom (unchanged; ClinVar where available)
            dz = r.loc[r["up_pos"].isin(dz_set), "up_pos"]
            if len(dz):
                ax.plot(dz, np.full(len(dz), 0.012 * ymax), "|", color=DISEASE_COLOR,
                        markersize=5, mew=0.8, zorder=4)
            # labelled sites: disagreement > floor AND in the elbow tail (top 2.5-5%)
            thr, pct = elbow_tail_threshold(r["disagreement"].values)
            lab_r = r[(r["disagreement"] > DISAGREE_FLOOR) & (r["disagreement"] >= thr)]
            lab_r = lab_r.sort_values("disagreement", ascending=False).head(MAX_LABELS)
            # faint marker of the effective labelling floor
            ax.axhline(max(thr, DISAGREE_FLOOR), color="#999", lw=0.6, ls=":", zorder=1)
            ax.scatter(lab_r["up_pos"], lab_r["disagreement"], color="black", s=22, zorder=5)
            for _, row in lab_r.iterrows():
                pos = int(row["up_pos"])
                wt = seq[pos - 1] if seq and 1 <= pos <= len(seq) else ""
                is_dz = pos in dz_set
                ax.annotate(f"{wt}{pos}", (pos, row["disagreement"]), fontsize=6.5,
                            ha="center", va="bottom", xytext=(0, 3),
                            textcoords="offset points",
                            fontweight="bold" if is_dz else "normal")
            ax.set_ylim(0, ymax)
            ax.set_ylabel(f"disagreement\n[{dms_id.replace(uni + '_', '')}]", fontsize=8)
            print(f"    {dms_id:44s} label>{max(thr, DISAGREE_FLOOR):.2f} "
                  f"(tail {pct:.1f}%)  n_labelled={len(lab_r)}")
        axes[-1].set_xlabel(f"{uni} residue position (UniProt numbering)")
        # legend is built per protein from the features actually present
        legend = []
        for lig in sorted(ligand_pos):
            legend.append(Patch(color=lig_color[lig], label=f"{lig}-binding residue"))
        if active_pos:
            legend.append(plt.Line2D([0], [0], color=ACTIVE_COLOR, ls=(0, (2, 2)),
                                     lw=1, label="active-site residue"))
        for i, (name, s0, e0) in enumerate(subdomains):
            legend.append(Patch(color=sub_color[i], alpha=0.4, label=f"{name} ({s0}-{e0})"))
        legend += [
            Patch(color=DISORDER_COLOR, alpha=0.3, label="disordered region"),
            plt.Line2D([0], [0], marker="|", color=DISEASE_COLOR, lw=0,
                       label=f"disease variant ({dz_src}); bold label"),
            plt.Line2D([0], [0], marker="o", color="black", lw=0,
                       label="labelled site (>0.5 & top 2.5-5%)"),
        ]
        ncol = min(4, max(2, (len(legend) + 1) // 2))
        axes[0].legend(handles=legend, fontsize=7.5, ncol=ncol, loc="upper center",
                       bbox_to_anchor=(0.5, 1.34))
        axes[0].set_title(f"{names.get(uni, uni)} [{uni}] — per-site model disagreement vs functional map",
                          fontsize=11, pad=26)
        plt.tight_layout()
        out = os.path.join(OUTDIR, f"fig_{cname}_structure_map.png")
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"{uni:18s} -> {os.path.basename(out)}  ({len(assays)} assays, {g['up_pos'].nunique()} sites)")


if __name__ == "__main__":
    main()
