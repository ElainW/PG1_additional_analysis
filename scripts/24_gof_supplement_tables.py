"""
Write down the gain-of-function (GoF) supplement variants used by
23_structure_maps_all.py — i.e. the UniProt natural-variant disease annotations
whose description matches the GoF keyword list (GOF_KW) — as a separate table per
protein, plus a combined summary.

These are the variants a loss-of-function-dominated ClinVar query misses (GCK
hyperinsulinism/HHF3, CARD11 BENTA, KRAS activating), added back to the disease
annotation. Source: data/uniprot/raw_json/<UniProt_ID>.json.

Outputs (in data/uniprot/additional annotations/):
  gof_supplement_<UniProt_ID>.csv   one per protein that has >=1 GoF match
  gof_supplement_summary.csv        all ClinVar proteins (n matches + positions)
"""
import os
import json
import pandas as pd
import common as C
from importlib import import_module

sm = import_module("23_structure_maps_all")   # reuse GOF_KW, CLINVAR_TSV, UNIPROT_DIR
OUTDIR = os.path.join(C.ANALYSIS, "data", "uniprot", "additional annotations")


def gene_name(uni):
    p = os.path.join(sm.UNIPROT_DIR, uni + ".json")
    genes = json.load(open(p)).get("genes", []) if os.path.exists(p) else []
    return genes[0].get("geneName", {}).get("value", "") if genes else ""


def gof_variants(uni):
    """One row per (position, wt, alt) GoF-matched natural variant, with the
    matched keyword(s) and the full UniProt description."""
    path = os.path.join(sm.UNIPROT_DIR, uni + ".json")
    if not os.path.exists(path):
        return pd.DataFrame()
    e = json.load(open(path))
    rows = []
    for f in e.get("features", []):
        if f["type"] != "Natural variant":
            continue
        desc = f.get("description", "") or ""
        hits = [k for k in sm.GOF_KW if k in desc.lower()]
        if not hits:
            continue
        alt = f.get("alternativeSequence", {})
        wt = alt.get("originalSequence", "")
        pos = f["location"]["start"]["value"]
        for a in (alt.get("alternativeSequences", []) or [""]):
            rows.append({
                "UniProt_ID": uni, "gene": gene_name(uni), "position": pos,
                "wt": wt, "alt": a, "variant": f"{wt}{pos}{a}",
                "matched_keywords": "; ".join(hits), "description": desc,
            })
    return pd.DataFrame(rows).sort_values("position") if rows else pd.DataFrame()


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    summary = []
    for uni in sm.CLINVAR_TSV:
        gv = gof_variants(uni)
        positions = sorted(int(p) for p in gv["position"].unique()) if not gv.empty else []
        summary.append({"UniProt_ID": uni, "gene": gene_name(uni),
                        "n_gof_variants": len(gv), "n_gof_positions": len(positions),
                        "positions": "; ".join(map(str, positions))})
        if not gv.empty:
            out = os.path.join(OUTDIR, f"gof_supplement_{uni}.csv")
            gv.to_csv(out, index=False)
            print(f"{uni:13s} {len(gv):2d} GoF variants at {positions} -> {os.path.basename(out)}")
        else:
            print(f"{uni:13s}  0 GoF variants")
    pd.DataFrame(summary).to_csv(os.path.join(OUTDIR, "gof_supplement_summary.csv"), index=False)
    print(f"\nSaved summary + per-protein tables to '{OUTDIR}'")


if __name__ == "__main__":
    main()
