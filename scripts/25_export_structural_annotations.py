"""
Write down, per protein, the structural / functional annotations that
23_structure_maps_all.py overlays on each disagreement map, so the exact feature
set used by the figures is persisted as a plain-text reference.

For every candidate protein (results/generalization_candidate_proteins.csv) we
record four annotation categories, all in UniProt numbering (== up_pos):
  * binding residues, listed separately per ligand (UniProt "Binding site")
  * catalytic active-site residues (UniProt "Active site")   [reported alongside]
  * subdomain / domain ranges (UniProt "Region" / "Domain", excl. disorder)
  * disordered region residue ranges (UniProt "Region" = Disordered)
  * disease-variant positions, using the SAME source 23 uses for its rug:
      - the 10 proteins with a data/clinvar/ export -> ClinVar P/LP (WT-matched)
        + the UniProt gain-of-function supplement;
      - every other protein          -> UniProt natural-variant disease flag
        (as already materialised in the sites CSV's `disease` column).

The parsing helpers are imported verbatim from 23_structure_maps_all.py so the
written-down annotations cannot drift from the figures.

Outputs (in data/uniprot/additional_annotations/):
  structural_annotations_all.csv   long format, one row per (protein, category, name)
  structural_annotations.md        human-readable, one section per protein
"""
import os
import importlib.util
import pandas as pd
import common as C

# --- import 23's parsing helpers (module name can't start with a digit) -------
_S23 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "23_structure_maps_all.py")
_spec = importlib.util.spec_from_file_location("structure_maps_all", _S23)
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)

OUTDIR = os.path.join(C.ANALYSIS, "data", "uniprot", "additional_annotations")


def fmt_ranges(positions):
    """Compact a set/list of ints into 'a-b, c, e-f'."""
    ps = sorted(set(int(p) for p in positions))
    if not ps:
        return ""
    out, s, prev = [], ps[0], ps[0]
    for p in ps[1:]:
        if p == prev + 1:
            prev = p
        else:
            out.append(f"{s}-{prev}" if s != prev else f"{s}")
            s = prev = p
    out.append(f"{s}-{prev}" if s != prev else f"{s}")
    return ", ".join(out)


def disease_positions(uni, seq, sites):
    """Reproduce 23's dz_set: ClinVar P/LP (+GoF supplement) where an export
    exists, otherwise the UniProt disease flag carried in the sites CSV."""
    cv = S.clinvar_positions(uni, seq)
    if cv is not None:
        supp = S.uniprot_disease_supplement(uni)
        src = "ClinVar P/LP + UniProt GoF" if supp else "ClinVar P/LP"
        return (cv | supp), src
    g = sites[sites["UniProt_ID"] == uni]
    return set(g.loc[g["disease"], "up_pos"].astype(int)), "UniProt natural-variant"


def disordered_ranges(uni):
    """Raw UniProt disordered Region spans (start, end)."""
    path = os.path.join(S.UNIPROT_DIR, uni + ".json")
    if not os.path.exists(path):
        return []
    import json
    e = json.load(open(path))
    out = []
    for f in e.get("features", []):
        if f["type"] == "Region" and "disorder" in (f.get("description", "") or "").lower():
            out.append((f["location"]["start"]["value"], f["location"]["end"]["value"]))
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    cand = pd.read_csv(S.CANDID)
    genes = {r["UniProt_ID"]: str(r["protein"]).split(" (")[0].strip()
             for _, r in cand.iterrows()}
    sites = pd.read_csv(S.SITES)

    rows = []          # long-format CSV rows
    md = ["# Structural / functional annotations per protein",
          "",
          "UniProt numbering (matches `up_pos` in "
          "`generalization_functional_disagreement_sites.csv`). Feature source: "
          "UniProt `raw_json/<UniProt_ID>.json`, parsed exactly as in "
          "`23_structure_maps_all.py`. Disease-variant source noted per protein.",
          ""]

    for uni in cand["UniProt_ID"]:
        gene = genes.get(uni, uni)
        seq = S.uniprot_seq(uni)
        ligand_pos, active_pos, subdomains = S.structural_features(uni)
        dis = disordered_ranges(uni)
        dz, dz_src = disease_positions(uni, seq, sites)

        md.append(f"## {gene} ({uni})")
        # ---- binding residues, per ligand ----
        if ligand_pos:
            for lig in sorted(ligand_pos):
                pos = ligand_pos[lig]
                rows.append([uni, gene, "binding", lig, len(pos), fmt_ranges(pos)])
                md.append(f"- **{lig}-binding** ({len(pos)}): {fmt_ranges(pos)}")
        else:
            md.append("- **binding**: none annotated")
        # ---- active-site residues ----
        if active_pos:
            rows.append([uni, gene, "active_site", "Active site",
                         len(active_pos), fmt_ranges(active_pos)])
            md.append(f"- **active site** ({len(active_pos)}): {fmt_ranges(active_pos)}")
        # ---- subdomains / domains ----
        if subdomains:
            for name, s0, e0 in subdomains:
                rows.append([uni, gene, "subdomain", name, e0 - s0 + 1, f"{s0}-{e0}"])
                md.append(f"- **subdomain** — {name}: {s0}-{e0}")
        else:
            md.append("- **subdomains**: none annotated")
        # ---- disordered regions ----
        if dis:
            for s0, e0 in dis:
                rows.append([uni, gene, "disordered", "Disordered", e0 - s0 + 1, f"{s0}-{e0}"])
            md.append("- **disordered**: " + "; ".join(f"{s0}-{e0}" for s0, e0 in dis))
        else:
            md.append("- **disordered**: none annotated")
        # ---- disease variants ----
        rows.append([uni, gene, "disease", dz_src, len(dz), fmt_ranges(dz)])
        md.append(f"- **disease variants** ({len(dz)}, {dz_src}): "
                  + (fmt_ranges(dz) or "none"))
        md.append("")

        print(f"{uni:16s} ligands={len(ligand_pos):2d} active={len(active_pos)} "
              f"subdom={len(subdomains)} disorder={len(dis)} disease={len(dz)}")

    df = pd.DataFrame(rows, columns=["UniProt_ID", "gene", "category", "name",
                                     "n_residues", "residues"])
    csv_out = os.path.join(OUTDIR, "structural_annotations_all.csv")
    md_out = os.path.join(OUTDIR, "structural_annotations.md")
    df.to_csv(csv_out, index=False)
    open(md_out, "w").write("\n".join(md))
    print(f"\nwrote {csv_out}  ({len(df)} rows)")
    print(f"wrote {md_out}")


if __name__ == "__main__":
    main()
