"""
Candidate proteins for testing whether the HXK4_HUMAN high- vs low-agreement
findings generalize — i.e. the multi-assay proteins (>=2 DMS assays), the same
set used by the functional-disagreement generalization (15).

For each we assemble the selection-relevant metadata (n_assays, taxon, seq_len,
selection types, multi-phenotype flag, MSA_Neff_L) and the available RCSB PDB
structures pulled from the cached UniProt cross-references (data/uniprot/*.json),
so a real structure can be chosen for a 3-D analysis like 17 (which used 1V4S for
GCK). ProteinGym also ships its own <UniProt_ID>.pdb (column proteingym_pdb).

Outputs:
  results/generalization_candidate_proteins.csv   (all PDB options per protein)
  figures/fig_generalization_candidates.png       (visual selection sheet)
"""
import os
import re
import json
import urllib.request
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import common as C

UNIPROT_DIR = os.path.join(C.ANALYSIS, "data", "uniprot")
LIGAND_CACHE = os.path.join(C.ANALYSIS, "data", "pdb_ligands_cache.json")

# crystallization additives / ions / glycans to ignore when deciding "holo"
LIGAND_BLACKLIST = {
    "HOH", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "PE4", "PE3",
    "MPD", "TRS", "EPE", "MES", "BME", "DMS", "DMSO", "IPA", "ACT", "ACY", "FMT",
    "CIT", "FLC", "TLA", "MLA", "MLI", "SCN", "AZI", "NO3", "CO3", "SO4", "PO4",
    "NH4", "IMD", "BOG", "LDA", "MRD", "BU3", "BTB", "TAR", "OGA", "144",
    "CAC", "P4G", "PGO", "PG5", "PG6", "PG0", "XPE", "7PE", "12P", "15P", "33O",
    "DIO", "PDO", "MG8", "C8E", "DMU", "SDS", "OCT", "D10", "9FO", "AE3", "AE4",
    "NAG", "BMA", "MAN", "FUC", "NDG", "GAL", "SIA", "GLA", "XYP", "SUC", "TRE",
    "NA", "K", "CL", "BR", "IOD", "F", "MG", "CA", "ZN", "MN", "FE", "FE2",
    "NI", "CO", "CU", "CU1", "CD", "HG", "SR", "CS", "RB", "LI", "BA", "PB",
}
LIGAND_MW_MIN = 130.0     # excludes lone ions and small buffers


def _load_cache():
    if os.path.exists(LIGAND_CACHE):
        return json.load(open(LIGAND_CACHE))
    return {}


def qualifying_ligands(raw):
    """Filter a cached raw [[comp_id, mw], ...] list to biologically-interesting
    ligands (drop additives/ions/glycans and anything below the mass floor)."""
    return [cid for cid, mw in raw
            if cid not in LIGAND_BLACKLIST and (mw or 0) >= LIGAND_MW_MIN]


def fetch_ligands(pdb_ids, cache):
    """Populate cache[pdb_id] = [[comp_id, mw], ...] (RAW, so blacklist changes
    take effect without re-fetching); batch-query RCSB GraphQL for uncached ids
    (skips silently if offline)."""
    need = [p for p in pdb_ids if p not in cache]
    for i in range(0, len(need), 50):
        chunk = need[i:i + 50]
        q = ('{entries(entry_ids:[%s]){rcsb_id nonpolymer_entities{nonpolymer_comp'
             '{chem_comp{id formula_weight}}}}}' % ",".join(f'"{p}"' for p in chunk))
        try:
            req = urllib.request.Request(
                "https://data.rcsb.org/graphql",
                data=json.dumps({"query": q}).encode(),
                headers={"Content-Type": "application/json"})
            data = json.load(urllib.request.urlopen(req, timeout=40))
            for e in (data.get("data", {}).get("entries") or []):
                comps = [[ent["nonpolymer_comp"]["chem_comp"]["id"],
                          ent["nonpolymer_comp"]["chem_comp"].get("formula_weight") or 0]
                         for ent in (e.get("nonpolymer_entities") or [])]
                cache[e["rcsb_id"]] = comps
            for p in chunk:                     # entries with no nonpolymer come back absent
                cache.setdefault(p, [])
        except Exception as ex:
            print(f"  [ligand fetch skipped for {len(chunk)} PDBs: {ex}]")
            break
    return cache
TAXON_COLOR = {"Human": "#4c72b0", "Eukaryote": "#55a868",
               "Prokaryote": "#c44e52", "Virus": "#8172b3"}


def protein_name(uni, fallback=""):
    """Human-readable label from the cached UniProt entry: 'GENE (Full name)',
    falling back to the reference molecule_name for sparse TrEMBL entries."""
    p = os.path.join(UNIPROT_DIR, uni + ".json")
    gene = name = None
    if os.path.exists(p):
        e = json.load(open(p))
        genes = e.get("genes", [])
        gene = genes[0].get("geneName", {}).get("value") if genes else None
        pdsc = e.get("proteinDescription", {})
        name = pdsc.get("recommendedName", {}).get("fullName", {}).get("value")
        if not name:
            subs = pdsc.get("submissionNames", [])
            name = subs[0].get("fullName", {}).get("value") if subs else None
    name = name or fallback or "?"
    return f"{gene} ({name})" if gene else name


def parse_coverage(chains, seq_len):
    """Union of residue ranges in a UniProt PDB 'Chains' string -> covered fraction."""
    covered = set()
    for m in re.finditer(r"(\d+)-(\d+)", chains or ""):
        covered.update(range(int(m.group(1)), int(m.group(2)) + 1))
    return len(covered) / seq_len if seq_len else 0.0


def pdb_options(uni, seq_len):
    p = os.path.join(UNIPROT_DIR, uni + ".json")
    if not os.path.exists(p):
        return []
    e = json.load(open(p))
    out = []
    for x in e.get("uniProtKBCrossReferences", []):
        if x.get("database") != "PDB":
            continue
        pr = {q["key"]: q["value"] for q in x.get("properties", [])}
        res = pr.get("Resolution", "")
        res_val = float(res.split()[0]) if res and res[0].isdigit() else np.nan
        out.append({"pdb": x["id"], "method": pr.get("Method", ""),
                    "resolution": res_val, "chains": pr.get("Chains", ""),
                    "coverage": round(parse_coverage(pr.get("Chains", ""), seq_len), 2)})
    return out


def best_pdb(options):
    """Prefer X-ray, maximize residue coverage first, then best resolution — so a
    near-full structure wins over a high-res peptide fragment."""
    xray = [o for o in options if o["method"] == "X-ray" and not np.isnan(o["resolution"])]
    pool = xray or [o for o in options if not np.isnan(o.get("resolution", np.nan))] or options
    if not pool:
        return None
    return sorted(pool, key=lambda o: (-o["coverage"], o.get("resolution", 99)))[0]


def main():
    ref = C.load_reference()
    multi = ref.groupby("UniProt_ID").filter(lambda g: len(g) >= 2)
    cache = _load_cache()
    rows = []
    for uni, g in multi.groupby("UniProt_ID"):
        seq_len = int(g["seq_len"].iloc[0])
        opts = pdb_options(uni, seq_len)
        bp = best_pdb(opts)
        xray_sorted = sorted([o for o in opts if o["method"] == "X-ray"],
                             key=lambda o: (-o["coverage"], o.get("resolution", 99)))
        # look up ligands for the top candidates (bounded), pick the best holo one
        cand = xray_sorted[:40]
        fetch_ligands([o["pdb"] for o in cand], cache)
        for o in cand:
            o["ligands"] = qualifying_ligands(cache.get(o["pdb"], []))
            o["holo"] = len(o["ligands"]) > 0
        holo = next((o for o in cand if o["holo"]), None)     # first = best coverage/res holo
        bp_holo = bool(qualifying_ligands(cache.get(bp["pdb"], []))) if bp else False
        top3 = sorted([o for o in opts if o["method"] == "X-ray"],
                      key=lambda o: (o["coverage"] < 0.7, o.get("resolution", 99)))[:3]
        rows.append({
            "UniProt_ID": uni,
            "protein": protein_name(uni, fallback=str(g["molecule_name"].iloc[0])),
            "n_assays": len(g),
            "taxon": g["taxon"].iloc[0],
            "seq_len": seq_len,
            "selection_types": ";".join(sorted(set(g["coarse_selection_type"]))),
            "multi_phenotype": g["coarse_selection_type"].nunique() > 1,
            "MSA_Neff_L": round(float(g["MSA_Neff_L"].iloc[0]), 2),
            "n_pdb": len(opts),
            "best_pdb": bp["pdb"] if bp else "",
            "best_pdb_res": bp["resolution"] if bp else np.nan,
            "best_pdb_coverage": bp["coverage"] if bp else np.nan,
            "best_pdb_has_ligand": bp_holo,
            "best_holo_pdb": holo["pdb"] if holo else "",
            "best_holo_res": holo["resolution"] if holo else np.nan,
            "best_holo_coverage": holo["coverage"] if holo else np.nan,
            "best_holo_ligands": ";".join(holo["ligands"]) if holo else "",
            "top_xray_pdbs": ", ".join(f"{o['pdb']}({o['resolution']}A,{int(o['coverage']*100)}%)" for o in top3),
            "proteingym_pdb": g["pdb_file"].iloc[0],
        })
    os.makedirs(os.path.dirname(LIGAND_CACHE), exist_ok=True)
    json.dump(cache, open(LIGAND_CACHE, "w"), indent=0, sort_keys=True)
    tab = pd.DataFrame(rows).sort_values(["multi_phenotype", "MSA_Neff_L"], ascending=[False, True])
    tab.to_csv(os.path.join(C.RESULTS, "generalization_candidate_proteins.csv"), index=False)
    print(f"{len(tab)} candidate proteins written.")
    with pd.option_context("display.max_colwidth", 30, "display.width", 200):
        print(tab[["UniProt_ID", "n_assays", "multi_phenotype", "MSA_Neff_L",
                   "best_pdb", "best_pdb_has_ligand", "best_holo_pdb", "best_holo_ligands"]].to_string(index=False))

    # ---------- visual selection sheet ----------
    t = tab.iloc[::-1].reset_index(drop=True)   # so highest interest ends up on top
    y = np.arange(len(t))
    neff = t["MSA_Neff_L"].clip(lower=0.03)     # for log display (0 -> 0.03)
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.barh(y, neff, color=[TAXON_COLOR[x] for x in t["taxon"]], height=0.6, zorder=3)
    ax.set_xscale("log")
    ax.axvspan(0.01, 1, color="red", alpha=0.06, zorder=0)          # too shallow
    ax.axvspan(1000, 1e5, color="orange", alpha=0.06, zorder=0)     # very deep
    ax.axvline(1, color="red", ls="--", lw=0.8); ax.axvline(1000, color="orange", ls="--", lw=0.8)
    ax.set_yticks(y)
    labels = [f"{r.UniProt_ID}{'  ★' if r.multi_phenotype else ''}" for r in t.itertuples()]
    ax.set_yticklabels(labels, fontsize=8)
    for tick, tax in zip(ax.get_yticklabels(), t["taxon"]):
        tick.set_color(TAXON_COLOR[tax])
    ax.set_xlim(0.02, 1e5)
    ax.set_xlabel("MSA_Neff_L  (log scale — shallow <1 shaded red, very deep >1000 shaded orange)")
    ax.set_title("Candidate proteins for generalizing the high- vs low-agreement findings\n"
                 "★ = multi-phenotype (>1 selection type); bar = MSA depth; text = assays · selection · best PDB",
                 fontsize=12)
    # annotate each bar with n_assays, selection types, best PDB
    for yi, r in zip(y, t.itertuples()):
        sel = r.selection_types.replace("OrganismalFitness", "OrgFit").replace("Expression", "Expr")
        if r.best_holo_pdb:                                   # prefer a ligand-bound structure
            pdb = f"holo {r.best_holo_pdb}({r.best_holo_ligands}) / apo-best {r.best_pdb}"
        elif r.best_pdb:
            pdb = f"{r.best_pdb}({r.best_pdb_res}Å,{int(r.best_pdb_coverage*100)}%) [no holo found]"
        else:
            pdb = "no PDB"
        ax.text(1.05e5, yi, f"{r.n_assays} assays · {sel} · {r.seq_len} aa · {pdb}",
                va="center", ha="left", fontsize=6.5)
    ax.set_xlim(0.02, 1e5)
    handles = [Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=9, label=k)
               for k, c in TAXON_COLOR.items()]
    ax.legend(handles=handles, title="taxon", fontsize=8, loc="lower right")
    plt.tight_layout()
    out = os.path.join(C.FIGURES, "fig_generalization_candidates.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")
    print("Saved results/generalization_candidate_proteins.csv")


if __name__ == "__main__":
    main()
