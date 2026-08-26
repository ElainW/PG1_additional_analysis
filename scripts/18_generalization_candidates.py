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
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import common as C

UNIPROT_DIR = os.path.join(C.ANALYSIS, "data", "uniprot")
TAXON_COLOR = {"Human": "#4c72b0", "Eukaryote": "#55a868",
               "Prokaryote": "#c44e52", "Virus": "#8172b3"}


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
    rows = []
    for uni, g in multi.groupby("UniProt_ID"):
        seq_len = int(g["seq_len"].iloc[0])
        opts = pdb_options(uni, seq_len)
        bp = best_pdb(opts)
        top3 = sorted([o for o in opts if o["method"] == "X-ray"],
                      key=lambda o: (o["coverage"] < 0.7, o.get("resolution", 99)))[:3]
        rows.append({
            "UniProt_ID": uni,
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
            "top_xray_pdbs": ", ".join(f"{o['pdb']}({o['resolution']}A,{int(o['coverage']*100)}%)" for o in top3),
            "proteingym_pdb": g["pdb_file"].iloc[0],
        })
    tab = pd.DataFrame(rows).sort_values(["multi_phenotype", "MSA_Neff_L"], ascending=[False, True])
    tab.to_csv(os.path.join(C.RESULTS, "generalization_candidate_proteins.csv"), index=False)
    print(f"{len(tab)} candidate proteins written.")
    print(tab[["UniProt_ID", "n_assays", "taxon", "multi_phenotype", "MSA_Neff_L",
               "best_pdb", "best_pdb_res", "best_pdb_coverage", "n_pdb"]].to_string(index=False))

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
        pdb = f"{r.best_pdb}({r.best_pdb_res}Å,{int(r.best_pdb_coverage*100)}%)" if r.best_pdb else "no PDB"
        ax.text(1.05e5, yi, f"{r.n_assays} assays · {sel} · {r.seq_len} aa · PDB {pdb}",
                va="center", ha="left", fontsize=7)
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
