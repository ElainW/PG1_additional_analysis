"""
Where do the high-disagreement sites land? -- descriptive characterization for
the paper (no enrichment test; see 55 for that).

For the five allosteric-map proteins (CYP2C9, GCK, KRAS, PTEN, SRC) we take the
high-disagreement sites exported by 27 (results/high_disagreement_sites/), collapse
them to unique residue positions, and label each by the region it sits in:

  terminal   -- beginning/end of the protein: within the first/last TERM_WIN
                residues, or inside a UniProt disordered segment that abuts a
                terminus (SRC 1-53, PTEN 352-403). These are the flexible tails.
  loop_turn  -- secondary structure is loop or turn (52 / UniProt SS), i.e. not
                helix or strand.
  functional -- binding/active/allosteric site or a disease residue, defined as the
                UNION of: the curated allosteric-map atlas (26/47, dir1), UniProt
                binding+-2 & active site (the pool's near_binding, 15/dir3), dir2
                ClinVar pathogenic/LP (parsed as in 23), and UniProt natural-variant
                disease (the pool's disease column). Sub-membership is reported too.
  core       -- everything else: structured (helix/strand), non-terminal, non-functional.

A site can satisfy several labels (a functional loop, a terminal disordered residue).
The headline table is a MUTUALLY EXCLUSIVE classification under the priority
functional > terminal > loop_turn > core, so the columns sum to 100% and read as
"of the high-disagreement sites, this fraction are functional, of the rest this many
are terminal, ...". The raw non-exclusive flags (pct_terminal_any, pct_loopturn_any,
pct_functional_any and each functional sub-type) are written alongside so nothing is
hidden by the priority choice.

Unit is the unique residue position (a residue flagged high in two assays counts
once); n_high_rows keeps the raw (assay x position) count for reference. Positions
are UniProt-numbered; these five assays are offset 0 so up_pos == SS/annotation pos.

Outputs -> results/high_disagreement_site_characterization/
  site_characterization_by_protein.csv   one row per protein: n sites + every
                                         percentage (exclusive classes + raw flags).
  site_characterization_sites.csv        one row per high-disagreement position with
                                         all its labels (audit trail).
"""
import os
import re
import json
import numpy as np
import pandas as pd
import common as C

TERM_WIN = 20   # first/last N residues counted as beginning/end of the protein

# gene, UniProt entry, allosteric-map file, tail file, ClinVar export
PROTEINS = [
    ("CYP2C9", "CP2C9_HUMAN", "CYP2C9_HUMAN.tsv", "cyp2c9_high_disagreement.tsv", "CYP2C9.tsv"),
    ("GCK",    "HXK4_HUMAN",  "GCK_HUMAN.tsv",    "gck_high_disagreement.tsv",    "GCK.tsv"),
    ("KRAS",   "RASK_HUMAN",  "RASK_HUMAN.tsv",   "kras_high_disagreement.tsv",   "KRAS.tsv"),
    ("PTEN",   "PTEN_HUMAN",  "PTEN_HUMAN.tsv",   "pten_high_disagreement.tsv",   "PTEN.tsv"),
    ("SRC",    "SRC_HUMAN",   "SRC_HUMAN.tsv",    "src_high_disagreement.tsv",    "SRC.tsv"),
]

ALLO_DIR = os.path.join(C.RESULTS, "allosteric_map")
CLINVAR_DIR = os.path.join(C.DATA_DIR, "clinvar")
UNIPROT_DIR = os.path.join(C.DATA_DIR, "uniprot", "raw_json")
ANNOT_ALL = os.path.join(C.DATA_DIR, "uniprot", "additional_annotations",
                         "structural_annotations_all.csv")
SS_ALL = os.path.join(C.RESULTS, "structure_features", "secondary_structure",
                      "secondary_structure_all.csv")
POOL_PATH = C.rp("generalization_functional_disagreement_sites.csv")
TAIL_DIR = C.rp("high_disagreement_sites")
OUTDIR = C.rp("high_disagreement_site_characterization")

# --- ClinVar parsing, as in 23 / 55 ---------------------------------------------
GOF_KW = ("gain-of-function", "gain of function", "constitutively activ", "activating",
          "increased activ", "hyperactiv", "increases activ", "increase in activ",
          "hhf", "hyperinsulin", "phhi", "benta")
_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
        "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
        "Tyr": "Y", "Val": "V"}
_PCHANGE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")


def clinvar_positions(clinvar_file, seq):
    df = pd.read_csv(os.path.join(CLINVAR_DIR, clinvar_file), sep="\t", index_col=False)
    pos = set()
    for nm in df["Name"]:
        m = _PCHANGE.search(str(nm))
        if not m:
            continue
        wt, p = _AA3.get(m.group(1), "?"), int(m.group(2))
        if seq and 1 <= p <= len(seq) and seq[p - 1] == wt:
            pos.add(p)
    return pos


def uniprot_gof_supplement(uni):
    path = os.path.join(UNIPROT_DIR, uni + ".json")
    if not os.path.exists(path):
        return set()
    e = json.load(open(path))
    return {f["location"]["start"]["value"] for f in e.get("features", [])
            if f["type"] == "Natural variant"
            and any(k in (f.get("description", "") or "").lower() for k in GOF_KW)}


def allosteric_positions(map_file):
    df = pd.read_csv(os.path.join(ALLO_DIR, map_file), sep="\t", comment="#")
    return set(df["position"].astype(int))


def terminal_disordered(upid, L):
    """Positions in a UniProt disordered segment that touches a terminus (start<=1 or
    end>=L) -- the flexible N/C tails (SRC 1-53, PTEN 352-403)."""
    ann = pd.read_csv(ANNOT_ALL)
    seg = ann[(ann["UniProt_ID"] == upid) & (ann["category"] == "disordered")]
    out = set()
    for r in seg["residues"]:
        for span in str(r).split(";"):
            span = span.strip()
            if "-" in span:
                a, b = (int(x) for x in span.split("-"))
            else:
                a = b = int(span)
            if a <= 1 or b >= L:
                out.update(range(a, b + 1))
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    pool = pd.read_csv(POOL_PATH)
    pool["up_pos"] = pool["up_pos"].astype(int)
    ref = C.drop_excluded_dms(C.load_reference())
    seqmap = ref.groupby("UniProt_ID")["target_seq"].first()

    ss = pd.read_csv(SS_ALL)
    ss_map = {(p, int(q)): s for p, q, s in zip(ss["protein"], ss["pos"], ss["ss"])}
    length = ss.groupby("protein")["pos"].max().to_dict()

    per_protein, per_site = [], []

    for gene, upid, map_file, tail_file, clinvar_file in PROTEINS:
        L = int(length[upid])
        u = pool[pool["UniProt_ID"] == upid]
        seq = seqmap.get(upid, "")

        allo = allosteric_positions(map_file)
        binding = set(u.loc[u["near_binding"] == True, "up_pos"])           # noqa: E712
        uni_disease = set(u.loc[u["disease"] == True, "up_pos"])            # noqa: E712
        clinvar = clinvar_positions(clinvar_file, seq) | uniprot_gof_supplement(upid)
        term_dis = terminal_disordered(upid, L)

        tail = pd.read_csv(os.path.join(TAIL_DIR, tail_file), sep="\t")
        tail["up_pos"] = tail["up_pos"].astype(int)
        n_rows = len(tail)
        positions = sorted(set(tail["up_pos"]))

        rows = []
        for p in positions:
            s = ss_map.get((upid, p), "NA")
            f_allo = p in allo
            f_bind = p in binding
            f_clin = p in clinvar
            f_udis = p in uni_disease
            functional = f_allo or f_bind or f_clin or f_udis
            terminal = (p <= TERM_WIN) or (p >= L - TERM_WIN + 1) or (p in term_dis)
            loopturn = s in ("loop", "turn")
            # mutually exclusive priority: functional > terminal > loop_turn > core
            if functional:
                cls = "functional"
            elif terminal:
                cls = "terminal"
            elif loopturn:
                cls = "loop_turn"
            else:
                cls = "core"
            rows.append(dict(protein=gene, UniProt_ID=upid, up_pos=p, ss=s, length=L,
                             is_functional=functional, is_terminal=terminal,
                             is_loop_turn=loopturn, exclusive_class=cls,
                             func_allosteric=f_allo, func_uniprot_binding=f_bind,
                             func_clinvar=f_clin, func_uniprot_disease=f_udis))
        d = pd.DataFrame(rows)
        per_site.append(d)
        n = len(d)

        def pct(mask):
            return round(100 * mask.sum() / n, 1) if n else np.nan

        per_protein.append(dict(
            protein=gene, UniProt_ID=upid, length=L,
            n_high_sites=n, n_high_rows=n_rows,
            # mutually exclusive classes (sum to 100%)
            pct_terminal=pct(d.exclusive_class == "terminal"),
            pct_loop_turn=pct(d.exclusive_class == "loop_turn"),
            pct_functional=pct(d.exclusive_class == "functional"),
            pct_core=pct(d.exclusive_class == "core"),
            # raw non-exclusive flags
            pct_terminal_any=pct(d.is_terminal),
            pct_loopturn_any=pct(d.is_loop_turn),
            pct_functional_any=pct(d.is_functional),
            # functional sub-types (non-exclusive)
            pct_func_allosteric=pct(d.func_allosteric),
            pct_func_uniprot_binding=pct(d.func_uniprot_binding),
            pct_func_clinvar=pct(d.func_clinvar),
            pct_func_uniprot_disease=pct(d.func_uniprot_disease),
        ))

    bp = pd.DataFrame(per_protein)
    sites = pd.concat(per_site, ignore_index=True)
    bp.to_csv(os.path.join(OUTDIR, "site_characterization_by_protein.csv"), index=False)
    sites.to_csv(os.path.join(OUTDIR, "site_characterization_sites.csv"), index=False)

    # simplified paper table: per protein, the mutually exclusive region breakdown
    simple = bp[["protein", "length", "n_high_sites",
                 "pct_functional", "pct_terminal", "pct_loop_turn", "pct_core"]].rename(
        columns={"n_high_sites": "n_high_disagreement_sites",
                 "pct_terminal": "pct_beginning_end"})
    simple.to_csv(os.path.join(OUTDIR, "site_characterization_simplified.csv"), index=False)

    pd.set_option("display.width", 220, "display.max_columns", 40)
    print(f"High-disagreement site characterization  (unit = unique residue; TERM_WIN={TERM_WIN})\n")
    print("Mutually exclusive region classes (priority functional > terminal > loop_turn > core):")
    print(bp[["protein", "length", "n_high_sites", "n_high_rows",
              "pct_functional", "pct_terminal", "pct_loop_turn", "pct_core"]].to_string(index=False))
    print("\nRaw non-exclusive flags + functional sub-types (%):")
    print(bp[["protein", "pct_functional_any", "pct_terminal_any", "pct_loopturn_any",
              "pct_func_allosteric", "pct_func_uniprot_binding", "pct_func_clinvar",
              "pct_func_uniprot_disease"]].to_string(index=False))
    tot = len(sites)
    print(f"\nPooled over {tot} sites — exclusive:  "
          f"functional {100*(sites.exclusive_class=='functional').mean():.0f}%  "
          f"terminal {100*(sites.exclusive_class=='terminal').mean():.0f}%  "
          f"loop_turn {100*(sites.exclusive_class=='loop_turn').mean():.0f}%  "
          f"core {100*(sites.exclusive_class=='core').mean():.0f}%")
    print(f"wrote 2 csv -> {OUTDIR}")


if __name__ == "__main__":
    main()
