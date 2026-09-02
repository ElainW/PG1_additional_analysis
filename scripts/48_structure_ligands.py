"""
Ligand / substrate inventory for the curated PDB structures.

For every selected structure (data/structures/pdb/) we read the HETATM records
directly from the coordinates, classify each non-water heteroatom group by its
FUNCTIONAL role, and record any bound protein partner chain (for the co-complexes).
The question answered per structure: is a functional ligand / substrate actually
bound, or is the site empty (apo / additives only)?

Role classes
  substrate     orthosteric substrate or substrate-analog in the active site
  cofactor      prosthetic group (heme, FAD, ...)
  nucleotide    ATP/GTP analog or GDP at the active/nucleotide site
  allosteric    allosteric modulator (GCK activators)
  inhibitor     active-site inhibitor / drug
  modification  covalent PTM modeled as a HETATM (phospho-Tyr, chromophore)
  metal         functional metal ion (Mg2+ coordinating a nucleotide)
  ion/additive  crystallization ion / cryoprotectant / lipid (not functional)

A structure is flagged FUNCTIONAL-LIGAND = yes if it carries any substrate,
cofactor, nucleotide, allosteric modulator, or inhibitor, or a bound protein partner.

Outputs (results/structure_features/):
  ligands_by_structure_long.csv   one row per (structure, ligand)
  ligands_by_structure.csv        one row per structure (wide summary)
  ligands_by_structure.md         human-readable table
"""
import os
import sys
from collections import defaultdict, OrderedDict

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

PDB_DIR = os.path.join(C.DATA_DIR, "structures", "pdb")
OUTDIR = os.path.join(C.RESULTS, "structure_features")
os.makedirs(OUTDIR, exist_ok=True)

# --- HET code -> (readable name, role) -------------------------------------
# roles: substrate, cofactor, nucleotide, allosteric, inhibitor, modification,
#        metal, ion/additive, lipid
HET = {
    "GLC": ("alpha-D-glucose", "substrate"),
    "MRK": ("GCK allosteric activator (MRK-compound)", "allosteric"),
    "0H4": ("GCK allosteric activator (0H4)", "allosteric"),
    "HEM": ("heme (protoporphyrin IX + Fe)", "cofactor"),
    "FLP": ("flurbiprofen (substrate/drug)", "substrate"),
    "TLA": ("L(+)-tartrate (substrate-analog)", "substrate"),
    "ANP": ("AMP-PNP (non-hydrolyzable ATP analog)", "nucleotide"),
    "MPZ": ("ATP-competitive inhibitor (purine)", "inhibitor"),
    "UA7": ("brodifacoum (anticoagulant inhibitor)", "inhibitor"),
    "GNP": ("GppNHp / GMPPNP (non-hydrolyzable GTP analog)", "nucleotide"),
    "GDP": ("GDP", "nucleotide"),
    "PTR": ("O-phosphotyrosine (regulatory PTM)", "modification"),
    "CRO": ("GFP-fusion chromophore", "modification"),
    "MG":  ("magnesium ion", "metal"),
    "NA":  ("sodium ion", "ion/additive"),
    "CL":  ("chloride ion", "ion/additive"),
    "SO4": ("sulfate", "ion/additive"),
    "PO4": ("phosphate", "ion/additive"),
    "GOL": ("glycerol (cryoprotectant)", "ion/additive"),
    "MPD": ("2-methyl-2,4-pentanediol (cryo)", "ion/additive"),
    "EDO": ("ethylene glycol (cryo)", "ion/additive"),
    "OLC": ("monoolein (LCP lipid)", "lipid"),
}
FUNCTIONAL_ROLES = {"substrate", "cofactor", "nucleotide", "allosteric", "inhibitor"}

# --- structures to inventory: (label, protein, relpath, target_chain, partners) ---
# partners: {chain: partner_name} for bound protein co-complexes
STRUCTS = [
    ("1V4T", "HXK4 / GCK",  "1V4T.pdb", "A", {}),
    ("1V4S", "HXK4 / GCK",  "1V4S.pdb", "A", {}),
    ("3VEV*", "HXK4 / GCK", "3VEV.pdb", "A", {}),   # * not in curated set
    ("1R9O", "CP2C9",       "1R9O.pdb", "A", {}),
    ("1D5R", "PTEN",        "1D5R.pdb", "A", {}),
    ("2SRC", "SRC",         "2SRC.pdb", "A", {}),
    ("1Y57", "SRC",         "1Y57.pdb", "A", {}),
    ("6WVH", "VKOR1",       "6WVH.pdb", "A", {}),
    ("6VJJ", "RASK / KRAS", "KRAS/6VJJ.pdb", "A", {"B": "RAF1 RBD (effector)"}),
    ("5O2S", "RASK / KRAS", "KRAS/5O2S.pdb", "A", {"B": "DARPin K27 (GDP-selective binder)"}),
    ("1NVW", "RASK / H-Ras","KRAS/1NVW.pdb", "R", {"S": "SOS1 (GEF)"}),
]


# curated per-structure notes (biological caveats not derivable from HETATM alone)
NOTES = {
    "1V4T": "apo super-open (inactive) — active site empty; sulfate/Na only",
    "1NVW": "mixed-state: GNP+Mg on chain Q (allosteric Ras); target SOS-substrate "
            "chain R is nucleotide-free (PO4 only)",
    "6WVH": "CRO chromophore + OLC lipid belong to the GFP crystallization fusion, "
            "not VKOR; VKORC1 = author residues 145-297",
    "5O2S": "4 KRAS copies in the AU, each GDP+Mg-bound; chain A analysed",
    "2SRC": "PTR = autoinhibitory phospho-Tyr (auth 527) modeled as a modified residue",
}


def parse_hetatm(path):
    """Return list of (resname, chain, resseq) for non-water HETATM groups (unique)."""
    seen = set()
    out = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("HETATM"):
                res = ln[17:20].strip()
                if res == "HOH":
                    continue
                chain = ln[21].strip()
                resseq = ln[22:27].strip()
                key = (res, chain, resseq)
                if key not in seen:
                    seen.add(key)
                    out.append((res, chain, resseq))
    return out


def fmt(counts):
    """{resname: n} -> 'RES xN; ...' sorted by descending count."""
    parts = []
    for res, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        parts.append(res + (f" x{n}" if n > 1 else ""))
    return "; ".join(parts) or "—"


def main():
    long_rows = []
    wide_rows = []
    for label, prot, rel, tchain, partners in STRUCTS:
        path = os.path.join(PDB_DIR, rel)
        groups = parse_hetatm(path)

        # count per (role, on_target) ; ions/additives counted file-wide
        target_roles = defaultdict(lambda: defaultdict(int))   # role -> {res: n} on target chain
        offtarget_func = defaultdict(int)                      # res -> n  (functional, other chains)
        additives = defaultdict(int)                           # res -> n  (ion/additive/lipid, any chain)
        for res, chain, _ in groups:
            name, role = HET.get(res, (res, "unknown"))
            on_target = (chain == tchain)
            long_rows.append(OrderedDict(
                structure=label, protein=prot, het_code=res, name=name,
                role=role, chain=chain, on_target_chain=on_target))
            if role in ("ion/additive", "lipid"):
                additives[res] += 1
            elif role in FUNCTIONAL_ROLES or role == "modification" or role == "metal":
                if on_target:
                    target_roles[role][res] += 1
                elif role in FUNCTIONAL_ROLES:
                    offtarget_func[res] += 1

        has_partner = bool(partners)
        functional = any(role in target_roles for role in FUNCTIONAL_ROLES) or has_partner

        def col(role):
            return fmt(target_roles.get(role, {}))

        partner_str = "; ".join(f"{v} (chain {k})" for k, v in partners.items()) or "—"

        note = NOTES.get(label, "")
        if offtarget_func:
            note = (note + "; " if note else "") + \
                   "off-target-chain: " + fmt(offtarget_func)

        wide_rows.append(OrderedDict(
            structure=label,
            protein=prot,
            functional_ligand_bound="yes" if functional else "NO (apo/additives only)",
            substrate=col("substrate"),
            cofactor=col("cofactor"),
            nucleotide=col("nucleotide"),
            allosteric=col("allosteric"),
            inhibitor=col("inhibitor"),
            metal_Mg=fmt(target_roles.get("metal", {})),
            bound_partner=partner_str,
            modification=col("modification"),
            ions_additives=fmt(additives),
            note=note or "—",
        ))

    long_df = pd.DataFrame(long_rows)
    wide_df = pd.DataFrame(wide_rows)
    long_df.to_csv(os.path.join(OUTDIR, "ligands_by_structure_long.csv"), index=False)
    wide_df.to_csv(os.path.join(OUTDIR, "ligands_by_structure.csv"), index=False)

    # markdown
    md = ["# Ligand / substrate inventory of the curated PDB structures",
          "",
          "Read directly from the `HETATM` records of each structure "
          "(`scripts/48_structure_ligands.py`). "
          "`functional_ligand_bound` = yes if a substrate, cofactor, nucleotide, "
          "allosteric modulator, inhibitor, **or** a bound protein partner is present; "
          "structures carrying only ions / cryoprotectants are flagged apo.",
          "",
          "`*` 3VEV is present in the directory but **not** part of the curated GCK selection.",
          ""]
    cols = ["structure", "protein", "functional_ligand_bound", "substrate",
            "cofactor", "nucleotide", "allosteric", "inhibitor", "metal_Mg",
            "bound_partner", "modification", "ions_additives", "note"]
    md.append("| " + " | ".join(cols) + " |")
    md.append("|" + "|".join(["---"] * len(cols)) + "|")
    for _, r in wide_df.iterrows():
        md.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    md.append("")
    with open(os.path.join(OUTDIR, "ligands_by_structure.md"), "w") as fh:
        fh.write("\n".join(md))

    print(wide_df.to_string(index=False))
    print(f"\nWrote:\n  {OUTDIR}/ligands_by_structure.csv"
          f"\n  {OUTDIR}/ligands_by_structure_long.csv"
          f"\n  {OUTDIR}/ligands_by_structure.md")


if __name__ == "__main__":
    main()
