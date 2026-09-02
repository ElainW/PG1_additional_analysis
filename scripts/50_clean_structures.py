"""
Clean the remaining curated structures for the structure-feature analyses:
remove solvent, strip crystallization additives, and extract the single target
chain, KEEPING functional ligands (substrate / cofactor / nucleotide / allosteric
modulator / inhibitor / functional metal) and inline modified residues.

6WVH (scripts/49) and the KRAS set (KRAS/*_*P.pdb) are already done; this handles
the six single-chain structures:

  1V4T  GCK super-open (apo)   keep: —                 strip: SO4, NA, waters
  1V4S  GCK closed             keep: GLC, MRK          strip: NA, waters
  1R9O  CYP2C9                 keep: HEM, FLP          strip: GOL, waters
  1D5R  PTEN                   keep: TLA               strip: waters
  2SRC  Src inactive           keep: ANP, PTR*         strip: waters
  1Y57  Src active             keep: MPZ               strip: SO4, waters
  * PTR = phospho-Tyr527, an inline modified residue (kept in place).

All six are already chain A only; we still filter to chain A defensively.
Author numbering is preserved. Atom serials are renumbered sequentially.

Output: data/structures/pdb/<ID>_cleaned.pdb  (processed structures -> stay in pdb/)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

PDB_DIR = os.path.join(C.DATA_DIR, "structures", "pdb")

# HET code -> role. KEEP_ROLES are functional and retained; everything else
# (ion/additive, lipid) plus HOH is stripped.
ROLE = {
    "GLC": "substrate", "FLP": "substrate", "TLA": "substrate",
    "HEM": "cofactor",
    "ANP": "nucleotide",
    "MRK": "allosteric", "0H4": "allosteric",
    "MPZ": "inhibitor", "UA7": "inhibitor",
    "PTR": "modification", "CRO": "modification",
    "MG": "metal",
    # additives / solvent (stripped)
    "NA": "ion/additive", "CL": "ion/additive", "SO4": "ion/additive",
    "PO4": "ion/additive", "GOL": "ion/additive", "MPD": "ion/additive",
    "EDO": "ion/additive", "OLC": "lipid",
}
KEEP_ROLES = {"substrate", "cofactor", "nucleotide", "allosteric",
              "inhibitor", "metal", "modification"}

# (pdb_id, target_chain)
STRUCTS = [
    ("1V4T", "A"), ("1V4S", "A"), ("1R9O", "A"),
    ("1D5R", "A"), ("2SRC", "A"), ("1Y57", "A"),
]


def clean(pdb_id, tchain):
    src = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    out = os.path.join(PDB_DIR, f"{pdb_id}_cleaned.pdb")

    preserved, coord_lines = [], []
    kept_het = {}          # resname -> count (for reporting)
    stripped_het = {}
    with open(src) as fh:
        for ln in fh:
            rec = ln[:6].strip()
            if rec in ("CRYST1", "SCALE1", "SCALE2", "SCALE3",
                       "ORIGX1", "ORIGX2", "ORIGX3"):
                preserved.append(ln)
            elif rec == "ATOM":
                if ln[21] == tchain:
                    coord_lines.append(ln)
            elif rec == "HETATM":
                res = ln[17:20].strip()
                if res == "HOH":
                    continue
                role = ROLE.get(res, "unknown")
                if ln[21] == tchain and role in KEEP_ROLES:
                    coord_lines.append(ln)
                    kept_het[res] = kept_het.get(res, 0) + 1
                else:
                    stripped_het[res] = stripped_het.get(res, 0) + 1
            elif rec == "TER":
                if ln[21] in (tchain, " "):
                    coord_lines.append(ln)

    # renumber serials sequentially
    out_lines = [
        f"REMARK   1 {pdb_id}_cleaned : chain {tchain} only; waters + crystallization\n",
        "REMARK   1 additives stripped; functional ligands / modified residues kept.\n",
        "REMARK   1 Author numbering preserved. scripts/50_clean_structures.py from "
        f"{pdb_id}.pdb.\n",
    ]
    out_lines.extend(preserved)
    serial = 0
    for ln in coord_lines:
        serial += 1
        out_lines.append(ln[:6] + f"{serial:>5}" + ln[11:])
    out_lines.append("END\n")
    with open(out, "w") as fh:
        fh.writelines(out_lines)

    def fmt(d):
        return ", ".join(f"{k}" for k in sorted(d)) or "—"
    print(f"{pdb_id}: kept ligands [{fmt(kept_het)}]  stripped [{fmt(stripped_het)}] + waters "
          f"-> {os.path.basename(out)}")


def main():
    for pdb_id, ch in STRUCTS:
        clean(pdb_id, ch)


if __name__ == "__main__":
    main()
