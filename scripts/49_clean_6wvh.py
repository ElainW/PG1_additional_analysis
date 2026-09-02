"""
Process 6WVH -> 6WVH_cleaned.pdb : strip the crystallization-fusion artifacts,
leaving VKORC1 + its functional ligand only.

6WVH is a fusion construct (see data/structures/pdb/README.md):
  author residues   2-144   = Aequorea GFP fusion, N-part   (REMOVE)
  author residues 145-297   = VKORC1_HUMAN (UniProt 3-155)  (KEEP)
  author residues 298-392   = Aequorea GFP fusion, C-part   (REMOVE)
  HETATM CRO (res 66)       = GFP chromophore               (REMOVE)
  HETATM OLC (res 402,403)  = monoolein / LCP lipid         (REMOVE)
  HETATM UA7 (res 401)      = brodifacoum (VKOR inhibitor)  (KEEP)
  HOH                       = waters                        (REMOVE)

Author numbering is preserved (UniProt = author - 142). Atom serials are
renumbered sequentially and a TER is written after each VKORC1 chain.

Two processed files are written (both -> stay in pdb/):
  6WVH_cleaned.pdb          both VKORC1 copies (chains A + B) + brodifacoum
  6WVH_cleaned_chainA.pdb   single VKORC1 copy (chain A only) + brodifacoum,
                            for per-residue SASA (the A/B crystal contact would
                            otherwise be scored as buried surface)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

PDB_DIR = os.path.join(C.DATA_DIR, "structures", "pdb")
SRC = os.path.join(PDB_DIR, "6WVH.pdb")
OUT = os.path.join(PDB_DIR, "6WVH_cleaned.pdb")
OUT_A = os.path.join(PDB_DIR, "6WVH_cleaned_chainA.pdb")

VKOR_LO, VKOR_HI = 145, 297          # author residue range of VKORC1
KEEP_HET = {"UA7"}                    # functional ligand to retain
CHAINS = ("A", "B")


def keep_line(ln):
    """True if this ATOM/HETATM belongs to VKORC1 or the retained ligand."""
    rec = ln[:6].strip()
    if rec not in ("ATOM", "HETATM"):
        return False
    chain = ln[21]
    if chain not in CHAINS:
        return False
    res = ln[17:20].strip()
    try:
        n = int(ln[22:26])
    except ValueError:
        return False
    if rec == "ATOM":
        return VKOR_LO <= n <= VKOR_HI          # VKORC1 protein only
    return res in KEEP_HET                        # UA7; drops CRO, OLC, HOH


def ter_record(serial, last_atom_line):
    """Build a TER line from the last ATOM line of a chain."""
    res = last_atom_line[17:20]
    chain = last_atom_line[21]
    resseq = last_atom_line[22:26]
    icode = last_atom_line[26]
    return f"TER   {serial:>5}      {res} {chain}{resseq}{icode}\n"


def write_pdb(path, kept, preserved, chains):
    """Emit a cleaned PDB for the given chain list (serials renumbered)."""
    chlabel = ",".join(chains)
    out = [
        "REMARK   1 6WVH_cleaned : VKORC1 (VKOR1_HUMAN, Q9BQB6) extracted from the\n",
        "REMARK   1 6WVH GFP-fusion crystal. Removed: GFP fusion protein (author\n",
        "REMARK   1 res 2-144 & 298-392), GFP chromophore CRO, LCP lipid OLC, waters.\n",
        f"REMARK   1 Kept: VKORC1 author res 145-297 (chain(s) {chlabel}) + brodifacoum UA7.\n",
        "REMARK   1 Author numbering preserved; UniProt position = author - 142.\n",
        "REMARK   1 Produced by scripts/49_clean_6wvh.py from 6WVH.pdb.\n",
    ]
    out.extend(preserved)
    serial = 0
    for chain in chains:
        atoms = kept[chain]["atoms"]
        hets = kept[chain]["het"]
        if not atoms and not hets:
            continue
        for ln in atoms:
            serial += 1
            out.append(ln[:6] + f"{serial:>5}" + ln[11:])
        if atoms:
            serial += 1
            out.append(ter_record(serial, atoms[-1]))
        for ln in hets:
            serial += 1
            out.append(ln[:6] + f"{serial:>5}" + ln[11:])
    out.append("END\n")
    with open(path, "w") as fh:
        fh.writelines(out)
    return serial


def main():
    # preserve unit-cell / coordinate-frame records so the file stays tool-valid
    preserved = []
    kept = {c: {"atoms": [], "het": []} for c in CHAINS}
    with open(SRC) as fh:
        for ln in fh:
            rec = ln[:6].strip()
            if rec in ("CRYST1", "SCALE1", "SCALE2", "SCALE3",
                       "ORIGX1", "ORIGX2", "ORIGX3"):
                preserved.append(ln)
                continue
            if not keep_line(ln):
                continue
            chain = ln[21]
            if rec == "ATOM":
                kept[chain]["atoms"].append(ln)
            else:
                kept[chain]["het"].append(ln)

    for chain in CHAINS:
        resA = sorted({int(l[22:26]) for l in kept[chain]["atoms"]})
        hets = sorted({l[17:20].strip() for l in kept[chain]["het"]})
        span = f"{resA[0]}-{resA[-1]}" if resA else "none"
        print(f"chain {chain}: VKORC1 residues {span} ({len(resA)} res), ligand {hets or '—'}")

    n_both = write_pdb(OUT, kept, preserved, CHAINS)
    n_a = write_pdb(OUT_A, kept, preserved, ("A",))
    print(f"Wrote {OUT}         (chains A+B, {n_both} records)")
    print(f"Wrote {OUT_A}  (chain A only, {n_a} records)")


if __name__ == "__main__":
    main()
