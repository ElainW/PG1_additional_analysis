"""
Static residue-residue contacts (GetContacts) for the cleaned structures.

Pipeline per structure:
  1. add hydrogens with PyMOL h_add  (required by GetContacts H-bond geometry)
  2. get_static_contacts.py --itypes all --sele "protein or resname PTR"
     (protein residues + the phospho-Tyr modified residue in 2SRC; other HETATM
      ligands are excluded -> residue-residue contacts only)

Outputs (results/structure_features/contacts/):
  hadd/<label>_H.pdb            hydrogen-added coordinates
  <label>_contacts.tsv         raw GetContacts output, all interaction types   [ALL structures]
  contacts_itype_summary.csv   per-structure count of each interaction type    [ALL structures]
  4types/<label>_4types.tsv    salt-bridge / pi-cation / sc-sc H-bond /
                               sc-bb H-bond contacts (these carry conformational
                               state specificity)                              [MULTI-STATE proteins]

Multi-state proteins (HXK4, RASK, SRC) get the 4-type tables that feed the
active/inactive residue classification (script 54).
"""
import os
import sys
import subprocess
from collections import Counter

import pandas as pd
import pymol2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

PDB_DIR = os.path.join(C.DATA_DIR, "structures", "pdb")
OUTDIR = os.path.join(C.RESULTS, "structure_features", "contacts")
HADD = os.path.join(OUTDIR, "hadd")
FOURD = os.path.join(OUTDIR, "4types")
for d in (OUTDIR, HADD, FOURD):
    os.makedirs(d, exist_ok=True)

GETCONTACTS = os.path.expanduser("~/tools/getcontacts/get_static_contacts.py")
SELE = "protein or resname PTR"
FOUR_TYPES = ["sb", "pc", "hbss", "hbsb"]   # salt bridge, pi-cation, sc-sc Hb, sc-bb Hb

# (label, protein, uniprot, state, relpath, multistate)
STRUCTS = [
    ("1V4T", "HXK4_HUMAN", "P35557", "super-open (inactive)", "1V4T_cleaned.pdb", True),
    ("1V4S", "HXK4_HUMAN", "P35557", "closed (active)",       "1V4S_cleaned.pdb", True),
    ("1R9O", "CP2C9_HUMAN", "P11712", "single",               "1R9O_cleaned.pdb", False),
    ("1D5R", "PTEN_HUMAN", "P60484", "single",                "1D5R_cleaned.pdb", False),
    ("2SRC", "SRC_HUMAN", "P12931", "inactive",               "2SRC_cleaned.pdb", True),
    ("1Y57", "SRC_HUMAN", "P12931", "active",                 "1Y57_cleaned.pdb", True),
    ("6WVH", "VKOR1_HUMAN", "Q9BQB6", "single",               "6WVH_cleaned_chainA.pdb", False),
    ("6VJJ", "RASK_HUMAN", "P01116", "ON/3P (active)",        "KRAS/6VJJ_ON_3P.pdb", True),
    ("5O2S", "RASK_HUMAN", "P01116", "OFF/2P (inactive)",     "KRAS/5O2S_OFF_2P.pdb", True),
    ("1NVW", "RASK_HUMAN", "P01116", "apo/0P",                "KRAS/1NVW_APO_0P.pdb", True),
]


def h_add(inp, out):
    with pymol2.PyMOL() as p:
        p.cmd.load(inp, "m")
        p.cmd.remove("solvent")
        p.cmd.h_add("m")
        p.cmd.save(out, "m")


def run_getcontacts(struct, out_tsv):
    subprocess.run([sys.executable, GETCONTACTS, "--structure", struct,
                    "--output", out_tsv, "--itypes", "all", "--sele", SELE],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_contacts(tsv):
    """Yield dicts for each contact line (first two atoms)."""
    rows = []
    with open(tsv) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            itype, a1, a2 = f[1], f[2], f[3]
            c1, r1, n1, at1 = a1.split(":")
            c2, r2, n2, at2 = a2.split(":")
            rows.append(dict(itype=itype,
                             chain1=c1, res1=r1, resnum1=int(n1), atom1=at1,
                             chain2=c2, res2=r2, resnum2=int(n2), atom2=at2))
    return rows


def main():
    summ = []
    for label, prot, uni, state, rel, multistate in STRUCTS:
        src = os.path.join(PDB_DIR, rel)
        hpdb = os.path.join(HADD, f"{label}_H.pdb")
        tsv = os.path.join(OUTDIR, f"{label}_contacts.tsv")
        h_add(src, hpdb)
        run_getcontacts(hpdb, tsv)
        rows = parse_contacts(tsv)
        counts = Counter(r["itype"] for r in rows)
        summ.append(dict(structure=label, protein=prot, state=state,
                         total=len(rows), **{t: counts.get(t, 0)
                         for t in ["sb", "pc", "ps", "ts", "hbss", "hbsb",
                                   "hbbb", "hp", "vdw"]}))

        if multistate:
            four = [r for r in rows if r["itype"] in FOUR_TYPES]
            df4 = pd.DataFrame(four)
            df4.to_csv(os.path.join(FOURD, f"{label}_4types.tsv"),
                       sep="\t", index=False)
        n4 = sum(counts.get(t, 0) for t in FOUR_TYPES)
        print(f"{label:5s} {prot:12s} {state:22s} total={len(rows):5d}  "
              f"sb={counts.get('sb',0):3d} pc={counts.get('pc',0):3d} "
              f"hbss={counts.get('hbss',0):3d} hbsb={counts.get('hbsb',0):3d}"
              f"{'   [4types written]' if multistate else ''}")

    pd.DataFrame(summ).to_csv(
        os.path.join(OUTDIR, "contacts_itype_summary.csv"), index=False)
    print(f"\nWrote raw *_contacts.tsv (all), contacts_itype_summary.csv, and "
          f"4types/*_4types.tsv (multi-state) to\n  {OUTDIR}")


if __name__ == "__main__":
    main()
