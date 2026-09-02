"""
Relative SASA (rSASA) for the cleaned structures.

Runs the FreeSASA CLI with the specified parameters on each single-chain,
protein-only cleaned structure:

    freesasa <pdb> -n 20 --format=rsa --radii=naccess

FreeSASA reads ATOM records only (functional HETATM ligands in the cleaned files
are ignored) -> protein-only surface, as instructed. For KRAS the partner-removed
single-chain per-state files (KRAS/*_*P.pdb) are used; for VKOR the single-copy
6WVH_cleaned_chainA.pdb; for the rest the cleaned chain-A files.

Outputs (results/structure_features/rsasa/):
  <label>.rsa            raw NACCESS-format FreeSASA output (per structure)
  rsasa_all.csv          parsed long table, all structures, per residue
                         (REL = relative SASA %, the burial metric)
"""
import os
import re
import sys
import shutil
import subprocess

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

PDB_DIR = os.path.join(C.DATA_DIR, "structures", "pdb")
OUTDIR = os.path.join(C.RESULTS, "structure_features", "rsasa")
os.makedirs(OUTDIR, exist_ok=True)

FREESASA = shutil.which("freesasa") or os.path.join(
    os.environ.get("CONDA_PREFIX", ""), "bin", "freesasa")

# (label, protein, uniprot, state, relative_path)
STRUCTS = [
    ("1V4T", "HXK4_HUMAN", "P35557", "super-open (inactive)", "1V4T_cleaned.pdb"),
    ("1V4S", "HXK4_HUMAN", "P35557", "closed (active)",       "1V4S_cleaned.pdb"),
    ("1R9O", "CP2C9_HUMAN", "P11712", "single",               "1R9O_cleaned.pdb"),
    ("1D5R", "PTEN_HUMAN", "P60484", "single",                "1D5R_cleaned.pdb"),
    ("2SRC", "SRC_HUMAN", "P12931", "inactive",               "2SRC_cleaned.pdb"),
    ("1Y57", "SRC_HUMAN", "P12931", "active",                 "1Y57_cleaned.pdb"),
    ("6WVH", "VKOR1_HUMAN", "Q9BQB6", "single",               "6WVH_cleaned_chainA.pdb"),
    ("6VJJ", "RASK_HUMAN", "P01116", "ON/3P (active)",        "KRAS/6VJJ_ON_3P.pdb"),
    ("5O2S", "RASK_HUMAN", "P01116", "OFF/2P (inactive)",     "KRAS/5O2S_OFF_2P.pdb"),
    ("1NVW", "RASK_HUMAN", "P01116", "apo/0P",                "KRAS/1NVW_APO_0P.pdb"),
]

# RES line columns of a FreeSASA .rsa (NACCESS layout):
#   RES <res> <chain> <num>  All-abs All-rel  Side-abs Side-rel  Main-abs Main-rel
#   Nonpolar-abs Nonpolar-rel  Polar-abs Polar-rel
RES_RE = re.compile(
    r"^RES\s+(?P<res>\S+)\s+(?P<chain>\S)\s*(?P<num>[-0-9]+)\s+"
    r"(?P<all_abs>[-0-9.]+)\s+(?P<all_rel>[-0-9.]+)\s+"
    r"(?P<side_abs>[-0-9.]+)\s+(?P<side_rel>[-0-9.]+)\s+"
    r"(?P<main_abs>[-0-9.]+)\s+(?P<main_rel>[-0-9.]+)\s+"
    r"(?P<np_abs>[-0-9.]+)\s+(?P<np_rel>[-0-9.]+)\s+"
    r"(?P<p_abs>[-0-9.]+)\s+(?P<p_rel>[-0-9.]+)")


def run_freesasa(path, out_rsa):
    with open(out_rsa, "w") as fh:
        subprocess.run([FREESASA, path, "-n", "20", "--format=rsa",
                        "--radii=naccess"], stdout=fh, check=True)


def parse_rsa(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            m = RES_RE.match(ln)
            if not m:
                continue
            d = m.groupdict()
            rows.append(dict(
                resname=d["res"], chain=d["chain"], resnum=int(d["num"]),
                rel_all=float(d["all_rel"]), rel_side=float(d["side_rel"]),
                rel_main=float(d["main_rel"]),
                abs_all=float(d["all_abs"]), abs_side=float(d["side_abs"]),
            ))
    return rows


def main():
    if not os.path.exists(FREESASA):
        sys.exit(f"freesasa CLI not found at {FREESASA}; activate env 'structfeat'.")
    allrows = []
    for label, prot, uni, state, rel in STRUCTS:
        path = os.path.join(PDB_DIR, rel)
        out_rsa = os.path.join(OUTDIR, f"{label}.rsa")
        run_freesasa(path, out_rsa)
        rows = parse_rsa(out_rsa)
        for r in rows:
            r.update(structure=label, protein=prot, uniprot=uni, state=state)
            allrows.append(r)
        print(f"{label:5s} {prot:12s} {state:20s} -> {len(rows)} residues  ({os.path.basename(out_rsa)})")

    df = pd.DataFrame(allrows, columns=[
        "structure", "protein", "uniprot", "state", "chain", "resnum",
        "resname", "rel_all", "rel_side", "rel_main", "abs_all", "abs_side"])
    df.to_csv(os.path.join(OUTDIR, "rsasa_all.csv"), index=False)
    print(f"\nWrote {OUTDIR}/rsasa_all.csv  ({len(df)} residue rows across {df.structure.nunique()} structures)")
    # quick buried-fraction sanity (REL<25% commonly 'buried')
    for label, g in df.groupby("structure", sort=False):
        buried = (g.rel_all < 25).mean() * 100
        print(f"  {label}: {buried:.0f}% residues buried (rel_all<25%)")


if __name__ == "__main__":
    main()
