"""
3-D redo of the "disagreement elevated near binding/active-site residues"
generalization test (15 does this in 1-D via UniProt Binding/Active site
annotations +-2; pooled 10/13 proteins, p=5e-13).

Here we replace the 1-D annotation window with a real 3-D heavy-atom distance
to the protein's bound functional ligand(s), read from the curated experimental
structures in data/structures/pdb/. Only proteins with a curated structure
carrying a genuine functional ligand (substrate, cofactor, nucleotide,
allosteric modulator, inhibitor -- per
results/structure_features/ligands_by_structure.md) are testable; that is 6 of
the 23 multi-assay proteins (vs 13 testable for the 1-D annotation test):
  HXK4_HUMAN (1V4S: GLC, MRK), CP2C9_HUMAN (1R9O: FLP, HEM),
  PTEN_HUMAN (1D5R: TLA), SRC_HUMAN (2SRC: ANP),
  VKOR1_HUMAN (6WVH chain A: UA7), RASK_HUMAN (KRAS 6VJJ ON state: GNP, MG).
Excluded as non-functional per the ligand inventory: ions/cryoprotectants
(NA, GOL), the VKOR1 GFP-fusion chromophore/lipid (CRO, OLC), and SRC's PTR
(a modified residue, not a small-molecule ligand).

UniProt<->PDB-author numbering offset, per protein, is determined EMPIRICALLY
at runtime: for every candidate offset in [-150, 150] we count how often the
PDB residue identity at (auth + offset) equals the wild-type amino acid the DMS
assay itself encodes at that position (parsed from the assay's own `mutant`
column, the most authoritative ground truth -- these are literally the residues
disagreement is being computed on), and take the best-matching offset,
requiring >=95% agreement over >=20 overlapping positions. This catches
numbering-convention drift instead of trusting any hardcoded/documented
value. It disagrees with the top-level PDB README for two of the six proteins
(HXK4/1V4S: README says auth+1, actual assay-verified offset is +0; SRC/2SRC:
README says auth+2, actual is +3, matching 1Y57's documented offset) --
printed as a warning, not silently overridden.

For each testable protein we compute, per residue (UniProt numbering), the
minimum heavy-atom distance to the nearest functional-ligand heavy atom, then
join to the per-site disagreement table (15's
generalization_functional_disagreement_sites.csv, z_disagree = within-assay
z-scored disagreement). Two tests, run per protein and pooled, mirroring 15's
per_protein()/pooled structure exactly so the 1-D vs 3-D result is directly
comparable:
  (a) continuous: Spearman(z_disagree, dist_to_ligand)  [H1: negative -- closer
      to ligand -> more disagreement]
  (b) binary "near-ligand" (dist <= threshold, tested at 5/8/10 A) vs rest:
      one-sided Mann-Whitney, H1 near > rest z_disagree (same direction/test as
      15's near_binding vs rest).

Outputs (results/functional_disagreement_3d/):
  <protein>_residue_distances.csv   per-residue dist_to_ligand (all resolved AAs)
  site_distances.csv                per-site (protein,DMS_id,up_pos) joined to
                                     dist_to_ligand + disagreement + z_disagree
  per_protein_tests.csv             per-protein Spearman + MW @ 3 thresholds
  pooled_tests.csv                  pooled Spearman + MW @ 3 thresholds
"""
import os
import re
import glob
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one
import common as C

warnings.simplefilter("ignore")
PDBDIR = os.path.join(C.ANALYSIS, "data", "structures", "pdb")
OUTDIR = os.path.join(C.RESULTS, "functional_disagreement_3d")
THRESHOLDS = [5.0, 8.0, 10.0]
MIN_MATCH_FRAC = 0.95
MIN_MATCH_N = 20

# protein -> (pdb path, chain, ligand resnames, documented offset [only a
# starting hint / cross-check -- the actual offset used is verified below])
PROTEINS = {
    "HXK4_HUMAN": dict(
        pdb=os.path.join(PDBDIR, "1V4S.pdb"), chain="A", ligands=["GLC", "MRK"],
        documented_offset=1),
    "CP2C9_HUMAN": dict(
        pdb=os.path.join(PDBDIR, "1R9O.pdb"), chain="A", ligands=["FLP", "HEM"],
        documented_offset=0),
    "PTEN_HUMAN": dict(
        pdb=os.path.join(PDBDIR, "1D5R.pdb"), chain="A", ligands=["TLA"],
        documented_offset=0),
    "SRC_HUMAN": dict(
        pdb=os.path.join(PDBDIR, "2SRC.pdb"), chain="A", ligands=["ANP"],
        documented_offset=2),
    "VKOR1_HUMAN": dict(
        pdb=os.path.join(PDBDIR, "6WVH_cleaned_chainA.pdb"), chain="A", ligands=["UA7"],
        documented_offset=-142),
    "RASK_HUMAN": dict(
        pdb=os.path.join(PDBDIR, "KRAS", "6VJJ_ON_3P.pdb"), chain="A", ligands=["GNP", "MG"],
        documented_offset=0),
}


def one(resname):
    try:
        return index_to_one(three_to_index(resname))
    except Exception:
        return "?"


def dms_wt_map(prot):
    """pos -> wild-type amino acid, parsed from the DMS assay's own `mutant`
    column (X123Y). Ground truth for what residue disagreement is scored on."""
    files = sorted(glob.glob(os.path.join(C.SCORES_DIR, f"{prot}_*.csv")))
    d = pd.read_csv(files[0], usecols=["mutant"])
    wt = {}
    for mu in d["mutant"].astype(str):
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", mu)
        if m:
            wt[int(m.group(2))] = m.group(1)
    return wt


def verify_offset(prot, pdb_resname, documented_offset):
    """Scan offsets, score by agreement with the DMS assay's own WT residues."""
    wt = dms_wt_map(prot)
    best = None
    for off in range(-150, 151):
        matches = sum(1 for auth, aa in pdb_resname.items() if wt.get(auth + off) == aa)
        total = sum(1 for auth in pdb_resname if (auth + off) in wt)
        if total >= MIN_MATCH_N:
            frac = matches / total
            if best is None or frac > best[1]:
                best = (off, frac, matches, total)
    assert best is not None and best[1] >= MIN_MATCH_FRAC, (
        f"{prot}: could not verify a numbering offset (best={best})")
    off, frac, matches, total = best
    flag = "" if off == documented_offset else "  ** differs from PDB README **"
    print(f"{prot:14s} offset verified: {off:+4d}  ({matches}/{total} = {frac:.1%} match; "
          f"documented={documented_offset:+d}){flag}")
    return off


def load_and_distances(prot, spec):
    struct = PDBParser(QUIET=True).get_structure(prot, spec["pdb"])[0][spec["chain"]]
    heavy, resname = {}, {}
    ligand_atoms = []
    for r in struct:
        het = r.id[0]
        coords = np.array([a.coord for a in r if a.element != "H"])
        if het == " ":
            n = r.id[1]
            resname[n] = one(r.resname)
            heavy[n] = coords
        elif r.resname in spec["ligands"]:
            ligand_atoms.append(coords)
    ligand_atoms = np.vstack(ligand_atoms) if ligand_atoms else np.empty((0, 3))
    assert len(ligand_atoms), f"{prot}: no ligand atoms found for {spec['ligands']}"

    offset = verify_offset(prot, resname, spec["documented_offset"])

    rows = []
    for n, coords in heavy.items():
        d = float(cdist(coords, ligand_atoms).min())
        rows.append({"auth_pos": n, "up_pos": n + offset,
                     "pdb_aa": resname[n], "dist_to_ligand": d})
    return pd.DataFrame(rows), offset


def per_protein_test(df):
    d = df.dropna(subset=["z_disagree", "dist_to_ligand"])
    rho, p_rho = spearmanr(d["z_disagree"], d["dist_to_ligand"])
    out = {"n_sites": len(d), "spearman_zdisagree_vs_dist": rho, "spearman_p": p_rho}
    for thr in THRESHOLDS:
        near = d[d["dist_to_ligand"] <= thr]["z_disagree"]
        far = d[d["dist_to_ligand"] > thr]["z_disagree"]
        if len(near) >= 5 and len(far) >= 5:
            p = mannwhitneyu(near, far, alternative="greater").pvalue
            out[f"n_near_{thr:g}A"] = len(near)
            out[f"n_far_{thr:g}A"] = len(far)
            out[f"mean_z_near_{thr:g}A"] = near.mean()
            out[f"mean_z_far_{thr:g}A"] = far.mean()
            out[f"delta_z_{thr:g}A"] = near.mean() - far.mean()
            out[f"mw_p_near_gt_far_{thr:g}A"] = p
        else:
            for k in (f"n_near_{thr:g}A", f"n_far_{thr:g}A", f"mean_z_near_{thr:g}A",
                      f"mean_z_far_{thr:g}A", f"delta_z_{thr:g}A", f"mw_p_near_gt_far_{thr:g}A"):
                out[k] = np.nan
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    sites = pd.read_csv(C.rp("generalization_functional_disagreement_sites.csv"))

    all_site_rows = []
    per_protein_rows = []
    for prot, spec in PROTEINS.items():
        dist, offset = load_and_distances(prot, spec)
        dist.to_csv(os.path.join(OUTDIR, f"{prot}_residue_distances.csv"), index=False)

        g = sites[sites["UniProt_ID"] == prot]
        m = g.merge(dist[["up_pos", "dist_to_ligand", "pdb_aa"]], on="up_pos", how="inner")
        all_site_rows.append(m)

        stats = per_protein_test(m)
        stats["UniProt_ID"] = prot
        stats["n_assays"] = g["DMS_id"].nunique()
        stats["ligands"] = ";".join(spec["ligands"])
        stats["offset_used"] = offset
        per_protein_rows.append(stats)
        print(f"  -> n_sites={stats['n_sites']:4d}  "
              f"rho(z,dist)={stats['spearman_zdisagree_vs_dist']:+.3f} "
              f"(p={stats['spearman_p']:.3g})  "
              f"8A: delta_z={stats.get('delta_z_8A', float('nan')):+.3f} "
              f"p={stats.get('mw_p_near_gt_far_8A', float('nan')):.3g}\n")

    all_sites = pd.concat(all_site_rows, ignore_index=True)
    all_sites.to_csv(os.path.join(OUTDIR, "site_distances.csv"), index=False)

    per_prot_df = pd.DataFrame(per_protein_rows)
    cols = ["UniProt_ID", "n_assays", "ligands", "offset_used", "n_sites",
            "spearman_zdisagree_vs_dist", "spearman_p"]
    for thr in THRESHOLDS:
        cols += [f"n_near_{thr:g}A", f"n_far_{thr:g}A", f"mean_z_near_{thr:g}A",
                 f"mean_z_far_{thr:g}A", f"delta_z_{thr:g}A", f"mw_p_near_gt_far_{thr:g}A"]
    per_prot_df = per_prot_df[cols]
    per_prot_df.to_csv(os.path.join(OUTDIR, "per_protein_tests.csv"), index=False)

    # pooled test (mirrors 15's pooled test: pool z_disagree across all testable
    # proteins' sites, near vs far by 3-D distance)
    d = all_sites.dropna(subset=["z_disagree", "dist_to_ligand"])
    rho, p_rho = spearmanr(d["z_disagree"], d["dist_to_ligand"])
    pooled = {"n_proteins": all_sites["UniProt_ID"].nunique(), "n_sites": len(d),
              "spearman_zdisagree_vs_dist": rho, "spearman_p": p_rho}
    for thr in THRESHOLDS:
        near = d[d["dist_to_ligand"] <= thr]["z_disagree"]
        far = d[d["dist_to_ligand"] > thr]["z_disagree"]
        p = mannwhitneyu(near, far, alternative="greater").pvalue
        pooled[f"n_near_{thr:g}A"] = len(near)
        pooled[f"n_far_{thr:g}A"] = len(far)
        pooled[f"mean_z_near_{thr:g}A"] = near.mean()
        pooled[f"mean_z_far_{thr:g}A"] = far.mean()
        pooled[f"delta_z_{thr:g}A"] = near.mean() - far.mean()
        pooled[f"mw_p_near_gt_far_{thr:g}A"] = p
    pd.DataFrame([pooled]).to_csv(os.path.join(OUTDIR, "pooled_tests.csv"), index=False)

    n_pos_delta = (per_prot_df["delta_z_8A"] > 0).sum()
    n_sig = (per_prot_df["mw_p_near_gt_far_8A"] < 0.05).sum()
    n_test = per_prot_df["delta_z_8A"].notna().sum()
    print("=== POOLED (3-D, <=8A threshold) ===")
    print(f"pooled Spearman(z_disagree, dist_to_ligand) = {rho:+.3f} (p={p_rho:.3g}, n={len(d)})")
    print(f"pooled delta_z (near-far) = {pooled['delta_z_8A']:+.3f}, "
          f"MW p (near>far) = {pooled['mw_p_near_gt_far_8A']:.3g}")
    print(f"per-protein: higher disagreement near ligand in {n_pos_delta}/{n_test} proteins; "
          f"significant (p<0.05) in {n_sig}/{n_test}")
    for thr in (5.0, 10.0):
        pn = (per_prot_df[f"delta_z_{thr:g}A"] > 0).sum()
        ps = (per_prot_df[f"mw_p_near_gt_far_{thr:g}A"] < 0.05).sum()
        pt = per_prot_df[f"delta_z_{thr:g}A"].notna().sum()
        print(f"  @ {thr:g}A: positive in {pn}/{pt}, significant in {ps}/{pt}; "
              f"pooled delta_z={pooled[f'delta_z_{thr:g}A']:+.3f}, "
              f"p={pooled[f'mw_p_near_gt_far_{thr:g}A']:.3g}")

    print(f"\nWrote outputs to {OUTDIR}")


if __name__ == "__main__":
    main()
