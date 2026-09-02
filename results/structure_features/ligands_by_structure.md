# Ligand / substrate inventory of the curated PDB structures

Read directly from the `HETATM` records of each structure (`scripts/48_structure_ligands.py`). `functional_ligand_bound` = yes if a substrate, cofactor, nucleotide, allosteric modulator, inhibitor, **or** a bound protein partner is present; structures carrying only ions / cryoprotectants are flagged apo.

`*` 3VEV is present in the directory but **not** part of the curated GCK selection.

| structure | protein | functional_ligand_bound | substrate | cofactor | nucleotide | allosteric | inhibitor | metal_Mg | bound_partner | modification | ions_additives | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1V4T | HXK4 / GCK | NO (apo/additives only) | — | — | — | — | — | — | — | — | SO4 x2; NA | apo super-open (inactive) — active site empty; sulfate/Na only |
| 1V4S | HXK4 / GCK | yes | GLC | — | — | MRK | — | — | — | — | NA | — |
| 3VEV* | HXK4 / GCK | yes | GLC | — | — | 0H4 | — | — | — | — | NA | — |
| 1R9O | CP2C9 | yes | FLP | HEM | — | — | — | — | — | — | GOL | — |
| 1D5R | PTEN | yes | TLA | — | — | — | — | — | — | — | — | — |
| 2SRC | SRC | yes | — | — | ANP | — | — | — | — | PTR | — | PTR = autoinhibitory phospho-Tyr (auth 527) modeled as a modified residue |
| 1Y57 | SRC | yes | — | — | — | — | MPZ | — | — | — | SO4 x5 | — |
| 6WVH | VKOR1 | yes | — | — | — | — | UA7 | — | — | CRO | OLC x2 | CRO chromophore + OLC lipid belong to the GFP crystallization fusion, not VKOR; VKORC1 = author residues 145-297; off-target-chain: UA7 |
| 6VJJ | RASK / KRAS | yes | — | — | GNP | — | — | MG | RAF1 RBD (effector) (chain B) | — | CL x3; EDO x2; MPD | — |
| 5O2S | RASK / KRAS | yes | — | — | GDP | — | — | MG | DARPin K27 (GDP-selective binder) (chain B) | — | — | 4 KRAS copies in the AU, each GDP+Mg-bound; chain A analysed; off-target-chain: GDP x3 |
| 1NVW | RASK / H-Ras | yes | — | — | — | — | — | — | SOS1 (GEF) (chain S) | — | PO4 x7 | mixed-state: GNP+Mg on chain Q (allosteric Ras); target SOS-substrate chain R is nucleotide-free (PO4 only); off-target-chain: GNP |
