# Structure features

Structure-derived features for the DMS-benchmark proteins, computed on the curated
PDB structures in [`../../data/structures/pdb/`](../../data/structures/pdb/) (see that
directory's README for structure selection, cleaning, and numbering offsets).

Environment: conda env **`structfeat`** (freesasa CLI 2.1.2, pymol-open-source,
vmd-python + GetContacts, biopython). Scripts `48`–`54` in [`../../scripts/`](../../scripts/).

## Inputs used per analysis

| Protein | rSASA / contacts input | notes |
|---|---|---|
| HXK4 (GCK) | `1V4S_cleaned`, `1V4T_cleaned` | active=closed 1V4S, inactive=super-open 1V4T |
| CP2C9 | `1R9O_cleaned` | single state |
| PTEN | `1D5R_cleaned` | single state |
| SRC | `1Y57_cleaned`, `2SRC_cleaned` | active=1Y57, inactive=2SRC |
| VKOR1 | `6WVH_cleaned_chainA` | single chain (avoids A↔B crystal contact) |
| RASK (KRAS) | `KRAS/6VJJ_ON_3P`, `5O2S_OFF_2P`, `1NVW_APO_0P` | active=ON, inactive=OFF, apo=3rd state |

## `ligands_by_structure.*`  (script 48)
Ligand / substrate inventory per structure (functional ligand vs apo, chain-aware).

## `rsasa/`  (script 51)
Relative SASA via `freesasa <pdb> -n 20 --format=rsa --radii=naccess` (Lee-Richards,
20 slices, NACCESS radii), protein-only (ATOM records; functional HETATM ignored).
- `<label>.rsa` — raw NACCESS-format output per structure
- `rsasa_all.csv` — parsed long table; `rel_all`/`rel_side` = relative SASA (% burial)

## `secondary_structure/`  (script 52)
Secondary structure from UniProt features (Helix→helix, Beta strand→strand,
Turn→turn, else loop), in UniProt numbering, per UniProt ID.
- `<UNIPROT>_<PROT>_ss.tsv` — per-residue pos/aa/ss
- `secondary_structure_all.csv`, `secondary_structure_summary.csv`

## `contacts/`  (scripts 53–54)
Static residue-residue contacts via GetContacts after PyMOL `h_add`, `--itypes all`,
`--sele "protein or resname PTR"` (protein residues + the phospho-Tyr in 2SRC).
- `<label>_contacts.tsv` — raw GetContacts, all interaction types (**all** structures)
- `contacts_itype_summary.csv` — per-structure interaction-type counts
- `hadd/<label>_H.pdb` — hydrogen-added coordinates fed to GetContacts
- `4types/<label>_4types.tsv` — salt bridge (`sb`), pi-cation (`pc`), sc–sc H-bond
  (`hbss`), sc–bb H-bond (`hbsb`) — the state-specific types (**multi-state** proteins)

### `contacts/classification/`  (script 54)
For the multi-conformation proteins (HXK4, RASK, SRC), each contacting residue is
labeled by comparing active vs inactive state 4-type contacts, **in UniProt numbering**:
`active` / `inactive` / `nonspecific` (identical partners both states) / `swapping`
(contacts in both states but to different partners).
- `<PROTEIN>_residue_classes.csv` — per-residue label + partner sets
- `classification_summary.csv` — counts per class per protein
- `RASK_apo_contacts_note.csv` — KRAS apo (1NVW) contacts, the 3rd state

⚠️ **SRC numbering:** 2SRC and 1Y57 share author numbering; both map to canonical
P12931 at **+3** (2SRC's DBREF says +2, which would misalign the states by one). The
classification uses +3 for both.
