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

---

# Results & interpretation

## rSASA — burial (`rsasa/rsasa_all.csv`)

Fraction of residues buried (`rel_all` < 25%), per structure:

| structure | protein / state | % buried |
|---|---|---:|
| 1V4T | GCK super-open (inactive) | 52 |
| 1V4S | GCK closed (active) | **59** |
| 1R9O | CYP2C9 | 55 |
| 1D5R | PTEN | 50 |
| 2SRC | Src inactive | 54 |
| 1Y57 | Src active | 48 |
| 6WVH | VKOR1 | 41 |
| 6VJJ | KRAS ON | 50 |
| 5O2S | KRAS OFF | 48 |
| 1NVW | KRAS apo | 42 |

- **GCK closes on glucose:** the closed active state (1V4S, 59% buried) buries more
  surface than the super-open apo state (1V4T, 52%) — the domain closure signature,
  consistent with the 158–180 loop ordering only in the closed state.
- ⚠️ **VKOR1 is a membrane protein.** rSASA is computed on the isolated protein with no
  membrane, so transmembrane-helix surface reads as "exposed" — its 41% buried is a
  lower bound and TM residues' rSASA should not be read as true solvent exposure.

## Secondary structure (`secondary_structure/secondary_structure_summary.csv`)

Per-residue helix / strand / turn / loop counts (UniProt numbering), fractions:

| protein | len | %helix | %strand | %turn | %loop | fold note |
|---|---:|---:|---:|---:|---:|---|
| HXK4 (GCK) | 465 | 44 | 27 | 4 | 25 | α/β kinase |
| RASK (KRAS) | 189 | 40 | 31 | 2 | 28 | G-domain α/β |
| CP2C9 | 490 | 50 | 14 | 4 | 32 | helix-rich P450 |
| PTEN | 403 | 24 | 36 | 1 | 39 | strand-rich phosphatase + C2 |
| SRC | 536 | 28 | 25 | 4 | 43 | SH3 (β) + SH2 + bilobal kinase |
| VKOR1 | 163 | **65** | 2 | 0 | 33 | all-α transmembrane bundle |

Fold-consistent: VKOR1 is almost entirely helical (TM helices), PTEN is strand-rich,
the P450 (CYP2C9) is helix-dominated.

## Contacts (`contacts/contacts_itype_summary.csv`)

Interaction counts scale with size; the four state-specific types (salt bridge `sb`,
π-cation `pc`, sc–sc H-bond `hbss`, sc–bb H-bond `hbsb`) are what the classification
uses. State-pair comparison of these:

| protein | state | sb | pc | hbss | hbsb |
|---|---|---:|---:|---:|---:|
| KRAS | ON (6VJJ) | **27** | 7 | 36 | 54 |
| KRAS | OFF (5O2S) | 10 | 3 | 28 | 56 |
| GCK | closed (1V4S) | 55 | 14 | 80 | 140 |
| GCK | super-open (1V4T) | 51 | 19 | 82 | 129 |
| Src | inactive (2SRC) | 33 | 19 | 65 | **199** |
| Src | active (1Y57) | 24 | 16 | 60 | 163 |

- **KRAS ON has ~2.7× the salt bridges of OFF** (27 vs 10) — GTP-analog + M²⁺ + engaged
  switch I/II form extra ionic contacts lost on GDP.
- **Src inactive has more sc–bb H-bonds** (199 vs 163) — the compact autoinhibited
  assembly (SH3–SH2 clamp onto the kinase) adds regulatory contacts.

## Active/inactive residue classification (`contacts/classification/`)

`classification_summary.csv`:

| protein | active-only | inactive-only | nonspecific | swapping |
|---|---:|---:|---:|---:|
| HXK4 (1V4S vs 1V4T) | 62 | 43 | 70 | 103 |
| SRC (1Y57 vs 2SRC) | 32 | 55 | **105** | 101 |
| RASK (6VJJ vs 5O2S) | 17 | 9 | 40 | 38 |

How to read a residue's row (`<PROTEIN>_residue_classes.csv`): `nonspecific` = same
contact partners in both states (conserved scaffold); `swapping` = contacts in both
states but to different partners (conformational rewiring); `active`/`inactive` = only
engaged in that one state.

- **KRAS — a clean switch.** Switch I (30–38) is **100% state-specific** (no nonspecific
  residues); Switch II (60–76) is mostly swapping; the P-loop / G-domain core carries the
  nonspecific (constant) contacts. Textbook GDP↔GTP behaviour. The apo/0P state (1NVW) is
  reported separately in `RASK_apo_contacts_note.csv` as a third state.
- **Src — autoinhibition is visible per residue.** More inactive-only than active-only
  residues (55 vs 32): the autoinhibited state has extra regulatory contacts. **Y419**
  (activation loop) is `swapping`; **Y530** (C-terminal regulatory tail) is `swapping`
  with **4 inactive-state contacts vs 1 active** — the pTyr530→SH2 clamp that locks the
  inactive kinase.
- **GCK — large domain motion.** High swapping (103) + active-only (62) reflects the big
  super-open→closed hinge motion; the 70 nonspecific residues are the conserved core.

### Caveats
- Contacts are **residue–residue only** (`--sele "protein or resname PTR"`); functional
  HETATM ligands (glucose, heme, nucleotide, inhibitors) are excluded, so ligand-mediated
  contacts are not in these tables.
- The classification is a **two-state** comparison per protein (active vs inactive); for
  KRAS the apo state is provided as a separate reference rather than folded into the 4-way
  labels.
