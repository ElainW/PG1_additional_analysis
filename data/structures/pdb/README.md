# Experimental (PDB) structures — curated per DMS-benchmark protein

Hand-picked crystal structures for the structure-informed variant-effect analyses,
one entry (or one **per conformational state**, for switches) per protein. Each pick
follows the DMS/allostery paper that generated the deep-mutational-scanning data, or,
where no structure was named, the selection logic in
[`../../../papers/pdb_selection_framework.md`](../../../papers/pdb_selection_framework.md).

These are the **experimental** coordinates. AlphaFold2 models for the same proteins
live alongside in [`../af2/`](../af2/) (one per UniProt ID, apo, no ligands).

> **Status.** This README records the curated *selection* and the structural features
> read directly from the coordinates (method, resolution, chain, modeled span, ligands,
> numbering offset). The deeper structure-feature analyses (per-residue distance-to-site,
> orthosteric/allosteric contact classification, state-to-state contact swapping) are
> **pending new instructions** and are not run here.

## Selected structures

| Protein (UniProt) | State | PDB · chain | Method / Res. | Key ligands (role) | Modeled span (UniProt) | DMS / source paper |
|---|---|---|---|---|---|---|
| **HXK4_HUMAN** — glucokinase / GCK (P35557) | super-open (inactive, apo) | **1V4T** · A | XRAY 3.40 Å | SO₄, Na⁺ *(apo — no substrate)* | ~16–462; loop **158–180 disordered** | Gersing et al. (GCK activity + abundance) |
| | open (intermediate) | *— not assigned —* | | | | structure uncertain in the literature |
| | closed (active) | **1V4S** · A | XRAY 2.30 Å | **GLC** (glucose, orthosteric), **MRK** (allosteric activator), Na⁺ | 15–462; **no gaps** | Gersing et al. |
| **RASK_HUMAN** — KRAS | ON / OFF / apo | see [`KRAS/`](KRAS/) | — | 6VJJ (RAF1), 5O2S (DARPin K27), 1NVW (SOS) | — | Weng, Faure & Lehner 2024 |
| **CP2C9_HUMAN** — cytochrome P450 2C9 (P11712) | single state | **1R9O** · A | XRAY 2.00 Å | **HEM** (heme, prosthetic cofactor), **FLP** (flurbiprofen, substrate/orthosteric), GOL | ~26–490; gaps 38–42, 214–220 | Amorosi et al. (CYP2C9 abundance + activity) |
| **PTEN_HUMAN** — PTEN phosphatase (P60484) | single state | **1D5R** · A | XRAY 2.10 Å | **TLA** (L-tartrate, substrate-analog at catalytic site) | 14–351; **282–312 disordered** (C-tail) | Matreyek; Mighell |
| **SRC_HUMAN** — c-Src kinase (P12931) | inactive / autoinhibited | **2SRC** · A | XRAY 1.50 Å | **ANP** (AMP-PNP, ATP analog at active site), **PTR** (phospho-Tyr, autoinhibitory C-terminal tail — the classic Tyr527) | ~86–535; no backbone gaps (pTyr modeled as PTR at auth 527) | Ahler; Chakraborty; Nguyen — Src allosteric landscape (Beltran & Lehner) |
| | active | **1Y57** · A | XRAY 1.91 Å | **MPZ** (ATP-competitive inhibitor), SO₄ | ~85–536; **no gaps** | " |
| **VKOR1_HUMAN** — VKORC1 (Q9BQB6) | active-site inhibited | **6WVH** · A,B | XRAY 1.99 Å | **UA7** (brodifacoum, anticoagulant inhibitor), OLC (LCP lipid), CRO *(fusion chromophore)* | VKORC1 = UNP **3–155** (fully modeled) | Chiasson et al. 2020 |
| | *(homolog reference)* | *4NV5 — not human* | | | Synechococcus VKOR homolog | cited in the DMS paper as the pre-human-structure reference |

State is matched to the phenotype the DMS reports on, and for a conformational **switch**
(GCK T/R allosteric enzyme; KRAS GDP/GTP GTPase; Src active/inactive kinase) one canonical
structure is kept **per state** so state-to-state contact changes can be compared — see the
framework doc.

## Numbering offsets (read before mapping DMS positions)

DMS/UniProt positions are **not** the author residue numbers in every PDB. `UniProt = auth + offset`:

| PDB | offset | note |
|---|---|---|
| 1V4T, 1V4S | **+1** | GCK auth 12 = UniProt 13 |
| 1R9O, 1D5R | **0** | 1:1 |
| 2SRC, 1Y57 | **+3** | Both Src crystals share the **same author numbering** (verified: 450 matched positions, only auth 527 differs as PTR vs TYR). Canonical P12931 alignment is **+3** for both (auth 416→Y419 activation loop, 527→Y530 regulatory tail). ⚠️ 2SRC's own `DBREF` says +2 (it cites a 535-aa reference) — using that would misalign the two states by one residue, so **use +3 for both** when comparing Src states. |
| 6WVH | **−142** | fusion construct (below) |
| KRAS set | see [`KRAS/README.md`](KRAS/README.md) | H-Ras vs KRAS, chain caveats |

## Special cases

- **GCK "open" state is unfilled.** GCK cycles super-open → open → closed. We have the two
  endpoints (super-open **1V4T**, closed **1V4S**); the intermediate *open* conformation has
  no agreed reference structure and is left unassigned. The switch signature is already visible
  in coverage alone: the mobile connecting region **158–180 (UniProt)** is **disordered in the
  super-open apo 1V4T** but **fully ordered in the glucose+activator-bound closed 1V4S**.
  *(3VEV, a glucose + activator-0H4 closed-like crystal, is present in this directory but is
  **excluded** from the curated GCK selection; it is referenced only by legacy scripts
  `18_compare_3vev_1v4s.py` / `20_native_3vev_vs_1v4s.py`.)*

- **6WVH is a fusion construct.** VKORC1_HUMAN (Q9BQB6, UniProt 3–155) occupies **author
  residues 145–297** only; it is flanked by an *Aequorea victoria* GFP-derived crystallization
  partner (UniProt A0A059PIQ0, author 1–144 and 298–383) — hence the `CRO` chromophore HETATM,
  which is **part of the fusion, not a VKOR ligand**. The VKORC1 portion is fully modeled and
  brodifacoum (`UA7`) sits in the VKOR active site. Two copies are present (chains A, B).
  A processed **`6WVH_cleaned.pdb`** (script `49_clean_6wvh.py`) keeps only VKORC1 (author
  145–297, chains A+B) + brodifacoum and drops the GFP fusion, `CRO`, `OLC`, and waters — use it
  for isolated-VKOR analyses; author numbering is preserved (UniProt = author − 142). The two
  VKORC1 copies pack against each other, so for **per-residue rSASA** use the single-chain
  **`6WVH_cleaned_chainA.pdb`** (chain A only) — otherwise the A↔B crystal contact is scored as
  buried surface. Use the two-chain file only where both copies are wanted.

- **PTEN — mapping the active site.** 1D5R carries no natural substrate; use the bound
  **`TLA` (L(+)-tartrate)**, the substrate-analog sitting in the catalytic pocket, to define the
  active site — active-site residues = those within contact/short heavy-atom distance of TLA
  (the same functional-site-by-bound-molecule logic used for the other enzymes).

- **VKOR 4NV5 is not human.** The VKOR DMS paper (Chiasson et al. 2020) built on **4NV5**, a
  *Synechococcus* bacterial VKOR homolog — the only VKOR structure available at the time. The
  human structure **6WVH** postdates the paper and was selected here; keep 4NV5 in mind only as
  a distant homolog reference (~sequence identity to human VKORC1 is low).

- **1NVW is mixed-state — use the right Ras chain.** See [`KRAS/README.md`](KRAS/README.md):
  chain R (nucleotide-free, catalytic/Cdc25 site) is the SOS exchange substrate; chain Q
  (GNP-bound, allosteric/Rem site) is a different interface.

- **KRAS — which coordinate file for which analysis.** Two file types live in [`KRAS/`](KRAS/):
  the **processed single-chain** `KRAS/*_*P.pdb` (`1NVW_APO_0P`, `5O2S_OFF_2P`, `6VJJ_ON_3P` —
  partner removed, waters stripped, superposed) and the **full original** co-complexes
  `KRAS/{6VJJ,5O2S,1NVW}.pdb`.
  - **Relative SASA (rSASA):** use the **processed `*_*P.pdb`** files, so burial is computed for
    KRAS *in isolation* — with the partner chain present, interface residues would be scored as
    buried by the partner rather than by KRAS itself.
  - **Binding-site distance mapping:** use the **full original** `.pdb` (partner chain retained),
    since the functional site is defined as residues within cutoff of the bound partner and that
    partner must be present.

## Files

**Originals (RCSB downloads):**

| File | Contents |
|---|---|
| `1V4T.pdb`, `1V4S.pdb` | GCK super-open (apo) and closed (glucose+activator) |
| `1R9O.pdb` | CYP2C9 + heme + flurbiprofen |
| `1D5R.pdb` | PTEN + tartrate |
| `2SRC.pdb`, `1Y57.pdb` | c-Src inactive (AMP-PNP + pTyr tail) and active (inhibitor) |
| `6WVH.pdb` | Human VKORC1–GFP fusion + brodifacoum |
| `3VEV.pdb` | GCK glucose + activator (closed-like) — **not** in the curated set; legacy scripts only |

**Processed (cleaned for analysis — waters + crystallization additives stripped, single target chain, functional ligands / modified residues kept):**

| File | Contents |
|---|---|
| `1V4T_cleaned.pdb` | GCK super-open, chain A, apo (no functional ligand) |
| `1V4S_cleaned.pdb` | GCK closed, chain A + GLC + MRK |
| `1R9O_cleaned.pdb` | CYP2C9, chain A + HEM + FLP |
| `1D5R_cleaned.pdb` | PTEN, chain A + TLA |
| `2SRC_cleaned.pdb` | c-Src inactive, chain A + ANP + PTR (pTyr527) |
| `1Y57_cleaned.pdb` | c-Src active, chain A + MPZ |
| `6WVH_cleaned.pdb` | VKORC1 (auth 145–297, chains A+B) + brodifacoum |
| `6WVH_cleaned_chainA.pdb` | VKORC1 single copy (chain A) + brodifacoum — for per-residue rSASA |
| `KRAS/*_*P.pdb` | KRAS ON/OFF/apo single-chain per-state set (own README) |

Cleaners: `scripts/50_clean_structures.py` (the six single-chain structures),
`scripts/49_clean_6wvh.py` (VKORC1 fusion), `scripts/visualize_kras_states.pml` (KRAS).

**Other:**

| File | Contents |
|---|---|
| `KRAS/` | KRAS ON/OFF/apo per-state set + PyMOL script (own README) |
| `../kras_rascore/` | rascore conformational-state tables for KRAS (0P/2P/3P) |

## Provenance

Original RCSB coordinates (legacy PDB format), fetched from
`https://files.rcsb.org/download/<ID>.pdb`. UniProt mappings and residue ranges above are
read from each file's `DBREF` records; modeled spans, internal gaps, and ligand lists are read
directly from the `ATOM`/`HETATM`/`REMARK 465` records of these coordinates.
