# KRAS structures — nucleotide-state comparison (ON / OFF / apo)

Curated KRAS/RAS structures for the effector/partner-binding analysis, one per
nucleotide state, plus the PyMOL script that superposes and visualizes them.
Rationale and the general selection logic live in
[`../../../../papers/pdb_selection_framework.md`](../../../../papers/pdb_selection_framework.md);
this file records the concrete KRAS choices and the per-state coordinate files.

## The three states

| State | Structure (chain) | Nucleotide | Bound partner | Method / Res. | Modified file |
|---|---|---|---|---|---|
| **ON (3P)** — active | 6VJJ chain A (KRAS4b) | GppNHp (GNP) + Mg²⁺ | RAF1 RBD (effector) | XRAY 1.4 Å | `6VJJ_ON_3P.pdb` |
| **OFF (2P)** — inactive | 5O2S chain E (KRAS) | GDP + Mg²⁺ | DARPin K27 (GDP-selective binder) | XRAY 3.22 Å | `5O2S_OFF_2P.pdb` |
| **apo (0P)** — nucleotide-free | 1NVW chain R (H-Ras) | none | SOS1 (GEF, catalytic site) | XRAY 2.7 Å | `1NVW_APO_0P.pdb` |

State is matched to the partner: effectors (RAF1) bind GTP-active KRAS; DARPin K27
recognizes inactive GDP-bound KRAS; SOS captures the nucleotide-free exchange
intermediate.

## Decisions recorded

- **Use co-complexes, one per state (bound + unbound).** The binding phenotype defines
  the functional site as residues < 5 Å from the partner, so each state is represented by
  a partner co-complex rather than a bare apo model. Characterizing KRAS as a switch needs
  multiple conformations of the same protein across states, in both bound and unbound
  forms — see the framework doc.

- **1NVW is a mixed-state structure — use chain R, not chain Q.** It carries two Ras copies
  on one SOS: chain **R** = nucleotide-free (0P) at the **catalytic (Cdc25) site** — the
  canonical SOS–Ras exchange interface — and chain **Q** = GNP/3P at the distal
  **allosteric (Rem) site**. Verified from coordinates: chain R contacts 40 Cdc25-domain
  residues (vs 2 Rem) including the catalytic hairpin 929–940, with switch I + II fully
  engaged. The entry-level ligand list shows GNP only because of chain Q; the 0P assignment
  is correct **for chain R**. Compute the SOS interface against chain R.

- **Structures follow Chenchun Weng, Faure & Lehner 2024 (*Nature*).** We use the KRAS
  partner co-complexes cited in that paper (here 6VJJ, 5O2S, 1NVW), not the apo entries
  enumerated in [`../kras_rascore/KRAS_0P.tsv`](../kras_rascore/KRAS_0P.tsv).

## Files

| File | Contents |
|---|---|
| `6VJJ.pdb`, `5O2S.pdb`, `1NVW.pdb` | Original RCSB coordinates (legacy PDB format) |
| `6VJJ_ON_3P.pdb`, `5O2S_OFF_2P.pdb`, `1NVW_APO_0P.pdb` | Modified: single Ras chain, waters removed, superposed on 6VJJ chain A (G-domain core); nucleotide + Mg²⁺ retained where present |
| `visualize_kras_states.pml` | PyMOL script — loads the originals, aligns, colors switch I/II per state, and re-saves the modified files. Canonical copy: [`../../../../scripts/visualize_kras_states.pml`](../../../../scripts/visualize_kras_states.pml) |

## Reproduce

```
pymol visualize_kras_states.pml     # or, in PyMOL:  @visualize_kras_states.pml
```

Run from this directory or the repo root; the script resolves the KRAS path either way.
Superposition detail: G-domain core fits at ~0.5 Å; switch I spreads to ~11 Å in the apo
(SOS-opened) state, and the ON↔OFF difference is largest in switch II (~3.9 Å).
