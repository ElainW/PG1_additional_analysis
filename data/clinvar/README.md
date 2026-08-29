# ClinVar disease-variant exports

TSV exports downloaded from the ClinVar website on **28 Aug 2026**, used by
`scripts/23_structure_maps_all.py` to annotate disease-variant positions on the
per-protein structure maps. Only **human** multi-assay proteins are covered.

## Base query (all files)
Search the gene symbol, then filter:
- **Germline classification** = Pathogenic + Likely pathogenic
- **Molecular consequence** = missense
- **Variant type** = single nucleotide
- (**Review status** — see per-file column below)

Script `23` parses the residue position from the canonical `p.(Xxx###Yyy)` change
in the `Name` column and keeps only variants whose wild-type matches our assay-
isoform sequence (drops splice-isoform mismatches, e.g. KRAS D153V on KRAS4A).

## Per-file review-status filter
| file | gene | UniProt entry | review-status filter | note |
|---|---|---|---|---|
| `GCK.tsv` | GCK | HXK4_HUMAN | practice guideline + expert panel | |
| `PTEN.tsv` | PTEN | PTEN_HUMAN | practice guideline + expert panel | |
| `TP53.tsv` | TP53 | P53_HUMAN | practice guideline + expert panel | |
| `SRC.tsv` | SRC | SRC_HUMAN | none | very few variants |
| `CYP2C9.tsv` | CYP2C9 | CP2C9_HUMAN | none | very few variants |
| `KCNE1.tsv` | KCNE1 | KCNE1_HUMAN | none | very few variants |
| `CARD11.tsv` | CARD11 | CAR11_HUMAN | none | very few variants |
| `VKORC1.tsv` | VKORC1 | VKOR1_HUMAN | none | very few variants |
| `KRAS.tsv` | KRAS | RASK_HUMAN | none | very few variants |
| `NPC1.tsv` | NPC1 | NPC1_HUMAN | ≥ multiple submitters | nothing more stringent available |

## Not included
- **HLA-A** (Q53Z42_HUMAN) and **SLC22A1** (S22A1_HUMAN): no qualifying ClinVar entry.
- Non-human proteins (BLAT_ECOLX, CCDB_ECOLI, DYR_ECOLI, F7YBW8_MESOW, HSP82_YEAST,
  KCNJ2_MOUSE, OXDA_RHOTO, RL40A_YEAST, RPC1_LAMBD, SPG1_STRSG, SPIKE_SARS2,
  A0A2Z5U3Z0_9INFA): omitted — they keep the UniProt natural-variant flag on the map.

## Gain-of-function supplement
ClinVar pathogenic queries are dominated by loss-of-function conditions and can
miss gain-of-function disease variants. Script `23` therefore adds UniProt natural-
variant positions whose disease **description** (case-insensitive substring) matches
any of these gain-of-function search keywords (`GOF_KW` in `23_structure_maps_all.py`):

    gain-of-function      hyperactiv          hhf
    gain of function      increases activ     hyperinsulin
    constitutively activ  increase in activ   phhi
    activating            increased activ     benta

This recovers e.g. GCK hyperinsulinism/HHF3 (positions 65, 91, 99, 214, 442, 455,
456 — incl. W99C), CARD11 BENTA (123, 134), and KRAS activating variants (12, 14,
22, 34, 58, 60) that a diabetes- or LoF-only ClinVar query excludes. Matches found
per protein: GCK 7, CARD11 2, KRAS 7; 0 for PTEN, TP53, SRC, CYP2C9, KCNE1, VKORC1,
NPC1.
