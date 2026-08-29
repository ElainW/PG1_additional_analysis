# additional_annotations

Curated / derived UniProt annotations for the generalization candidate proteins
(`results/generalization_candidate_proteins.csv`), kept alongside the raw UniProt
entries in `../raw_json/`. Everything here is in **UniProt numbering**, which
matches the `up_pos` column of
`results/generalization_functional_disagreement_sites.csv` and the residue axis
of the `23_structure_maps_all.py` figures.

## Structural / functional annotations

Written by `scripts/25_export_structural_annotations.py`, which imports the
feature-parsing helpers from `scripts/23_structure_maps_all.py` verbatim so these
records cannot drift from what the structure-map figures overlay.

- **`structural_annotations.md`** — human-readable, one section per protein.
- **`structural_annotations_all.csv`** — long format, one row per
  `(UniProt_ID, category, name)` with columns
  `UniProt_ID, gene, category, name, n_residues, residues`
  (`residues` compacted as `a-b, c, e-f`).

Categories:

| category      | source (UniProt feature)                    | `name` holds            |
|---------------|---------------------------------------------|-------------------------|
| `binding`     | `Binding site`, grouped per ligand          | ligand name (ATP, substrate, heme, FAD, GTP, ...) |
| `active_site` | `Active site`                               | `Active site`           |
| `subdomain`   | `Region` / `Domain` (excl. disorder)        | region/domain name      |
| `disordered`  | `Region` described as *Disordered*          | `Disordered`            |
| `disease`     | see disease-source rule below               | disease source          |

`subdomain` applies the same legibility filters as the figures: whole-protein
domains (≥85% of length) and short `Region` motifs (<15 aa) are dropped. Metal
ions (e.g. `Mg(2+)`) appear as their own `binding` ligand because UniProt records
them as `Binding site` entries. `disordered` ranges are the raw UniProt spans
(untrimmed).

### Disease-variant source

Matches the disease rug in `23_structure_maps_all.py`:

- The 10 human proteins with a `data/clinvar/` export → **ClinVar pathogenic /
  likely-pathogenic** missense (WT-matched to the assay isoform) **plus the
  UniProt gain-of-function supplement** below.
- Every other protein → the **UniProt natural-variant disease flag** as
  materialised in the sites CSV's `disease` column.

## Gain-of-function supplement

ClinVar pathogenic queries are dominated by loss-of-function conditions and miss
gain-of-function disease variants (e.g. GCK hyperinsulinism/HHF3, CARD11 BENTA,
KRAS activating). These files add UniProt natural-variant positions whose disease
description carries GoF language (see `GOF_KW` in `23_structure_maps_all.py`). Written by
`scripts/24_gof_supplement_tables.py`.

- **`gof_supplement_<UniProt_ID>.csv`** — per-protein GoF variants
  (`position, wt, alt, variant, matched_keywords, description`).
- **`gof_supplement_summary.csv`** — per-protein counts and positions.

Only proteins with a `data/clinvar/` export are supplemented; a protein with 0
GoF variants still appears in the summary with an empty position list.

## Regenerating

```bash
cd analysis/scripts
python3 24_gof_supplement_tables.py            # gof_supplement_*.csv
python3 25_export_structural_annotations.py    # structural_annotations.{md,csv}
```

Both read the UniProt entries in `../raw_json/`; rerun them if those change.
