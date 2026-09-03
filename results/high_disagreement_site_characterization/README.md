# High-disagreement site characterization

Descriptive characterization of *where* the high-disagreement sites land, for the
five proteins with a curated allosteric-map atlas (CYP2C9, GCK, KRAS, PTEN, SRC).
Produced by [`../../scripts/56_high_disagreement_site_characterization.py`](../../scripts/56_high_disagreement_site_characterization.py)
(base `pg2` env). This is purely descriptive — the enrichment test lives in
[`../high_disagreement_allosteric_enrichment/`](../high_disagreement_allosteric_enrichment/)
(script `55`), which finds no functional enrichment.

## Method

High-disagreement sites are the tail exported by `27`
([`../high_disagreement_sites/`](../high_disagreement_sites/)), collapsed to unique
residue positions (a residue flagged in two assays counts once; `n_high_rows` keeps
the raw assay×position count). Positions are UniProt-numbered; these five assays are
offset 0, so `up_pos` == the secondary-structure / annotation position.

Each site is labelled by the region it sits in:

| Label | Definition | Source |
|---|---|---|
| **functional** | binding / active / allosteric site **or** disease residue — the UNION of the four sub-flags below | dir1 + dir3 + dir2 |
| **beginning/end** (`terminal`) | first/last `TERM_WIN`=20 residues, or inside a UniProt disordered segment abutting a terminus (SRC 1–53, PTEN 352–403) | length + `structural_annotations_all.csv` |
| **loop/turn** | UniProt secondary structure is loop or turn (not helix/strand) | `52` → `secondary_structure_all.csv` |
| **core** | structured (helix/strand), non-terminal, non-functional (the remainder) | — |

Functional sub-flags (non-exclusive):
- `func_allosteric` — curated allosteric-map atlas (`26`/`47`, dir1)
- `func_uniprot_binding` — UniProt binding±2 & active site (the pool's `near_binding`, `15`, dir3)
- `func_clinvar` — dir2 ClinVar pathogenic/LP, parsed as in `23` (WT-matched) ∪ UniProt GoF supplement
- `func_uniprot_disease` — UniProt natural-variant disease (the pool's `disease` column)

The headline percentages are a **mutually exclusive** classification under the
priority `functional > terminal > loop_turn > core`, so they sum to 100% per protein.
The raw non-exclusive flags are kept alongside so nothing is hidden by the priority
choice (e.g. loop/turn is 42–69% of sites when counted non-exclusively).

## Files

### `site_characterization_simplified.csv`
The paper table: one row per protein — `length`, `n_high_disagreement_sites`, and the
four mutually exclusive percentages (`pct_functional`, `pct_beginning_end`,
`pct_loop_turn`, `pct_core`).

### `site_characterization_by_protein.csv`
Full per-protein table: the exclusive classes plus the raw non-exclusive flags
(`pct_functional_any`, `pct_terminal_any`, `pct_loopturn_any`) and every functional
sub-type percentage.

### `site_characterization_sites.csv`
Per-site audit trail: one row per high-disagreement residue with its `ss`, `length`,
each boolean label, the assigned `exclusive_class`, and the four functional sub-flags.

## Notes

- **GCK is the outlier (82.6% functional)** and is driven almost entirely by UniProt
  disease annotation (`func_uniprot_disease` 78%): GCK carries hundreds of MODY /
  hyperinsulinism natural variants, so nearly every residue overlaps a disease
  position. Its binding/allosteric share is only ~22%.
- Counted non-exclusively, loops/turns dominate every protein (42–69%), the clearest
  single structural signal — most high-disagreement residues are in loops/turns
  whatever else they also are.
- `TERM_WIN` (terminal window, default 20) and the priority order are set at the top
  of script `56`; change there and re-run to regenerate all three files.
