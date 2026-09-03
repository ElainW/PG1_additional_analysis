# High-disagreement × functional-region enrichment

Does the high-disagreement tail (from `27`, [`../high_disagreement_sites/`](../high_disagreement_sites/))
land in functional regions more than the rest of the assayed sites, for the five
proteins with a curated allosteric-map atlas (CYP2C9, GCK, KRAS, PTEN, SRC)?
Produced by [`../../scripts/55_high_disagreement_allosteric_enrichment.py`](../../scripts/55_high_disagreement_allosteric_enrichment.py).
Exploratory; the descriptive paper table is in
[`../high_disagreement_site_characterization/`](../high_disagreement_site_characterization/) (`56`).

## Method
Universe = every assayed (DMS_id, up_pos) in the disagreement pool (`15`). A site is
"high" if it is in the `27` tail. Functional region is defined five ways and each is
tested per protein by a one-sided Fisher exact (H1: high sites MORE often functional)
plus a pooled 2×2, at (assay×position) and position-collapsed granularity:

- **ALLOSTERIC** — any position in the allosteric-map atlas (dir1).
- **ALLOSTERIC_ONLY** — atlas positions whose `site_class` contains "allosteric"
  (drops orthosteric / active-site / binding-interface residues; CYP2C9 has none).
- **UNIPROT** — UniProt binding±2 & active site (pool `near_binding`, dir3).
- **CLINVAR** — dir2 ClinVar P/LP (parsed as in `23`) ∪ UniProt GoF supplement.
- **UNION** — ALLOSTERIC ∪ UNIPROT ∪ CLINVAR.

## Finding: no enrichment

Across every definition the pooled odds ratio is ≤1 or non-significant — the
high-disagreement tail is **not** enriched in functional regions, and under the
broader definitions it is mildly depleted (pooled, assay×position):

| Definition | % high in region | % rest in region | OR | p |
|---|---|---|---|---|
| ALLOSTERIC (whole atlas) | 24.3% | 26.0% | 0.91 | 0.71 |
| ALLOSTERIC_ONLY | 9.6% | 7.3% | 1.34 | 0.20 |
| UNIPROT binding/active | 3.7% | 8.3% | 0.42 | 0.99 |
| ClinVar P/LP + GoF | 11.0% | 11.4% | 0.96 | 0.60 |
| Union of all three | 28.7% | 35.3% | 0.74 | 0.95 |

### Allosteric-class sites only — per protein (positive + negative)
Restricting to the allosteric-labelled sites gives the strongest signal but it does
not generalise: only **SRC is (borderline) positive**, GCK trends, and KRAS/PTEN go
the other way, so the pooled test is not significant (OR 1.34, p=0.20).

| Protein | allosteric sites | % high in them | % rest | OR | p |
|---|---|---|---|---|---|
| SRC | 99 | 29.2% | 14.2% | 2.49 | **0.049** (positive, borderline) |
| GCK | 46 | 10.8% | 4.7% | 2.44 | 0.11 (trend, ns) |
| KRAS | 82 | 13.3% | 23.3% | 0.51 | 0.90 (depleted) |
| PTEN | 40 | 0.0% | 5.5% | 0.00 | 1.00 (depleted) |
| CYP2C9 | 0 | — | — | — | (no allosteric-class sites) |

The SRC result is consistent with the per-protein z-score test in
[`../generalization_functional_disagreement_per_protein.csv`](../generalization_functional_disagreement_per_protein.csv),
where SRC is the only protein with significantly elevated disagreement at functional
sites. On its own it is marginal and does not survive as a general effect.

The one robust signal is **continuous**, not tail-based: mean disagreement magnitude
is modestly higher across allosteric-functional sites in all five proteins
(Mann-Whitney p = 0.0005–0.012; see `magnitude_by_protein.csv`).

## Files
- `enrichment_by_protein.csv` — one row per (protein × definition × granularity): 2×2 counts, % high/rest functional, odds ratio + 95% CI, one-sided Fisher p.
- `enrichment_pooled.csv` — pooled 2×2 per (definition × granularity).
- `magnitude_by_protein.csv` — per protein, mean disagreement in allosteric-functional vs non-functional sites + Mann-Whitney p.
