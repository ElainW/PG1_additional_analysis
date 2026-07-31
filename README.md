# ProteinGym v1 — model predictions analyzed *by position*

Scripts and outputs for: (1) does each model's Spearman vs `DMS_score` differ by
**taxon** and by **model_type**, and (2) how the **ranking of mutations at each
mutation site** differs across models and model types.

Data: 217 DMS substitution assays in `O2/20250613/*.csv`, 97 models (columns
`Site_Independent … Protriever`), metadata in
`ProteinGym/reference_files/DMS_substitutions.csv`, model config in
`ProteinGym/config.json`.

## How to run
```bash
cd analysis/scripts
python3 00_check_directionality.py   # directionality sanity check
python3 01_compute_spearman.py       # per-assay Spearman -> results/spearman_long.csv
python3 02_spearman_by_taxa.py       # Analysis 1 (taxa / model_type) + figs 1-4
python3 03_persite_ranking.py        # Analysis 2 (within-site ranking) + figs 5-8
```
Requires `pandas, numpy, scipy, matplotlib, seaborn`.

## Directionality (important)
`config.json` lists a per-model `directionality` (+1 = higher score ⇒ higher
fitness; -1 = higher score ⇒ lower fitness). **The score columns in
`O2/20250613/*.csv` are already oriented so that higher = fitter for every
model.** Evidence (`00_check_directionality.py`, `results/directionality_check.csv`):
all 8 `directionality=-1` models correlate *positively* with `DMS_score`
(EVE +0.46, RSALOR +0.48, DeepSequence +0.44, Wavenet +0.39, Unirep +0.18…).
If the columns were raw model output these would be negative. So the config
directionality was applied upstream and we **do not** re-apply it — matching
ProteinGym's own `performance_DMS_benchmarks.py`, which calls
`spearmanr(DMS_score, score)` directly.

## Methods
**Analysis 1 — Spearman by taxon.** Per assay × model, Spearman of the model
column vs `DMS_score` over all variants (`spearman_long.csv`). Aggregated to
model×taxon and model_type×taxon means. Per-model and per-model_type
Kruskal–Wallis tests ask whether the per-assay Spearman distribution differs
across the four taxa; p-values BH-FDR corrected.

**Analysis 2 — within-site mutation ranking.** Single substitutions only; group
variants by residue position (the "site"); within each site (≥3 substitutions)
convert each model's score to a normalized rank in [0,1]. Pearson correlation of
these within-site ranks, pooled over variants and averaged across assays
(Fisher-z), gives (a) a model×model agreement matrix and (b) each model's
within-site agreement with the experiment's own within-site ranking. This
removes the easy across-position signal and isolates *which* substitution a
model prefers at a position.

## Key findings
**Taxon.** Pooling all models, mean Spearman is similar for Human (0.42),
Eukaryote (0.43), Prokaryote (0.41) but markedly lower for **Virus (0.35)**.
38/97 models show a statistically significant taxon effect (BH q<0.05). Small,
sequence-only models are the most taxon-sensitive (ESM2-35M/150M, ESMC, small
SaProt have the largest Human↔other gaps and worst Virus performance), whereas
MSA and Structure&MSA models are the most taxon-robust.

**Model type.** Ranking by mean Spearman is stable across taxa:
Structure&MSA (0.51) > Function-annotation/ESM3 (0.48) ≈ Single-seq&Structure
(0.46) > MSA (0.44) > Structure (0.37) > Single-sequence (0.36). The biggest
taxon swing is for sequence-only and single-seq&structure types on Virus.

**Within-site ranking.** Models agree far less on *within-site* rankings than on
overall scores. Model types cluster: MSA↔Structure&MSA agree most (0.70–0.73),
while pure **Structure** models are the outliers (≈0.37–0.44 agreement with
sequence/MSA types) — they capture a different, structure-driven local
preference. Recovering the experiment's within-site ranking is hardest overall
(top model ≈0.35): best are structure-aware models with MSA/sequence context
(AIDO.Protein-RAG, VenusREM, ProSST-2048, ESM3, S3F); worst are autoregressive
single-sequence LMs (xTrimoPGLM-CLM, ProtGPT2, Unirep, CARP). Within-site DMS
agreement is again lowest on Virus for most types.

## Outputs
`results/`
- `directionality_check.csv` — per-model mean Spearman + config directionality.
- `spearman_long.csv` — per assay × model Spearman (+ taxon, model_type, n).
- `spearman_model_by_taxon.csv`, `spearman_modeltype_by_taxon.csv` — means.
- `kruskal_model_across_taxa.csv`, `kruskal_modeltype_across_taxa.csv` — tests.
- `persite_model_agreement_matrix.csv` — 98×98 within-site agreement (incl. DMS).
- `persite_modeltype_agreement_matrix.csv` — model_type×model_type agreement.
- `persite_dms_agreement_long.csv` / `_by_model.csv` / `_by_modeltype_taxon.csv`.

`figures/`
- `fig1_model_by_taxon_heatmap.png` — model × taxon mean Spearman.
- `fig2_modeltype_by_taxon_box.png` — per-assay Spearman by type × taxon.
- `fig3_modeltype_by_taxon_heatmap.png` — model_type × taxon mean.
- `fig4_taxon_overall_box.png` — overall Spearman by taxon.
- `fig5_persite_model_agreement_clustermap.png` — clustered model×model within-site agreement.
- `fig6_persite_modeltype_agreement_heatmap.png` — model_type×model_type within-site agreement.
- `fig7_persite_dms_agreement_by_model.png` — each model's within-site agreement with DMS.
- `fig8_persite_dms_agreement_modeltype_taxon.png` — within-site DMS agreement by type × taxon.
