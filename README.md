# ProteinGym v1 — model predictions analyzed *by position*

Scripts, results, and findings for a by-position study of the ProteinGym v1
zero-shot substitution models: how per-model performance varies by **taxon** and
**model_type**, how models rank mutations **within a site**, where they
**disagree**, and whether that disagreement tracks functional/structural context.
Deep-dive protein: **HXK4_HUMAN** (glucokinase / GCK).

Data (not tracked in this repo): 217 DMS substitution assays in
`../O2/20250613/*.csv`, 97 models (columns `Site_Independent … Protriever`),
metadata `../ProteinGym/reference_files/DMS_substitutions.csv`, config
`../ProteinGym/config.json`. UniProt entries for the multi-assay proteins are
cached in `data/uniprot/` (fetched by scripts 08/15).

Requires `pandas, numpy, scipy, matplotlib, seaborn` (plus `biopython` for the
3-D analysis `17`). Figures are **not** tracked (regenerable); see `.gitignore`.

## Directionality (read first)
`config.json` lists a per-model `directionality` (+1 = higher score ⇒ higher
fitness; −1 ⇒ lower). **The score columns in `O2/20250613/*.csv` are already
oriented so higher = fitter for every model** — verified in `00`: all 8
`directionality=−1` models correlate *positively* with `DMS_score` (EVE +0.46,
RSALOR +0.48, DeepSequence +0.44…). So we never re-apply directionality, matching
ProteinGym's own `performance_DMS_benchmarks.py`.

## Model subsets used
Two curated subsets exist for different jobs: the **45 family-best** set
deduplicates the 97 so redundant same-family variants don't inflate agreement
statistics (used wherever we pool over many models), while the **23-model
subset** is a smaller hand-picked representative panel meant to be *read* in a
single figure (per-model bars, labelled heatmaps) without 97 rows of clutter.
They are not nested and use slightly different `model_type` labels, so numbers
are only comparable within the same subset.
- **97 models** — full config list (analyses 01–03).
- **45 family-best** (`05`) — one model per family (parameter-size / single-vs-
  ensemble variants collapsed), keeping the highest overall-Spearman member.
  Used for disagreement work (06, 08, 11, 15, 16).
- **23-model subset** (`04`) — a hand-curated representative set (user-supplied
  `model_dict`, with its own `model_type` labels). Used for within-site
  Spearman-vs-experiment work (09, 12, 13, 14).

## How to run
```bash
cd scripts
python3 00_check_directionality.py            # directionality verification
python3 01_compute_spearman.py                # per-assay Spearman -> spearman_long.csv
python3 02_spearman_by_taxa.py                # taxa / model_type + figs 1-4
python3 03_persite_ranking.py                 # within-site agreement matrices + figs 5-8
python3 04_subset_figs.py                     # Fig5_subset + Fig6 by selection_type (23 models)
python3 05_fig6_familybest.py                 # Fig6 from 45 family-best models
python3 06_hxk4_site_disagreement.py          # HXK4 per-site disagreement + overview figs
python3 07_w99_deepdive.py                    # W99 rankings (45-model view)
python3 08_gck_structure_map.py               # disagreement vs GCK functional map (UniProt)
python3 09_w99_subset_figs.py                 # W99 rankings + per-model-vs-DMS stat (23 models)
python3 10_assay_overlap.py                   # multi-assay proteins sharing coordinates
python3 11_hxk4_lowagreement_batch.py         # 06/09 batched over all low-agreement HXK4 sites
python3 12_abundance_physchem.py              # abundance vs hydrophobicity + charge
python3 13_agreement_vs_accuracy.py           # among-model agreement vs accuracy-vs-experiment
python3 14_quartile_by_agreement.py           # 'winners' collapse, high- vs low-agreement sites
python3 15_functional_disagreement_generalization.py   # generalization across multi-assay proteins
python3 16_collapse_modeltype_withinsite.py   # collapse by model_type; global vs within-site split
python3 17_gck_3d_distances.py                # GCK disagreement vs 3-D ligand distance (PDB 1V4S)
python3 18_generalization_candidates.py       # multi-assay candidate proteins + RCSB PDB options
```

## Reproduce
`make all` runs scripts `00 → 16` in order (each reads the previous ones'
`results/`), regenerating every CSV and figure and logging timestamps + output to
`RUN.log`. This is how to confirm the numbers quoted below stay in sync with the
code. Requires the `data/` UniProt cache (tracked) for `08`/`15`; no network
needed. `make clean` removes generated `results/`/`figures/` (leaves `data/`).

## Key findings

**Taxon & model_type (02).** Pooled over models, mean Spearman is similar for
Human (0.42) / Eukaryote (0.43) / Prokaryote (0.41) but lower for **Virus
(0.35)**; 38/97 models show a significant taxon effect (BH q<0.05), the most
taxon-sensitive being small sequence-only models. Model-type ranking is stable:
Structure&MSA (0.51) > ESM3/Function (0.48) ≈ Single-seq&Structure (0.46) >
MSA (0.44) > Structure (0.37) > Single-sequence (0.36).

**Within-site ranking agreement (03–05).** Models agree far less on *within-site*
rankings than on overall scores. Types cluster: MSA↔Structure&MSA agree most
(0.70–0.73); pure **Structure** models are the outliers (≈0.37–0.44 with
sequence/MSA types). Recovering the experiment's within-site order is hard (best
model ≈0.35): best are structure-aware+context models (AIDO.Protein-RAG,
VenusREM, ProSST-2048, ESM3, S3F); worst are autoregressive single-sequence LMs.

**HXK4 per-site disagreement (06–07, 11).** Sites where models most disagree on
within-site ranking overlap heavily between the two GCK assays (17/20 of the
top-20 shared) → disagreement is a property of the *position*, not the phenotype.
Pure **Structure** models are the biggest dissenters. Curated low-agreement
sites: T49, R63, W99, H156, Y214, Q287, R345, Y349, R377, V401, M402 (terminal
tails L2/D4/R5/L460–463 excluded).

**Functional context on GCK (08, 17).** In 1-D, disagreement is significantly
**higher at binding-site residues** than elsewhere (activity p=0.005, abundance
p=0.018). But in real 3-D (PDB 1V4S, glucose + activator; `17`) high-disagreement
residues are **not** closer to the ligand pockets (Spearman ≈ −0.04 to −0.06, n.s.)
and **not** spatially clustered (permutation p≈0.88) — so the 1-D "binding-site"
signal reflects sequence-annotation overlap, not a physical pocket neighbourhood.

**Multi-assay coordinate overlap (10).** Of 186 proteins, 24 have ≥2 assays; **23
have ≥1 assay pair sharing mutated residues** (all on one identical `target_seq`).
Exception: F7YBW8_MESOW. P53_HUMAN_Kotler_2018 dropped (different `target_seq`).
Many pairs deliberately span selection types (abundance vs activity vs binding).

**Abundance ≠ physicochemistry (12).** A hydrophobicity+charge baseline predicts
the abundance assay only weakly (Spearman 0.22, R²=0.08); controlling for those
features removes only ~13% of model performance (mean raw 0.37 → partial 0.32).
Being "more physicochemical" does not make a model better.

**Consensus predicts accuracy (13).** Among-model agreement at a site positively
predicts model-vs-experiment accuracy (Spearman 0.25 activity / 0.38 abundance);
accuracy roughly doubles from low- to high-agreement sites — but with a hard
ceiling (~0.3–0.4) and much scatter, so consensus is a useful, imperfect
confidence signal.

**"Winners" collapse (14, 16).** Model ranking-vs-experiment collapses to ~0 on
the top fitness quartile (the classic "poor prioritizer of winners"), in **both**
high- and low-agreement sites — agreement doesn't rescue it. Structure-only
models do **not** collapse worse (slightly better on the activity tail).
Redefining "winners" as best-substitution-*within-site* rather than
globally-most-fit roughly quadruples retained signal (top-25% ≈0.12–0.38),
showing the collapse is largely an artifact of the global definition mixing in
across-position signal.

**Generalization of the disagreement→function link (15).** Across all multi-assay
proteins, only the binding-site part generalizes: disagreement is elevated near
**binding/active-site residues** in 10/13 proteins (pooled p=5×10⁻¹³). The
**disease-variant** and **disordered-region** patterns do **not** generalize
(disordered regions trend the *opposite* way — lower disagreement — in P53 and
HSP82). Disease-annotation saturation (HXK4 412/465, P53 333/393 positions) makes
a combined "functional" flag misleading; the three signals must be tested apart.

**Disease-mutation anecdotes (from UniProt P35557).** W99 → **W99C / HHF3**
(gain-of-function, activating) — the assay agrees (W99 is an activating position),
but conservation-based models rank it backwards. R377 → **R377H / MODY2**
(loss-of-function) — non-functional in both assays.

## Caveats (carried across analyses)
- **Range restriction** attenuates every quartile/tail Spearman; read absolute
  tail values as "≈0", trust the *comparisons* (between splits, groups, types).
- **"Winners" within-site metric is pooled across sites**, so it reintroduces
  across-position ranking — an easier question than global gain-of-function.
- **`Structure` = only 3 models** in the family-best set (ESM-IF1, MIF,
  ProteinMPNN); those lines are noisy.
- **"Near binding" is sequence-proximity ±2**, not 3D contact.
- **Disease/disordered generalization** rests on few proteins (3 and 5) and on
  UniProt annotations; GCK's "disordered" T49/R63 were user-asserted, not
  UniProt-annotated (GCK has no UniProt disordered region).
- Per-site disagreement uses sites with **≥5 substitutions** (06/11; ≥3 in 03).
- **"Disagreement at functional sites" is an annotation-level, not a structural,
  claim.** The 1-D binding-site enrichment (08, p=0.005/0.018) does **not** survive
  in 3-D (`17`): high-disagreement residues are neither closer to the glucose/
  activator pockets (Spearman ≈ −0.05, n.s.) nor spatially clustered (permutation
  p≈0.88). So the signal reflects overlap with UniProt binding-site *annotations*,
  not a physical pocket neighbourhood — interpret the functional-site result
  cautiously and, where a structure exists, prefer a real 3-D distance test.
- **Candidate-protein `best_pdb` (18) optimizes coverage then resolution, not
  ligand presence.** For a ligand-distance analysis like `17`, pick a *liganded*
  structure yourself (e.g. GCK 1V4S with glucose+activator, not the apo 4DCH the
  heuristic returns); `top_xray_pdbs` lists alternates with coverage %.
- **`MSA_Neff_L` is not good for some assays — both extremes hurt MSA-dependent
  models.** Effective MSA depth per length varies enormously across the 217 assays
  (median 30.5, range 0 → 12,395). **Too shallow** (36/217 are `Low` /
  `MSA_Neff_L < 1`, e.g. CAR11 0.04, P53 0.39, viral RBDs ≈1): the alignment
  carries little position-specific signal, so MSA-based models (EVmutation, EVE,
  GEMME, MSA-Transformer) and Structure&MSA models are starved of information.
  **Too deep** (Neff_L in the thousands): very large/over-diverse alignments can
  fold in distant homologs and paralogs that dilute the couplings for *this*
  sequence, again degrading specificity — so the depth→MSA-model-quality relation
  is non-monotonic, best at an intermediate depth. Either way, any **model_type
  comparison involving MSA-dependent types is confounded by MSA depth** (02, 03,
  05, 15, 16) — check per-assay `MSA_Neff_L` (both tails) before attributing an
  MSA-model score to the model. Relevant low-Neff multi-assay proteins in `10`:
  CAR11 0.04, P53 0.39, NPC1 0.67, SPIKE 1.06.

## Outputs (`results/`)
Per-assay & taxa: `directionality_check.csv`, `spearman_long.csv`,
`spearman_model_by_taxon.csv`, `spearman_modeltype_by_taxon.csv`,
`kruskal_{model,modeltype}_across_taxa.csv`.
Within-site agreement: `persite_model_agreement_matrix.csv`,
`persite_modeltype_agreement_matrix.csv`, `persite_dms_agreement_*.csv`,
`persite_modeltype_agreement_by_selectiontype*.csv`, `family_best_selection.csv`.
HXK4 site work: `hxk4_*_site_disagreement.csv`, `site_W99_*_rank_matrix.csv`,
`gck_site_functional_annotation_*.csv`, `gck_disagreement_vs_function_summary.csv`,
`gck_3d_residue_distances.csv`, `gck_3d_summary.csv`.
Overlap: `assay_position_overlap_{proteins,pairs,selected_assays}.csv`.
Physchem / accuracy / collapse: `abundance_physchem_model_table.csv`,
`hxk4_agreement_vs_accuracy_*.csv`, `hxk4_quartile_by_agreement*.csv`,
`collapse_modeltype_withinsite.csv`.
Generalization: `generalization_functional_disagreement_{sites,per_protein}.csv`,
`generalization_candidate_proteins.csv` (24 multi-assay proteins + RCSB PDB
options, for choosing a protein/structure to extend the analysis — `18`).

Figures (regenerable, not tracked) land in `figures/`, with the HXK4 low-agreement
batch under `figures/site_disagreement/HXK4_HUMAN/`.

## Version control
Repo scoped to `analysis/` only; remote
`https://github.com/ElainW/PG1_additional_analysis.git`. Figures and large
glucokinase media assets (`*.gif/*.mp4/*.pse`) are git-ignored.
