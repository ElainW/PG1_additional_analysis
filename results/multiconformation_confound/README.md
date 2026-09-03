# Multiple-conformation confound at high-disagreement sites

Script `57_multiconformation_confound.py`. Tests whether within-site model
disagreement and structure-vs-sequence model accuracy at the multi-conformation
proteins (HXK4/GCK, RASK/KRAS, SRC) are explained by structure-only models
ingesting a single conformation.

## Candidate sites (step 1)

Candidate = residue in a long-range side-chain<->side-chain contact classified
**active-only** or **swapping** between the two states (script 54). 353 candidate residues across the three proteins.

| protein | active | swapping | inactive | nonspecific |
|---|---:|---:|---:|---:|
| HXK4_HUMAN | 62 | 103 | 43 | 70 |
| RASK_HUMAN | 17 | 38 | 9 | 40 |
| SRC_HUMAN | 32 | 101 | 55 | 105 |

## Step 2 — disagreement elevation

Median within-site disagreement (z per assay), candidate vs rest (Mann-Whitney, H1 candidate greater):

| protein | n cand | n rest | med z cand | med z rest | p (cand>rest) | rank-biserial |
|---|---:|---:|---:|---:|---:|---:|
| HXK4_HUMAN | 330 | 594 | -0.201 | -0.072 | 0.934 | -0.060 |
| RASK_HUMAN | 108 | 250 | 0.143 | -0.348 | 0.000543 | 0.218 |
| SRC_HUMAN | 169 | 503 | 0.094 | -0.199 | 0.00151 | 0.152 |
| POOLED | 607 | 1347 | -0.050 | -0.181 | 0.0268 | 0.054 |

## Step 3 — structure-only vs sequence-ingesting models

Structure-only = model_type 'Structure' (MIF, ESM-IF1, ProteinMPNN); sequence-ingesting =
every other model_type (single-sequence, MSA, sequence+structure, sequence+
structure+function). Per site we take each group's mean within-site Spearman
(vs DMS) and the seq-minus-structure gap.

| protein | n cand | cand struct | cand seq | cand gap med | base gap med | paired p (seq>struct) | interaction p (gap cand>base) |
|---|---:|---:|---:|---:|---:|---:|---:|
| HXK4_HUMAN | 330 | 0.198 | 0.250 | 0.047 | 0.051 | 5.48e-09 | 0.531 |
| RASK_HUMAN | 108 | 0.271 | 0.332 | 0.084 | 0.126 | 3.65e-05 | 0.997 |
| SRC_HUMAN | 173 | 0.005 | 0.048 | 0.041 | 0.058 | 0.000686 | 0.897 |
| POOLED | 611 | 0.156 | 0.207 | 0.051 | 0.069 | 1.91e-14 | 0.979 |

### Does the structure deficit track disagreement? (confound signature)

Per-site Spearman(z_disagree, seq-struct gap) and the gap split at z=0. A real
single-conformation confound predicts the structure deficit GROWS with
disagreement (positive rho, hi>lo).

| protein | n sites | rho(disagree, gap) | p | gap med hi-z | gap med lo-z | p (hi>lo) |
|---|---:|---:|---:|---:|---:|---:|
| HXK4_HUMAN | 924 | -0.124 | 0.000166 | 0.027 | 0.071 | 0.997 |
| RASK_HUMAN | 358 | -0.253 | 1.22e-06 | 0.040 | 0.135 | 1 |
| SRC_HUMAN | 672 | -0.141 | 0.000242 | 0.025 | 0.090 | 0.999 |
| POOLED | 1954 | -0.154 | 6.62e-12 | 0.027 | 0.089 | 1 |

## Findings

**Step 2 — partial support.** Active/swapping candidate sites carry elevated
within-site disagreement in the two clean two-state switches — RASK/KRAS
(median z +0.14 vs -0.35, p=5e-4) and SRC (+0.09 vs -0.20, p=1.5e-3) — but
NOT in HXK4/GCK (p=0.93), whose whole-domain glucose closure rewires so many
residues that the contact classes do not separate. The per-class breakdown
(`disagreement_by_class.csv`) localises the signal to the **swapping** class
(RASK swapping median z +0.30, SRC +0.18); **active-only** sites are not
elevated in any protein.

**Step 3 — sequence models win, but not via the hypothesised mechanism.**
Sequence-ingesting models beat structure-only models at candidate sites in all
three proteins (pooled paired Wilcoxon p=2e-14; at SRC the structure-only
models have ~zero within-site skill). But the advantage is NOT candidate-
specific: the interaction test is null, and the per-site seq-minus-structure
gap *negatively* tracks disagreement (pooled rho=-0.15, p=7e-12; gap is
smaller, not larger, at high-disagreement sites). So the structure deficit is a
general property of these proteins, not a targeted blind spot at conformation-
changing residues — high disagreement there reflects genuine difficulty that
hobbles all model classes.

**Switch residues in the high-disagreement tail are functional (allosteric).**
Of the candidate residues that reach the strict high-disagreement tail (27),
**swapping = 6/9 (67%)** are functional and 5 of those 6 carry an *allosteric*
annotation (GCK 99, 287; KRAS 65, 68; SRC 387), whereas **active-only = 2/5
(40%)** are functional and only via disease, none allosteric/binding. Where
models disagree AND the residue rewires its contacts between states, it tends
to be a real allosteric/switch position (`switch_high_disagreement_functional.csv`,
`high_disagreement_candidate_counts.csv`).

## Output files (by step)

**Step 1 — candidate enumeration**
- `candidate_sites.csv` — one row per contacting residue of HXK4/RASK/SRC (UniProt pos, resname, `classification` = active/inactive/nonspecific/swapping from script 54, and `is_candidate` = class in {active, swapping}). The candidate list every later step draws from.

**Step 2 — disagreement elevation**
- `disagreement_by_class.csv` — per protein x contact-class summary of z_disagree (n, median, mean, fraction above the assay mean, fraction in the assay's top quartile). The descriptive table behind 'the signal is in the swapping class'.
- `disagreement_tests.csv` — the step-2 test: per protein + pooled, candidate vs rest medians and one-sided Mann-Whitney (H1 candidate > rest) with rank-biserial effect size.

**Step 3 — structure-only vs sequence-ingesting models**
- `model_class_persite_spearman.csv` — per (protein, assay, site): the mean within-site Spearman(model, DMS) of the structure-only models and of the sequence-ingesting models, their gap, the site's contact class, is_candidate, and z_disagree. The per-site table both step-3 tests aggregate.
- `model_class_tests.csv` — step 3a: per protein + pooled, candidate-site structure vs sequence means, the paired Wilcoxon (H1 seq > struct), and the interaction test (is the gap larger at candidate than baseline sites).
- `gap_vs_disagreement_tests.csv` — step 3b: per protein + pooled, Spearman(z_disagree, gap) and the gap split at z=0 (hi vs lo disagreement) — tests whether the structure deficit GROWS with disagreement (the confound signature). It does not.

**Step 4 — candidates in the strict high-disagreement tail** — step-1 candidates intersected with the elbow-thresholded tail (27) and labelled functional via 56 (needs those upstream outputs present):
- `high_disagreement_candidate_counts.csv` — per protein, active-only and swapping counts that reach the strict high-disagreement tail (script 27), plus the other (non-candidate) high-disagreement residues; unique residues, no double-count across assays.
- `switch_high_disagreement_functional.csv` — the active-only/swapping residues in the high-disagreement tail, each labelled functional / functional-source by joining script 56's `site_characterization_sites.csv` (allosteric / UniProt-binding / ClinVar / disease).
