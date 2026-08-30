# Planned analyses

## Expand the allosteric map beyond Src

**Plan.** Generate the allosteric map for **HXK4_HUMAN, RASK_HUMAN, CP2C9_HUMAN,
PTEN_HUMAN, VKOR1_HUMAN** — in addition to **SRC_HUMAN**, already in
`results/allosteric_map/` (`SRC_HUMAN.tsv`). These proteins were selected because
(a) there are existing allostery-mapping papers to source/curate site annotations
from, (b) they have a large number of DMS mutations, and (c) their MSA has a
suitable Neff/L (`MSA_Neff_L`).

**Structure / conformation selection (PDB IDs).** Per protein, the conformational
states and the PDB structures collected so far (`?` = state exists but structure not
yet selected):

| Protein | Conformations | PDB IDs collected |
| --- | --- | --- |
| HXK4_HUMAN | 3 states: super-open, open, closed | super-open `1V4T`, open `?`, closed `1V4S` |
| RASK_HUMAN | ≥2 (need ≥2 structures) | active (RAF1+GTP-bound) `6VJJ`; second `?` |
| SRC_HUMAN | 3 states (only 2 collected) | inactive `2SRC`, active `1Y57`; third state `?` |
| CP2C9_HUMAN | — | `1R9O` |
| PTEN_HUMAN | — | `1D5R` |
| VKOR1_HUMAN | — | `6WVH` (Claude-selected, extrapolated from other literature structures) |

**Open gaps to resolve before mapping.** (1) HXK4 **open**-state PDB; (2) RASK
**second** structure (e.g. an inactive/GDP-bound state to pair with active `6VJJ`);
(3) SRC **third** conformational state (a third PDB beyond `2SRC`/`1Y57`). Note VKOR1
`6WVH` is an extrapolated pick and should be double-checked against the literature.

**Ties to:** `results/allosteric_map/` (per-protein atlases),
"Conservation vs allosteric / active sites (Src)" section below (the Src template
this generalizes), `28_msa_conservation.py` (`MSA_Neff_L` for MSA-suitability screen).

## Multiple-conformation confound at high-disagreement sites

**Scope.** For every mapped protein above with **more than one conformation** —
**HXK4_HUMAN** (3 states: super-open/open/closed), **RASK_HUMAN**, **SRC_HUMAN**
(3 states) — test whether
within-site model disagreement is driven by the models not accounting for
**multiple conformations**.

**Motivation / hypothesis.** These proteins have different active states associated
with different conformations. **None** of our structure-based models here considers
multiple conformations, so residues that behave differently across active vs
inactive states should be systematically mishandled by structure-only models.
The signature to look for: high-disagreement residues involved in **long-range
side-chain↔side-chain interactions** — salt bridges, pi-cation interactions,
side-chain–side-chain bonds, side-chain–backbone H-bonds — that are present in the
**active-only** conformation and in **swapping** residues (interactions that change
partner/state between conformations). Sequence models (or any model that ingests
sequence — sequence-only, or multimodal sequence+structure / sequence+structure+
functional-annotation; consider whether an MSA is required) are expected to
circumvent this problem because they are not tied to a single input conformation.

**Steps.**
1. Enumerate the candidate sites: residues in long-range side-chain↔side-chain
   interactions that are active-only or swapping across the conformations (per
   protein: HXK4 3-state, RASK, SRC 3-state).
2. For all such sites, check whether we see **some level** of within-site model
   disagreement — it need not be in the absolute-highest tier, but should be
   somewhat elevated relative to baseline.
3. Where such disagreement is present, test whether **sequence models or
   multimodal models that include sequence** (sequence+structure, or
   sequence+structure+functional-annotation) **outperform structure-based
   models** at those sites.

**Ties to:** `results/allosteric_map/` (per-protein conformation/interaction
annotations), `28_msa_conservation.py` + H2 findings (model-class within-site
accuracy by class), `16_collapse_modeltype_withinsite.py` (model-type structure),
`results/qc_directionality/directionality_check.csv` (`model_type` classes).

## Sequence/structure features as a trust signal (or direct predictor) for model calls

**Question.** Can sequence and structural features tell us **how much to trust** the
existing DL/ML models — or can we use the sequence/structural features **directly**
(instead of the models) to predict mutation effects, especially at low-agreement
sites?

**Structural features to quantify.**
- **Intramolecular interactions & dynamics:** long-range side-chain↔side-chain
  contacts, classified by conformational role — **active-only**, **inactive-only**,
  **swapping**, and **static** (see the Src kinase allosteric-landscape paper's
  quantification of structural features and contacts).
- **rSASA** (relative solvent-accessible surface area): which low-agreement residues
  are in the **core (rSASA < 0.25)** vs on the **surface**? Are the core ones
  involved in key interactions for proper folding?
- **Secondary-structure type** annotated from UniProt (helix, beta sheet, turn,
  loop) or from prior publications.
- **Distance from the catalytic site** and **distance to the nucleotide**.

**Sequence features to quantify.**
- Change in **net charge**, change in **hydrophobicity**.
- **Mutation type** — WT and mutant amino-acid identity.

**Analysis.** Relate these features to within-site model disagreement / agreement:
(a) do they explain where models are untrustworthy (trust signal), and (b) how far a
feature-only predictor gets on its own vs the DL/ML models.

**Ties to:** the two sections above (allosteric maps + multiple-conformation sites),
Src kinase allosteric-landscape paper (structural-feature/contact quantification),
`17_gck_3d_distances.py` (catalytic/ligand distance template),
`28_msa_conservation.py` (conservation feature, per-site disagreement).

## R377H (and model-undercalled GCK LoF) via the GCK allostery map

**Motivation.** Script `22` showed models systematically *under-call* LoF disease
variants (majority ranked in the tolerated top-half), R377H (MODY2) being the
anecdote: model within-site percentile ≈0.89 (predicted benign) despite being a
clinical loss-of-function, and R377 is **not** an annotated ATP/glucose-binding
residue → likely a functional/allosteric or kinetic LoF that conservation-based
models miss. We have an allostery map for GCK, so test this directly.

**Hypothesis.** The GCK LoF variants that models under-call are enriched for
**activity-specific (stability-independent), allosteric/distal** effects — high on
the activity–abundance residual map, and located away from the orthosteric site —
whereas the LoF variants models *do* catch are stability-disrupting / near the
active site.

**Allostery map source (corrected).** The activity–abundance residual (Liao &
Lehner 2026): activity below its LOESS expectation given abundance = function loss
independent of stability. **Their GCK residue map is from DIFFERENT DMS assays than
our pipeline** — the paper used ~8,255 GCK variants (abundance from VAMP-seq, one of
their nine VAMP-seq proteins; activity from a yeast-based functional assay), whereas
`O2/20250613` has Gersing_2022_activity (8,570) + Gersing_2023_abundance (8,396).
The differing counts show these are not the same processed datasets. So do NOT
assume our two Gersing assays reproduce their map. Options: (a) obtain the paper's
GCK residual map / source data (Supplementary Table 2 + their code), or (b) if we
want a map on OUR data, re-derive the activity–abundance residual from the two
Gersing assays ourselves and treat it as an independent map (label it as such, not
"the Liao & Lehner map"), first checking assay/variant overlap. Pair either with
1V4S orthosteric glucose/ATP distances from script `17`.

**Steps.**
1. Either (a) obtain the paper's GCK residual map, or (b) compute our own
   per-variant activity–abundance residual from the two Gersing assays (LOESS of
   activity vs abundance; residual = vertical deviation) — reproducing the paper's
   *method*, on *our* data, and labelling it as an independent map.
2. For each model-undercalled MODY2/LoF variant (from `22`), record: activity–
   abundance residual, min side-chain heavy-atom distance to glucose/ATP (1V4S),
   model within-site percentile, and abundance effect.
3. Test: do model-undercalled LoF variants have (a) larger activity–abundance
   residuals (more activity-specific) and (b) larger distance-to-orthosteric-site
   (more distal/allosteric) than model-caught LoF variants?
4. Place R377H on this map as the worked example; check whether it is a high-
   residual, distal (allosteric) variant.
5. Extend to the other model-undercalled LoF sets (PTEN, KCNE1) where paired
   activity+abundance assays exist, to see if "under-call = allosteric/kinetic LoF"
   generalizes.

**Ties to:** `22_disease_gof_lof.py` (undercalled LoF set), `17_gck_3d_distances.py`
(1V4S distances), `papers/pdb_selection_framework.md` (1V4S = orthosteric-ligand
structure), Liao & Lehner 2026 (activity–abundance residual method).

## Conservation vs allosteric / active sites (Src)

**Motivation.** We now have redundancy-reduced per-column conservation (`28`,
theta-reweighted) and a curated Src allosteric-site atlas
(`results/allosteric_map/SRC_HUMAN.tsv`, from Beltran, Naqvi, Faure & Lehner 2026).
Test whether evolutionary conservation distinguishes **active-site** from
**allosteric** positions, and whether within-site model disagreement relates to
conservation *differently* across these mechanistic classes — i.e. whether
conservation captures orthosteric/active-site importance but misses allosteric,
distal function (a candidate explanation for where MSA-based models fail).

**Data / join.** The atlas `position` is UniProt P12931 numbering = our `up_pos`
for `SRC_HUMAN`. Join to `results/MSA_conservation/msa_conservation_by_site.csv`
on (`UniProt_ID=SRC_HUMAN`, `up_pos`) for the three Src assays (Ahler_2019,
Nguyen_2022, Chakraborty_2023_binding-DAS_25uM). Atlas fields: `site_class`
(active_site=7, allosteric=26, allosteric_second_shell=15, other/blank),
`direction` (inhibitory / activating), `pockets`, `region`.

**Hypotheses.**
1. Active-site residues are the most conserved (highest weighted `conservation`
   and within-assay `conservation_percentile`); allosteric second-shell
   intermediate; bare allosteric sites *less* conserved despite being functional.
2. Activating (GoF, ddGa<0) sites are less conserved than inhibitory (LoF) sites.
3. High-disagreement Src sites (`26`, `27`; `per_protein_conservation/src_high_disagreement_conservation.tsv`)
   are enriched for **allosteric** rather than active-site class, and tend to be
   conserved-but-allosteric or variable — where sequence/MSA models lack the signal.

**Steps.**
1. Join the Src atlas to per-assay conservation and to the disagreement sites.
2. Compare `conservation` / `conservation_percentile` across `site_class`
   (Kruskal–Wallis + pairwise Mann–Whitney; active_site vs allosteric the key test).
3. Split by `direction` (inhibitory vs activating) and test conservation difference.
4. Cross with per-site disagreement: does Spearman(disagreement, conservation)
   differ by `site_class`? Are the labelled high-disagreement Src sites allosteric,
   and are they high- or low-conservation?
5. If a signal holds, frame it as: conservation encodes active-site constraint but
   under-represents allosteric/distal function → mechanistic account of MSA-model
   blind spots. Extend to other proteins if/when more allosteric maps are curated.

**Caveat.** Only Src has an allosteric atlas so far (~48 annotated positions, one
protein) — treat as a case study, not a generalization; small per-class n limits
power (7 active-site residues).

**Ties to:** `28_msa_conservation.py` (conservation),
`26_src_allosteric_sites.py` + `results/allosteric_map/SRC_HUMAN.tsv` (atlas),
`results/MSA_conservation/per_protein_conservation/src_high_disagreement_conservation.tsv` (labelled sites),
Beltran, Naqvi, Faure & Lehner 2026 (Src allosteric landscape).

## Mechanism of the per-assay disagreement ↔ conservation sign

**Motivation.** Script `28` shows the sign of Spearman(within-site disagreement,
conservation) is protein-specific (−0.59 to +0.42): in some proteins (P53, HSP82,
VKOR1) models disagree most at the **most-conserved** sites, in others (NPC1, PTEN,
RL40A, RPC1, CCDB) at the **least-conserved**. Why? Two hypotheses.

**H1 (DMS-score spread) — DONE, rejected. H2 (model-class performance) — DONE,
supported.** The full findings, numbers, and the threshold caveat now live in the
README (§"Site conservation vs disagreement", the two-hypotheses subsection);
outputs are `results/MSA_conservation/H1/H1_*` and `H2/H2_*` from script `34`. In brief:
H1's premise (conserved ⇒ wider within-site DMS spread ⇒ disagreement) fails — the
opposite holds, high-disagreement sites carry a *compressed* within-site DMS signal
(relfrac ≈46 vs ≈62 at the low-disagreement sites, robust to the conservation
split); **note** the "conserved-site baseline" is threshold-sensitive (≈68–69 at
50/50, ≈59–62 at 80/20, ≈46–52 at 90/10), so always state the split. H2 holds under
the strict split: at P53's conserved high-disagreement sites seq+MSA models beat
structure-only (0.144 vs 0.040), whereas at PTEN/HXK4's variable sites all classes
are uninformative (0.069 vs 0.035). What **remains** is the generalization beyond
P53/PTEN/HXK4 (the Steps below).

**Steps (remaining generalization work).**
1. Compute per-site, per-model within-site rank accuracy for all assays; summarise
   by `model_type` and by conservation tertile (extend the P53/PTEN probe to all 53).
2. Decompose the disagreement metric itself: is high disagreement driven by
   *between-class* divergence (MSA vs single-seq vs structure) or *within-class*
   spread? Correlate the between-class component with conservation per assay.
3. Test H2 directly: per assay, does (mean MSA-class accuracy − mean single-seq
   accuracy) *increase with conservation* more in positive-sign proteins than
   negative-sign ones? And is disagreement at low-conservation sites explained by
   uniformly low accuracy (all classes near 0)?
4. Relate the per-assay sign to MSA depth (`MSA_Neff_L`) and taxon — is
   "disagree-at-variable" the low-information / shallow-MSA regime?

**Ties to:** `28_msa_conservation.py` (conservation, per-assay ρ, high/low tails),
`results/qc_directionality/directionality_check.csv` (model_type classes),
`O2/20250613/*.csv` (per-model scores + `DMS_score`), `16_collapse_modeltype_withinsite.py`
(model-type within-site structure).

## Synonymous-variant baseline for the DMS percentile/spread plots (script 37)

**Motivation.** `37_conservation_dms_panels.py` draws, per assay (PTEN, CP2C9, HXK4,
P53), a 2x2 figure whose bottom row is residue-level: (3) residue-mean DMS-score
percentile and (4) per-residue percentile **spread** = max−min of the substitution-
level DMS percentiles, across the 4 categories (shared low-agreement, shared
high-agreement, highly conserved, all residues). The spread panel currently has **no
neutral reference**: the WT is a single value with no spread, and the processed
`O2/20250613/*.csv` score files are **missense-only** (0 synonymous variants,
`X###X`), so there is nothing to anchor "how wide is a no-effect residue's spread".

**Plan.** Obtain the **synonymous** variants for these assays (from the papers'
raw/supplementary data — Amorosi 2021, Mighell 2018, Matreyek 2021, Gersing
2022/2023, Kotler 2018, Giacomelli 2018) and use them as the neutral baseline:
- Panel 3: mark the synonymous-variant DMS distribution (its spread, not just the
  single WT line) on the percentile axis — the empirical "no amino-acid change" band.
- Panel 4: compute the synonymous variants' own percentile spread as the reference
  for how much per-residue divergence is expected under no functional effect
  (measurement noise), and overlay it.
- **Do NOT plot the "highly conserved" category on these synonymous-baseline plots**
  — synonymous (silent) variants have no amino-acid change, so a conservation-based
  residue category is not meaningful for them; restrict to shared low-agreement,
  shared high-agreement, and all residues (plus the synonymous baseline series).

**Implementation notes.** Keep synonymous variants **out** of the missense residue
means / percentile ranking (exclude `wt==mut`) so they don't skew the categories;
carry them as a separate baseline series. Decide encoding once the data is in
(`X###X` rows appended to the score CSVs vs a separate synonymous file).

**Ties to:** `37_conservation_dms_panels.py` (the merged 4-panel plots +
`dms_percentile_summary.csv`), `O2/20250613/*.csv` (current missense-only scores),
the "H1 (DMS-score spread)" subsection above (same residue-level DMS-spread question,
per-category).
