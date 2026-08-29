# Planned analyses

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
3. High-disagreement Src sites (`26`, `27`; `src_high_disagreement_conservation.tsv`)
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
`results/MSA_conservation/src_high_disagreement_conservation.tsv` (labelled sites),
Beltran, Naqvi, Faure & Lehner 2026 (Src allosteric landscape).

## Mechanism of the per-assay disagreement ↔ conservation sign

**Motivation.** Script `28` shows the sign of Spearman(within-site disagreement,
conservation) is protein-specific (−0.59 to +0.42): in some proteins (P53, HSP82,
VKOR1) models disagree most at the **most-conserved** sites, in others (NPC1, PTEN,
RL40A, RPC1, CCDB) at the **least-conserved**. Why? Two hypotheses.

**H1 — experimental (DMS-score spread).** At more-conserved positions the DMS
scores of the ~19 substitutions are spread over a larger fraction of the assay's
score range (more functional signal), which could drive within-site disagreement.
*Status: tested, rejected — and the true pattern is the opposite.* Within-assay
Spearman(conservation, per-position DMS-score std) = **+0.20** (positive in 41/48
assays) — the assay-wide premise holds — but Spearman(disagreement, DMS-spread) =
**−0.08** (positive in only 12/48). Zoomed into P53's three Giacomelli assays under the strict split (script `34`;
high = percentile ≥80, low = ≤20 or NaN;
`results/MSA_conservation/H1_P53_dms_relfrac_by_site_p80_20.csv`), with each site's
DMS *relative fraction* = max−min of its substitutions' DMS-percentiles within the
assay: the **high-disagreement** selected sites are **compressed** (mean relfrac
≈47) vs the **low-disagreement** sites (≈60–63) and the assay-wide high-conservation
baseline (≈59–62); robust at 90/10 (≈45 vs ≈60). Mechanism: at a critical conserved
residue nearly every substitution is similarly deleterious → narrow within-site DMS
range → the fine within-site order is noise → models disagree. Disagreement tracks
*compressed* within-site signal, not the larger spread H1 posited. (ProteinGym
assays are pre-selected for a decent overall DMS range, so this is a per-site, not
per-assay, compression.)

**H2 — model-class performance.** Models disagree most at conserved sites *when the
MSA / Structure&MSA classes do much better there than sequence/structure-only
models* (a real accuracy gap → divergent rankings), and disagree most at
low-conservation sites *when all classes do poorly* (little information). Classes
per `results/qc_directionality/directionality_check.csv:model_type`.
*Status: supported under the strict conservation split (script `34`; high =
percentile ≥80, low = ≤20 or NaN; stable at 90/10). Per-site Spearman by class,
`results/MSA_conservation/H2_persite_spearman_by_class_p80_20.csv`.*
- **Disagree-at-conserved (P53), high-conservation shared sites — a real class
  gap.** sequence+MSA per-site Spearman **0.144** vs structure-only **0.040**
  (seq+MSA > structure in 18/27 site-assay pairs). Conserved high-disagreement
  sites carry signal that sequence/MSA models partly read and structure-only
  models don't → the classes genuinely diverge.
- **Disagree-at-variable (PTEN, HXK4), low-conservation shared sites — everyone
  loses.** sequence+MSA **0.069** vs structure-only **0.035**, ~78% of sites below
  0.2 (HXK4 W99, percentile 13.8, the archetype). No class is informative.
- So the two regimes are mechanistically distinct: disagree-at-conserved =
  compressed within-site signal that sequence/MSA models exploit unevenly vs
  structure-only; disagree-at-variable = a low-information regime for every class.
  Threshold sensitivity (80/20 vs 90/10) and the full per-assay generalization
  (beyond P53/PTEN/HXK4) remain to complete.

**Steps.**
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
