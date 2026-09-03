# Planned analyses

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

**Ties to:** the multiple-conformation section above, `results/allosteric_map/`
(finished per-protein atlases),
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
