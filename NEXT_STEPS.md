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
