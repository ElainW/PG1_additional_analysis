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

**Allostery map source.** The activity–abundance residual (Liao & Lehner 2026):
activity below its LOESS expectation given abundance = function loss independent of
stability. We can reproduce it directly from the two GCK assays we already have
(Gersing_2022 activity + Gersing_2023 abundance; ~8,255 variants), or use the
paper's map if available. Pair with 1V4S distances (orthosteric glucose/ATP sites)
from script `17`.

**Steps.**
1. Compute per-variant activity–abundance residual for GCK (LOESS of activity vs
   abundance; residual = vertical deviation), reproducing the paper's method.
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
