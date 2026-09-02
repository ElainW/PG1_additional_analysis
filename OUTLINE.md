# Paper outline — ICBINB-BIO 2026 submission

**Working title:** *Good filters, poor optimizers: why ProteinGym leaderboard
rankings don't tell you which model to trust for protein optimization*

**Venue.** ICBINB-BIO 2026 (negative / unexpected-failure results in AI for
biology). Full paper **8 pages** (excl. refs/appendix), tiny paper 4 pages;
workshop LaTeX template; OpenReview, **double-blind**. Required elements:
(1) biological problem + metrics, (2) modeling approach + theoretical basis,
(3) quantified negative/unexpected results, (4) rigorous failure-cause
diagnostics + boundary conditions + actionable takeaways. Bonus awards for
surprising negative results and pedagogical value.

**Target:** full 8-page paper (the venue rewards depth of failure analysis).

**Double-blind housekeeping:** don't print the repo URL or author-identifying
paths inside the PDF; if a reviewer code link is needed, use an anonymized
supplementary link.

---

## Spine (repeat in abstract, end of §1, and §6)

> When you must pick **one** ProteinGym model to score mutations for optimizing a
> *new* protein — no site-saturation DMS, no solved structure, no time to
> fine-tune — the standard recipe (highest **Spearman** or **NDCG@10** in the
> assay category) barely predicts which model is actually useful, because it
> measures **global** rank correlation while optimization lives in **within-site**
> ranking of substitutions at **functionally important** positions.

Spinner et al. showed that a battery of *in silico* sequence, structure, and
biochemical/biophysical descriptors do not correlate well with model performance
or enzyme function. We take a different route. Empirically we only broaden their
setting — more protein families and more of the high-performing fitness predictors
benchmarked in ProteinGym (217 DMS assays × 97 zero-shot models) — rather than
claiming a large generalization. And crucially we do **not** use those in-silico
descriptors: instead we propose **new sequence- and structure-based metrics** that
explain *where and why* the models fail (their failure modes) and thereby inform
what to do instead.

---

## Figure budget (main text 3–4; rest to appendix)

> **Main-text figures:** **F0** (method flow), **F3** (winners collapse), **F6**
> (co-localization 2-panel — *may instead go to appendix as a table if space is
> tight*), and **F20** (mechanism panel assembled from selected F10 + F12 panels).
> That is 3–4 depending on F6. **4th-slot candidates** if a slot frees up: **F1**
> (intro schematic), **F16** (functional-enrichment stakes — needs regeneration),
> **F9** (DMS-range explainer), **F11** (conservation null). Everything else →
> appendix. F0b may fold into F0 (one method panel) or drop to appendix.

### Confirmed main-text figures (in main-text order)

| # | Figure | Section | Source |
|---|--------|---------|--------|
| F0 | **Within-site disagreement computation flow** (rank ~19 subs per model → pairwise Spearman every model pair → 1 − mean = disagreement) — the basis of all analyses | §2 | `render_within_site_flow.py` → `figures/within_site_disagreement_flow.pdf` (exists) |
| F3 | **Winners collapse by model_type — pooled across all 216 assays** (each assay equal weight ⇒ weighted by #assays per category), global vs within-site split; ρ vs fitness subset all→top-50/25/10%. Global panel collapses to ≈0 on the top decile; within-site panel far gentler (retains ≈0.15–0.27); Structure in bold (mid-pack within-site, among the least-bad on the global tail) | §3.3 | `16` `collapse_modeltype_withinsite` (exists; fig `figures/persite_analysis/fig_collapse_modeltype_pooled.png`) |
| F6 | **Co-localization forest — 2-panel excerpt:** P53_HUMAN (obs overlap ≫ null band, strongly significant) vs VKOR1_HUMAN (obs ≈ null, n.s.) — the co-localization contrast in one figure. *Space permitting; if tight, move to appendix and present the co-localization stats as a compact table instead.* | §4.1 | `33` `fig_colocalization_forest.png` (excerpt these two panels) |
| F20 | **Mechanism panel — selected from F10 + F12:** within-site DMS relfrac (H1) + per-site class accuracy by conservation regime (H2) — compressed within-site signal at conserved high-disagreement sites that seq+MSA models read but structure-only don't. Full sets stay in appendix (F10, F12) | §4.3–§4.4 | `34` H1/H2 + `figures/h2_persite_spearman_density/` (select panels — see to-do) |

### Confirmed appendix figures (in section order)

| # | Figure | Section | Source |
|---|--------|---------|--------|
| F18 | Within-site rank agreement between **model types** (MSA↔Structure&MSA highest ≈0.70–0.73; pure Structure the outliers ≈0.37–0.44) | §2 | `03` `persite_modeltype_agreement_matrix.csv` (exists as heatmap) |
| F19 | Within-site rank agreement matrix, **23-model subset** (labelled per-model heatmap) | §2 | `04` `figures/persite_analysis/Fig5_subset.pdf` (exists) |
| F17 | Per-protein disagreement structure map exemplar (GCK) with binding/disease/allostery overlay | §3.1 | `23` `fig_<gene>_structure_map.png`; allostery `26` |
| F2 | **Aggregate** per-site (within-site) Spearman vs overall Spearman, across assays × models — within-site sits far below the diagonal (much harder than global) | §3.2 | `figures/persite_vs_global/aggregate_persite_vs_global.pdf` (exists) |
| F2b | **Example within-site breakdown outlier:** P53_HUMAN_Giacomelli_2018_WT_Nutlin — overall ρ≈0.4–0.5 yet per-site Spearman ≈0 for *all* models (Pearson 0.02, ρ=−0.17, n=23) | §3.2 | `figures/persite_vs_global/P53_HUMAN_Giacomelli_2018_WT_Nutlin.*` (exists) |
| F22 | **Winners-collapse by model_type — high-disagreement sites vs all ≥5-sub sites** (50 multi-assay assays with labelled high-disagreement sites; 45 family-best). Adds the high-vs-all contrast that F3 (pooled, all sites) doesn't carry. 2 panels (global / within-site split), model_type lines. **Global panel** shows both site sets (solid = high-disagreement, dashed = all): both collapse to ≈0 on the top decile (high starts lower, ~0.34 vs 0.42; Structure holds up best on the global tail). **Within-site panel** shows **all sites only** (dashed) — deliberately omits the high-disagreement lines because there are too few high-disagreement sites for a reliable within-site quartile pool; the all-sites within-site collapse is gentle (0.42→0.19). *Caveat: within-site split is pooled across positions (reintroduces across-position signal), not the per-site accuracy of §4.2/F21.* | §3.3 | `38` `figures/persite_analysis/fig_collapse_modeltype_highdisagreement.png` (exists) |
| F8 | **HXK4 winners collapse, high- vs low-agreement sites** (activity/abundance; 23-model subset): ρ vs fitness subset all→top-50/25/10% — collapses to ≈0 on the high-fitness tail in **both** agreement classes (agreement doesn't rescue the collapse). *Complements F3 rather than duplicating it: different dataset (HXK4 only, 23-model subset) and a different axis (high- vs low-**agreement** sites), where F3 is model_type × global-vs-within-site pooled across assays. Decision: F3 stays the main-text winners-collapse figure, F8 stays in the appendix.* | §3.3 | `14` `hxk4_quartile_by_agreement*` (exists) |
| F21 | **Both metric-selected models degrade at high-disagreement sites (H2):** per-site accuracy of the category-top-Spearman and category-top-NDCG@10 models at assay-specific high-disagreement sites (2 panels; whole-assay ρ≈0.36 → per-site ρ≈0.09) | §4.2 | `34` `figures/h2_persite_spearman_density/fig_assayspecific_category_best_{spearman,ndcg}_panel.png` (exist) |
| F4 | **fig9 — do better-ranked models win at high- vs low-disagreement sites** (45 family-best): within-site Spearman vs category-rank tier (Top1/3/5/10/All) × {category Spearman, NDCG@10} × {high, low} disagreement, by assay category. Tiers ≈flat (Top1≈All) ⇒ category rank barely predicts within-site accuracy | §3.4 & §4.2 | `29` `fig9_*` (exists). **Dataset:** 45 family-best models; rank from full per-category leaderboard (all category assays); within-site sites from `27`/`30` (52-assay / 23-multi-assay-protein set, no Stability); pooled over all high-disagreement sites (A/B not split) |
| F5 | Spearman vs NDCG@10 rank–rank scatter — broadly correlated across the full ranking, but the rank-1 model differs by metric in all 5 categories (see T2) | §3.5 | `29` leaderboard |
| F6b | Full **high**-disagreement co-localization forest, all 20 proteins (obs vs 95% null band, red=BH q<0.05) | §4.1 | `33` `fig_colocalization_forest.png` (exists) |
| F6c | Full **low**-disagreement co-localization forest, all 20 proteins (control — the easy tail co-localizes too, but less) | §4.1 | `33` `fig_colocalization_forest_low.png` (exists) |
| F7 | Co-localization summary: pooled fraction + bootstrap CI, high vs low | §4.1 | `32` `colocalization_pooled_ci.csv` |
| F10 | Per-site within-site DMS relfrac, high- vs low-disagreement sites | §4.3 | `34` H1 relfrac tables |
| F12 | Per-site class accuracy (seq/MSA vs structure) by conservation regime — P53 as case, non-generalization | §4.4 | `34` H2 tables |
| F13 | Physchem baseline + partial-correlation control (abundance) | §4.5 | `12` `fig_abundance_physchem.png` |

### Others — optional 4th-slot candidates & conditional (in section order)

| # | Figure | Section | Status | Source |
|---|--------|---------|--------|--------|
| F1 | Schematic: optimization-as-scoring-function + global-vs-within-site axes | §1–§2 | 4th-slot candidate | new (schematic; only new one to draw) |
| F0b | GCK per-site within-site disagreement **distribution**, both assays, high/low tails + 0.5 boundary + median — the metric on a real protein | §2 | fold into F0, else appx | `31` `figures/site_disagreement/distribution/HXK4_HUMAN.*` (exists) |
| F16 | Functional enrichment: disagreement near binding/active/**allosteric** + disease residues, per protein | §3.1 | 4th-slot candidate — **REGENERATE** | `15` re-run with updated ClinVar/UniProt annotations (`23`/`25`) + allosteric map (`26`, expanding); current 10/13 p=5e-13 render is stale |
| F9 | Within-site DMS range & |mean−WT| vs disagreement, pooled across assays | §4.3 | 4th-slot candidate | new (from `O2/20250613` + disagreement tables) |
| F11 | Conservation null: pooled ρ(disagreement, conservation) unweighted vs weighted; per-assay ρ spread (−0.59…+0.42) | §4.4 | 4th-slot candidate | `28` `conservation_disagreement_spearman_by_assay.csv` |

**Regenerability / work needed:** most figures already exist on disk
(F0/F0b/F18/F19 and the co-localization/collapse/fig9 renders — `31`/`03`/`04`/
`33`/`16`/`29`). Still to do (**later, not now**): **F1** (new schematic),
**F9** (new cut), **F16** (re-run `15` with updated ClinVar/UniProt + allosteric
annotations), **F20** (assemble selected panels from F10/F12 +
`h2_persite_spearman_density/`).

---

## Planned tables

| # | Table | Section | Source |
|---|-------|---------|--------|
| T1 | **All assays selected for this study** — 52 DMS substitution assays across 23 multi-assay proteins (the controlled-comparison set of §2). Columns: `UniProt_ID`, `DMS_id`, `coarse_selection_type`, `taxon`, `selected_pdb`, `MSA_Neff_L`, `MSA_perc_cov`, `seq_len`, `DMS_wt_score`. Spans 4 assay categories (Activity / Binding / Expression / OrganismalFitness), 4 taxa, and the full Neff/L range (0.0 → 337 — what the §5 MSA-depth boundary draws on). *Coverage caveats:* `selected_pdb` populated only for the 6 proteins with a curated experimental structure (HXK4/RASK/CP2C9/PTEN/SRC/VKOR1 — per-state for switches); `DMS_wt_score` populated for the 11 assays present in `DMS_wt_score.csv` (P53/PTEN/RASK/CP2C9/HXK4). **Footnote — deliberate exclusion:** `P53_HUMAN_Kotler_2018` is *not* in this set because it sits on a **different `target_seq`** than the three P53 Giacomelli assays, so its positions are not in the same coordinate frame; the selection (`10`) requires every kept assay for a protein to share one identical `target_seq` (enforced by assertion), which the paired within-site comparison depends on. **Kotler is now excluded from *all* cross-assay analyses** (the R72P `target_seq` difference means it is not on the shared P53 coordinate frame) via `common.EXCLUDE_DMS`, filtered at the master-site-table source (`15`): the co-localization / shared high-disagreement (`32`/`33`), MSA-conservation (`28`, `34`–`45`), and P53 structure-map (`23`) panels are all built from the three Giacomelli assays only. It survives only in **non-pooled, global/single-assay** contexts (e.g. the 217-assay Spearman leaderboard; its own `DMS_wt_score` of −2.4) where it is never mixed with the Giacomelli coordinate frame | appendix | `results/assay_selection/selected_assays_annotated.csv` (exists; built by `36_selected_assay_table.py` joining `assay_position_overlap_selected_assays.csv` + `generalization_candidate_proteins.csv` + PDB README selected structures + `data/DMS_wt_score.csv`) |
| T2 | **Per-category "recipe" leaderboard — top models.** For each of the 5 `coarse_selection_type` categories, the rank-1 model by **mean Spearman** and the rank-1 by **mean NDCG@10** (45 family-best models; 10 rows), annotated with **model class** (`model_type`). Columns: `coarse_selection_type`, `top_by`, `model_type`, `family`, `model`, `n_assays`, `mean_spearman`, `mean_ndcg10`, `spearman_rank`, `ndcg_rank`. **Key result: the top-Spearman and top-NDCG model differ in all 5/5 categories** — e.g. Expression top-Spearman VenusREM is only NDCG-rank **28**; Activity top-Spearman AIDO.Protein-RAG-16B is NDCG-rank **11** — so which model the recipe "selects" depends on the metric (sharpens §3.4/§3.5: the two metrics broadly correlate across the full ranking yet pick different winners). Top-Spearman picks are almost all **Structure & MSA** (AIDO/VenusREM; Stability→ProSST, Single seq & Structure); top-NDCG picks are mixed (ESM3 seq+struct+func, Protriever/TranceptEVE MSA). Rank is over **all** assays per category (incl. Stability, 66 assays). | main | `results/persite_analysis/family_best_category_leaderboard_top.csv` (top cut, `37`) ← full 225-row `family_best_category_leaderboard.csv` (appendix, `29`) |
| T3 | **The 23-model representative subset** — the hand-curated panel used in the within-site / labelled-heatmap figures (F0b/F18/F19/F5-family, W99, `09`/`12`/`13`/`14`). Columns: `model`, `model_type`. Composition: MSA 8, Structure & MSA 5, Single sequence 4, Single seq & Structure 3, Structure 2, seq+struct+func 1. Documents which models the 23-model figures show (distinct from the 45 family-best set used for leaderboard/collapse work). | appendix | `results/persite_analysis/subset_23_models.csv` (exists; `39`, from `04_subset_figs.SUBSET` + directionality QC) |

---

## 1. Introduction & biological problem  *(element #1)*

- **Goal:** use the ProteinGym (PG) models for protein optimization — sample amino
  acids at a given position from a PG model's *pseudo-likelihood*, whether within
  a Gibbs-sampling framework (Spinner et al.) or other
  search-based algorithms (McCarter et al., arXiv 2026), to build up multi-mutant
  sequences. Cheaper than one-mutation-at-a-time lab loops, respects epistasis,
  avoids stacking deleterious mutations; ML-designed multi-mutants tend to work.
- **Realistic setting:** a *new* protein, no DMS / structure / fine-tuning budget
  → pick a model off the shelf and trust it.
- **Recipe under test:** ProteinGym → 97 models → pick best by Spearman or
  NDCG@10 in the relevant assay category. *Does that solve model selection?*
- **Thesis (the "I can't believe it's not better"):** No — the selection metric
  is computed on a quantity (global rank correlation) that isn't the quantity
  optimization needs (within-site ranking of substitutions at functionally
  important positions).
- **The broader premise — ProteinGym is treated as the gold standard, but the
  leaderboard doesn't transfer to design.** The field defaults to "higher
  ProteinGym Spearman/NDCG@10 ⇒ better model," yet **higher ProteinGym performance
  does not necessarily correlate with a better success rate in actual protein
  design/optimization.** We argue this gap stems from **biases baked into the
  current evaluation metrics** — they reward the wrong axis (global vs within-site,
  §2) and are scored against DMS labels of **limited experimental resolution**,
  where measurement noise/detection bias caps the achievable correlation
  (§3.2, §4.3, §5). So the disconnect is a property of the *evaluation*, not just
  of any one model.
- **Position vs Spinner et al.:** cite as precedent, but we don't over-generalize
  from it — we broaden the empirical setting (more protein families, more
  high-performing PG-benchmarked models) and, distinctively, focus on **proteins
  measured by multiple types of DMS assays** (the same residues scored under
  different phenotypes — abundance vs activity vs binding), which lets us separate
  model error from phenotype. And rather than reusing their in-silico descriptors,
  we propose **new sequence/structure metrics** to explain the failure modes.
- **Stakes (why the hard sites are the ones that matter):** the sites where model
  selection fails — the high-disagreement sites — are **enriched for functional
  importance** (binding / active / allosteric residues; §3.1). In a real
  optimization task you are *working on exactly these positions*, so the recipe
  fails precisely where the task lives. This is what makes the negative result
  consequential rather than academic.
- **Goal / contribution (what a better evaluation buys):** the constructive aim
  behind the negative result is to **mitigate — ideally resolve — the biases in
  current computational-model evaluation**, which in turn lets the field:
  1. **think critically about the experimental assays themselves** (e.g. DMS) —
     *how much to spend* on them and *how reliable* they are for training and
     evaluating computational predictors (limited resolution ⇒ a noisy yardstick);
  2. **guide method choice by goal** — a suite of better, task-aligned metrics
     (within-site ranking, winners, site-level triage) so users pick a model for
     *their* optimization objective rather than for leaderboard rank.
- **Figs:** F1 (schematic + two axes; appx / 4th-slot candidate).

## 2. Setup, metrics & approach  *(elements #1–#2)*

- Data: 217 DMS substitution assays × 97 zero-shot models; directionality
  normalization (scores already higher = fitter; verified in `00`).
- Two evaluation axes: **global** (Spearman / NDCG@10, what the leaderboard uses)
  vs **within-site** (rank the ~19 substitutions at a fixed position — what "which
  mutation here?" asks). Optimization lives on the second.
- **Within-site disagreement = 1 − mean pairwise Spearman** across models
  (decision rule: agree → any model works; disagree → need to know which is right).
  Computed per site: rank each model's scores over the ~19 substitutions → Spearman
  between every model pair → mean pairwise Spearman = consensus → disagreement =
  1 − that. **This pipeline is the basis of every analysis in the paper (F0).**
- Controlled comparison: multi-assay proteins, shared coordinates, different
  `coarse_selection_type` → separates *model wrong* from *model measures a
  different phenotype*. The selected set is **52 assays across 23 proteins**
  (tabulated in **T1**, appendix).
- **Figs:** F0 (within-site disagreement computation flow — the foundational
  method figure; optionally combined with F0b, the GCK disagreement distribution,
  into one method panel), F1 (shared with §1). Appendix support: F18 (model-type
  within-site agreement — Structure models are the outliers), F19 (23-model
  within-site agreement matrix).

## 3. The recipe fails  *(element #3 — core, ordered as an argument)*

**3.1 Why these sites matter: high-disagreement sites are enriched for functional
importance.** The stakes, established before the failures — the reason the negative
result is consequential. In an actual optimization task you work on the functionally
important residues (binding / active / allosteric), and those are exactly where
models disagree, so the recipe fails where the task lives:
- **Binding/active-site enrichment generalizes.** Across multi-assay proteins,
  within-site disagreement is elevated near **binding/active-site residues** in
  **10/13 proteins** (pooled p=**5×10⁻¹³**, `15`); on GCK specifically the 1-D
  binding-site enrichment holds in both assays (activity p=0.005, abundance
  p=0.018, `08`).
- **Allostery (being expanded).** Curated allosteric-site atlases (currently
  SRC_HUMAN, `26`; HXK4/RASK/CP2C9/PTEN/VKOR1 in progress per NEXT_STEPS) let us
  test whether disagreement also concentrates at **allosteric/distal** functional
  sites that conservation-based models miss. *[User to add.]*
- **Honest scope caveats (state up front):**
  - Functional annotations (ClinVar P/LP + UniProt binding/GoF supplement) are
    high-quality for only a **few well-characterized human proteins with decent
    MSAs** — GCK, PTEN, TP53, SRC, CYP2C9, KRAS, etc. — so the enrichment is a
    curated-protein result, not a proteome-wide claim.
  - **Not all functional signals generalize:** the disease-variant and
    disordered-region patterns do **not** (disordered trends *opposite* — lower
    disagreement — in P53/HSP82); disease-annotation saturation (HXK4 412/465,
    P53 333/393) makes a combined "functional" flag misleading → test binding,
    disease, disorder **apart**.
  - The binding enrichment is **annotation-level, not 3-D**: on GCK it does *not*
    survive as a physical-pocket-proximity effect in real 3-D (`17`), though a few
    hotspots (R63, Y214, Q287) do contact ligand. Frame as "enriched at
    functionally *annotated* residues."
- **Figs:** F16 (appx, 4th-slot candidate — cross-protein binding enrichment),
  F17 (appx — GCK structure-map exemplar with allostery overlay).

**3.2 Global rank correlation mostly tracks within-site — until it doesn't.**
For most models on most assays, global and within-site Spearman are correlated,
so the leaderboard isn't uniformly wrong. The failure is **assay-specific and
catastrophic**: on some assays within-site signal collapses even for globally
strong models — e.g. **P53_HUMAN Giacomelli** (within-site Spearman breaks down
completely). You can't tell from a global number whether you're on a well-behaved
or a broken-within-site assay, and the new-protein setting gives no way to check.
- **Figs:** F2 (aggregate global-vs-within-site, appx) + F2b (P53 WT_Nutlin
  outlier where within-site collapses for all models, appx).

**3.3 The winners collapse.** Ranking power vanishes on the top fitness quartile
(ESM3 0.49→0.04→−0.01 at top-25%/10%; reproduces the CP2C9 cytochrome case;
generalizes across 209 assays). Models are good *filters of deleterious variants,
poor prioritizers of winners* — and optimization is exactly the winner regime.
**Two things F3 makes visible** (pooled over all 216 assays): (i) the collapse is a
property of the **global** fitness definition — under a **within-site** winners
split it is far gentler (retains ≈0.15–0.27 vs ≈0 globally), so much of the collapse
is the global definition mixing in across-position signal (cf. §4); (ii) it holds
across **model types**, and Structure-only models (bold) do **not** collapse worse
(among the least-bad on the global tail).
- **Figs:** F3 (main — winners collapse by model_type, global vs within-site,
  pooled across all assays); F22 (appx — the high-disagreement-vs-all-sites contrast
  on the multi-assay subset); F8 (appx — HXK4 collapse, high- vs low-agreement
  sites: agreement doesn't rescue).

**3.4 Headline: category-level selection barely helps.** Better-ranked models win
at high-disagreement sites only weakly (ρ(rank, within-site) ≈ **−0.05 to
−0.15**); the **per-site oracle (0.56–0.73) dwarfs every category tier**. Picking
the top model by Spearman/NDCG in the assay category does **not** solve it;
site-specific choice would help far more.
- **Figs:** F4 (fig9 — category-rank tiers ≈flat at high- vs low-disagreement
  sites; appx). The per-site oracle number is quoted in text.

**3.5 Spearman and NDCG@10 disagree on the winner — and neither is the right
axis.** Across the full per-category ranking the two metrics broadly correlate, but
they select **different rank-1 models in all 5/5 categories** (T2; e.g. Expression's
top-Spearman VenusREM is NDCG-rank **28**, Activity's top-Spearman AIDO is NDCG-rank
**11**). So the recipe is ambiguous at the selection step — "best by Spearman" and
"best by NDCG@10" name different models — and, more fundamentally, **neither metric
captures within-site or winner ranking** (§3.2–3.4). Whichever you pick, you're
selecting on the wrong axis. (Top-Spearman picks are almost all Structure & MSA;
top-NDCG picks are mixed — see T2.)
- **Figs:** F5 (appx — rank–rank scatter: broadly correlated but top picks differ).
  **Table:** T2 (top Spearman/NDCG model per category, with model class).

## 4. Failure-cause diagnostics  *(element #4)*

Why do models disagree, and is the disagreement exploitable for selection? We
established *why we care* in §3.1 (these sites are functionally important). Here we
split high-disagreement sites into two categories, then ask what property (if any)
tells you which model to trust.

**4.1 Category A — shared high-disagreement sites (recur across ≥2 assays).**
*(Also establishes disagreement is site-intrinsic, not phenotypic noise — the
recurrence across a protein's assays is the evidence.)*
Significantly enriched: **131/238 (55%)** union sites shared by ≥2 assays
(cluster-bootstrap 95% CI 41–70%; 18/20 proteins BH q<0.05; large fold-enrichment,
e.g. P53 obs 36 vs null [0,4]); holds at lowest Neff/L → not a shallow-MSA
artifact. **Control:** low-disagreement tail also recurs but **less** — 110/256
(43%), CI 34–54%, 16/20 significant → hard sites recur *more* (not "any extreme
repeats"). Cause: intrinsically ill-posed positions; disagreement here is *not*
selectable-away.
- **Figs:** F6 (main — 2-panel excerpt, P53 significant vs VKOR1 n.s.); F6b (appx
  — full high-disagreement forest, 20 proteins), F6c (appx — full low-disagreement
  forest control), F7 (appx — pooled fraction + bootstrap CI).

**4.2 Category B — assay-specific high-disagreement sites (one assay only).**
Two candidate causes: (i) genuine model-performance spread, (ii) assay-specific
measurement issues (→ §5 DMS caveat). **Spine result:** where (i) holds, *which*
model wins is still poorly predicted by leaderboard rank (pooled ρ ≈ −0.05…−0.15,
`29`); per-site oracle far exceeds any tier → the recipe can't identify the right
chooser even where model choice matters.
- **Both metric-selected top models degrade here (tested in H2).** The
  category-top-**Spearman** *and* category-top-**NDCG@10** models both fall from a
  whole-assay median ρ ≈ **0.36** to a per-site median ρ ≈ **0.09** at these
  assay-specific high-disagreement sites (`34` / `results/MSA_conservation/H2/
  h2_assayspecific_persite_spearman_by_site{,_ndcg}.csv`, baseline in
  `h2_assayspecific_rank_tiers.csv`) — so selecting by either metric doesn't
  protect you at exactly the sites that matter. Consistent with fig9/F4's
  high-vs-low-disagreement tiers.
- **Caveat / to-run:** `29` is **pooled over all high-disagreement sites, not
  split A vs B**. The sharper claim (rank↔usefulness fails *specifically at
  assay-specific sites*) needs a re-run of `29` restricted to Category-B sites
  (partition by the `32` co-localization output). Flag as planned; don't
  over-claim from the pooled number.
- **Figs:** F4 (fig9 — the ≈flat rank tiers at high-disagreement sites; appx,
  shared with §3.4); F21 (appx — H2 per-site accuracy of the Spearman- and
  NDCG-selected models at assay-specific high-disagreement sites).

**4.3 What intrinsic site property *does* explain assay-specific disagreement:
the DMS-score structure at the position.** Across the many assays with
assay-specific high-disagreement sites (not just P53), models disagree where there
is little within-site signal to rank:
- **Tolerated sites** — substitutions' DMS scores sit close to the **WT DMS
  score**; mutations barely move fitness, so the fine ordering is noise and models
  diverge.
- **Hard-to-distinguish sites** — the **range** of DMS scores across the site's
  substitutions is **smaller than the assay average** (compressed within-site
  spread; cf. H1 relfrac ≈46 at high- vs ≈62 at low-disagreement sites). Even
  uniformly deleterious sites count: substitutions are near-indistinguishable, so
  ranking is ill-posed.
- **Takeaway:** disagreement at these sites signals the *measurement's resolving
  power at that position*, not (mainly) model quality — which is why no model and
  no leaderboard rank reliably wins there (§4.2). Ties to §5 DMS caveat: low
  within-site range is where assay noise dominates.
- **Figs:** F20 (main — selected relfrac panels, shared with §4.4); F9 (appx,
  4th-slot candidate — range & |mean−WT| vs disagreement, pooled), F10 (appx —
  full relfrac high vs low).

**4.4 Conservation does *not* cleanly explain who wins (null result).** We tested
whether MSA-based, redundancy-reduced conservation (`28`) predicts where
sequence/MSA models beat structure-only models. It does not hold up:
- Global link near-zero after redundancy reduction (ρ 0.18 → **0.06**; raw signal
  a homolog-redundancy artifact); per-assay sign spans **−0.59 to +0.42**.
- **No clear signal that high conservation by itself lets sequence/MSA models rank
  better within-site.** The P53 class gap (seq+MSA > structure at conserved sites,
  0.139 vs 0.027) is real *for P53* but does **not** generalize into "conserved ⇒
  seq/MSA wins."
- ICBINB-worthy negative: the intuitive feature (conservation) is *not* the
  diagnostic; the DMS-score structure (4.3) is.
- **Figs:** F20 (main — selected class-accuracy panels, shared with §4.3); F11
  (appx, 4th-slot candidate — conservation null), F12 (appx — full P53 class case;
  see also H2 panels to-do).

**4.5 Ruling out other trivial reductions.** Physicochemistry (hydro+charge ~13%)
and 3-D pocket proximity (1-D binding-site signal fails in real 3-D, Pearson-0.996
structure-robust) — disagreement isn't reducible to a single simple feature.
- **Figs:** F13 (appx).

## 5. Boundary conditions & actionable takeaways  *(element #4)*

**Boundary conditions:**
1. **We score one site at a time from WT, not sequentially conditioned mutations.**
   Our data are mostly **single** substitutions, so a
   within-site ranking is evaluated against the wild-type background; the real
   optimization mutates one site per step but each step is **conditioned on the
   mutations already fixed** (the Gibbs/search loop of §1). Single-substitution DMS
   cannot capture that conditioning — within-site conclusions transfer only if
   epistasis is mild.
2. **Improving the model selection metric does not, by itself, resolve
   optimization.** Even
   with a well-chosen model, the **samplers themselves fail to elicit sequences
   that are high pseudo-likelihood under their own models** (a separate search/
   inference failure, not a scoring failure) — this paper does not address that.
3. **We are not addressing how models select *which* sites to mutate** — only how
   they rank substitutions at a given site. Site selection is out of scope.
4. MSA-dependent models confounded by depth (both tails of `MSA_Neff_L` hurt).
5. Analyses skew to **well-studied human proteins**.
6. **Range restriction** attenuates every tail Spearman — read comparisons, not
   absolute tail values.
7. **DMS treated as ground truth but isn't a clean fitness readout** —
   detection/selection bias (missing-not-at-random dropouts) + heterogeneous assay
   quality (dynamic range, replicate noise, coverage). Inseparable from §4.3: a
   compressed within-site range can be real site biology **or** an assay
   resolving-power limit, and we can't always tell which.

**Actionable takeaways:**
1. **Don't select a single model by category-level Spearman/NDCG@10 for
   optimization.** Global rank mostly tracks within-site but breaks
   catastrophically on some assays (§3.2), collapses on the winners tail (§3.3),
   and barely predicts which model is right at disagreement sites (§4.2); the
   per-site oracle dwarfs every tier.
2. **Choose the metric from the goal, not a one-number leaderboard.** Report a
   *suite* of task-aligned metrics — **within-site ranking** and **top-quartile
   ("winners") accuracy** for optimization, in addition to global Spearman/NDCG@10
   for the deleterious-filter use case — and pick by objective. Spearman and
   NDCG@10 are redundant to each other and don't stand in for the optimization
   axes. (Constructive counterpart to the §1 premise: better ProteinGym numbers
   don't transfer to design because the reported metric doesn't match the goal.)
3. **Treat low within-site DMS range as a red flag that the position is
   unrankable.** Where substitutions cluster near WT (tolerated) or span a
   below-average range (hard to distinguish), disagreement reflects the assay's
   resolving power, not model quality — no model and no better metric will fix it.
   Route these sites to a **higher-resolution assay**, not more replicates of the
   same low-resolution one.
4. **Don't expect conservation to tell you which model to trust** — high
   conservation ⇒ sequence/MSA wins does **not** generalize (§4.4); a per-protein
   case, not a selection rule.
5. **Use model disagreement as site-level triage, with realistic expectations.**
   Agree → any model; disagree → first check within-site DMS range (takeaway 3)
   before asking *which* model, because much disagreement is unrankable rather
   than resolvable. Consensus↔accuracy is real but capped (~0.3–0.4).
6. **Use disagreement to assess assays and budget experiments, not just models.**
   The same signals run in reverse: positions where models disagree *and*
   within-site DMS range is compressed flag where the **assay** (not the model) is
   uninformative. This gives a data-driven read on **how reliable a DMS is as a
   training/evaluation yardstick** and **where extra experimental spend is (or
   isn't) worth it** — invest in higher-resolution readouts at the unrankable,
   functionally important sites rather than blanket re-measurement. (Realizes §1
   goal 1: think critically about the experimental assays themselves.)
7. **The ceiling worth chasing is the per-site oracle** — site-adaptive model
   choice and better within-site-resolved assays, not a better global leaderboard.
- **Figs:** none (MSA_Neff_L boundary is a textual point).

## 6. Conclusion / what to eval on next

Restate the gap (Spinner et al., generalized + mechanized) and the constructive
turn: optimization-relevant in silico metrics (within-site + winners), site-level
triage by within-site DMS range, and honest accounting of DMS measurement limits.

## Appendix

Selected-assay table (**T1**); directionality validation (`00`); model subsets
(23 / 45 / 97); co-localization
stats (`32`/`33`); DMS-score-structure detail (`34` H1); conservation null +
P53 class case (`28`/`34` H2); physchem control (`12`); 3-D robustness
(`17`/`18`/`20`); GoF/LoF (`22`); per-protein disagreement maps (`23`/`27`);
`MSA_Neff_L` boundary (`02`/`28`).

---

## Open to-run analyses before submission

- **Split `29` by Category A vs B** (shared vs assay-specific high-disagreement
  sites) so §4.2 can claim the rank↔usefulness failure *specifically* at
  assay-specific sites, not just pooled.
- **F9 (new):** pooled within-site DMS range & |mean−WT| vs disagreement across
  assays — the generalization of the tolerance/compression explanation beyond P53.
- **(Optional) clean degradation plot (§4.2):** the "both metric-selected models
  degrade at high-disagreement sites" result is already in H2 (F21) + fig9 (F4), but
  a single purpose-built figure — whole-assay baseline ρ vs per-site ρ at
  high-disagreement sites, Spearman-selected and NDCG-selected side by side — would
  state it most cleanly. Data ready in `h2_assayspecific_rank_tiers.csv` +
  `h2_assayspecific_persite_spearman_by_site{,_ndcg}.csv`. Make only if a crisper
  figure is wanted.
- **F16 regeneration (later, §3.1):** re-run `15` (functional enrichment) with the
  updated ClinVar/UniProt annotations (`23`/`25`) and the expanded allosteric map
  (`26`), then refresh the 10/13 · p=5e-13 numbers. *Deferred — do not run now.*
- **F20 assembly (later, §4.3–§4.4):** compose the main mechanism figure from
  selected F10 (H1 relfrac) + F12 (H2 class-accuracy) panels — overlaps the
  `h2_persite_spearman_density/` panel-selection to-do below.
- **Pick panels from `figures/h2_persite_spearman_density/` (H2 / §4.4):** choose
  which per-site-Spearman-by-model-class panels to include (main vs appx) — e.g.
  `fig_shared_highcons_persite_spearman_panel.png` (seq+MSA > structure at conserved
  P53 sites) vs `fig_shared_lowcons_persite_spearman_panel.png` (all classes ≈0 at
  variable sites), plus `fig_p53_highcons_classwise.png`,
  `fig_assayspecific_modelclass_best_panel.png`,
  `fig_assayspecific_rank_tiers*.png`. These are the visual evidence for the H2
  class-gap / conservation-null story (§4.4) and the assay-specific rank-tier point
  (§4.2). *[User to select.]*
- **Recombine `figures/conservation_dms_panels/` into a 3-regime figure:** compose
  a single panel showing the **three conservation levels of high-disagreement
  sites** — conserved (e.g. P53 Giacomelli), variable (e.g. PTEN, HXK4), and
  intermediate/mixed (e.g. CP2C9) — from the existing per-assay panels
  (`P53_*`, `PTEN_*`, `HXK4_*`, `CP2C9_*`). Illustrates that high disagreement
  occurs across the conservation spectrum (ties §4.3 tolerance/compression + §4.4
  conservation null). Candidate appendix figure for §4.4.
- Confirm which figs are regenerable vs need new scripts (F2, F9 appear to need
  new cuts).
- **Allostery expansion (§3.1):** extend the allosteric-site atlas beyond
  SRC_HUMAN (HXK4/RASK/CP2C9/PTEN/VKOR1, per NEXT_STEPS) and test whether
  high-disagreement sites concentrate at allosteric/distal functional residues —
  strengthens the "functionally important sites" stakes. *[User adding info.]*
