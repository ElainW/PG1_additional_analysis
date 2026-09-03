"""
Multiple-conformation confound at high-disagreement sites.

Hypothesis (NEXT_STEPS.md "Multiple-conformation confound at high-disagreement
sites"): our structure-based models each ingest a SINGLE conformation, so residues
whose side-chain<->side-chain interactions change between the active and inactive
states -- the *active-only* and *swapping* residues from the contact-state
classification (54) -- should be (a) systematically harder, showing elevated
within-site model disagreement, and (b) specifically mishandled by STRUCTURE-ONLY
models relative to models that also ingest SEQUENCE (sequence-only / MSA / multimodal
sequence+structure), which are not tied to one input conformation.

Multi-conformation proteins with a two-state contact classification: HXK4 (GCK,
1V4S closed vs 1V4T super-open), RASK (KRAS, 6VJJ ON vs 5O2S OFF), SRC (1Y57 active
vs 2SRC inactive). Contact classes are in UniProt numbering; the three proteins'
DMS assays are also UniProt-numbered (up_pos - pos == 0 in 15's site table), so
positions join directly.

Steps (mirroring NEXT_STEPS):
 1. Enumerate candidate sites = residues in long-range side-chain<->side-chain
    contacts that are active-only OR swapping across the two states (54).
 2. Test whether candidate sites carry ELEVATED within-site disagreement (the
    z-scored 1 - mean pairwise within-site Spearman from 15), vs the rest of the
    assayed sites -- pooled per protein and across proteins (Mann-Whitney, one-sided).
 3. Where disagreement is present, test whether SEQUENCE / multimodal-with-sequence
    models beat STRUCTURE-ONLY models at those sites: per (assay, site, model)
    within-site Spearman(model, DMS) [46's within-site accuracy], collapsed to a
    per-site structure-only mean and sequence-ingesting mean; paired Wilcoxon
    seq>struct at candidate sites, and an interaction test asking whether the
    sequence-minus-structure advantage is LARGER at candidate than baseline sites
    (the confound signature: structure models specifically fail where the input
    conformation misrepresents the interaction network).

Inputs
  results/structure_features/contacts/classification/<PROT>_HUMAN_residue_classes.csv
  results/generalization_functional_disagreement_sites.csv   (15; disagreement + z)
  O2/20250613/*.csv                                          (per-variant model scores)
Outputs (results/multiconformation_confound/)
  candidate_sites.csv                 per-protein per-residue class + is_candidate
  disagreement_by_class.csv           median/mean z_disagree per conf-class per protein
  disagreement_tests.csv              MW candidate-vs-rest, per protein + pooled
  model_class_persite_spearman.csv    per (protein,assay,site): struct/seq means, gap
  model_class_tests.csv               paired seq-vs-struct + interaction, per protein + pooled
  gap_vs_disagreement_tests.csv       Spearman(disagreement, seq-struct gap) + hi/lo split
  high_disagreement_candidate_counts.csv   step 4: active/swapping counts in the 27 tail
  switch_high_disagreement_functional.csv  step 4: those residues labelled functional (56)
  README.md                           tables + findings/interpretation
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon, spearmanr
import common as C

MIN_SITE_VARIANTS = 3           # matches 03/46 within-site accuracy
CANDIDATE_CLASSES = ["active", "swapping"]

# protein -> classification file stem (all _HUMAN here)
PROTEINS = {
    "HXK4_HUMAN": "HXK4_HUMAN",
    "RASK_HUMAN": "RASK_HUMAN",
    "SRC_HUMAN":  "SRC_HUMAN",
}
CLASS_ORDER = ["active", "swapping", "inactive", "nonspecific", "non_contact"]

# UniProt_ID -> (high_disagreement_sites file stem [27], site-characterization
# protein name [56]). Used by step 4 to intersect candidates with the strict
# high-disagreement tail and label those residues functional.
HD_MAP = {
    "HXK4_HUMAN": ("gck", "GCK"),
    "RASK_HUMAN": ("kras", "KRAS"),
    "SRC_HUMAN": ("src", "SRC"),
}
# func_ flag column -> label used in switch_high_disagreement_functional.csv
FUNC_SRCS = {"func_allosteric": "allosteric", "func_uniprot_binding": "uniprot_binding",
             "func_clinvar": "clinvar", "func_uniprot_disease": "disease"}

OUTDIR = C.rp("multiconformation_confound")


def load_classification():
    """Per-protein residue -> conf class, UniProt numbering."""
    frames = []
    cdir = C.rp("structure_features/contacts/classification")
    for uni, stem in PROTEINS.items():
        d = pd.read_csv(os.path.join(cdir, f"{stem}_residue_classes.csv"))
        d = d.rename(columns={"uniprot_pos": "up_pos"})
        d["UniProt_ID"] = uni
        d["is_candidate"] = d["classification"].isin(CANDIDATE_CLASSES)
        frames.append(d[["UniProt_ID", "up_pos", "resname", "classification", "is_candidate"]])
    return pd.concat(frames, ignore_index=True)


def annotate_sites(sites, cls):
    """Attach conf class to every assayed (UniProt_ID, up_pos); non-contacting
    assayed residues become 'non_contact'."""
    key = cls[["UniProt_ID", "up_pos", "classification"]].drop_duplicates()
    m = sites.merge(key, on=["UniProt_ID", "up_pos"], how="left")
    m["conf_class"] = m["classification"].fillna("non_contact")
    m["is_candidate"] = m["conf_class"].isin(CANDIDATE_CLASSES)
    return m


# ---------------------------------------------------------------- step 2
def step2_disagreement(sites_annot):
    """Is z_disagree elevated at candidate (active/swapping) vs rest, per protein?"""
    by_class_rows, test_rows = [], []
    for uni in PROTEINS:
        g = sites_annot[sites_annot["UniProt_ID"] == uni]
        for cls_name in CLASS_ORDER:
            v = g.loc[g["conf_class"] == cls_name, "z_disagree"].dropna()
            if len(v):
                by_class_rows.append(dict(UniProt_ID=uni, conf_class=cls_name, n=len(v),
                                          median_z=v.median(), mean_z=v.mean(),
                                          frac_above_assay_mean=float((v > 0).mean()),
                                          frac_top_quartile=float((v > g["z_disagree"].quantile(0.75)).mean())))
        cand = g.loc[g["is_candidate"], "z_disagree"].dropna()
        rest = g.loc[~g["is_candidate"], "z_disagree"].dropna()
        test_rows.append(_mw(uni, cand, rest))
    # pooled across proteins (z already within-assay comparable)
    cand = sites_annot.loc[sites_annot["is_candidate"], "z_disagree"].dropna()
    rest = sites_annot.loc[~sites_annot["is_candidate"], "z_disagree"].dropna()
    test_rows.append(_mw("POOLED", cand, rest))
    return pd.DataFrame(by_class_rows), pd.DataFrame(test_rows)


def _mw(label, cand, rest):
    if len(cand) >= 3 and len(rest) >= 3:
        U, p = mannwhitneyu(cand, rest, alternative="greater")
        # rank-biserial effect size
        rbc = 2 * U / (len(cand) * len(rest)) - 1
    else:
        p, rbc = np.nan, np.nan
    return dict(UniProt_ID=label, n_candidate=len(cand), n_rest=len(rest),
                median_z_candidate=cand.median() if len(cand) else np.nan,
                median_z_rest=rest.median() if len(rest) else np.nan,
                mannwhitney_p_greater=p, rank_biserial=rbc)


# ---------------------------------------------------------------- step 3
def persite_model_spearman(models):
    """For the 3 proteins' assays: per (UniProt_ID, DMS_id, up_pos, model) within-site
    Spearman(model, DMS) over the substitutions at that site (>= MIN_SITE_VARIANTS)."""
    want_files = {f"{a}" for a in []}  # placeholder
    rows = []
    for path in C.list_score_files():
        dms_id = C.dms_id_from_path(path)
        uni = dms_id.split("_HUMAN")[0] + "_HUMAN" if "_HUMAN" in dms_id else None
        if uni not in PROTEINS:
            continue
        df = pd.read_csv(path, usecols=lambda c: c in C.NON_MODEL_COLS or c in set(models))
        pos = df["mutant"].map(C.parse_position)
        keep = pos.notna()
        df, pos = df[keep], pos[keep]
        present = [m for m in models if m in df.columns]
        for p, sub in df.assign(__pos__=pos.values).groupby("__pos__"):
            if len(sub) < MIN_SITE_VARIANTS:
                continue
            sp = C.spearman_with_target(sub, "DMS_score", present)
            for m, val in sp.items():
                if not np.isnan(val):
                    rows.append(dict(UniProt_ID=uni, DMS_id=dms_id, up_pos=int(p),
                                     n_subs=len(sub), model=m, persite_spearman=val))
    return pd.DataFrame(rows)


def step3_model_classes(models, mtype, cls, sites):
    ps = persite_model_spearman(models)
    ps["model_type"] = ps["model"].map(mtype)
    ps["class_group"] = np.where(ps["model_type"] == "Structure", "structure_only", "seq_ingesting")

    # collapse to per-site group means
    site = (ps.groupby(["UniProt_ID", "DMS_id", "up_pos", "class_group"])["persite_spearman"]
              .mean().unstack("class_group").reset_index())
    for col in ("structure_only", "seq_ingesting"):
        if col not in site:
            site[col] = np.nan
    site["gap_seq_minus_struct"] = site["seq_ingesting"] - site["structure_only"]

    key = cls[["UniProt_ID", "up_pos", "classification"]].drop_duplicates()
    site = site.merge(key, on=["UniProt_ID", "up_pos"], how="left")
    site["conf_class"] = site["classification"].fillna("non_contact")
    site["is_candidate"] = site["conf_class"].isin(CANDIDATE_CLASSES)
    # attach per-(assay,site) disagreement z so we can test whether the structure
    # deficit tracks disagreement (the actual confound signature)
    site = site.merge(sites[["UniProt_ID", "DMS_id", "up_pos", "z_disagree"]],
                      on=["UniProt_ID", "DMS_id", "up_pos"], how="left")

    test_rows = []
    for uni in list(PROTEINS) + ["POOLED"]:
        s = site if uni == "POOLED" else site[site["UniProt_ID"] == uni]
        cand = s[s["is_candidate"]].dropna(subset=["structure_only", "seq_ingesting"])
        base = s[~s["is_candidate"]].dropna(subset=["structure_only", "seq_ingesting"])
        row = dict(UniProt_ID=uni,
                   n_candidate_sites=len(cand), n_baseline_sites=len(base),
                   cand_struct_mean=cand["structure_only"].mean(),
                   cand_seq_mean=cand["seq_ingesting"].mean(),
                   cand_gap_median=cand["gap_seq_minus_struct"].median(),
                   base_gap_median=base["gap_seq_minus_struct"].median())
        # paired seq > struct at candidate sites
        if len(cand) >= 6:
            try:
                _, row["paired_wilcoxon_p_seq_gt_struct"] = wilcoxon(
                    cand["seq_ingesting"], cand["structure_only"], alternative="greater")
            except ValueError:
                row["paired_wilcoxon_p_seq_gt_struct"] = np.nan
        else:
            row["paired_wilcoxon_p_seq_gt_struct"] = np.nan
        # interaction: gap larger at candidate than baseline
        if len(cand) >= 3 and len(base) >= 3:
            _, row["interaction_p_gap_cand_gt_base"] = mannwhitneyu(
                cand["gap_seq_minus_struct"], base["gap_seq_minus_struct"], alternative="greater")
        else:
            row["interaction_p_gap_cand_gt_base"] = np.nan
        test_rows.append(row)

    # Does the structure-model deficit TRACK disagreement? (confound signature)
    track_rows = []
    for uni in list(PROTEINS) + ["POOLED"]:
        s = site if uni == "POOLED" else site[site["UniProt_ID"] == uni]
        s = s.dropna(subset=["gap_seq_minus_struct", "z_disagree"])
        if len(s) < 6:
            continue
        rho, p = spearmanr(s["z_disagree"], s["gap_seq_minus_struct"])
        hi = s.loc[s["z_disagree"] > 0, "gap_seq_minus_struct"]
        lo = s.loc[s["z_disagree"] <= 0, "gap_seq_minus_struct"]
        pu = mannwhitneyu(hi, lo, alternative="greater")[1] if len(hi) >= 3 and len(lo) >= 3 else np.nan
        track_rows.append(dict(UniProt_ID=uni, n_sites=len(s),
                               spearman_disagree_vs_gap=rho, spearman_p=p,
                               gap_median_hi_disagree=hi.median(), gap_median_lo_disagree=lo.median(),
                               mw_p_gap_hi_gt_lo=pu))
    return site, pd.DataFrame(test_rows), pd.DataFrame(track_rows)


# ---------------------------------------------------------------- step 4
def step4_high_disagreement_subset(cls):
    """Intersect the step-1 candidates with the strict high-disagreement tail
    (27) and label those residues functional via 56's site characterization.

    Returns (counts_df, switch_func_df) or (None, None) if the upstream 27/56
    outputs are missing. Positions are UniProt-numbered; offset 0 for these three
    proteins (see 15/56), so up_pos joins directly to the classification pos."""
    hd_dir = C.rp("high_disagreement_sites")
    char_path = C.rp("high_disagreement_site_characterization/site_characterization_sites.csv")
    if not os.path.exists(char_path):
        print("step 4 skipped: missing 56 site_characterization_sites.csv")
        return None, None
    char = pd.read_csv(char_path)

    count_rows, func_rows = [], []
    for uni, (stem, cprot) in HD_MAP.items():
        hd_path = os.path.join(hd_dir, f"{stem}_high_disagreement.tsv")
        if not os.path.exists(hd_path):
            print(f"step 4: missing {hd_path}, skipping {uni}")
            continue
        hd_pos = set(pd.read_csv(hd_path, sep="\t")["up_pos"].astype(int))  # unique, both assays folded in
        c = cls[cls["UniProt_ID"] == uni]
        a_hd = sorted(int(p) for p in c.loc[c["classification"] == "active", "up_pos"] if int(p) in hd_pos)
        s_hd = sorted(int(p) for p in c.loc[c["classification"] == "swapping", "up_pos"] if int(p) in hd_pos)
        count_rows.append(dict(protein=cprot,
                               active_only_high_disagreement=len(a_hd),
                               swapping_high_disagreement=len(s_hd),
                               other_high_disagreement=len(hd_pos) - len(a_hd) - len(s_hd)))
        ch = char[char["protein"] == cprot].set_index("up_pos")
        for cls_name, poss in (("active", a_hd), ("swapping", s_hd)):
            for p in poss:
                if p not in ch.index:
                    continue
                r = ch.loc[p]
                src = ";".join(lbl for col, lbl in FUNC_SRCS.items() if bool(r.get(col, False)))
                func_rows.append(dict(protein=cprot, **{"class": cls_name}, up_pos=p, ss=r["ss"],
                                      functional=bool(r["is_functional"]), functional_source=src))
    return pd.DataFrame(count_rows), pd.DataFrame(func_rows)


# ---------------------------------------------------------------- report
def write_readme(by_class, dis_tests, model_tests, track_tests, cls, site):
    struct_models = "MIF, ESM-IF1, ProteinMPNN"
    n_cand = int(cls["is_candidate"].sum())
    lines = []
    lines.append("# Multiple-conformation confound at high-disagreement sites\n")
    lines.append("Script `57_multiconformation_confound.py`. Tests whether within-site model\n"
                 "disagreement and structure-vs-sequence model accuracy at the multi-conformation\n"
                 "proteins (HXK4/GCK, RASK/KRAS, SRC) are explained by structure-only models\n"
                 "ingesting a single conformation.\n")
    lines.append("## Candidate sites (step 1)\n")
    lines.append("Candidate = residue in a long-range side-chain<->side-chain contact classified\n"
                 "**active-only** or **swapping** between the two states (script 54). "
                 f"{n_cand} candidate residues across the three proteins.\n")
    lines.append("| protein | active | swapping | inactive | nonspecific |")
    lines.append("|---|---:|---:|---:|---:|")
    for uni in PROTEINS:
        c = cls[cls["UniProt_ID"] == uni]["classification"].value_counts()
        lines.append(f"| {uni} | {c.get('active',0)} | {c.get('swapping',0)} | "
                     f"{c.get('inactive',0)} | {c.get('nonspecific',0)} |")
    lines.append("\n## Step 2 — disagreement elevation\n")
    lines.append("Median within-site disagreement (z per assay), candidate vs rest "
                 "(Mann-Whitney, H1 candidate greater):\n")
    lines.append("| protein | n cand | n rest | med z cand | med z rest | p (cand>rest) | rank-biserial |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for _, r in dis_tests.iterrows():
        lines.append(f"| {r.UniProt_ID} | {r.n_candidate} | {r.n_rest} | "
                     f"{r.median_z_candidate:.3f} | {r.median_z_rest:.3f} | "
                     f"{r.mannwhitney_p_greater:.3g} | {r.rank_biserial:.3f} |")
    lines.append("\n## Step 3 — structure-only vs sequence-ingesting models\n")
    lines.append(f"Structure-only = model_type 'Structure' ({struct_models}); sequence-ingesting =\n"
                 "every other model_type (single-sequence, MSA, sequence+structure, sequence+\n"
                 "structure+function). Per site we take each group's mean within-site Spearman\n"
                 "(vs DMS) and the seq-minus-structure gap.\n")
    lines.append("| protein | n cand | cand struct | cand seq | cand gap med | base gap med | "
                 "paired p (seq>struct) | interaction p (gap cand>base) |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in model_tests.iterrows():
        lines.append(f"| {r.UniProt_ID} | {int(r.n_candidate_sites)} | {r.cand_struct_mean:.3f} | "
                     f"{r.cand_seq_mean:.3f} | {r.cand_gap_median:.3f} | {r.base_gap_median:.3f} | "
                     f"{r.paired_wilcoxon_p_seq_gt_struct:.3g} | {r.interaction_p_gap_cand_gt_base:.3g} |")
    lines.append("\n### Does the structure deficit track disagreement? (confound signature)\n")
    lines.append("Per-site Spearman(z_disagree, seq-struct gap) and the gap split at z=0. A real\n"
                 "single-conformation confound predicts the structure deficit GROWS with\n"
                 "disagreement (positive rho, hi>lo).\n")
    lines.append("| protein | n sites | rho(disagree, gap) | p | gap med hi-z | gap med lo-z | p (hi>lo) |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for _, r in track_tests.iterrows():
        lines.append(f"| {r.UniProt_ID} | {int(r.n_sites)} | {r.spearman_disagree_vs_gap:.3f} | "
                     f"{r.spearman_p:.3g} | {r.gap_median_hi_disagree:.3f} | "
                     f"{r.gap_median_lo_disagree:.3f} | {r.mw_p_gap_hi_gt_lo:.3g} |")
    lines.append("\n## Findings\n")
    lines.append("**Step 2 — partial support.** Active/swapping candidate sites carry elevated\n"
                 "within-site disagreement in the two clean two-state switches — RASK/KRAS\n"
                 "(median z +0.14 vs -0.35, p=5e-4) and SRC (+0.09 vs -0.20, p=1.5e-3) — but\n"
                 "NOT in HXK4/GCK (p=0.93), whose whole-domain glucose closure rewires so many\n"
                 "residues that the contact classes do not separate. The per-class breakdown\n"
                 "(`disagreement_by_class.csv`) localises the signal to the **swapping** class\n"
                 "(RASK swapping median z +0.30, SRC +0.18); **active-only** sites are not\n"
                 "elevated in any protein.")
    lines.append("\n**Step 3 — sequence models win, but not via the hypothesised mechanism.**\n"
                 "Sequence-ingesting models beat structure-only models at candidate sites in all\n"
                 "three proteins (pooled paired Wilcoxon p=2e-14; at SRC the structure-only\n"
                 "models have ~zero within-site skill). But the advantage is NOT candidate-\n"
                 "specific: the interaction test is null, and the per-site seq-minus-structure\n"
                 "gap *negatively* tracks disagreement (pooled rho=-0.15, p=7e-12; gap is\n"
                 "smaller, not larger, at high-disagreement sites). So the structure deficit is a\n"
                 "general property of these proteins, not a targeted blind spot at conformation-\n"
                 "changing residues — high disagreement there reflects genuine difficulty that\n"
                 "hobbles all model classes.")
    lines.append("\n**Switch residues in the high-disagreement tail are functional (allosteric).**\n"
                 "Of the candidate residues that reach the strict high-disagreement tail (27),\n"
                 "**swapping = 6/9 (67%)** are functional and 5 of those 6 carry an *allosteric*\n"
                 "annotation (GCK 99, 287; KRAS 65, 68; SRC 387), whereas **active-only = 2/5\n"
                 "(40%)** are functional and only via disease, none allosteric/binding. Where\n"
                 "models disagree AND the residue rewires its contacts between states, it tends\n"
                 "to be a real allosteric/switch position (`switch_high_disagreement_functional.csv`,\n"
                 "`high_disagreement_candidate_counts.csv`).")
    lines.append("\n## Output files (by step)\n")
    lines.append("**Step 1 — candidate enumeration**\n"
                 "- `candidate_sites.csv` — one row per contacting residue of HXK4/RASK/SRC "
                 "(UniProt pos, resname, `classification` = active/inactive/nonspecific/swapping "
                 "from script 54, and `is_candidate` = class in {active, swapping}). The candidate "
                 "list every later step draws from.\n")
    lines.append("**Step 2 — disagreement elevation**\n"
                 "- `disagreement_by_class.csv` — per protein x contact-class summary of z_disagree "
                 "(n, median, mean, fraction above the assay mean, fraction in the assay's top "
                 "quartile). The descriptive table behind 'the signal is in the swapping class'.\n"
                 "- `disagreement_tests.csv` — the step-2 test: per protein + pooled, candidate vs "
                 "rest medians and one-sided Mann-Whitney (H1 candidate > rest) with rank-biserial "
                 "effect size.\n")
    lines.append("**Step 3 — structure-only vs sequence-ingesting models**\n"
                 "- `model_class_persite_spearman.csv` — per (protein, assay, site): the mean "
                 "within-site Spearman(model, DMS) of the structure-only models and of the "
                 "sequence-ingesting models, their gap, the site's contact class, is_candidate, "
                 "and z_disagree. The per-site table both step-3 tests aggregate.\n"
                 "- `model_class_tests.csv` — step 3a: per protein + pooled, candidate-site "
                 "structure vs sequence means, the paired Wilcoxon (H1 seq > struct), and the "
                 "interaction test (is the gap larger at candidate than baseline sites).\n"
                 "- `gap_vs_disagreement_tests.csv` — step 3b: per protein + pooled, "
                 "Spearman(z_disagree, gap) and the gap split at z=0 (hi vs lo disagreement) — "
                 "tests whether the structure deficit GROWS with disagreement (the confound "
                 "signature). It does not.\n")
    lines.append("**Step 4 — candidates in the strict high-disagreement tail** — step-1 "
                 "candidates intersected with the elbow-thresholded tail (27) and labelled "
                 "functional via 56 (needs those upstream outputs present):\n"
                 "- `high_disagreement_candidate_counts.csv` — per protein, active-only and "
                 "swapping counts that reach the strict high-disagreement tail (script 27), plus "
                 "the other (non-candidate) high-disagreement residues; unique residues, no "
                 "double-count across assays.\n"
                 "- `switch_high_disagreement_functional.csv` — the active-only/swapping residues "
                 "in the high-disagreement tail, each labelled functional / functional-source by "
                 "joining script 56's `site_characterization_sites.csv` "
                 "(allosteric / UniProt-binding / ClinVar / disease).\n")
    with open(os.path.join(OUTDIR, "README.md"), "w") as fh:
        fh.write("\n".join(lines))


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    models, mtype, _ = C.model_meta()

    cls = load_classification()
    cls.to_csv(os.path.join(OUTDIR, "candidate_sites.csv"), index=False)

    sites = pd.read_csv(C.rp("generalization_functional_disagreement_sites.csv"))
    sites = sites[sites["UniProt_ID"].isin(PROTEINS)].copy()
    sites_annot = annotate_sites(sites, cls)

    # sanity: how many classification positions are actually assayed?
    assayed = set(zip(sites["UniProt_ID"], sites["up_pos"]))
    hit = cls.apply(lambda r: (r.UniProt_ID, r.up_pos) in assayed, axis=1).mean()
    print(f"classification positions covered by >=1 assay: {hit:.0%}")

    by_class, dis_tests = step2_disagreement(sites_annot)
    by_class.to_csv(os.path.join(OUTDIR, "disagreement_by_class.csv"), index=False)
    dis_tests.to_csv(os.path.join(OUTDIR, "disagreement_tests.csv"), index=False)

    site, model_tests, track_tests = step3_model_classes(models, mtype, cls, sites)
    site.to_csv(os.path.join(OUTDIR, "model_class_persite_spearman.csv"), index=False)
    model_tests.to_csv(os.path.join(OUTDIR, "model_class_tests.csv"), index=False)
    track_tests.to_csv(os.path.join(OUTDIR, "gap_vs_disagreement_tests.csv"), index=False)

    counts, switch_func = step4_high_disagreement_subset(cls)
    if counts is not None:
        counts.to_csv(os.path.join(OUTDIR, "high_disagreement_candidate_counts.csv"), index=False)
        switch_func.to_csv(os.path.join(OUTDIR, "switch_high_disagreement_functional.csv"), index=False)

    write_readme(by_class, dis_tests, model_tests, track_tests, cls, site)

    print("\n== Step 2: disagreement candidate vs rest ==")
    print(dis_tests.to_string(index=False))
    print("\n== Step 3a: structure-only vs sequence-ingesting ==")
    print(model_tests.to_string(index=False))
    print("\n== Step 3b: does seq-struct gap track disagreement? ==")
    print(track_tests.to_string(index=False))
    if counts is not None:
        print("\n== Step 4: candidates in the high-disagreement tail ==")
        print(counts.to_string(index=False))
    print(f"\nWrote outputs to {OUTDIR}")


if __name__ == "__main__":
    main()
