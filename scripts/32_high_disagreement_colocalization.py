"""
Do the high-disagreement sites (the black dots exported by 27) co-localize across
assays of the *same* protein?

For each UniProt_ID we take the labelled high-disagreement sites of every assay
(results/high_disagreement_sites/<common>_high_disagreement.tsv), form the UNION of
positions across all its assays, and count how many of those positions are "shared" =
labelled high-disagreement in >= 2 of the assays. For proteins with >2 assays a
position counts as shared if ANY two assays both flag it (that is exactly "labelled in
>=2 assays"). We annotate each assay with its MSA Neff/L, since a natural worry is that
disagreement at a site is just an artefact of a shallow alignment -- if the same sites
light up across assays that measure different phenotypes, that worry is answered.

Excluded per instruction:
  SPG1_STRSG   -- 2nd assay has very few sites
  F7YBW8_MESOW -- too few positions total across its 2 assays (empty tsv here)
  NPC1_HUMAN   -- the two assays barely overlap in coordinates

Significance: each assay gets a permutation p (N=20000), per-assay p's are combined per
protein by Fisher's method, and the protein-level p's are BH-FDR corrected across proteins
(q_colocalization). Effect size is reported as fold-enrichment (obs/exp) plus a
permutation-null 95% band, both per assay and rolled up per protein (protein band =
convolution of the per-assay nulls). The headline pooled fraction gets a cluster-bootstrap
95% CI (resampling proteins, not sites).

Outputs -> results/high_disagreement_colocalization/
  colocalization_summary.csv     one row per protein (union n, shared n, fraction; rollup
                                 obs/exp totals, null band, fold-enrichment; Fisher p + BH q; Neff/L)
  colocalization_shared_sites.csv one row per shared position (which assays flag it)
  colocalization_per_assay.csv    one row per assay (n sites, obs overlap, null band, fold, p_perm, Neff/L)
  colocalization_pooled_ci.csv    one row: pooled co-localized fraction + cluster-bootstrap 95% CI
"""
import os
import sys
import glob
from collections import Counter
import numpy as np
import pandas as pd
import common as C

N_PERM = 20000
N_BOOT = 10000          # cluster-bootstrap resamples for the pooled-fraction CI
RNG = np.random.default_rng(0)

# Which disagreement tail to test: "high" (default, from 27) or "low" (from 30).
# The permutation null is identical -- only the flagged site set differs; the
# testable-position pool (generalization_functional_disagreement_sites.csv) is shared.
MODE = sys.argv[1] if len(sys.argv) > 1 else "high"
assert MODE in ("high", "low"), "usage: 32_high_disagreement_colocalization.py [high|low]"

INDIR = C.rp(f"{MODE}_disagreement_sites")
OUTDIR = C.rp(f"{MODE}_disagreement_colocalization")
META = C.REF_PATH

EXCLUDE = {"SPG1_STRSG", "F7YBW8_MESOW", "NPC1_HUMAN"}


def perm_pvalue(pool, k, other_union, obs):
    """Per-assay permutation test of co-localization. Draw k positions at random
    (no replacement) from this assay's testable `pool`, keeping the other assays'
    real flagged sets (`other_union`) fixed; statistic = # drawn positions that
    fall in other_union. p = P(perm overlap >= observed).
    Returns (p, expected_overlap, null_lo95, null_hi95) where [lo, hi] are the
    2.5/97.5 percentiles of the permutation null -- a reference band to read `obs`
    against, free of the 1/(N_PERM+1) p-value floor. `hits` (or None) is the full
    permutation null vector, returned so the caller can convolve nulls across a
    protein's assays for a protein-level band."""
    pool = np.asarray(sorted(pool))
    if k == 0 or len(other_union) == 0:
        return np.nan, 0.0, np.nan, np.nan, None
    others = np.array([p in other_union for p in pool])  # boolean mask over pool
    hits = np.empty(N_PERM, dtype=int)
    for i in range(N_PERM):
        idx = RNG.choice(len(pool), size=k, replace=False)
        hits[i] = others[idx].sum()
    p = (1 + (hits >= obs).sum()) / (N_PERM + 1)
    lo, hi = np.percentile(hits, [2.5, 97.5])
    return p, float(hits.mean()), float(lo), float(hi), hits


def fold_enrichment(obs, exp):
    """obs/exp, guarding exp==0: NaN when obs is also 0, +inf otherwise."""
    if exp > 0:
        return obs / exp
    return float("nan") if obs == 0 else float("inf")


def bh_qvalues(pvals):
    """Benjamini-Hochberg FDR q-values. NaN p-values pass through as NaN and are
    excluded from the count m. Valid on discrete permutation p-values (BH only
    needs each p to be a valid p-value); on ties/discreteness it is conservative."""
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    mask = ~np.isnan(p)
    pv = p[mask]
    m = pv.size
    if m == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    q_ranked = ranked * m / np.arange(1, m + 1)
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]  # enforce monotonicity
    q_ranked = np.minimum(q_ranked, 1.0)
    out = np.empty(m)
    out[order] = q_ranked
    q[mask] = out
    return q


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    meta = pd.read_csv(META)
    neff = dict(zip(meta["DMS_id"], meta["MSA_Neff_L"]))
    neffcat = dict(zip(meta["DMS_id"], meta["MSA_Neff_L_category"]))

    sites = pd.read_csv(C.rp("generalization_functional_disagreement_sites.csv"))
    pools = {d: set(g["up_pos"].astype(int))          # testable positions per assay
             for d, g in sites.groupby("DMS_id")}

    summary, shared_rows, assay_rows = [], [], []

    for f in sorted(glob.glob(os.path.join(INDIR, "*.tsv"))):
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if len(df) == 0:
            continue
        uni = df["UniProt_ID"].iloc[0]
        if uni in EXCLUDE:
            continue

        assays = {dms: set(g["up_pos"].astype(int))
                  for dms, g in df.groupby("DMS_id")}
        union = set().union(*assays.values())

        cnt = Counter(p for s in assays.values() for p in s)
        shared = {p for p, c in cnt.items() if c >= 2}

        # per-assay permutation p-value of co-localization with the OTHER assays
        pvals = {}
        obs_tot, exp_tot = 0, 0.0          # protein-level co-occurrence totals
        null_sum = np.zeros(N_PERM)        # convolved permutation null across assays
        for d, s in sorted(assays.items()):
            other_union = set().union(*[assays[o] for o in assays if o != d])
            obs = len(s & other_union)
            p, exp, lo, hi, hits = perm_pvalue(pools.get(d, s), len(s), other_union, obs)
            pvals[d] = p
            fold = fold_enrichment(obs, exp)
            obs_tot += obs
            exp_tot += exp
            if hits is not None:
                null_sum += hits
            assay_rows.append({
                "UniProt_ID": uni, "DMS_id": d,
                "n_high_disagree_sites": len(s),
                "pool_positions": len(pools.get(d, s)),
                "obs_overlap_with_others": obs,
                "exp_overlap_perm": round(exp, 2),
                "null_lo95": (round(lo, 2) if lo == lo else np.nan),
                "null_hi95": (round(hi, 2) if hi == hi else np.nan),
                "fold_enrichment": (round(fold, 2) if np.isfinite(fold) else fold),
                "p_perm": (round(p, 5) if p == p else np.nan),
                "Neff_L": round(neff.get(d, float("nan")), 2),
                "Neff_L_category": neffcat.get(d, ""),
            })

        # combine the assay p-values into one protein-level p (Fisher's method)
        ps = [p for p in pvals.values() if p == p]
        if ps:
            from scipy.stats import combine_pvalues
            p_prot = combine_pvalues(ps, method="fisher")[1]
        else:
            p_prot = float("nan")

        # protein-level rollup of the per-assay effect: pooled co-occurrence count
        # (each shared pair counted from both assays' side), its permutation-null
        # 95% band (convolution of the per-assay nulls), and fold = obs_tot/exp_tot.
        prot_fold = fold_enrichment(obs_tot, exp_tot)
        band_lo, band_hi = np.percentile(null_sum, [2.5, 97.5])

        min_neff = min(neff.get(d, float("nan")) for d in assays)
        summary.append({
            "UniProt_ID": uni,
            "n_assays": len(assays),
            "n_union_sites": len(union),
            "n_shared_sites": len(shared),
            "frac_colocalized": round(len(shared) / len(union), 3) if union else 0.0,
            "obs_overlap_total": obs_tot,
            "exp_overlap_total": round(exp_tot, 2),
            "null_lo95": round(float(band_lo), 2),
            "null_hi95": round(float(band_hi), 2),
            "fold_enrichment": (round(prot_fold, 2) if np.isfinite(prot_fold) else prot_fold),
            "p_colocalization": p_prot,          # numeric (Fisher); NaN if no p's
            "min_Neff_L": round(min_neff, 2),
            "Neff_L_category": neffcat.get(next(iter(assays)), ""),
        })

        wt = dict(zip(df["up_pos"].astype(int), df["residue"]))
        for p in sorted(shared):
            flagging = sorted(d for d, s in assays.items() if p in s)
            shared_rows.append({
                "UniProt_ID": uni, "up_pos": p, "residue": wt.get(p, ""),
                "n_assays_flagging": len(flagging),
                "assays_flagging": ";".join(flagging),
            })

    sm = pd.DataFrame(summary).sort_values("frac_colocalized", ascending=False)

    # BH-FDR across proteins on the protein-level (Fisher) p. Family = one test per
    # protein; the permutation origin of the underlying p's is not an obstacle to BH.
    sm["q_colocalization"] = bh_qvalues(sm["p_colocalization"].to_numpy())
    sm.to_csv(os.path.join(OUTDIR, "colocalization_summary.csv"), index=False)
    pd.DataFrame(shared_rows).to_csv(
        os.path.join(OUTDIR, "colocalization_shared_sites.csv"), index=False)

    # BH-FDR across ALL per-assay tests (the family shown in the per-assay forest
    # plot, script 33). q_perm is the multiple-testing-aware significance per assay.
    pa = pd.DataFrame(assay_rows)
    pa.insert(pa.columns.get_loc("p_perm") + 1,
              "q_perm", bh_qvalues(pa["p_perm"].to_numpy()))
    pa.to_csv(os.path.join(OUTDIR, "colocalization_per_assay.csv"), index=False)

    tot_u, tot_s = sm["n_union_sites"].sum(), sm["n_shared_sites"].sum()

    # Cluster-bootstrap 95% CI on the pooled co-localized fraction: resample PROTEINS
    # (not sites) with replacement so within-protein site correlation is respected.
    u = sm["n_union_sites"].to_numpy()
    s = sm["n_shared_sites"].to_numpy()
    nprot = len(u)
    boot = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = RNG.integers(0, nprot, nprot)
        boot[b] = s[idx].sum() / u[idx].sum()
    ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
    pooled = tot_s / tot_u
    pd.DataFrame([{
        "n_proteins": nprot,
        "pooled_union_sites": int(tot_u),
        "pooled_colocalized_sites": int(tot_s),
        "pooled_frac_colocalized": round(pooled, 4),
        "ci95_lo": round(ci_lo, 4),
        "ci95_hi": round(ci_hi, 4),
        "n_boot": N_BOOT,
    }]).to_csv(os.path.join(OUTDIR, "colocalization_pooled_ci.csv"), index=False)

    print(sm.to_string(index=False))
    print(f"\nProteins: {len(sm)}   with >=1 co-localized site: {(sm.n_shared_sites>0).sum()}")
    print(f"Significant (BH q<0.05): {(sm.q_colocalization < 0.05).sum()}")
    print(f"Pooled union sites: {tot_u}   co-localized (>=2 assays): {tot_s}  "
          f"= {pooled:.1%}  (95% CI {ci_lo:.1%}-{ci_hi:.1%}, cluster bootstrap over proteins)")
    print(f"Median per-protein fraction: {sm.frac_colocalized.median():.2f}")
    print(f"-> wrote 4 tables to {OUTDIR}")


if __name__ == "__main__":
    main()
