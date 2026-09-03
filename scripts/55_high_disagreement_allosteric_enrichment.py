"""
Are high-disagreement sites enriched in the curated functional regions of the
allosteric-map proteins?

For the five proteins that have a hand-curated orthosteric/allosteric residue
atlas (results/allosteric_map/: CYP2C9, GCK, KRAS, PTEN, SRC) we ask whether the
model-disagreement tail exported by 27 (results/high_disagreement_sites/) falls
in those "functional regions" more than the rest of the assayed sites would
predict.

Universe / testable pool: results/generalization_functional_disagreement_sites.csv
(15) -- every assayed (DMS_id, up_pos) with a disagreement score. Unit of the
primary test is (assay x position), matching how 27 exports the tail and 32
tests colocalization; a position-level collapse (high if flagged in ANY of its
assays) is reported alongside as robustness.

Functional region is defined from the three source directories the question
names, kept separate so each is testable on its own and then combined:
  ALLOSTERIC -- position in results/allosteric_map/<gene>.tsv (dir1). The curated
                orthosteric/allosteric atlas; the definition the question asks about.
  ALLOSTERIC_ONLY -- the subset of the atlas whose site_class contains "allosteric"
                (allosteric / allosteric_second_shell / allosteric_activator_site /
                binding_interface;allosteric), i.e. drops the orthosteric / active-
                site / binding-interface residues. CYP2C9 has none of these.
  UNIPROT    -- UniProt binding/active-site +-2 (dir3). Taken from the pool's own
                `near_binding` column (15 builds it from the UniProt Binding/Active
                site/Site/Metal features, window +-2). Structural only, no disease.
  CLINVAR    -- ClinVar pathogenic/likely-pathogenic missense positions (dir2),
                parsed here EXACTLY as 23 does: the p.(3-letter) change in the
                Name column, kept only where the wild-type matches the assay-isoform
                sequence, UNION the UniProt gain-of-function supplement (GoF disease
                variants a LoF-dominated ClinVar query misses). NB this is a
                different source from the pool's `disease` column, which 15 builds
                from UniProt natural-variant disease keywords, NOT from dir2.
  UNION       -- ALLOSTERIC OR UNIPROT OR CLINVAR (all three source dirs folded in).

Test: per protein a one-sided Fisher exact (H1: high-disagreement sites are MORE
often functional than the rest), odds ratio + 95% CI, plus a pooled 2x2 across
the five proteins. As a continuous companion we also Mann-Whitney the raw
disagreement magnitude of allosteric-functional vs non-functional assayed sites
(H1: functional greater), which probes the bulk of the distribution rather than
the thresholded tail.

Naming note: the allosteric-map and high-disagreement files are keyed by gene
symbol (CYP2C9, GCK) while the ProteinGym pool / UniProt json / ClinVar map are
keyed by UniProt entry name (CP2C9_HUMAN, HXK4_HUMAN). PROTEINS carries both.
ClinVar positions are matched against the assay-isoform sequence and used as
up_pos directly (these five assays are UniProt-canonical, offset 0), following 23.

Outputs -> results/high_disagreement_allosteric_enrichment/
  enrichment_by_protein.csv   one row per (protein x definition): 2x2 counts,
                              % of high sites functional, % of rest functional,
                              odds ratio + CI, one-sided Fisher p.
  enrichment_pooled.csv       one row per (definition x granularity): pooled 2x2
                              + OR + p, at assay-position and position-collapsed.
  magnitude_by_protein.csv    per protein: mean disagreement in allosteric-
                              functional vs non-functional sites + Mann-Whitney p.
"""
import os
import re
import json
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from scipy.stats.contingency import odds_ratio
import common as C

# gene (map/tail key), UniProt entry (pool/json/clinvar key), allosteric-map file,
# high-disagreement tail file, ClinVar export file (dir2).
PROTEINS = [
    ("CYP2C9", "CP2C9_HUMAN", "CYP2C9_HUMAN.tsv", "cyp2c9_high_disagreement.tsv", "CYP2C9.tsv"),
    ("GCK",    "HXK4_HUMAN",  "GCK_HUMAN.tsv",    "gck_high_disagreement.tsv",    "GCK.tsv"),
    ("KRAS",   "RASK_HUMAN",  "RASK_HUMAN.tsv",   "kras_high_disagreement.tsv",   "KRAS.tsv"),
    ("PTEN",   "PTEN_HUMAN",  "PTEN_HUMAN.tsv",   "pten_high_disagreement.tsv",   "PTEN.tsv"),
    ("SRC",    "SRC_HUMAN",   "SRC_HUMAN.tsv",    "src_high_disagreement.tsv",    "SRC.tsv"),
]
DEFS = ["ALLOSTERIC", "ALLOSTERIC_ONLY", "UNIPROT", "CLINVAR", "UNION"]

ALLO_DIR = os.path.join(C.RESULTS, "allosteric_map")
CLINVAR_DIR = os.path.join(C.DATA_DIR, "clinvar")
UNIPROT_DIR = os.path.join(C.DATA_DIR, "uniprot", "raw_json")
POOL_PATH = C.rp("generalization_functional_disagreement_sites.csv")
TAIL_DIR = C.rp("high_disagreement_sites")
OUTDIR = C.rp("high_disagreement_allosteric_enrichment")

# --- ClinVar parsing, replicated verbatim from 23_structure_maps_all.py ---------
GOF_KW = ("gain-of-function", "gain of function", "constitutively activ", "activating",
          "increased activ", "hyperactiv", "increases activ", "increase in activ",
          "hhf", "hyperinsulin", "phhi", "benta")
_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
        "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
        "Tyr": "Y", "Val": "V"}
_PCHANGE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")


def clinvar_positions(clinvar_file, seq):
    """ClinVar P/LP missense positions: p.(3-letter) change in the Name column,
    kept only where the wild-type matches the assay-isoform sequence."""
    df = pd.read_csv(os.path.join(CLINVAR_DIR, clinvar_file), sep="\t", index_col=False)
    pos = set()
    for nm in df["Name"]:
        m = _PCHANGE.search(str(nm))
        if not m:
            continue
        wt, p = _AA3.get(m.group(1), "?"), int(m.group(2))
        if seq and 1 <= p <= len(seq) and seq[p - 1] == wt:
            pos.add(p)
    return pos


def uniprot_gof_supplement(uni):
    """UniProt natural-variant positions whose disease description is gain-of-function
    (the GoF variants a LoF-dominated ClinVar query misses)."""
    path = os.path.join(UNIPROT_DIR, uni + ".json")
    if not os.path.exists(path):
        return set()
    e = json.load(open(path))
    pos = set()
    for f in e.get("features", []):
        if f["type"] == "Natural variant":
            d = (f.get("description", "") or "").lower()
            if any(k in d for k in GOF_KW):
                pos.add(f["location"]["start"]["value"])
    return pos


def allosteric_positions(map_file):
    """Set of residue positions curated in an allosteric-map tsv (skips # header)."""
    df = pd.read_csv(os.path.join(ALLO_DIR, map_file), sep="\t", comment="#")
    return set(df["position"].astype(int))


def allosteric_only_positions(map_file):
    """Positions whose site_class contains 'allosteric' -- the allosteric subset of
    the atlas (drops orthosteric / active-site / binding-interface residues)."""
    df = pd.read_csv(os.path.join(ALLO_DIR, map_file), sep="\t", comment="#")
    sc = df["site_class"].fillna("").str.lower()
    return set(df.loc[sc.str.contains("allosteric"), "position"].astype(int))


def fisher_row(a, b, c, d):
    """One-sided (greater) Fisher on [[a,b],[c,d]] with OR 95% CI. a/b = high in/out
    of functional; c/d = rest in/out of functional."""
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    try:
        ci = odds_ratio([[a, b], [c, d]]).confidence_interval(0.95)
        lo, hi = ci.low, ci.high
    except Exception:
        lo = hi = np.nan
    return dict(
        high_func=a, high_nonfunc=b, rest_func=c, rest_nonfunc=d,
        pct_high_func=100 * a / (a + b) if (a + b) else np.nan,
        pct_rest_func=100 * c / (c + d) if (c + d) else np.nan,
        odds_ratio=orr, or_ci_lo=lo, or_ci_hi=hi, p_greater=p,
    )


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    pool = pd.read_csv(POOL_PATH)
    pool["up_pos"] = pool["up_pos"].astype(int)
    # assay-isoform sequence per protein (for the ClinVar wild-type match), as in 15/23
    ref = C.drop_excluded_dms(C.load_reference())
    seqmap = ref.groupby("UniProt_ID")["target_seq"].first()

    by_protein = []          # per (protein, definition), assay-position unit
    mag_rows = []            # continuous disagreement magnitude (allosteric def)
    pooled = {d: {"ap": np.zeros((2, 2)), "pos": np.zeros((2, 2))} for d in DEFS}

    for gene, upid, map_file, tail_file, clinvar_file in PROTEINS:
        u = pool[pool["UniProt_ID"] == upid].copy()
        if u.empty:
            raise SystemExit(f"{upid} not in pool {POOL_PATH}")
        seq = seqmap.get(upid, "")

        allo = allosteric_positions(map_file)
        allo_only = allosteric_only_positions(map_file)
        uniprot = set(u.loc[u["near_binding"] == True, "up_pos"])           # noqa: E712
        clinvar = clinvar_positions(clinvar_file, seq) | uniprot_gof_supplement(upid)
        defs = {
            "ALLOSTERIC": allo,
            "ALLOSTERIC_ONLY": allo_only,
            "UNIPROT": uniprot,
            "CLINVAR": clinvar,
            "UNION": allo | uniprot | clinvar,
        }

        tail = pd.read_csv(os.path.join(TAIL_DIR, tail_file), sep="\t")
        tail["up_pos"] = tail["up_pos"].astype(int)
        tail_keys = set(zip(tail["DMS_id"], tail["up_pos"]))       # assay x position
        tail_pos = set(tail["up_pos"])                             # collapsed

        u["is_high"] = [(dms, p) in tail_keys for dms, p in zip(u["DMS_id"], u["up_pos"])]
        missing = tail_keys - set(zip(u["DMS_id"], u["up_pos"]))
        if missing:
            print(f"WARN {gene}: {len(missing)} tail (assay,pos) not in pool")

        pos_all = set(u["up_pos"])
        pos_high = tail_pos & pos_all

        for dname, F in defs.items():
            # assay x position 2x2
            f = u["up_pos"].isin(F)
            a = int((u["is_high"] & f).sum())
            b = int((u["is_high"] & ~f).sum())
            c = int((~u["is_high"] & f).sum())
            d = int((~u["is_high"] & ~f).sum())
            row = fisher_row(a, b, c, d)
            row.update(protein=gene, UniProt_ID=upid, definition=dname,
                       unit="assay_position", n_universe=len(u), n_func=int(f.sum()),
                       n_high=int(u["is_high"].sum()))
            by_protein.append(row)
            pooled[dname]["ap"] += [[a, b], [c, d]]

            # position-collapsed 2x2
            fp = pos_all & F
            a2 = len(pos_high & fp)
            b2 = len(pos_high - fp)
            c2 = len((pos_all & fp) - pos_high)
            d2 = len((pos_all - fp) - pos_high)
            pooled[dname]["pos"] += [[a2, b2], [c2, d2]]

        # continuous magnitude, allosteric definition
        fa = u["up_pos"].isin(allo)
        x = u.loc[fa, "disagreement"].values
        y = u.loc[~fa, "disagreement"].values
        mwu_p = mannwhitneyu(x, y, alternative="greater").pvalue if len(x) and len(y) else np.nan
        mag_rows.append(dict(protein=gene, UniProt_ID=upid,
                             n_func=len(x), n_nonfunc=len(y),
                             mean_disagree_func=float(np.mean(x)) if len(x) else np.nan,
                             mean_disagree_nonfunc=float(np.mean(y)) if len(y) else np.nan,
                             mannwhitney_p_greater=mwu_p))

    by_protein = pd.DataFrame(by_protein)
    mag = pd.DataFrame(mag_rows)

    pooled_rows = []
    for dname in DEFS:
        for gran, key in (("assay_position", "ap"), ("position", "pos")):
            (a, b), (c, d) = pooled[dname][key]
            r = fisher_row(int(a), int(b), int(c), int(d))
            r.update(definition=dname, unit=gran)
            pooled_rows.append(r)
    pooled_df = pd.DataFrame(pooled_rows)

    front = ["protein", "UniProt_ID", "definition", "unit", "n_universe", "n_func", "n_high"]
    by_protein = by_protein[front + [c for c in by_protein.columns if c not in front]]
    pfront = ["definition", "unit"]
    pooled_df = pooled_df[pfront + [c for c in pooled_df.columns if c not in pfront]]

    by_protein.to_csv(os.path.join(OUTDIR, "enrichment_by_protein.csv"), index=False)
    pooled_df.to_csv(os.path.join(OUTDIR, "enrichment_pooled.csv"), index=False)
    mag.to_csv(os.path.join(OUTDIR, "magnitude_by_protein.csv"), index=False)

    pd.set_option("display.width", 220, "display.max_columns", 40)
    print("\n== Per-protein enrichment, ALLOSTERIC definition (assay x position) ==")
    sub = by_protein[(by_protein.definition == "ALLOSTERIC") & (by_protein.unit == "assay_position")]
    print(sub[["protein", "n_universe", "n_func", "n_high", "high_func",
               "pct_high_func", "pct_rest_func", "odds_ratio", "p_greater"]].to_string(index=False))
    print("\n== Pooled across the 5 proteins (each source dir, and their union) ==")
    print(pooled_df[["definition", "unit", "high_func", "high_nonfunc", "pct_high_func",
                     "pct_rest_func", "odds_ratio", "or_ci_lo", "or_ci_hi", "p_greater"]].to_string(index=False))
    print("\n== Continuous disagreement magnitude (allosteric-functional vs rest) ==")
    print(mag.to_string(index=False))
    print(f"\nwrote 3 csv -> {OUTDIR}")


if __name__ == "__main__":
    main()
