"""
Shared utilities for the ProteinGym-v1 by-position analysis.

Paths, config loading, model->model_type map, DMS_id->taxon map, mutant parsing,
and a fast vectorized Spearman helper (rank + Pearson).

IMPORTANT — directionality note
--------------------------------
config.json records a `directionality` field per model (raw model output:
+1 => higher score = higher fitness, -1 => higher score = lower fitness).

Empirically (see 00_check_directionality.py) the score columns in
O2/20250613/*.csv are ALREADY oriented so that higher score = higher predicted
fitness for every model -- all 8 directionality=-1 models (EVE, DeepSequence,
RSALOR, ...) correlate strongly *positively* with DMS_score. The config
directionality was applied upstream when these merged files were produced.

Therefore we do NOT re-apply directionality here: scores are used as-is so that
"higher = fitter" holds uniformly across models. This matches ProteinGym's own
performance_DMS_benchmarks.py, which calls spearmanr(DMS_score, score) directly.
"""
import os
import re
import json
import glob
import numpy as np
import pandas as pd

PG2 = os.path.expanduser("~/Desktop/MarksLab/PG2")
SCORES_DIR = os.path.join(PG2, "O2", "20250613")
CONFIG_PATH = os.path.join(PG2, "ProteinGym", "config.json")
REF_PATH = os.path.join(PG2, "ProteinGym", "reference_files", "DMS_substitutions.csv")

ANALYSIS = os.path.join(PG2, "analysis")
RESULTS = os.path.join(ANALYSIS, "results")
FIGURES = os.path.join(ANALYSIS, "figures")
DATA_DIR = os.path.join(ANALYSIS, "data")
WT_DMS_PATH = os.path.join(DATA_DIR, "DMS_wt_score.csv")

NON_MODEL_COLS = {"mutant", "mutated_sequence", "DMS_score", "DMS_score_bin"}

# Stable display order for the taxa
TAXA = ["Human", "Eukaryote", "Prokaryote", "Virus"]

# ---------------------------------------------------------------------------
# Assays excluded from every cross-assay / conservation / disagreement-site
# analysis. P53_HUMAN_Kotler_2018 sits on a DIFFERENT target_seq than the three
# P53 Giacomelli assays -- they are identical except codon 72 (Giacomelli R72,
# Kotler P72, the TP53 R72P polymorphism) -- so its positions are NOT in the same
# coordinate frame and it must not be pooled with the Giacomelli assays when
# counting shared high-disagreement sites, correlating disagreement vs MSA
# conservation, or building the P53 structure/colocalization panels. Filtered at
# the master-site-table source (15) so it propagates to every downstream reader
# (23/25/27/28/30/31/32/33/34-45); the hardcoded P53 assay lists in 40/46 drop it
# directly. Kotler may still be scored in single-assay/global-217 contexts (e.g.
# 01 spearman), where it is not pooled with Giacomelli.
EXCLUDE_DMS = {"P53_HUMAN_Kotler_2018"}


def drop_excluded_dms(df, col="DMS_id"):
    """Remove EXCLUDE_DMS rows from a frame keyed on a DMS_id column."""
    return df[~df[col].isin(EXCLUDE_DMS)].copy()

# Canonical colour per assay category (coarse_selection_type). Matches the fig9
# barplot's default-cycle order (categories sorted: Activity, Binding, Expression,
# OrganismalFitness, Stability -> tab10 C0..C4) so every by-category figure agrees.
CATEGORY_COLORS = {
    "Activity":          "#1f77b4",   # C0 blue
    "Binding":           "#ff7f0e",   # C1 orange
    "Expression":        "#2ca02c",   # C2 green
    "OrganismalFitness": "#d62728",   # C3 red
    "Stability":         "#9467bd",   # C4 purple
}

# Stable display order for model types (most -> least represented)
MODEL_TYPE_ORDER = [
    "MSA",
    "Single sequence",
    "Structure",
    "Single sequence & Structure",
    "Structure & MSA",
    "Single sequence, Structure & Function annotations",
]


# ---------------------------------------------------------------------------
# Results layout. Outputs are foldered by analysis under results/. rp() maps a
# base filename to its subfolder so writers and readers stay in sync; files with
# no rule (generalization_*, and the high_disagreement_sites/ & MSA_conservation/
# subdirs) resolve to results/ root, where they already live.
# ---------------------------------------------------------------------------
_RESULTS_SUBDIR_RULES = [
    ("spearman_",            "spearman_taxa"),
    ("kruskal_",             "spearman_taxa"),
    ("persite_",             "persite_analysis"),
    ("family_best",          "persite_analysis"),
    ("hxk4_",                "hxk4_pilot_analysis"),
    ("site_",                "HXK4_HUMAN_w99_deepdive"),
    ("gck_",                 "gck_3d"),
    ("assay_position_",      "assay_selection"),
    ("abundance_physchem",   "abundance_physchem"),
    ("collapse_",            "collapse_withinsite"),
    ("disease_gof_lof",      "disease_gof_lof"),
    ("directionality_check", "qc_directionality"),
]


def results_subdir(name):
    """Subfolder (under results/) an output filename belongs to, '' for root."""
    for prefix, sub in _RESULTS_SUBDIR_RULES:
        if name.startswith(prefix):
            return sub
    return ""


def rp(name):
    """Resolve an output filename to its full path under results/, creating the
    subfolder if needed. Filenames with no folder rule resolve to results/ root."""
    folder = os.path.join(RESULTS, results_subdir(name))
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, name)


def load_config():
    with open(CONFIG_PATH) as f:
        cfg = json.load(f)
    return cfg["model_list_zero_shot_substitutions_DMS"]


def model_meta():
    """Return (models list, dict model->model_type, dict model->directionality)."""
    cfg = load_config()
    models = list(cfg.keys())
    mtype = {m: cfg[m]["model_type"] for m in models}
    direction = {m: cfg[m]["directionality"] for m in models}
    return models, mtype, direction


def load_reference():
    ref = pd.read_csv(REF_PATH)
    return ref


def dms_to_taxon():
    ref = load_reference()
    return dict(zip(ref["DMS_id"], ref["taxon"]))


def load_wt_dms():
    """{DMS_id: wt_score} — per-assay wild-type DMS score, from the
    version-controlled data/DMS_wt_score.csv (single source of truth)."""
    wt = pd.read_csv(WT_DMS_PATH)
    return dict(zip(wt["DMS_id"], wt["wt_score"].astype(float)))


def list_score_files():
    return sorted(glob.glob(os.path.join(SCORES_DIR, "*.csv")))


def dms_id_from_path(path):
    return os.path.basename(path)[:-4]


_POS_RE = re.compile(r"^[A-Z](\d+)[A-Z]$")


def parse_position(mutant):
    """Return integer residue position for a single substitution like 'I291A',
    or np.nan for multi-mutants / unparseable strings."""
    m = _POS_RE.match(str(mutant))
    return int(m.group(1)) if m else np.nan


def spearman_with_target(df, target_col, score_cols):
    """Vectorized Spearman of each score column against target_col.

    Computes rank of every column then Pearson correlation of ranks with the
    target's ranks, column by column (handles NaNs per column independently).
    Returns dict score_col -> spearman (np.nan if degenerate).
    """
    out = {}
    t = df[target_col]
    t_rank_full = t.rank()
    for c in score_cols:
        s = df[c]
        mask = t.notna() & s.notna()
        n = int(mask.sum())
        if n < 3:
            out[c] = np.nan
            continue
        # rank within the masked subset for correctness
        tr = t[mask].rank()
        sr = s[mask].rank()
        if tr.std(ddof=0) == 0 or sr.std(ddof=0) == 0:
            out[c] = np.nan
            continue
        out[c] = float(np.corrcoef(tr.values, sr.values)[0, 1])
    return out


def fisher_mean(corrs):
    """Average correlations via Fisher z-transform, ignoring NaN and clipping."""
    a = np.asarray(corrs, dtype=float)
    a = a[~np.isnan(a)]
    if a.size == 0:
        return np.nan
    a = np.clip(a, -0.999999, 0.999999)
    z = np.arctanh(a)
    return float(np.tanh(z.mean()))
