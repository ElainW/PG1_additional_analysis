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

NON_MODEL_COLS = {"mutant", "mutated_sequence", "DMS_score", "DMS_score_bin"}

# Stable display order for the taxa
TAXA = ["Human", "Eukaryote", "Prokaryote", "Virus"]

# Stable display order for model types (most -> least represented)
MODEL_TYPE_ORDER = [
    "MSA",
    "Single sequence",
    "Structure",
    "Single sequence & Structure",
    "Structure & MSA",
    "Single sequence, Structure & Function annotations",
]


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
