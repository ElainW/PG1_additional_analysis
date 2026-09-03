# Scripts

Analysis scripts for the ProteinGym by-position / structure-feature work. Two
environments are used.

## 1. Base analysis env (scripts `00`–`47`)

The original by-position analyses (`make all`, driven by the [`Makefile`](../Makefile))
run under a plain Python 3 with the scientific stack:

```bash
mamba create -n pg2 -c conda-forge python=3.11 \
    pandas numpy scipy matplotlib seaborn biopython
conda activate pg2
make all          # runs scripts 00..34 in dependency order
```

These read the tracked `data/` cache and write to `results/` / `figures/`.

## 2. Structure-feature env `structfeat` (scripts `48`–`54`)

The structure-feature pipeline needs external structural-biology tools (FreeSASA
CLI, PyMOL, VMD-python + GetContacts). Set up a dedicated env once:

```bash
# core tools (pymol h_add, vmd-python for getcontacts, freesasa python, io)
mamba create -n structfeat -c conda-forge python=3.10 \
    freesasa pymol-open-source vmd-python numpy pandas biopython
conda activate structfeat
```

### FreeSASA command-line binary (built from source)

conda-forge ships only the FreeSASA *Python* package, not the `freesasa` CLI that
script `51` calls. Build the CLI into the env:

```bash
cd /tmp
curl -sL -o freesasa-2.1.2.zip \
    https://github.com/mittinatten/freesasa/releases/download/2.1.2/freesasa-2.1.2.zip
unzip -q freesasa-2.1.2.zip && cd freesasa-2.1.2
./configure --disable-xml --disable-json --prefix=$CONDA_PREFIX
make && make install
freesasa --version      # -> FreeSASA 2.1.2
```

(The freesasa Python API reproduces the identical calculation if you prefer not to
build the CLI, but script `51` invokes the binary.)

### GetContacts

Not on PyPI — clone the repo and install its Python deps into the env:

```bash
git clone https://github.com/getcontacts/getcontacts.git ~/tools/getcontacts
pip install cython scipy matplotlib scikit-learn networkx expiringdict
python ~/tools/getcontacts/get_static_contacts.py --help    # smoke test
```

Script `53` expects `get_static_contacts.py` at `~/tools/getcontacts/`; if you
clone elsewhere, edit `GETCONTACTS` at the top of `53_getcontacts.py`.

### Run order (structure features)

Run from the repo root with `structfeat` active. `48`–`50` prepare structures
(no external tools); `51`–`54` are the analyses.

```bash
conda activate structfeat
python scripts/48_structure_ligands.py          # ligand/substrate inventory
python scripts/49_clean_6wvh.py                 # VKORC1-only cleaned + single chain
python scripts/50_clean_structures.py           # clean the other structures
python scripts/51_rsasa.py                       # FreeSASA rSASA  -> rsasa/
python scripts/52_secondary_structure.py         # UniProt SS      -> secondary_structure/ (needs network)
python scripts/53_getcontacts.py                 # h_add + getcontacts -> contacts/
python scripts/54_contact_state_classification.py  # active/inactive residue classes
```

Inputs are the curated structures in
[`../data/structures/pdb/`](../data/structures/pdb/) (see that README for selection,
cleaning, and numbering offsets); outputs go to
[`../results/structure_features/`](../results/structure_features/) (see its README for
what each file holds and how to read the results).

## 3. High-disagreement × functional-region analyses (scripts `55`–`56`)

Run under the base `pg2` env; they read tracked outputs only (the disagreement
pool from `15`, the tail from `27`, the allosteric maps from `26`/`47`, ClinVar
from `data/clinvar/`, and the secondary structure from `52`), so no structural
tools are needed.

```bash
conda activate pg2
python scripts/55_high_disagreement_allosteric_enrichment.py   # enrichment test (exploratory)
python scripts/56_high_disagreement_site_characterization.py   # paper table: where the tail lands
```

- `55` — for the five allosteric-map proteins (CYP2C9, GCK, KRAS, PTEN, SRC),
  tests whether the high-disagreement tail is enriched in functional regions
  (allosteric map / UniProt binding±2 / dir2 ClinVar / their union) vs the rest
  of the assayed sites (one-sided Fisher, per protein + pooled). Exploratory —
  finds no enrichment; kept for the record. -> `results/high_disagreement_allosteric_enrichment/`
- `56` — descriptive characterization for the paper: per protein, the number of
  high-disagreement sites and the fraction landing in terminal (beginning/end),
  loop/turn, functional, or structured-core regions (mutually exclusive, plus raw
  non-exclusive flags). -> `results/high_disagreement_site_characterization/`
