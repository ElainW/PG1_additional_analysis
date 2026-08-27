# Choosing a PDB structure for allosteric / distance analyses

How the Lehner group selected structures in three allosteric-landscape papers, and
a procedure to extrapolate that choice to other proteins with many PDBs. This is
the rationale behind the `literature_pdb` / `recommended_pdb` / `conformational_switch`
columns in `results/generalization_candidate_proteins.csv` (script `18`).

## What the papers did

**Liao & Lehner 2026, *Nat Commun* — "Allostery is a widespread cause of LOF variant pathogenicity."**
The core metric is the **minimal side-chain heavy-atom distance from each residue to the
bound orthosteric ligand**, and orthosteric/allosteric residues are defined by **van der
Waals contact** with a bound substrate/modulator (ChimeraX VDW-overlap). The structure must
therefore be **holo with the substrate (or a substrate-analog) in the active site**:
- GCK **1V4S** — glucose (GLC, orthosteric) + allosteric activator MRK
- PTEN **1D5R** — L(+)-tartrate (TLA), a substrate-analog at the catalytic site
- PTP1B 1PTY (phosphotyrosine), CASP1 2HBQ, CHK1 2E9N (inhibitor 76A), IDH1 6IOO (NADP)

Regulator-type caveat: for GCK they annotated allosteric residues from the literature
(not crystal contacts) "as its allosteric regulator (MRK) is also an activator," and they
**excluded prothrombin because its regulator is an activator, not an inhibitor.**

**Beltran & Lehner 2026, *Sci Adv* — "The allosteric landscape of the Src kinase."**
Src is modeled as a three-state switch (unfolded / folded-inactive / folded-active), so they
used **one structure per conformational state**:
- **2SRC** — inactive/autoinhibited state, with **AMP-PNP** (ATP analog) bound. Used for
  active-site distances (to the nucleotide, catalytic D388, substrate-positioning P428) and rSASA.
- **1Y57** — active state.
- `getcontacts` was run on **both** structures (after PyMOL `h_add`) to classify every residue
  as **active-only / inactive-only / swapping / nonspecific** contact. The *state comparison*
  is what reveals conformational allostery.

**Weng, Faure & Lehner 2024, *Nature* — KRAS energetic/allosteric landscape.**
The phenotype is **effector/partner binding**, so the "functional site" is the binding
interface, defined as residues **within 5 Å of the bound partner** — hence a **co-complex**
structure per partner (Fig 4/5 legends):
> "Binding interface is defined by indicated binding partner distance < 5 Å. PDB IDs:
> **RAF1, 6VJJ**; PIK3CG, 1HE8; RALGDS, 1LFD; SOS1, 1NVW; DARPin K27, 5O2S; DARPin K55, 5O2T."

State matters here too: effectors (RAF1 → 6VJJ) bind **GTP-active** KRAS, while "the DARPin
K27 specifically recognizes **inactive GDP-bound KRAS**" (Fig 5b panels `KRAS-GTP/RAF1` vs
`KRAS-GDP/DARPin K27`). So **6VJJ = the KRAS–RAF1 effector complex** (GTP/active state).

## The selection logic, generalized (priority order)

1. **Functional-site match.** Pick a structure where the site the phenotype reports on is
   *defined by a bound molecule*: substrate/substrate-analog for catalytic activity
   (GCK glucose, PTEN tartrate, DHFR folate), nucleotide for kinases/GTPases (AMP-PNP,
   GTP-analog), effector/peptide for a **binding** phenotype (KRAS–RAF1, MHC–peptide). This
   is stronger than "any ligand" — it must be the **orthosteric** ligand or the **partner**.
2. **State.** If the protein is a conformational switch (kinase active/inactive; GTPase
   GDP/GTP; transporter in/out; channel open/closed; allosteric enzyme T/R), pick **one
   canonical structure per relevant state** and compare their contacts (active-only /
   inactive-only / swapping). Match the state to the phenotype (GTP for effector binding,
   GDP for a GDP-state binder).
3. **Coverage + quality.** Must cover the mutated region; prefer good resolution, correct
   species/ortholog, and the biological oligomeric assembly.
4. **Canonical / reference.** Prefer the field-standard structure used in prior mechanistic
   work (2SRC, 1Y57, 1V4S, 1D5R are the reference structures), which also aids comparability.

Distance metric follows the site type: **min side-chain heavy-atom distance to the ligand**
(enzymes) or **distance/contact to the bound partner** (binding phenotypes).

## Applying it to the candidate set (`18`)

The CSV encodes this as: `conformational_switch` (name-keyword hint: kinase/GTPase/channel/
transporter/ATPase), `literature_pdb` (the four paper choices, ground truth), and
`recommended_pdb` = `literature_pdb` → else `best_holo_pdb` → else `best_pdb`.

| protein | switch? | functional site | recommended |
|---|---|---|---|
| GCK (HXK4) | allosteric enzyme | glucose (+activator) | **1V4S** (lit) |
| PTEN | single state | substrate-analog at catalytic site | **1D5R** (lit) |
| SRC | kinase → 2 states | ATP-analog (inactive) + active | **2SRC;1Y57** (lit) |
| KRAS (RASK) | GTPase → GDP/GTP | effector complex, <5 Å to partner | **6VJJ** (lit, RAF1) |
| DHFR (DYR) | single | folate + NADP | 8UW0 |
| OXDA / CP2C9 | single | FAD / heme+substrate | 1C0P / 1R9O |
| KCNE1 / KCNJ2 / SLC22A1 / NPC1 | channel/transporter → states | (state pair) | best_holo (state-blind) |

**Limitations of the automated columns.** `best_holo_pdb` flags *any* qualifying ligand — it
cannot tell an orthosteric substrate from a surface fragment (e.g. PTEN's functional tartrate
TLA is on the additive blacklist, so `best_holo_pdb` is blank and only `literature_pdb`/
`recommended_pdb` gives 1D5R), and it has **no state dimension** (for switches the automated
pick returns a single structure; SRC/KRAS need a *pair* or a state-matched complex). Treat
`recommended_pdb` as the curated pick and `best_holo_pdb`/`best_pdb` as heuristic fallbacks.
