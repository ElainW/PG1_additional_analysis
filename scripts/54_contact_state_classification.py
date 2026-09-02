"""
Active/inactive residue classification from state-specific contacts.

For each protein with multiple conformations (HXK4, RASK, SRC) we compare the
state-specific contacts (the 4 types from script 53: salt bridge, pi-cation,
sc-sc H-bond, sc-bb H-bond) between the active and inactive structures and label
each contacting residue (following Beltran & Lehner / Weng, Faure & Lehner):

  active       contacts only in the active state
  inactive     contacts only in the inactive state
  nonspecific  identical contact partner set in both states
  swapping     contacts in both states but to (some) different partners

A "contact" is an unordered residue pair. Residue numbers are converted to
UniProt numbering FIRST (per-structure offset), which is essential for SRC where
2SRC = auth+2 but 1Y57 = auth+3 (the same author number is a different residue).

State pairs:
  HXK4  active=1V4S (closed)   inactive=1V4T (super-open)
  SRC   active=1Y57            inactive=2SRC
  RASK  active=6VJJ (ON/GTP)   inactive=5O2S (OFF/GDP)     apo=1NVW (3rd state, reported separately)

Outputs (results/structure_features/contacts/classification/):
  <PROTEIN>_residue_classes.csv    per-residue label + partner sets (UniProt numbering)
  classification_summary.csv       counts per class per protein
  RASK_apo_contacts_note.csv       apo (1NVW) contact partners, for the 3rd-state view
"""
import os
import sys
from collections import defaultdict

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

CONTACTS = os.path.join(C.RESULTS, "structure_features", "contacts")
FOURD = os.path.join(CONTACTS, "4types")
OUTDIR = os.path.join(CONTACTS, "classification")
os.makedirs(OUTDIR, exist_ok=True)

# UniProt = author + OFFSET (canonical UniProt).
# NB: 2SRC and 1Y57 share the SAME author numbering; canonical P12931 alignment is
# +3 for BOTH (author 416->Y419 activation loop, 527->Y530 regulatory tail). The
# 2SRC DBREF cites a 535-aa reference (+2), which would misalign the two Src states
# by one residue, so we use +3 for both.
OFFSET = {"1V4S": 1, "1V4T": 1, "2SRC": 3, "1Y57": 3, "6VJJ": 0, "5O2S": 0, "1NVW": 0}

PAIRS = [
    ("HXK4_HUMAN", "1V4S", "1V4T"),
    ("SRC_HUMAN",  "1Y57", "2SRC"),
    ("RASK_HUMAN", "6VJJ", "5O2S"),
]


def load_contacts(label):
    """Return (partners, resname) in UniProt numbering.
    partners[u] = set of UniProt partner positions u contacts (4-type)."""
    off = OFFSET[label]
    partners = defaultdict(set)
    resname = {}
    df = pd.read_csv(os.path.join(FOURD, f"{label}_4types.tsv"), sep="\t")
    for _, r in df.iterrows():
        u1 = int(r.resnum1) + off
        u2 = int(r.resnum2) + off
        if u1 == u2:
            continue
        partners[u1].add(u2)
        partners[u2].add(u1)
        resname[u1] = r.res1
        resname[u2] = r.res2
    return partners, resname


def classify(active_p, inactive_p):
    labels = {}
    for u in set(active_p) | set(inactive_p):
        A = active_p.get(u, set())
        I = inactive_p.get(u, set())
        if A and not I:
            labels[u] = "active"
        elif I and not A:
            labels[u] = "inactive"
        elif A == I:
            labels[u] = "nonspecific"
        else:
            labels[u] = "swapping"
    return labels


def fmt(s):
    return ",".join(str(x) for x in sorted(s)) if s else ""


def main():
    summ = []
    for prot, active, inactive in PAIRS:
        ap, arn = load_contacts(active)
        ip, irn = load_contacts(inactive)
        resname = {**irn, **arn}
        labels = classify(ap, ip)

        rows = []
        for u in sorted(labels):
            A, I = ap.get(u, set()), ip.get(u, set())
            rows.append(dict(
                uniprot_pos=u, resname=resname.get(u, "?"), classification=labels[u],
                n_active=len(A), n_inactive=len(I), n_shared=len(A & I),
                active_partners=fmt(A), inactive_partners=fmt(I)))
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(OUTDIR, f"{prot}_residue_classes.csv"), index=False)

        vc = df.classification.value_counts().to_dict()
        summ.append(dict(protein=prot, active_struct=active, inactive_struct=inactive,
                         n_residues=len(df),
                         active=vc.get("active", 0), inactive=vc.get("inactive", 0),
                         nonspecific=vc.get("nonspecific", 0),
                         swapping=vc.get("swapping", 0)))
        print(f"{prot:12s} active={active} inactive={inactive}  "
              f"residues={len(df):3d} | active={vc.get('active',0):3d} "
              f"inactive={vc.get('inactive',0):3d} nonspecific={vc.get('nonspecific',0):3d} "
              f"swapping={vc.get('swapping',0):3d}")

    pd.DataFrame(summ).to_csv(os.path.join(OUTDIR, "classification_summary.csv"), index=False)

    # KRAS apo (1NVW) as a 3rd-state reference: list its contacting residues
    ap_apo, rn_apo = load_contacts("1NVW")
    apo_rows = [dict(uniprot_pos=u, resname=rn_apo.get(u, "?"),
                     apo_partners=fmt(ap_apo[u])) for u in sorted(ap_apo)]
    pd.DataFrame(apo_rows).to_csv(
        os.path.join(OUTDIR, "RASK_apo_contacts_note.csv"), index=False)
    print(f"\nRASK apo (1NVW) 3rd-state contacts: {len(apo_rows)} residues "
          f"(RASK_apo_contacts_note.csv)")
    print(f"Wrote per-protein *_residue_classes.csv + classification_summary.csv to\n  {OUTDIR}")


if __name__ == "__main__":
    main()
