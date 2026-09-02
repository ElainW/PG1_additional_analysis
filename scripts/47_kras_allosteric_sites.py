"""Assemble the KRAS allosteric-site table by curating figures of
Weng, Faure, Escobedo & Lehner 2023, 'The energetic and allosteric landscape
for KRAS inhibition', Nature 626, 643-652 (2023) (doi:10.1038/s41586-023-06954-0).

Numbering = RASK_HUMAN / UniProt P01116 (KRAS4B, aa 1-188; the reference
sequence used by the Weng deep-mutagenesis study, verified: G12, D33, E37,
D38, Y40, T58, G60, R68, F78, P110, F141, I163 ...).
Only residue NUMBERS are collected (mutation identities are ignored), one row
per residue aggregated across the sources below.
Writes results/allosteric_map/RASK_HUMAN.tsv .

Sources curated (see 26_src_allosteric_sites.py for the sibling Src table):
 - Ostrem et al. Nature 2013 switch-II inhibitory allosteric pocket (G12C site)
 - Fig 2d,e  : RAF1 binding-interface residues most important for binding
 - Fig 3b    : 18 major allosteric sites for RAF1 (10 nucleotide-pocket + 8 novel)
 - Fig 3f / main text (p6): the four allosterically-active surface pockets
 - Fig 5a & Extended Data Fig 5 : binding-interface contact residues, 6 partners
 - Fig 6a,b  : major allosteric sites across all six interaction partners
 - Extended Data Fig 2a : additional KRAS beta-sheet allosteric residues
"""
import json, os
from collections import defaultdict

import common as C
ANALYSIS = C.ANALYSIS
SEQ = json.load(open(os.path.join(ANALYSIS, "data/uniprot/raw_json/RASK_HUMAN.json")))["sequence"]["value"]
OUT = os.path.join(ANALYSIS, "results/allosteric_map/RASK_HUMAN.tsv")

# accumulator: pos -> dict of sets
R = defaultdict(lambda: {"class": set(), "dir": set(), "region": set(),
                         "pockets": set(), "evidence": set(), "notes": set()})

def add(positions, evidence, cls=None, direction=None, region=None, pocket=None, note=None):
    for p in positions:
        r = R[p]
        r["evidence"].add(evidence)
        if cls: r["class"].add(cls)
        if direction: r["dir"].add(direction)
        if region: r["region"].add(region)
        if pocket: r["pockets"].add(pocket)
        if note: r["notes"].add(note)

INH, ACT = "inhibitory", "activating"
IFACE, ALLO = "binding_interface", "allosteric"

# ---- KRAS beta-sheet strands / shells (Extended Data Fig. 2a): used as region tags ----
STRANDS = {  # (start,end) inclusive : label
    (35, 41): "beta-strand 2 (1st shell)",
    (54, 60): "beta-strand 3 (2nd shell)",
    (4, 10):  "beta-strand 1 (3rd shell)",
    (77, 83): "beta-strand 4 (4th shell)",
    (111, 117): "beta-strand 5 (5th shell)",
    (140, 143): "beta-strand 6 (6th shell)",
}
def strand_region(p):
    for (a, b), lab in STRANDS.items():
        if a <= p <= b:
            return lab
    return None
def functional_region(p):
    # coarse functional regions from Fig 1i annotation of KRAS
    r = []
    if 10 <= p <= 17: r.append("P-loop")
    if 25 <= p <= 40: r.append("switch I")
    if 57 <= p <= 76: r.append("switch II")
    if p > 166:       r.append("hypervariable region (HVR)")
    return r

# ---- Ostrem et al. Nature 2013: switch-II inhibitory allosteric pocket (KRAS-G12C) ----
OSTREM = [7, 9, 12, 58, 60, 68, 72, 78, 99, 100]
for p in OSTREM:
    add([p], "Ostrem2013", cls=ALLO, direction=INH, pocket="switchII pocket (Ostrem2013)",
        note="switch-II inhibitory allosteric pocket (KRAS-G12C covalent-inhibitor site, Ostrem et al. Nature 2013)")

# ---- Fig 2d,e: RAF1 binding-interface residues most important for binding ----
add([33, 36, 37, 38, 40, 41], "Fig2de", cls=IFACE, direction=INH,
    region="RAF1-binding interface",
    note="interface residue most important for RAF1 binding (charged E37/D38, hydrophobic I36/Y40, plus D33/R41)")

# ---- Fig 5a & Extended Data Fig 5: binding-interface contact residues (6 partners) ----
# union of KRAS residues contacting SOS1, RALGDS, RAF1, PIK3CG, DARPin K55, DARPin K27
INTERFACE = [1, 3, 5, 17, 18, 20, 21, 22, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 36,
             37, 38, 39, 40, 41, 42, 43, 44, 45, 50, 52, 54, 56, 57, 58, 59, 60, 61, 63,
             64, 65, 66, 67, 68, 69, 70, 71, 73, 74, 95, 102, 103, 105, 149, 153]
for p in INTERFACE:
    add([p], "Fig5a/EDFig5", cls=IFACE, direction=INH,
        region="binding interface (6 KRAS effectors/binders)",
        note="binding-interface contact residue for one or more of SOS1/RALGDS/RAF1/PIK3CG/DARPin-K55/DARPin-K27")

# ---- Fig 3b: 18 major allosteric sites for RAF1 binding ----
# 10 in the nucleotide (GTP) binding pocket (blue) ...
F3B_GTP = [15, 16, 17, 18, 28, 32, 34, 35, 57, 60]
for p in F3B_GTP:
    add([p], "Fig3b", cls=ALLO, direction=INH, pocket="nucleotide (GTP) pocket",
        region="major allosteric site for RAF1 (nucleotide-binding pocket)",
        note="major allosteric site for RAF1 binding, located in the physiological nucleotide-binding pocket")
# ... 8 additional novel major allosteric sites (yellow; also listed in text p2) ...
F3B_NOVEL = [7, 10, 54, 58, 59, 110, 141, 163]
for p in F3B_NOVEL:
    add([p], "Fig3b", cls=ALLO, direction=INH,
        region="major allosteric site for RAF1 (novel)",
        note="additional major allosteric site (V7,G10,D54,T58,A59,P110,F141,I163)")
# extra location notes for the novel sites (text p2)
add([54], "text_p2", note="D54 adjacent to the binding interface")
add([58], "text_p2", note="T58 connects the binding interface to the nucleotide-binding pocket")
add([59], "text_p2", note="A59 connects the binding interface to the nucleotide-binding pocket")

# ---- Fig 6a,b: major allosteric sites across all six interaction partners ----
# (blue = nucleotide-pocket, yellow = novel) enumerated by reading every labelled
#  blue/yellow point across the six Fig 6a panels (RAF1/PIK3CG/RALGDS/SOS1/K27/K55);
#  this recovers Y96, Q99 and I163, which the Fig 6b panel headers omit.
F6 = [10, 15, 16, 17, 18, 20, 28, 31, 32, 34, 35, 37, 54, 55, 57, 58, 59, 60, 63, 68,
      69, 71, 72, 75, 76, 78, 83, 96, 99, 110, 112, 117, 156, 163]
for p in F6:
    add([p], "Fig6", cls=ALLO, direction=INH,
        region="major allosteric site (across 6 binding partners)",
        note="major allosteric site for one or more of the six KRAS interaction partners")
# gamma-phosphate contacts (solid boxes; text: K16 & G60 contact the gamma-phosphate)
add([16, 34, 35, 60], "Fig6", note="gamma-phosphate-contacting residue")
# nucleotide-contacting residues (dashed boxes in Fig 6b)
add([15, 17, 18, 28, 32, 57, 59, 117], "Fig6", note="nucleotide-contacting residue")
# specificity switches: allosterically INCREASE binding to a specific partner
add([16, 55, 60, 156], "Fig6_specificity", direction=ACT,
    note="allosterically inhibits most partners but INCREASES binding to DARPin K27")
add([76], "Fig6_specificity", direction=ACT,
    note="E76 allosterically increases binding to DARPin K27 (little effect on others)")
add([54], "Fig6_specificity", direction=ACT,
    note="D54 inhibits four interactions but retains/enhances binding to PIK3CG and RALGDS")
add([71, 72], "Fig6_specificity",
    note="Y71/M72 specifically inhibit binding to DARPin K55")

# ---- Extended Data Fig 2a: additional KRAS beta-sheet allosteric residues ----
# residues with binding-energy changes not otherwise reported, with shell/strand text
ED2A = [55, 6, 4, 7, 8, 9, 10]
for p in ED2A:
    add([p], "EDFig2a", cls=ALLO, direction=INH,
        note="beta-sheet allosteric residue (Extended Data Fig. 2a)")

# ---- Fig 3f / main text (p6): the four allosterically-active surface pockets ----
# pocket-forming residues (structural annotation -> 'pockets' column, like the Src table);
# the allosterically-active subset per pocket is n=10/9/6/9 (Extended Data Fig. 2d-g).
POCKETS = {
    "P1 switch-I/II pocket":       [5, 6, 7, 39, 54, 55, 56, 70, 71, 72, 73, 74, 75],
    "P2 switch-II (sotorasib) pocket": [61, 62, 63, 64, 65, 90, 91, 92, 93, 94],
    "P3 C-terminal distal pocket": [97, 101, 107, 108, 109, 110, 111, 136, 137, 138, 139,
                                    140, 161, 162, 163, 164, 165, 166],
    "P4 behind effector-binding loop": [17, 21, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
                                        34, 35, 36, 37, 38, 39, 40, 57],
}
for name, lst in POCKETS.items():
    for p in lst:
        add([p], "Fig3f/EDFig2dg", direction=INH, pocket=name,
            note="residue lining an allosterically-active KRAS surface pocket")

# ---- attach region tags (beta-strand shells + functional regions) to every collected residue ----
for p in list(R):
    s = strand_region(p)
    if s: R[p]["region"].add(s)
    for fr in functional_region(p):
        R[p]["region"].add(fr)

# ---------- aggregate & write ----------
def join_set(s):
    return " | ".join(sorted(s))
def join_dir(s):
    return ";".join(x for x in (INH, ACT) if x in s)
def join_class(s):
    return ";".join(x for x in (IFACE, ALLO) if x in s)

cols = ["UniProt_ID", "position", "wt_aa", "residue", "site_class", "direction",
        "region", "pockets", "evidence", "notes"]
rows = []
for p in sorted(R):
    r = R[p]
    wt = SEQ[p-1] if 1 <= p <= len(SEQ) else "?"
    rows.append([
        "RASK_HUMAN", p, wt, f"{wt}{p}", join_class(r["class"]), join_dir(r["dir"]),
        join_set(r["region"]), join_set(r["pockets"]),
        ";".join(sorted(r["evidence"])), join_set(r["notes"]),
    ])

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write("# KRAS allosteric-site atlas curated from Weng, Faure, Escobedo & Lehner 2023,\n")
    f.write("# 'The energetic and allosteric landscape for KRAS inhibition', Nature 626, 643-652.\n")
    f.write("# Numbering: RASK_HUMAN / UniProt P01116 (KRAS4B). Residue numbers only.\n")
    f.write("# site_class: binding_interface (orthosteric effector contact) and/or allosteric.\n")
    f.write("# direction: inhibitory (mutations reduce binding) vs activating (specificity gain).\n")
    f.write("# pockets: structural pocket membership (Ostrem switch-II, GTP pocket, P1-P4 surface pockets).\n")
    f.write("# NOTE: position 153 is E153 in the reference sequence (paper labels it D153).\n")
    f.write("\t".join(cols) + "\n")
    for row in rows:
        f.write("\t".join(str(x) for x in row) + "\n")

print(f"wrote {OUT}  ({len(rows)} residues)")
na = sum(1 for r in rows if "activating" in r[5])
ni = sum(1 for r in rows if "inhibitory" in r[5])
nif = sum(1 for r in rows if IFACE in r[4])
nal = sum(1 for r in rows if ALLO in r[4])
print(f"  site_class: binding_interface={nif} allosteric={nal}")
print(f"  direction: inhibitory={ni} activating(any)={na}")
