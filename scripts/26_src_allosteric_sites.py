"""Assemble the Src allosteric-site table by parsing Beltran et al. 2026
(Sci Adv, 'The allosteric landscape of the Src kinase') + its supplement.
Numbering = SRC_HUMAN / UniProt P12931 (verified: Y419, Y530, D389, K298 ...).
Writes results/allosteric_map/SRC_HUMAN.tsv (one row per residue, aggregated).
"""
import json, os
from collections import defaultdict

import common as C
ANALYSIS = C.ANALYSIS
SEQ = json.load(open(os.path.join(ANALYSIS, "data/uniprot/raw_json/SRC_HUMAN.json")))["sequence"]["value"]
OUT = os.path.join(ANALYSIS, "results/allosteric_map/SRC_HUMAN.tsv")

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
AS, A2, AL = "active_site", "allosteric_second_shell", "allosteric"

# ---- Fig 3F: previously described allosteric network connecting ATP & substrate sites (ref 52) ----
f3f = {298:AS,313:AL,317:AL,328:A2,387:A2,391:AS,407:AS,408:A2,409:A2,
       412:A2,430:A2,431:A2,457:A2,465:AL}
for p,c in f3f.items():
    add([p], "Fig3F", cls=c, direction=INH, region="ATP-substrate allosteric network (ref 52)")

# ---- Fig 3H: activation loop (per-residue stars read from heatmap; cross-checked
#      against Fig3F for overlapping residues) ----
f3h = {406:A2,407:AS,408:A2,409:A2,410:A2,411:A2,412:A2,419:A2,421:AL,424:AL,
       427:A2,428:AS,429:A2,430:A2,431:A2,432:A2,433:AL,435:AL}
for p,c in f3h.items():
    add([p], "Fig3H", cls=c, direction=INH, region="activation loop")

# ---- Fig 3I: kinase spines (R spine + C spine ddGa) ----
# R spine: strong inactivating regardless of abundance
rspine = {317:AL,328:A2,387:A2,408:A2}
for p,c in rspine.items():
    add([p], "Fig3I_Rspine", cls=c, direction=INH, region="R-spine (regulatory spine)")
# C spine: only ATP-contacting residues (A296,V284,L396) enriched for inactivating (red star)
add([284,296,396], "Fig3I_Cspine", cls=AS, direction=INH,
    region="C-spine (catalytic spine); direct ATP contact")
add([349,454], "Fig3I_Cspine", region="C-spine (catalytic spine)",
    note="C-spine; small ddGa, effect on activity ~ explained by fold stability")

# ---- Fig S4 / p6: 11 major activating allosteric sites (text layer of supplement p8) ----
fs4 = [368,378,381,443,444,445,446,458,506,528,530]
for p in fs4:
    add([p], "FigS4", cls=AL, direction=ACT, region="C-lobe (major activating allosteric site)")
# αF pocket subset explicitly named as activating hotspot (p6)
add([381,443,444,445,506], "p6_activating", region="alphaF pocket (activating hotspot)",
    note="alphaF pocket, previously defined activating hotspot (ref 44)")

# ---- p6: inhibitory allosteric mutations, named residues ----
add([310,313,314,317,320], "p6_inhibitory", direction=INH, region="helix alphaC",
    note="helix alphaC inhibitory cluster; E313 salt bridge with K298 in active state")
add([419], "p6_inhibitory", direction=INH, region="activation loop",
    note="Y419 autophosphorylation site; locks Src active when phosphorylated")

# ---- p8: electrostatic switch network (dynamic/'swapping' contacts, extremely detrimental) ----
sw = {407:"inactive D407-K298 -> breaks on activation",
      298:"switch: D407-K298 (inactive) -> E313-K298 (active)",
      313:"switch: E313-R412 (inactive) -> E313-K298 (active)",
      412:"switch: E313-R412 (inactive) -> R412-Y419 (active)",
      389:"inactive D389-Y419 (catalytic Asp of HRD)",
      419:"switch: D389-Y419 (inactive) -> R412-Y419 (active)"}
for p,n in sw.items():
    add([p], "p8_dynamic_contacts", direction=INH,
        region="electrostatic switch network", note=n)

# ---- p10: dynamic C-terminal regulatory tail (full-length vs KD) ----
# All tail residues (>526) are designated 'allosteric': they lie far from the
# ligand/substrate sites yet still modulate kinase activity.
add([530], "p10_Cterm_tail", cls=AL, direction=ACT, region="C-terminal regulatory tail",
    note="inhibitory phosphosite; mutations more activating in full-length Src (SH2 release)")
add([527,528,531], "p10_Cterm_tail", cls=AL, direction=ACT, region="C-terminal regulatory tail",
    note="adjacent to Y530; stronger activating effects in full-length Src")
add([529,532,533,535], "p10_Cterm_tail", cls=AL, direction=INH, region="C-terminal regulatory tail",
    note="stronger inhibitory effects in full-length Src (autoinhibitory tail motif)")
add([395,397], "p10_FL_cluster", cls=AL, direction=ACT, region="C-spine (catalytic spine)",
    note="C-spine I395/V397; part of C-lobe cluster more activating in full-length Src")

# ---- Fig 5C: known Src inhibitory pockets (pocket-forming residues) ----
P13 = [284,285,296,297,298,326,339,340,341,342,396,406,407,408]
P7  = [453,456,457,460,462,463,464,465,466,467,484,485,486,487,488,495,498,502]
for p in P13:
    add([p], "Fig5C", direction=INH, pocket="P13(ATP)",
        note="orthosteric ATP-binding pocket (competitive-inhibitor site)")
for p in P7:
    add([p], "Fig5C", direction=INH, pocket="P7(G-loop)",
        note="P7 G-loop pocket, targeted by small-molecule inhibitors")
# ---- Fig 5E: novel allosteric Src pockets (all inhibitory-enriched) ----
P4 = [316,317,319,320,325,379,380,381,383,384,385,386,387,408,411]
P2 = [280,391,419,425,426,427,428,429,430,431,464,465,471]
P1 = [317,320,321,322,323,324,325,326,375,376,377,379,380,403,404,405,406,408]
P15= [421,422,423,424,425,426,427,428,429,470,471,472,473,474,475,476]
for name,lst in [("P4",P4),("P2",P2),("P1",P1),("P15",P15)]:
    for p in lst:
        add([p], "Fig5E", direction=INH, pocket=name,
            note="novel allosteric pocket enriched for inhibitory mutations")

# C-terminal regulatory tail: every residue beyond the kinase domain (>526) is
# allosteric by definition here (distal to ligand/substrate, still modulates activity).
for p in list(R):
    if p > 526:
        R[p]["class"].add(AL)

# ---------- aggregate & write ----------
def pick_class(s):
    for c in (AS, A2, AL):          # priority: active > 2nd shell > allosteric
        if c in s: return c
    return ""
def join_dir(s):
    return ";".join(x for x in (INH, ACT) if x in s)

cols = ["UniProt_ID","position","wt_aa","residue","site_class","direction",
        "region","pockets","evidence","notes"]
rows = []
for p in sorted(R):
    r = R[p]
    wt = SEQ[p-1] if 1 <= p <= len(SEQ) else "?"
    direction = join_dir(r["dir"])
    # residues only known via inhibitory pockets inherit that direction (already INH)
    rows.append([
        "SRC_HUMAN", p, wt, f"{wt}{p}", pick_class(r["class"]), direction,
        " | ".join(sorted(r["region"])), ";".join(sorted(r["pockets"])),
        ";".join(sorted(r["evidence"])), " | ".join(sorted(r["notes"])),
    ])

with open(OUT, "w") as f:
    f.write("# Src allosteric-site atlas parsed from Beltran, Naqvi, Faure & Lehner 2026,\n")
    f.write("# 'The allosteric landscape of the Src kinase', Sci Adv 12, eaea2726.\n")
    f.write("# Numbering: SRC_HUMAN / UniProt P12931. site_class from Fig3 stars\n")
    f.write("# (active_site=red, allosteric_second_shell=salmon, allosteric=yellow).\n")
    f.write("# direction: inhibitory (loss of activity, ddGa>0) vs activating (ddGa<0).\n")
    f.write("\t".join(cols) + "\n")
    for row in rows:
        f.write("\t".join(str(x) for x in row) + "\n")

print(f"wrote {OUT}  ({len(rows)} residues)")
na = sum(1 for r in rows if r[5]=="activating"); ni = sum(1 for r in rows if r[5]=="inhibitory")
nb = sum(1 for r in rows if ";" in r[5])
print(f"  inhibitory={ni} activating={na} both={nb}")
