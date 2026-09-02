# ============================================================
# KRAS nucleotide-state comparison: ON (3P) vs OFF (2P) vs apo (0P)
#   6VJJ chain A = KRAS-GppNHp (GNP)  + RAF1  -> ACTIVE / ON  (3P)
#   5O2S chain E = KRAS-GDP           + K27   -> INACTIVE/OFF (2P)
#   1NVW chain R = H-Ras nucleotide-free (SOS catalytic site) -> apo (0P)
#
# Run in PyMOL:  @visualize_kras_states.pml
# ============================================================

reinitialize
bg_color white
set ray_shadows, 0
set cartoon_transparency, 0.0

# --- fetch (legacy PDB guarantees chain IDs A/E/R) ---
fetch 6vjj, type=pdb, async=0
fetch 5o2s, type=pdb, async=0
fetch 1nvw, type=pdb, async=0

# --- put OFF and apo into 6VJJ's frame, fitting on the G-domain core ---
#     (align auto-rejects the moving switch loops as outliers)
align 5o2s and chain E and polymer, 6vjj and chain A and polymer
align 1nvw and chain R and polymer, 6vjj and chain A and polymer

# --- isolate the three Ras copies of interest (PDB id in the name) ---
create 6VJJ_ON_3P,  6vjj and chain A
create 5O2S_OFF_2P, 5o2s and chain E
create 1NVW_APO_0P, 1nvw and chain R
delete 6vjj
delete 5o2s
delete 1nvw

# --- strip all waters ---
remove solvent

hide everything
show cartoon, 6VJJ_ON_3P or 5O2S_OFF_2P or 1NVW_APO_0P
color grey80,  6VJJ_ON_3P
color skyblue, 5O2S_OFF_2P
color salmon,  1NVW_APO_0P

# --- highlight switch I (30-38) and switch II (60-76) ---
color marine, 6VJJ_ON_3P  and resi 30-38
color orange, 6VJJ_ON_3P  and resi 60-76
color blue,   5O2S_OFF_2P and resi 30-38
color red,    5O2S_OFF_2P and resi 60-76
color purple, 1NVW_APO_0P and resi 30-38
color yellow, 1NVW_APO_0P and resi 60-76

# --- the definitive ON/OFF/apo marker: the nucleotide ---
show sticks, (6VJJ_ON_3P and resn GNP) or (5O2S_OFF_2P and resn GDP)
show spheres, (6VJJ_ON_3P or 5O2S_OFF_2P) and resn MG
set sphere_scale, 0.3, resn MG
color green, resn MG
util.cnc (6VJJ_ON_3P and resn GNP) or (5O2S_OFF_2P and resn GDP)

# --- label the switch regions once each ---
set label_size, 16
label 6VJJ_ON_3P and name CA and resi 35, '"switch I"'
label 6VJJ_ON_3P and name CA and resi 68, '"switch II"'

# --- view centred on the switches ---
orient 6VJJ_ON_3P and resi 30-76
zoom  6VJJ_ON_3P and resi 30-76, 5

# ------------------------------------------------------------
# Optional: toggle each state on its own to see the switch-I spread
#   disable 5O2S_OFF_2P
#   disable 1NVW_APO_0P
# Optional: bring in the partners to check binding-site overlap
#   (RAF1 = 6vjj chain B; DARPin K27 = 5o2s chain F; SOS = 1nvw chain S)
# ------------------------------------------------------------
