# Structural / functional annotations per protein

UniProt numbering (matches `up_pos` in `generalization_functional_disagreement_sites.csv`). Feature source: UniProt `raw_json/<UniProt_ID>.json`, parsed exactly as in `23_structure_maps_all.py`. Disease-variant source noted per protein.

## S (SPIKE_SARS2)
- **binding**: none annotated
- **subdomain** — BetaCoV S1-NTD: 14-303
- **subdomain** — BetaCoV S1-CTD: 334-527
- **subdomain** — Putative superantigen; may bind T-cell receptor alpha/TRAC: 280-301
- **subdomain** — Receptor-binding domain (RBD): 319-541
- **subdomain** — Receptor-binding motif; binding to human ACE2: 437-508
- **subdomain** — Fusion peptide 1: 816-837
- **subdomain** — Fusion peptide 2: 835-855
- **subdomain** — Heptad repeat 1: 920-970
- **subdomain** — Heptad repeat 2: 1163-1202
- **disordered**: 1-26; 67-80; 142-164; 173-185; 246-262; 319-335; 455-490; 519-541; 621-640; 673-687; 828-853
- **disease variants** (0, UniProt natural-variant): none

## KCNE1 (KCNE1_HUMAN)
- **binding**: none annotated
- **subdomain** — interaction with KCNQ1 C-terminus: 109-129
- **disordered**: none annotated
- **disease variants** (6, ClinVar P/LP): 1, 7, 52, 70, 73, 87

## SRC (SRC_HUMAN)
- **ATP-binding** (10): 276-284, 298
- **active site** (1): 389
- **subdomain** — SH3: 84-145
- **subdomain** — SH2: 151-248
- **subdomain** — Protein kinase: 270-523
- **disordered**: 1-53
- **disease variants** (2, ClinVar P/LP): 527, 531

## PTEN (PTEN_HUMAN)
- **binding**: none annotated
- **active site** (1): 124
- **subdomain** — Phosphatase tensin-type: 14-185
- **subdomain** — C2 tensin-type: 190-350
- **disordered**: 352-403
- **disease variants** (72, ClinVar P/LP): 1, 6, 15, 20, 22-24, 27, 31, 34-36, 38-39, 44, 47, 50, 53, 57, 61, 65-66, 70-71, 92-93, 95-96, 105, 107, 111-112, 119, 123-125, 127, 130-132, 134, 136, 141, 155, 157, 159, 165, 167, 170-171, 173-174, 177, 181, 188, 200, 204, 211, 241, 246-247, 251-252, 256, 267, 277, 298, 316, 320, 322, 335, 345

## GCK (HXK4_HUMAN)
- **ATP-binding** (19): 78-83, 228, 295-296, 332-336, 411-415
- **substrate-binding** (9): 151-152, 168-169, 204-205, 231, 256, 290
- **subdomain** — Hexokinase small subdomain: 67-203
- **subdomain** — Hexokinase large subdomain: 204-443
- **disordered**: none annotated
- **disease variants** (129, ClinVar P/LP + UniProt GoF): 19, 30, 33, 36-38, 40, 43, 49, 53-54, 56-57, 60, 62, 65, 72, 77-78, 91, 99, 108, 117, 122, 124, 130, 137, 146, 150-157, 160, 168, 174-175, 177-178, 182, 187-189, 191, 197, 201-203, 208-210, 214, 220-221, 224, 226-228, 230, 235, 246, 250, 255, 258-259, 261, 263-264, 268, 271, 274-275, 277, 294-295, 300, 304, 314, 317-318, 326, 337, 339-340, 344, 348, 352, 360, 363, 367, 371, 373-375, 377-389, 392, 394, 397, 410, 412, 414, 416, 419, 423-424, 436-437, 441-442, 449-450, 454-456

## Kcnj2 (KCNJ2_MOUSE)
- **binding**: none annotated
- **subdomain** — Polyphosphoinositide (PIP2)-binding: 181-208
- **disordered**: 383-428
- **disease variants** (0, UniProt natural-variant): none

## HLA-A (Q53Z42_HUMAN)
- **Zn(2+)-binding** (2): 221-222
- **subdomain** — Ig-like: 209-295
- **disordered**: 338-365
- **disease variants** (0, UniProt natural-variant): none

## VKORC1 (VKOR1_HUMAN)
- **(S)-warfarin-binding** (2): 80, 139
- **phylloquinone-binding** (2): 135, 139
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (5, ClinVar P/LP): 29, 45, 58, 77, 98

## RPL40A (RL40A_YEAST)
- **Zn(2+)-binding** (4): 96, 99, 110, 115
- **subdomain** — Ubiquitin-like: 1-76
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none

## SLC22A1 (S22A1_HUMAN)
- **binding**: none annotated
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none

## KRAS (RASK_HUMAN)
- **GTP-binding** (15): 13-18, 28-30, 116-117, 119, 145-147
- **Mg(2+)-binding** (2): 17, 57
- **subdomain** — Hypervariable region: 166-185
- **disordered**: none annotated
- **disease variants** (9, ClinVar P/LP + UniProt GoF): 5, 12, 14, 22, 34, 58, 60, 65, 156

## CYP2C9 (CP2C9_HUMAN)
- **heme-binding** (1): 435
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (1, ClinVar P/LP): 208

## ccdB (CCDB_ECOLI)
- **binding**: none annotated
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none

## DAO1 (OXDA_RHOTO)
- **(R)-lactate-binding** (4): 223, 238, 285, 335
- **FAD-binding** (15): 12, 15, 35-36, 47-48, 52, 54, 162, 285, 334-335, 337-339
- **anthranilate-binding** (5): 58, 223, 238, 285, 335
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none

## spg (SPG1_STRSG)
- **binding**: none annotated
- **subdomain** — 2 X 37 AA repeats: 104-215
- **subdomain** — 2 X 55 AA repeats: 228-352
- **subdomain** — 5 X 5 AA repeats of [DE]-D-A-K-K: 386-410
- **disordered**: 358-422
- **disease variants** (0, UniProt natural-variant): none

## CARD11 (CAR11_HUMAN)
- **binding**: none annotated
- **subdomain** — CARD: 18-110
- **subdomain** — PDZ: 667-755
- **subdomain** — Guanylate kinase-like: 973-1140
- **subdomain** — Linker: 111-128
- **subdomain** — Inhibitory domain (ID): 450-666
- **disordered**: 460-626
- **disease variants** (13, ClinVar P/LP + UniProt GoF): 30, 40, 43, 47-49, 72, 80, 123, 126, 134, 775, 833

## TP53 (P53_HUMAN)
- **Zn(2+)-binding** (4): 176, 179, 238, 242
- **subdomain** — Interaction with CCAR2: 1-320
- **subdomain** — Interaction with HRMT1L2: 1-83
- **subdomain** — Transcription activation (acidic): 1-44
- **subdomain** — Interaction with WWOX: 66-110
- **subdomain** — Interaction with HIPK1: 100-370
- **subdomain** — Required for interaction with ZNF385A: 100-300
- **subdomain** — Required for interaction with FBXO42: 113-236
- **subdomain** — Interaction with AXIN1: 116-292
- **subdomain** — Interaction with E4F1: 256-294
- **subdomain** — Interaction with CARM1: 300-393
- **subdomain** — Interaction with HIPK2: 319-360
- **subdomain** — Oligomerization: 325-356
- **subdomain** — Basic (repression of DNA-binding): 368-387
- **subdomain** — Interaction with MORN3: 374-393
- **disordered**: 50-96; 282-325; 351-393
- **disease variants** (37, ClinVar P/LP): 106, 110-111, 127, 132, 135, 138, 143, 152, 158, 163-164, 173, 175, 178-181, 190, 193-194, 197, 214, 220, 237, 245, 248, 254, 272-273, 277, 280, 283, 337, 342, 344, 347

## NPC1 (NPC1_HUMAN)
- **cholesterol-binding** (2): 41, 79
- **subdomain** — SSD: 620-785
- **subdomain** — Important for cholesterol binding and cholesterol transfer from NPC1 to liposomes: 175-205
- **disordered**: none annotated
- **disease variants** (82, ClinVar P/LP): 1, 63, 90, 92, 113, 169, 177, 248, 372, 375, 378, 404, 434, 471, 474, 501, 509, 512, 518, 543, 558, 612, 615, 634, 640, 652, 664, 667, 691, 703, 712, 724, 734, 764, 775, 789, 825, 842, 865, 874, 887, 910, 926-927, 934, 940, 944-945, 948, 950, 954, 958, 968, 978, 992, 1001, 1007, 1015, 1019, 1034-1036, 1054, 1061, 1087-1089, 1141-1142, 1145, 1151, 1156, 1167-1168, 1186, 1189, 1197, 1205, 1212-1213, 1240, 1244

## HSP82 (HSP82_YEAST)
- **ATP-binding** (16): 33, 37, 79, 84, 92, 98-100, 119-124, 171, 380
- **subdomain** — 5 X 5 AA repeats of [DE]-[DE]-[DE]-K-K; highly charged region: 221-263
- **disordered**: 219-262; 679-709
- **disease variants** (0, UniProt natural-variant): none

## Influenza hemagglutinin (A0A2Z5U3Z0_9INFA)
- **binding**: none annotated
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none

## cI (RPC1_LAMBD)
- **binding**: none annotated
- **subdomain** — HTH cro/C1-type: 19-77
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none

## Addiction module antidote protein, CopG/Arc/MetJ family (F7YBW8_MESOW)
- **binding**: none annotated
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none

## bla (BLAT_ECOLX)
- **substrate-binding** (3): 232-234
- **active site** (4): 68, 71, 128, 164
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none

## folA (DYR_ECOLI)
- **NADP(+)-binding** (21): 7, 13-19, 45-46, 63-64, 76, 95-102
- **substrate-binding** (5): 5, 27, 52, 57, 113
- **subdomains**: none annotated
- **disordered**: none annotated
- **disease variants** (0, UniProt natural-variant): none
