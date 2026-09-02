"""
Secondary-structure annotation from UniProt for each DMS-benchmark protein.

UniProt records carry experimentally-derived secondary-structure features
(type = 'Helix', 'Beta strand', 'Turn'), in UniProt sequence numbering. We fetch
the UniProtKB JSON, project those features onto every residue of the canonical
sequence, and label the remainder 'loop':

    helix   <- UniProt 'Helix'
    strand  <- UniProt 'Beta strand'   (beta sheet)
    turn    <- UniProt 'Turn'
    loop    <- everything else (incl. positions with no SS annotation)

Outputs (results/structure_features/secondary_structure/):
  <UNIPROT>_<PROTEIN>_ss.tsv   per-residue: pos, aa, ss
  secondary_structure_all.csv  combined long table
  secondary_structure_summary.csv  per-protein counts per SS class
"""
import os
import sys
import json
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

OUTDIR = os.path.join(C.RESULTS, "structure_features", "secondary_structure")
os.makedirs(OUTDIR, exist_ok=True)

PROTEINS = [
    ("HXK4_HUMAN", "P35557"),
    ("RASK_HUMAN", "P01116"),
    ("CP2C9_HUMAN", "P11712"),
    ("PTEN_HUMAN", "P60484"),
    ("SRC_HUMAN", "P12931"),
    ("VKOR1_HUMAN", "Q9BQB6"),
]

SS_MAP = {"Helix": "helix", "Beta strand": "strand", "Turn": "turn"}


def fetch(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.load(r)


def main():
    allrows = []
    summ = []
    for prot, acc in PROTEINS:
        data = fetch(acc)
        seq = data["sequence"]["value"]
        n = len(seq)
        ss = ["loop"] * n                     # default
        counts_feats = {"helix": 0, "strand": 0, "turn": 0}
        for feat in data.get("features", []):
            lab = SS_MAP.get(feat.get("type"))
            if lab is None:
                continue
            counts_feats[lab] += 1
            s = feat["location"]["start"]["value"]
            e = feat["location"]["end"]["value"]
            if s is None or e is None:
                continue
            for pos in range(int(s), int(e) + 1):
                if 1 <= pos <= n:
                    ss[pos - 1] = lab

        rows = [dict(uniprot=acc, protein=prot, pos=i + 1, aa=seq[i], ss=ss[i])
                for i in range(n)]
        # per-protein tsv
        pd.DataFrame(rows)[["pos", "aa", "ss"]].to_csv(
            os.path.join(OUTDIR, f"{acc}_{prot}_ss.tsv"), sep="\t", index=False)
        allrows.extend(rows)

        from collections import Counter
        c = Counter(r["ss"] for r in rows)
        summ.append(dict(protein=prot, uniprot=acc, length=n,
                         helix=c["helix"], strand=c["strand"], turn=c["turn"],
                         loop=c["loop"],
                         n_helix_segs=counts_feats["helix"],
                         n_strand_segs=counts_feats["strand"],
                         n_turn_segs=counts_feats["turn"]))
        print(f"{prot:12s} {acc}  len={n:4d}  "
              f"helix={c['helix']:3d} strand={c['strand']:3d} "
              f"turn={c['turn']:3d} loop={c['loop']:3d}  "
              f"({counts_feats['helix']}H/{counts_feats['strand']}E/{counts_feats['turn']}T segments)")

    pd.DataFrame(allrows).to_csv(
        os.path.join(OUTDIR, "secondary_structure_all.csv"), index=False)
    pd.DataFrame(summ).to_csv(
        os.path.join(OUTDIR, "secondary_structure_summary.csv"), index=False)
    print(f"\nWrote per-protein *_ss.tsv + secondary_structure_all.csv + _summary.csv to\n  {OUTDIR}")


if __name__ == "__main__":
    main()
