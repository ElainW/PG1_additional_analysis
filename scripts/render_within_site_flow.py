"""
Render the within-site disagreement flowchart (the mermaid diagram in README.md)
to a vector PDF. Not part of the numbered `make all` pipeline (no mermaid renderer
is available locally, so this reproduces the same diagram with matplotlib).

Output: figures/within_site_disagreement_flow.pdf
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import common as C

# colours matching the mermaid classDefs
PREP = dict(fc="#FDEBD0", ec="#E0A96D", tc="#5b4636")
COMP = dict(fc="#E6E0F5", ec="#8172B3", tc="#3D2F66")
RES = dict(fc="#D6EBD5", ec="#5FA05F", tc="#274D27")
LOOP_EC = "#8172B3"

CX, W, H = 5.0, 6.4, 1.15   # box centre-x, width, height


def box(ax, cy, title, sub, style):
    ax.add_patch(FancyBboxPatch((CX - W / 2, cy - H / 2), W, H,
                                boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=style["fc"], ec=style["ec"], lw=1.6, zorder=3))
    ax.text(CX, cy + 0.17, title, ha="center", va="center", fontsize=12.5,
            fontweight="bold", color=style["tc"], zorder=4)
    ax.text(CX, cy - 0.26, sub, ha="center", va="center", fontsize=9.5,
            style="italic", color=style["tc"], zorder=4)


def arrow(ax, y0, y1):
    ax.annotate("", xy=(CX, y1), xytext=(CX, y0),
                arrowprops=dict(arrowstyle="-|>", color="#888", lw=1.6,
                                shrinkA=0, shrinkB=0), zorder=2)


def main():
    fig, ax = plt.subplots(figsize=(7.2, 9.6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 13.2); ax.axis("off")

    yA, yB, yC, yD, yE, yF = 12.3, 10.7, 8.35, 6.85, 5.35, 3.35

    # dashed "repeat for each site" container around C, D, E
    top, bot = 9.55, 4.55
    ax.add_patch(Rectangle((CX - W / 2 - 0.45, bot), W + 0.9, top - bot,
                           fc="#F7F5FC", ec=LOOP_EC, lw=1.4, ls=(0, (6, 4)), zorder=1))
    ax.text(CX, top - 0.32, "repeat for each site", ha="center", va="center",
            fontsize=11, color=LOOP_EC, zorder=4)

    box(ax, yA, "Select proteins with >1 assay category", "disentangle phenotype vs position", PREP)
    box(ax, yB, "Select positions with ≥5 substitutions", "one site = one residue position", PREP)
    box(ax, yC, "within-site rank", "rank the ~19 subs, each of 23 models", COMP)
    box(ax, yD, "rank correlation", "Spearman, every model pair", COMP)
    box(ax, yE, "mean pairwise Spearman", "average over all pairs = consensus", COMP)
    box(ax, yF, "within-site disagreement",
        "= 1 − mean pairwise Spearman  (high = disagree)", RES)

    for y0, y1 in [(yA - H / 2, yB + H / 2), (yB - H / 2, yC + H / 2),
                   (yC - H / 2, yD + H / 2), (yD - H / 2, yE + H / 2),
                   (yE - H / 2, yF + H / 2)]:
        arrow(ax, y0, y1)

    # legend
    ly = 1.7
    for i, (lab, st) in enumerate([("inputs / prep", PREP), ("per-site computation", COMP),
                                   ("result", RES)]):
        x = 1.1 + i * 3.1
        ax.add_patch(Rectangle((x, ly), 0.32, 0.32, fc=st["fc"], ec=st["ec"], lw=1.4))
        ax.text(x + 0.45, ly + 0.16, lab, ha="left", va="center", fontsize=10, color="#444")

    fig.suptitle("Within-site model disagreement — computation", fontsize=13, y=0.965)
    out = os.path.join(C.FIGURES, "within_site_disagreement_flow.pdf")
    os.makedirs(C.FIGURES, exist_ok=True)
    plt.savefig(out, bbox_inches="tight")
    plt.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out} (+ .png)")


if __name__ == "__main__":
    main()
