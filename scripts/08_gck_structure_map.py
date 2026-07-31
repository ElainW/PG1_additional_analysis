"""
Map the HXK4_HUMAN (GCK / glucokinase) per-site model-disagreement onto known
functional regions from UniProt P35557 (data/uniprot_P35557_HXK4_HUMAN.json).

Functional annotations used:
  * Hexokinase small subdomain (67-203) / large subdomain (204-443)
  * ATP-binding residues and substrate(glucose)-binding residues
  * Known disease natural-variant positions (mostly MODY2 / hyperinsulinism)

We then ask: do the sites where models disagree most fall on functional residues?
We annotate every site, tabulate the top disagreement sites with their context,
and test whether binding-site / disease positions have different disagreement
than the rest. A 1-D sequence track figure overlays disagreement on the
functional map for both assays.

Outputs:
  results/gck_site_functional_annotation_<assay>.csv
  results/gck_disagreement_vs_function_summary.csv
  figures/fig_gck_structure_map.png
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
import common as C

ASSAYS = {
    "activity": "HXK4_HUMAN_Gersing_2022_activity",
    "abundance": "HXK4_HUMAN_Gersing_2023_abundance",
}
UNIPROT = os.path.join(C.ANALYSIS, "data", "uniprot_P35557_HXK4_HUMAN.json")


def load_features():
    d = json.load(open(UNIPROT))
    atp, substrate, disease = set(), set(), {}
    subdomain = {}
    for f in d["features"]:
        t = f["type"]
        s = f["location"]["start"]["value"]
        e = f["location"]["end"]["value"]
        if t == "Binding site":
            lig = f.get("ligand", {}).get("name", "")
            target = atp if lig == "ATP" else substrate
            for p in range(s, e + 1):
                target.add(p)
        elif t == "Region":
            desc = f.get("description", "")
            tag = "small" if "small" in desc else ("large" if "large" in desc else None)
            if tag:
                for p in range(s, e + 1):
                    subdomain[p] = tag
        elif t == "Natural variant":
            desc = f.get("description", "")
            if any(k in desc for k in ["MODY", "hyperinsulin", "PNDM", "diabetes", "HHF", "PHHI"]):
                disease[s] = desc
    return atp, substrate, disease, subdomain


def annotate(assay_key, dms_id, atp, substrate, disease, subdomain):
    res = pd.read_csv(os.path.join(C.RESULTS, f"hxk4_{dms_id}_site_disagreement.csv"))
    res["subdomain"] = res["position"].map(subdomain).fillna("other")
    res["ATP_site"] = res["position"].isin(atp)
    res["substrate_site"] = res["position"].isin(substrate)
    res["binding_site"] = res["ATP_site"] | res["substrate_site"]
    res["disease_variant"] = res["position"].isin(disease)
    res["functional"] = res["binding_site"] | res["disease_variant"]
    res.to_csv(os.path.join(C.RESULTS, f"gck_site_functional_annotation_{assay_key}.csv"), index=False)
    return res


def main():
    atp, substrate, disease, subdomain = load_features()
    print(f"ATP-binding residues: {sorted(atp)}")
    print(f"substrate(glucose)-binding residues: {sorted(substrate)}")
    print(f"disease-variant positions: {len(disease)}")

    annotated = {}
    summary_rows = []
    for key, dms_id in ASSAYS.items():
        res = annotate(key, dms_id, atp, substrate, disease, subdomain)
        annotated[key] = res
        # enrichment tests: disagreement at functional vs other
        for grp, mask in [("binding_site", res["binding_site"]),
                          ("disease_variant", res["disease_variant"]),
                          ("functional(any)", res["functional"])]:
            a = res.loc[mask, "disagreement"]
            b = res.loc[~mask, "disagreement"]
            if len(a) >= 3 and len(b) >= 3:
                U, p = mannwhitneyu(a, b, alternative="two-sided")
            else:
                p = np.nan
            summary_rows.append({
                "assay": key, "group": grp, "n_group": int(mask.sum()),
                "mean_disagree_group": round(a.mean(), 3),
                "mean_disagree_rest": round(b.mean(), 3),
                "mannwhitney_p": p,
            })

        # annotated view of the top-15 disagreed sites
        top = res.head(15).copy()
        top["site"] = top["wt_aa"] + top["position"].astype(str)
        top["context"] = top.apply(lambda r: ", ".join(filter(None, [
            "ATP-site" if r["ATP_site"] else "",
            "glucose-site" if r["substrate_site"] else "",
            "disease-variant" if r["disease_variant"] else "",
            f"{r['subdomain']} subdomain" if r["subdomain"] != "other" else "",
        ])) or "no annotated function", axis=1)
        print(f"\n==== {key} ({dms_id}): top-15 disagreed sites in functional context ====")
        print(top[["site", "n_subs", "mean_pairwise_spearman", "dms_consensus_spearman", "context"]]
              .to_string(index=False))

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(C.RESULTS, "gck_disagreement_vs_function_summary.csv"), index=False)
    print("\n=== disagreement at functional vs other residues ===")
    print(summary.to_string(index=False))

    # =================== FIGURE: 1-D sequence track ===================
    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    for ax, key in zip(axes, ASSAYS):
        res = annotated[key].sort_values("position")
        ax.plot(res["position"], res["disagreement"], color="#444", lw=0.8, zorder=2)
        ax.fill_between(res["position"], res["disagreement"], color="#bbb", alpha=0.4, zorder=1)
        # mark binding sites
        for p in sorted(atp):
            ax.axvline(p, color="#1f77b4", lw=0.8, alpha=0.6, zorder=0)
        for p in sorted(substrate):
            ax.axvline(p, color="#d62728", lw=0.8, alpha=0.6, zorder=0)
        # annotate top-12 disagreement sites
        top = res.sort_values("disagreement", ascending=False).head(12)
        ax.scatter(top["position"], top["disagreement"], color="black", s=25, zorder=3)
        for _, r in top.iterrows():
            tags = []
            if r["ATP_site"]: tags.append("ATP")
            if r["substrate_site"]: tags.append("glc")
            if r["disease_variant"]: tags.append("dz")
            lab = f"{r['wt_aa']}{int(r['position'])}" + (f"\n({'/'.join(tags)})" if tags else "")
            ax.annotate(lab, (r["position"], r["disagreement"]), fontsize=6.5,
                        ha="center", va="bottom", xytext=(0, 3), textcoords="offset points")
        ax.axvspan(67, 203, color="orange", alpha=0.06, zorder=0)
        ax.axvspan(204, 443, color="green", alpha=0.06, zorder=0)
        ax.set_ylabel(f"disagreement\n(1 - mean pairwise Spearman)\n[{key}]", fontsize=9)
        ax.set_ylim(0, max(0.85, res["disagreement"].max() + 0.08))
    axes[1].set_xlabel("GCK residue position")
    legend = [
        Patch(color="#1f77b4", label="ATP-binding residue"),
        Patch(color="#d62728", label="glucose-binding residue"),
        Patch(color="orange", alpha=0.3, label="small subdomain (67-203)"),
        Patch(color="green", alpha=0.3, label="large subdomain (204-443)"),
        plt.Line2D([0], [0], marker="o", color="black", lw=0, label="top-12 disagreed sites"),
    ]
    axes[0].legend(handles=legend, fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.28))
    axes[0].set_title("GCK (HXK4_HUMAN): per-site model disagreement vs functional map", fontsize=12, pad=28)
    plt.tight_layout()
    out = os.path.join(C.FIGURES, "fig_gck_structure_map.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
