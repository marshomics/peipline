#!/usr/bin/env python3
"""Assign every C71 sequence to the published four-class partition.

Wang et al. 2025 (Int J Biol Macromol, doi:10.1016/j.ijbiomac.2025.141813,
PDB 8JX4 / 8Z4F) show that the peptidase C71 catalytic centre carries a
conserved motif -- N197 C198 x D200 x x Q203 ... H233 ... D250, PeiW numbering --
but that two residues sitting 3-4 A from the His/Asp centre are *not* conserved,
and that they set the enzyme's activity. They divide the family on those two
positions:

    class I    V252 + C265     PeiW, PeiP: the active reference
    class II   V252 + V265     retains ~50% activity when engineered into PeiW-CD
    class III  T/S252 + I265   strongly reduced
    class IV   A252 + M/W265   strongly reduced; C265 -> bulky hydrophobic kills it

Every non-class-I combination they built into PeiW-CD lost catalytic activity.
So this is not a taxonomy of sequences, it is a prediction about function, made
from two alignment columns, by somebody who then went and measured it.

Three things this gives the pipeline that it could not generate for itself.

An external prior for k. Our active-site subgroups come from k-means on a
BLOSUM-embedded barcode, with k chosen by cosine silhouette. The literature says
four. `pei_class_vs_subgroup.tsv` reports the adjusted Rand index between them.
Agreement is evidence; disagreement is the more interesting result and needs to
be explained rather than smoothed over.

A non-circular partition for SDP calling, with one caveat that has to be handled
rather than ignored: V252 is two residues from the catalytic D250, so it sits
INSIDE the +/-5 triad-flank barcode. A pei_class partition is therefore partly
determined by a barcode column. sdp.py excludes both class positions from the
column set whenever this partition is used, and the exclusion is recorded.

A functional expectation for the screen: class I is the active form. A large
class III/IV population among triad-positive sequences would mean the triad
filter is retaining proteins that cannot cleave, which is a finding about the
filter, not about the biology.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import PALETTE, load_config, palette, read_fasta, savefig, set_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

CITATION = ("Wang et al. 2025, Int J Biol Macromol, "
            "doi:10.1016/j.ijbiomac.2025.141813 (PDB 8JX4, 8Z4F)")

# Only class I is the demonstrated-active combination.
ACTIVITY = {"I": "active (PeiW, PeiP)",
            "II": "~50% of PeiW-CD activity when engineered",
            "III": "strongly reduced",
            "IV": "strongly reduced (bulky hydrophobic at position 2)",
            "unassigned": "residue pair not among the four published classes"}


def classify(r1, r2, classes):
    for name, spec in classes.items():
        if r1 in spec["p1"] and r2 in spec["p2"]:
            return name
    return "unassigned"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--afa", required=True)
    ap.add_argument("--groove-json", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--ssn", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-agreement", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    pc = cfg["specificity"].get("pei_class") or {}
    fmts, dpi = tuple(cfg["plots"]["formats"]), cfg["plots"]["dpi"]
    set_style()
    os.makedirs(a.figdir, exist_ok=True)
    for p in (a.out, a.out_agreement):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    if not pc.get("enabled", False):
        pd.DataFrame(columns=["seq_id", "pei_class"]).to_csv(a.out, sep="\t", index=False)
        pd.DataFrame().to_csv(a.out_agreement, sep="\t", index=False)
        return

    gdef = json.load(open(a.groove_json))
    cc = gdef.get("class_positions") or {}
    if not cc:
        sys.exit("[pei_class] groove_definition.json has no class_positions; "
                 "run groove_map.py with specificity.pei_class.enabled")
    c1, c2 = cc["position_1"], cc["position_2"]
    print(f"[pei_class] {CITATION}", file=sys.stderr)
    print(f"[pei_class] position 1: {c1['residue']}{c1['auth']} -> match col "
          f"{c1['match_col']}", file=sys.stderr)
    print(f"[pei_class] position 2: {c2['residue']}{c2['auth']} -> match col "
          f"{c2['match_col']}", file=sys.stderr)

    # The reference structure must itself be class I, or the columns are wrong.
    ref = classify(c1["residue"], c2["residue"], pc["classes"])
    if ref != "I":
        sys.exit(f"[pei_class] the reference structure classifies as {ref}, not I. "
                 f"PeiW and PeiP are class I by construction, so the class columns "
                 f"are wrong. Check groove_map's residue-to-column mapping.")
    print("[pei_class] sanity: the reference structure is class I", file=sys.stderr)

    names, rows = zip(*read_fasta(a.afa))
    i1, i2 = c1["match_col"], c2["match_col"]
    r1 = [s[i1] if i1 < len(s) else "-" for s in rows]
    r2 = [s[i2] if i2 < len(s) else "-" for s in rows]
    cls = [classify(x, y, pc["classes"]) for x, y in zip(r1, r2)]

    wdf = pd.read_csv(a.weights, sep="\t").set_index("seq_id")
    df = pd.DataFrame({"seq_id": names, "res_position_1": r1, "res_position_2": r2,
                       "pei_class": cls})
    df["residue_pair"] = df["res_position_1"] + df["res_position_2"]
    df["expected_activity"] = df["pei_class"].map(ACTIVITY)
    df["weight"] = wdf["weight"].reindex(names).to_numpy()
    df["cluster"] = wdf["cluster"].reindex(names).to_numpy()
    df.to_csv(a.out, sep="\t", index=False)

    vc = df["pei_class"].value_counts()
    eff = df.groupby("pei_class")["weight"].sum()
    print("\n[pei_class] class distribution (raw / redundancy-weighted):", file=sys.stderr)
    for k in ["I", "II", "III", "IV", "unassigned"]:
        if k in vc.index:
            print(f"  class {k:10s} {vc[k]:>8,}  eff {eff[k]:>9.1f}   {ACTIVITY[k]}",
                  file=sys.stderr)

    frac_I = vc.get("I", 0) / len(df)
    if frac_I < 0.5:
        print(f"\n[pei_class] only {100*frac_I:.1f}% of triad-positive sequences are "
              f"class I. Classes III and IV lose activity in vitro, so a large "
              f"non-class-I population means the triad filter is retaining proteins "
              f"that probably cannot cleave. That is a statement about the filter.",
              file=sys.stderr)
    if vc.get("unassigned", 0) > 0.3 * len(df):
        print(f"\n[pei_class] {vc['unassigned']:,} sequences carry a residue pair "
              f"outside the four published classes. The classification was built on "
              f"the 30 proteins then in InterPro IPR022119; a screen of 350k genomes "
              f"is expected to exceed it. These are candidates for a fifth class, "
              f"not errors.", file=sys.stderr)

    # --- agreement with our own partitions -----------------------------------
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
    assign = pd.read_csv(a.assign, sep="\t", dtype={"seq_id": str}).set_index("seq_id")
    ssn = pd.read_csv(a.ssn, sep="\t", dtype={"seq_id": str}).set_index("seq_id")

    # dereplicate: an over-sequenced genus must not manufacture agreement
    derep = df.sort_values("seq_id").drop_duplicates("cluster")
    rows_a = []
    for label, series in (("active_site_subgroup", assign["subgroup"]),
                          ("ssn_cluster", ssn["ssn_cluster"])):
        s = derep["seq_id"].map(series)
        m = s.notna() & (derep["pei_class"] != "unassigned")
        if m.sum() < 20:
            continue
        rows_a.append({
            "partition": label,
            "n": int(m.sum()),
            "n_pei_classes": int(derep.loc[m, "pei_class"].nunique()),
            "n_partition_groups": int(s[m].nunique()),
            "adjusted_rand": adjusted_rand_score(derep.loc[m, "pei_class"], s[m]),
            "adjusted_mutual_info": adjusted_mutual_info_score(
                derep.loc[m, "pei_class"], s[m]),
        })
    agree = pd.DataFrame(rows_a)
    agree.to_csv(a.out_agreement, sep="\t", index=False)
    if len(agree):
        print("\n[pei_class] agreement with our own partitions "
              f"(on {len(derep):,} cluster representatives):", file=sys.stderr)
        print(agree.to_string(index=False), file=sys.stderr)
        ari = agree.loc[agree["partition"] == "active_site_subgroup", "adjusted_rand"]
        if len(ari) and float(ari.iloc[0]) < 0.2:
            print("\n[pei_class] the k-means active-site subgroups do NOT recover the "
                  "published classes. Either the barcode is reading something else, "
                  "or the classes are not the dominant axis of variation here. Say "
                  "which; do not average them together.", file=sys.stderr)

    # --- figure --------------------------------------------------------------
    order = [k for k in ["I", "II", "III", "IV", "unassigned"] if k in vc.index]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6),
                             gridspec_kw={"width_ratios": [1, 1.3]})
    ax = axes[0]
    ax.bar(order, [vc[k] for k in order], color=palette(len(order)), width=0.6)
    for i, k in enumerate(order):
        ax.text(i, vc[k], f"{vc[k]:,}", ha="center", va="bottom", fontsize=6)
    ax.set_yscale("log")
    ax.set_ylabel("Sequences")
    ax.set_title("A  Published C71 classes", loc="left")
    ax.tick_params(axis="x", rotation=30)

    ax = axes[1]
    if "subgroup" in assign.columns:
        ct = pd.crosstab(derep["pei_class"],
                         derep["seq_id"].map(assign["subgroup"]))
        ct = ct.reindex([k for k in order if k in ct.index])
        im = ax.imshow(ct.to_numpy(), cmap="magma", aspect="auto")
        ax.set_xticks(range(ct.shape[1]))
        ax.set_xticklabels([f"SG{c}" for c in ct.columns])
        ax.set_yticks(range(ct.shape[0]))
        ax.set_yticklabels(ct.index)
        for i in range(ct.shape[0]):
            for j in range(ct.shape[1]):
                v = ct.iat[i, j]
                ax.text(j, i, f"{v}", ha="center", va="center", fontsize=6,
                        color="w" if v < ct.to_numpy().max() * 0.6 else "k")
        fig.colorbar(im, ax=ax, fraction=0.046).set_label("cluster reps")
        r = agree.loc[agree["partition"] == "active_site_subgroup", "adjusted_rand"]
        ax.set_title(f"B  Class vs k-means subgroup"
                     + (f" (ARI = {float(r.iloc[0]):.3f})" if len(r) else ""),
                     loc="left")
    fig.tight_layout()
    savefig(fig, a.figdir, "28_pei_class", fmts, dpi)


if __name__ == "__main__":
    main()
