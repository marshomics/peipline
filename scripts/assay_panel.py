#!/usr/bin/env python3
"""Nominate the proteins to actually put in a tube.

Everything upstream is a hypothesis generator. This turns it into a purchase
order: N sequences chosen to span the specificity space, so that whichever way
the assay comes out, it discriminates.

Selection is max-min (farthest-point / Kennard-Stone) in a feature space built
from the things that could plausibly determine specificity:

  the replicated SDP residues (one-hot, weighted by how many partitions
  replicated them), the PMBR repeat count, the accessory binding module, and
  the host's literature P1 residue.

Max-min, not k-medoids: we want the corners of the space, not its modes. A
cluster centroid is the protein you can already predict; a corner is the one
that tests the model. PeiP and PeiW (or their closest relatives in the data) are
seeded first, because a panel with no characterised member cannot be calibrated.

Which assay to ask for, and why it is not the obvious one
--------------------------------------------------------
Two substrate formats exist for this family and they do not measure the same
thing.

  pNA dipeptides  Glu-gamma-Ala-pNA and relatives (Subedi et al. 2015). Cheap,
                  continuous, spectrophotometric. The p-nitroanilide leaving
                  group sits where the acyl acceptor belongs, so S1' is EMPTY.

  isopeptides     Glu-gamma-Ala-epsilon-Lys and relatives, sub1-6 (Wang et al.
                  2025). Carry a real lysine acceptor. Discontinuous.

They disagree about threonine at P1, and the disagreement is the most useful fact
in this file:

  Subedi:  Glu-gamma-Thr-pNA  -> no detectable activity, PeiW or PeiP.
           This substrate EXISTS. JPT synthesised it and both enzymes were
           assayed against it. An earlier version of this script asserted that no
           Thr-pNA had ever been made. That was wrong.
  Subedi:  M. ruminantium M1, Thr at P1 -> not lysed by either enzyme.
  Wang:    Glu-gamma-Thr-epsilon-Lys -> cleaved by PeiW, barely by PeiP.

The two reconcile if S1 and S1' are coupled: PeiW's S1 accommodates the extra
methyl of threonine only when a genuine epsilon-amino acceptor occupies S1'.
Hence the rule this script enforces:

    a P1 claim cannot be tested with a pNA substrate, unless the residue is
    Ala or Ser -- the two the pNA series is known to cleave.

  prediction_p1        which residue this protein's S1 pocket must read, AND the
                       substrate format that can test it.

  prediction_p1_prime  whether the acceptor arm matters. M. smithii PS has about
                       a quarter of its lysine replaced by ornithine, one
                       methylene shorter. NO published substrate tests this: the
                       pNA series has no acceptor at all and sub1-6 all carry Lys.

  assay_metal          both enzymes are dead after EDTA (<1% activity) and both
                       are rescued by Ca. Only PeiW is rescued by Mn, Mg, Ba or
                       Ni; for PeiP those give under 15%. Metal selectivity is
                       the sharpest measured difference between the two
                       characterised enzymes, and neither the groove nor the
                       four-class partition predicts it. So every panel member is
                       assayed against the cation series, not just Ca.

  prediction_binding   whether the protein can dock at all. Three PMB motifs are
                       required to bind an intact pseudomurein sacculus; two bind
                       only lysozyme-exposed murein fragments; one binds nothing
                       (Visweswaran et al. 2011). A soluble isopeptide never
                       engages the module, so a protein with fewer than three
                       motifs can cleave the substrate and still fail to lyse a
                       cell. Assay those on the isopeptide only: a whole-cell
                       negative would be a binding result reported as a catalytic
                       one.

  assay_ph             the module binds near its own pI -- completely at pH 9.0,
                       partially at 6.5, not at all at 4.0 -- and aggregates at
                       pH 7.0. Every published Pei lysis assay runs at 7.0-7.85.
                       PMB pIs span 3-10, so the pI is computed per protein.

This last point reopens something. PeiW cleaves the Glu-gamma-Thr-epsilon-Lys
isopeptide but cannot lyse Thr-walled M. ruminantium M1. If the failure were
catalytic the isopeptide would resist too, and it does not. So the failure may be
at the docking step, not the active site. The discriminating experiment is cheap:
assay the isolated PeiW catalytic domain against M1 sacculi alongside full-length
PeiW. Note the tension, though -- the PMB domain is thought to read NAG in the
glycan, and Thr-for-Ala is a substitution in the peptide, so a binding explanation
needs the two to be coupled somehow. Unresolved either way.

A protein's `pei_class` (Wang et al.'s four-class partition) is also carried
through, because class I is the only combination with demonstrated activity. A
panel member outside class I is a prediction that it will not cleave, which is
itself worth testing.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import PALETTE, load_config, palette, read_fasta, savefig, set_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY-"


def onehot(rows, cols, weights):
    X = np.zeros((len(rows), len(cols) * len(AA)), dtype=np.float32)
    for i, s in enumerate(rows):
        for k, c in enumerate(cols):
            j = AA.find(s[c]) if c < len(s) else -1
            if j >= 0:
                X[i, k * len(AA) + j] = weights[k]
    return X


def maxmin(X, n, seed_idx):
    """Kennard-Stone: repeatedly take the point farthest from everything chosen."""
    chosen = list(seed_idx)
    if not chosen:
        chosen = [int(np.argmax(np.linalg.norm(X - X.mean(0), axis=1)))]
    d = np.full(len(X), np.inf)
    for c in chosen:
        d = np.minimum(d, np.linalg.norm(X - X[c], axis=1))
    while len(chosen) < n:
        i = int(np.argmax(d))
        if not np.isfinite(d[i]) or d[i] == 0:
            break
        chosen.append(i)
        d = np.minimum(d, np.linalg.norm(X - X[i], axis=1))
    return chosen


# What has actually been synthesised and measured, so the panel asks for
# experiments that exist rather than experiments that sound plausible.
#
# Subedi et al. 2015, chromogenic series (S1' empty, pNA in the acceptor site):
#   EgammaA-pNA  active     both enzymes; the reference substrate
#   EA-pNA       active     the Glu-Ala bond need not be a gamma-isopeptide
#   EgammaS-pNA  active     Ser substitutes for Ala at P1
#   EgammaT-pNA  INACTIVE   Thr does not substitute -- with S1' empty
#   DbetaA-pNA   poor       Asp does not substitute for Glu at P2
#   A-pNA        poor       P2 Glu supplies most of the binding energy
#
# Wang et al. 2025, isopeptide series sub1-6 (S1' carries a real Lys):
#   both prefer Glu-gamma-Ala-epsilon-Lys;
#   PeiW ALSO cleaves Glu-gamma-Thr-epsilon-Lys, PeiP barely;
#   PeiP is weaker throughout (70-90% degradation vs near-complete for PeiW);
#   the PB (PMBR) repeats improve recognition of Glu-gamma-Thr/Ser and
#   Asp-beta-Ala, so the binding module is not a passive tether.
#
# The Thr rows are the reason `format` exists below. A pNA substrate cannot test
# a P1 residue that needs an occupied S1'.
P1_PREDICTION = {
    "Ala": ("S1 reads Ala. EgammaA-pNA is the reference substrate for both PeiW "
            "and PeiP [format: pNA]"),
    "Ser": ("S1 must accept Ser. EgammaS-pNA is cleaved by both enzymes, so the "
            "cheap continuous assay is valid here [format: pNA]"),
    "Thr": ("S1 must accept Thr. Do NOT use EgammaT-pNA: it exists, and neither "
            "enzyme cleaves it. Assay the Glu-gamma-Thr-epsilon-Lys ISOPEPTIDE, "
            "which PeiW does cleave. If this protein cuts the pNA form, that is "
            "new [format: isopeptide]"),
    "unsupported": ("wall chemistry asserted only in secondary sources; assay the "
                    "Ala / Ser / Thr isopeptide panel and let it decide "
                    "[format: isopeptide]"),
    "no_pseudomurein": "host has no pseudomurein sacculus; negative control",
    "unknown": ("no primary literature for this host; assay the sub1-6 isopeptide "
                "panel (Glu-gamma-Ala/Thr/Ser, Asp-beta-Ala) [format: isopeptide]"),
}
P1P_PREDICTION = {
    "Lys/Orn": ("host wall carries ~25% Orn at P1'. No published substrate tests "
                "this: the pNA series has no acceptor and sub1-6 all carry Lys. "
                "Synthesise the Ala-delta-Orn isopeptide"),
    "Lys": "",
    "unsupported": "",
    "no_pseudomurein": "",
    "unknown": "",
}

# Residues for which the cheap continuous pNA assay is a valid test of S1.
# Everything else needs an acceptor in S1'.
PNA_VALID_P1 = {"Ala", "Ser"}

# Subedi et al. 2015. PeiW is rescued from EDTA by all five; PeiP by Ca alone.
# The panel asks for the full series on every protein, because a second
# Ca-strict enzyme, or a Mg-tolerant one, would be the first structural handle on
# a site that no one has located.
METAL_SERIES = ["Ca", "Mn", "Mg", "Ba", "Ni"]
ASSAY_METAL = ("EDTA-treat, then titrate " + " / ".join(METAL_SERIES) +
               ". PeiW is rescued by all five, PeiP by Ca only (<15% for the "
               "rest). Report the rescue profile: it is the sharpest measured "
               "difference between the two characterised enzymes")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--afa", required=True)
    ap.add_argument("--sdp", required=True)
    ap.add_argument("--arch", required=True)
    ap.add_argument("--cellwall", required=True)
    ap.add_argument("--idmap", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    scfg = cfg["specificity"]
    fmts, dpi = tuple(cfg["plots"]["formats"]), cfg["plots"]["dpi"]
    set_style()
    os.makedirs(a.figdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    names, rows = zip(*read_fasta(a.afa))
    sdp = pd.read_csv(a.sdp, sep="\t")
    rep = sdp[sdp["replicated"]]
    if rep.empty:
        print("[panel] no replicated SDPs. Falling back to the groove columns, "
              "and saying so.", file=sys.stderr)
        cols = sorted(sdp["match_col"].unique())[:20]
        wts = np.ones(len(cols))
    else:
        cols = rep["match_col"].astype(int).tolist()
        wts = rep["n_partitions"].to_numpy(float)
    print(f"[panel] feature space: {len(cols)} SDP columns", file=sys.stderr)

    arch = pd.read_csv(a.arch, sep="\t", dtype={"seq_id": str}).set_index("seq_id")
    idmap = pd.read_csv(a.idmap, sep="\t", dtype=str).set_index("seq_id")
    cw = pd.read_csv(a.cellwall, sep="\t", dtype={"sample": str}).set_index("sample")
    assign = pd.read_csv(a.assign, sep="\t", dtype={"seq_id": str}).set_index("seq_id")
    wdf = pd.read_csv(a.weights, sep="\t").set_index("seq_id")

    # one representative per identity cluster: never nominate two near-identical
    # proteins, however far apart they are in SDP space
    keep = wdf.reset_index().sort_values("seq_id").drop_duplicates("cluster")["seq_id"]
    sel = [i for i, n in enumerate(names) if n in set(keep)]
    names_s = [names[i] for i in sel]
    rows_s = [rows[i] for i in sel]
    print(f"[panel] {len(names_s)} cluster representatives", file=sys.stderr)

    X = onehot(rows_s, cols, wts)

    meta = pd.DataFrame({"seq_id": names_s})
    meta["sample"] = meta["seq_id"].map(idmap["sample"])
    meta["protein_id"] = meta["seq_id"].map(idmap["protein_id"])
    meta["n_pmbr"] = meta["seq_id"].map(arch["n_pmbr"]).fillna(0).astype(int)
    meta["architecture_class"] = meta["seq_id"].map(arch["architecture_class"])
    for col in ("pmbr_binding_competent", "predicted_binding", "pmbr_pi",
                "pmbr_count_fragile", "assay_ph"):
        meta[col] = meta["seq_id"].map(arch[col]) if col in arch.columns else pd.NA
    meta["accessory"] = meta["seq_id"].map(arch["accessory_binding_domains"]).fillna("")
    meta["subgroup"] = meta["seq_id"].map(assign["subgroup"])
    meta["p1_residue"] = meta["sample"].map(cw["p1_residue"]).fillna("unknown")
    meta["p1_prime_residue"] = (meta["sample"].map(cw["p1_prime_residue"])
                                if "p1_prime_residue" in cw.columns else "unknown")
    meta["p1_prime_residue"] = meta["p1_prime_residue"].fillna("unknown")
    meta["p1_source"] = (meta["sample"].map(cw["p1_source"])
                         if "p1_source" in cw.columns else "none").fillna("none")
    meta["pathway_call"] = meta["sample"].map(cw["pathway_call"]).fillna("unknown")
    meta["species"] = meta["sample"].map(cw["species"])

    # categorical axes matter as much as the residues. Both sides of the
    # scissile bond are axes: P1 is what S1 reads, P1' is what S1' reads.
    cat = pd.get_dummies(meta[["n_pmbr", "architecture_class", "p1_residue",
                               "p1_prime_residue"]].astype(str)
                         ).to_numpy(dtype=np.float32)
    # scale so no single axis dominates the max-min geometry
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-9)
    Cn = cat / (np.linalg.norm(cat, axis=1, keepdims=True) + 1e-9)
    F = np.hstack([Xn, Cn])

    # seed with anything annotated as PeiP/PeiW-like, so the panel is calibratable
    seed = [i for i, s in enumerate(meta["species"].fillna(""))
            if "Methanothermobacter" in str(s)][:2]
    if seed:
        print(f"[panel] seeding with {len(seed)} Methanothermobacter protein(s), "
              f"the genus PeiW and PeiP came from", file=sys.stderr)
    else:
        print("[panel] no Methanothermobacter protein in the data; the panel has "
              "no characterised reference and the assay will need an external "
              "PeiP control", file=sys.stderr)

    pick = maxmin(F, int(scfg["panel_size"]), seed)
    out = meta.iloc[pick].copy()
    out["panel_rank"] = range(1, len(out) + 1)
    out["is_seed"] = [int(i in seed) for i in pick]
    out["sdp_residues"] = ["".join(rows_s[i][c] for c in cols) for i in pick]
    out["prediction_p1"] = out["p1_residue"].map(P1_PREDICTION).fillna(
        P1_PREDICTION["unknown"])
    out["prediction_p1_prime"] = out["p1_prime_residue"].map(P1P_PREDICTION).fillna("")
    # The substrate format is a hard consequence of the P1 call, not a comment.
    # Getting it wrong means assaying a protein against a substrate that cannot
    # answer the question, and reading the negative as biology.
    out["substrate_format"] = np.where(
        out["p1_residue"].isin(PNA_VALID_P1), "pNA (continuous)",
        "isopeptide (S1' must be occupied)")
    out["assay_metal"] = ASSAY_METAL
    # A soluble substrate never engages the binding module, so a protein with <3
    # PMB motifs can still cleave a pNA dipeptide while being unable to lyse a
    # cell. Keeping the two readouts separate is the whole point: a negative in a
    # whole-cell assay is uninterpretable without knowing whether the enzyme
    # could dock at all.
    # Three states, not two. A protein with NO PMB array is outside the rule, not
    # on the wrong side of it: PeiR (D3DZZ6) carries no PF09373 and lyses
    # M. ruminantium M1. Collapsing "no module" into "cannot dock" would have this
    # script tell you not to assay the one Pei with a rumen phenotype.
    competent = out["pmbr_binding_competent"].fillna(0).astype(int)
    applies = out["pmbr_rule_applies"].fillna(0).astype(int) if \
        "pmbr_rule_applies" in out.columns else (out["n_pmbr"] > 0).astype(int)
    out["can_lyse_cells"] = np.where(applies == 0, pd.NA, competent)
    out["prediction_binding"] = np.select(
        [applies == 0, competent == 1],
        ["no PMB array at all. The 3-motif rule has no jurisdiction: PeiR has no "
         "PMB motif and lyses M. ruminantium M1. Assay cells AND the soluble "
         "isopeptide, and find out what this one docks with",
         "docks on an intact sacculus; whole-cell lysis and soluble substrate "
         "both interpretable"],
        default="1-2 PMB motifs: the array is predicted unable to dock on an "
                "intact sacculus (Visweswaran 2011). Assay the soluble "
                "isopeptide. A whole-cell negative would test the module, not "
                "the active site")
    keep = ["panel_rank", "seq_id", "sample", "protein_id", "species",
            "subgroup", "n_pmbr", "pmbr_binding_competent", "pmbr_rule_applies",
            "pmbr_count_fragile", "pmbr_pi", "architecture_class", "accessory",
            "p1_residue", "p1_prime_residue", "p1_source", "pathway_call",
            "sdp_residues", "is_seed", "substrate_format", "prediction_p1",
            "prediction_p1_prime", "prediction_binding", "can_lyse_cells",
            "assay_ph", "assay_metal"]
    out = out[[c for c in keep if c in out.columns]]

    n_nodock = int(((applies == 1) & (competent == 0)).sum())
    n_nojur = int((applies == 0).sum())
    if n_nodock:
        print(f"[panel] {n_nodock} of {len(out)} picks carry 1-2 PMB motifs. Their "
              f"arrays are predicted unable to bind an intact sacculus; assay them "
              f"on the soluble isopeptide, or a whole-cell negative will be a "
              f"binding result reported as a catalytic one.", file=sys.stderr)
    if n_nojur:
        print(f"[panel] {n_nojur} picks carry no PMB motif at all. The 3-motif rule "
              f"says nothing about them. PeiR is in this class and lyses cells, so "
              f"assay them BOTH ways: a cell-lysis positive from a protein with no "
              f"PMB array is the most interesting result this panel can produce.",
              file=sys.stderr)
    n_frag = int(pd.to_numeric(out["pmbr_count_fragile"], errors="coerce")
                 .fillna(0).sum())
    if n_frag:
        print(f"[panel] {n_frag} picks have a repeat count that straddles the "
              f"3-motif cliff depending on the domain E-value. Count their motifs "
              f"by hand before ordering anything.", file=sys.stderr)

    n_iso = int((out["substrate_format"] != "pNA (continuous)").sum())
    if n_iso:
        print(f"[panel] {n_iso} of {len(out)} picks cannot be assayed with the "
              f"chromogenic pNA series: their predicted P1 is neither Ala nor "
              f"Ser, and Glu-gamma-Thr-pNA is not cleaved by either "
              f"characterised enzyme even though the Thr isopeptide is. Order "
              f"the isopeptides for these.", file=sys.stderr)

    n_unsup = int((out["p1_source"] == "disputed").sum())
    n_het = int((out["p1_source"] == "genus_heterogeneous").sum())
    if n_unsup or n_het:
        print(f"[panel] {n_unsup} panel members have a wall chemistry claimed only "
              f"in secondary sources; {n_het} sit in a genus whose characterised "
              f"species disagree. Both are testable, which is why they are in the "
              f"panel, but neither has a prior.", file=sys.stderr)
    out.to_csv(a.out, sep="\t", index=False)
    print("\n[panel] nominated:", file=sys.stderr)
    print(out.to_string(index=False), file=sys.stderr)

    # --- figure: where the panel sits ---------------------------------------
    from sklearn.decomposition import PCA
    if F.shape[1] > 2 and len(F) > 3:
        Z = PCA(n_components=2, random_state=1).fit_transform(F)
        fig, ax = plt.subplots(figsize=(4.2, 3.4))
        labs = meta["p1_residue"].fillna("unknown")
        uniq = sorted(labs.unique())
        cmap = dict(zip(uniq, palette(len(uniq))))
        for u in uniq:
            m = (labs == u).to_numpy()
            ax.scatter(Z[m, 0], Z[m, 1], s=6, color=cmap[u], linewidths=0,
                       alpha=0.6, label=f"{u} (n={int(m.sum())})", rasterized=True)
        ax.scatter(Z[pick, 0], Z[pick, 1], s=60, facecolors="none",
                   edgecolors="k", linewidths=1.0, label="panel", zorder=5)
        for r, i in enumerate(pick, 1):
            ax.annotate(str(r), (Z[i, 0], Z[i, 1]), fontsize=6, ha="center",
                        va="center", zorder=6)
        ax.set_xlabel("PC1 of SDP + architecture space")
        ax.set_ylabel("PC2")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=6)
        ax.set_title(f"Assay panel: {len(pick)} maximally separated proteins",
                     loc="left")
        fig.tight_layout()
        savefig(fig, a.figdir, "27_assay_panel", fmts, dpi)


if __name__ == "__main__":
    main()
