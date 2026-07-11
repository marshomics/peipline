#!/usr/bin/env python3
"""What PeiW and PeiP have actually been measured to cut, and what they lyse.

Single source: Schofield LR, Beattie AK, Tootill CM, Dey D, Ronimus RS (2015)
"Biochemical Characterisation of Phage Pseudomurein Endoisopeptidases PeiW and
PeiP Using Synthetic Peptides", Archaea 2015:828693. doi:10.1155/2015/828693.
(Author list verified against Crossref for this DOI; earlier drafts of this
pipeline miscited it as "Subedi et al." -- that attribution was wrong.)

Two tables, and they are the only experimental ground truth this pipeline has.

  SUBSTRATES   a chromogenic p-nitroanilide series that walks P1 and P2 one
               residue at a time.
  LYSIS_PANEL  eleven methanogen strains scored for susceptibility to purified
               recombinant PeiW and PeiP in an agarose plate lysate assay.

Why this module exists separately from cellwall_reference.py
------------------------------------------------------------
Kandler & Koenig say what the wall is made of. Schofield et al. say what the enzyme
cuts. Those are different claims from different experiments, and conflating them
is how a chemistry table quietly becomes a prediction. The lysis panel is used
here to *test* predictions derived from the chemistry, never to generate them.

The Thr problem, stated plainly
-------------------------------
Two papers disagree about threonine at P1, and the disagreement is informative
rather than embarrassing.

  Schofield et al. 2015   Glu-gamma-Thr-pNA: no detectable activity, either enzyme.
                       M. ruminantium M1, whose wall carries Thr at P1: not
                       lysed by either enzyme.
  Wang et al. 2025     the Glu-gamma-Thr-epsilon-Lys ISOPEPTIDE is cleaved by
                       PeiW (and barely by PeiP).

The substrates differ on the far side of the scissile bond. A pNA leaving group
is not a lysine: it is a chromophore sitting where the acyl acceptor belongs, so
the pNA series leaves S1' empty. The two results reconcile if PeiW's S1 pocket
tolerates the extra methyl of threonine only when S1' is occupied by a genuine
epsilon-amino acceptor -- that is, if S1 and S1' are energetically coupled.

That is a real, falsifiable claim, and it has a consequence for this pipeline:

    P1 preference cannot be read off the pNA series.

Any panel member nominated on the strength of its predicted P1 residue must be
assayed with an isopeptide, not a p-nitroanilide. `assay_panel.py` enforces this.

What the whole-cell result adds
-------------------------------
Cleaving a soluble isopeptide is not the same as lysing a sacculus. PeiW cleaves
the Thr isopeptide and still cannot lyse the Thr-walled organism. So either the
rate is too low to matter in vivo, or the four PMBR repeats fail to dock on a
Thr wall, or the wall is protected some other way. Wang et al. found the PB
repeats improve recognition of Glu-gamma-Thr, which makes the second explanation
the interesting one and leaves the question open.

The observation that closes the loop
-------------------------------------
Methanobrevibacter ruminantium M1 is the Thr-type organism. It resists PeiW and
PeiP. Its own prophage, Phi-mru, encodes an endoisopeptidase, PeiR, that lyses
it -- and PeiR shows little homology to PeiW or PeiP (Schofield et al. 2015, citing
Attwood et al.). A divergent enzyme for a divergent wall, in the one host where
we can check.

That is the prespecified hypothesis of the whole screen: **Pei sequence tracks
host wall chemistry**. The SDP, four-class and phylogenetic-regression analyses
downstream are tests of it, and they were designed before this panel was
consulted.
"""
from __future__ import annotations

CITATION = ("Schofield et al. 2015, Archaea 2015:828693 "
            "(doi:10.1155/2015/828693)")

# ---------------------------------------------------------------------------
# Synthetic p-nitroanilide substrates, agarose plate assay, 39 C, overnight,
# anaerobic. Scored by the presence and size of a yellow zone; the zone edge is
# diffuse, so the paper reports an ordering, not a rate. Reproduced as an
# ordering.
#
# `p2` / `p1` are the residues either side of the Glu-Ala bond; `linkage` is how
# they are joined; `acceptor` is what sits in S1'. Every one of these has pNA
# there, which is the point.
# ---------------------------------------------------------------------------
SUBSTRATES = {
    "EgammaA-pNA": {
        "formula": "Glu-gamma-Ala-pNA", "p2": "Glu", "p1": "Ala",
        "linkage": "gamma-isopeptide", "acceptor": "pNA",
        "peiW": "active", "peiP": "active", "rank": 1,
        "note": "best substrate for both; used for all kinetics",
    },
    "EA-pNA": {
        "formula": "H-Glu-Ala-pNA", "p2": "Glu", "p1": "Ala",
        "linkage": "alpha-peptide", "acceptor": "pNA",
        "peiW": "active", "peiP": "active", "rank": 2,
        "note": "the Glu-Ala bond need NOT be a gamma-isopeptide; only the "
                "Ala-Lys bond does",
    },
    "EgammaS-pNA": {
        "formula": "Glu-gamma-Ser-pNA", "p2": "Glu", "p1": "Ser",
        "linkage": "gamma-isopeptide", "acceptor": "pNA",
        "peiW": "active", "peiP": "active", "rank": 3,
        "note": "serine substitutes for alanine at P1, though neither enzyme's "
                "native host has Ser there",
    },
    "EgammaT-pNA": {
        "formula": "Glu-gamma-Thr-pNA", "p2": "Glu", "p1": "Thr",
        "linkage": "gamma-isopeptide", "acceptor": "pNA",
        "peiW": "inactive", "peiP": "inactive", "rank": None,
        "note": "NO detectable activity. This substrate exists and was tested "
                "(JPT Peptide Technologies). Cf. Wang et al. 2025, where the "
                "Thr ISOPEPTIDE is cleaved by PeiW: the difference is S1'",
    },
    "DbetaA-pNA": {
        "formula": "Asp-beta-Ala-pNA", "p2": "Asp", "p1": "Ala",
        "linkage": "beta-isopeptide", "acceptor": "pNA",
        "peiW": "poor", "peiP": "poor", "rank": 5,
        "note": "aspartate cannot substitute for glutamate at P2",
    },
    "A-pNA": {
        "formula": "Ala-pNA", "p2": None, "p1": "Ala",
        "linkage": None, "acceptor": "pNA",
        "peiW": "poor", "peiP": "poor", "rank": 6,
        "note": "too small to bind; the P2 glutamate contributes most of the "
                "binding energy",
    },
}

# Residues demonstrated to be accepted / rejected at P1 by the pNA series.
# Rejection here means "rejected when S1' is empty" -- see the module docstring.
P1_ACCEPTED_PNA = {"Ala", "Ser"}
P1_REJECTED_PNA = {"Thr"}
P2_ACCEPTED = {"Glu"}
P2_REJECTED = {"Asp"}

# ---------------------------------------------------------------------------
# Divalent metal. EDTA-treated enzyme retains <1% activity; a divalent cation
# restores it fully. The two enzymes differ, sharply, and nothing in the groove
# or in the four-class partition predicts the difference.
# ---------------------------------------------------------------------------
METAL = {
    "requirement": "absolute; <1% activity after 0.5 mM EDTA, assayed without metal",
    "peiW_order": ["Ca", "Mn", "Mg", "Ba", "Ni"],
    "peiW_note": "promiscuous: all five restore activity, in this order",
    "peiP_order": ["Ca"],
    "peiP_note": "strict: Mn, Mg, Ba and Ni each give <15% of the Ca activity",
    "thermostability": ("Ca also stabilises. PeiW at 80 C: 100% activity after "
                        "60 min with Ca, 50% lost after 5 min without. PeiP at "
                        "70 C: 50% after 30 min with Ca, after 8 min without."),
    "not_tested": ["Zn", "Fe"],   # precipitated in the assay
}

# Kinetics on EgammaA-pNA, 60 C, pH 7.85. Km in the mM range for both: the small
# substrate never engages the PMBR module, so these are catalytic-domain numbers.
KINETICS = {
    "PeiW": {"km_mM": 6.25, "km_sd": 0.44, "kcat_s": 5.73, "kcat_sd": 0.22},
    "PeiP": {"km_mM": 4.14, "km_sd": 0.14, "kcat_s": 1.18, "kcat_sd": 0.03},
}

# Inhibition, 1 mM, plate assay. Reported because it is why the family was called
# "transglutaminase-like" for a decade before there was a structure. Wang et al.
# 2025 place the fold in the papain superfamily; both statements are compatible,
# because the Cys-His-Asp triad is shared by clan CA and by the transglutaminases,
# and dansylcadaverine is an active-site-directed amine, not a fold probe.
INHIBITORS = {
    "N-ethylmaleimide": {"class": "cysteine protease", "peiW": 0.80, "peiP": 0.90},
    "dansylcadaverine": {"class": "transglutaminase", "peiW": 0.60, "peiP": 0.50},
    "cystamine": {"class": "cysteine protease", "peiW": None, "peiP": None,
                  "note": "100% inhibition only when the reductant is omitted"},
    "E64": {"class": "cysteine protease", "peiW": 0.10, "peiP": 0.30},
    "PMSF": {"class": "serine protease", "peiW": 0.10, "peiP": 0.20},
}

# ---------------------------------------------------------------------------
# Plate lysate panel. Purified recombinant PeiW / PeiP, 0.2 mg, on agarose
# containing washed methanogen cells; scored by zone of clearing after 16 h.
#
# HONESTY NOTE. The machine-readable full text of this paper renders italicised
# binomials as empty strings, so several species names in the lysis paragraph
# were lost. Every row below whose `species_confidence` is not "explicit" was
# reconstructed from the DSM number given in Materials and Methods, and the
# reconstruction has NOT been checked against the printed table. Rows marked
# `verified=False` are excluded from the automated check and are listed for a
# human to confirm against Table 4 of the paper.
#
# The rows that matter for the hypothesis are all `explicit`.
# ---------------------------------------------------------------------------
LYSIS_PANEL = [
    {
        "label": "M. thermautotrophicus dH", "dsm": "DSM 1053",
        "species": "s__Methanothermobacter thermautotrophicus",
        "species_confidence": "explicit", "verified": True,
        "wall": "pseudomurein", "peiW": "lysed", "peiP": "lysed",
        "note": "strain deltaH, named in the lysis sentence",
    },
    {
        "label": "Methanobrevibacter ruminantium M1", "dsm": "DSM 1093",
        "species": "s__Methanobrevibacter ruminantium",
        "species_confidence": "explicit", "verified": True,
        "wall": "pseudomurein", "peiW": "not_lysed", "peiP": "not_lysed",
        "note": "THE key row. Thr at P1 (Kandler & Koenig 1978, same strain M1). "
                "Resistant to both enzymes. Its own prophage Phi-mru encodes "
                "PeiR, which does lyse it and is not homologous to PeiW/PeiP",
    },
    {
        "label": "Methanobrevibacter sp. SM9", "dsm": None,
        "species": None,
        "species_confidence": "explicit", "verified": True,
        "wall": "pseudomurein", "peiW": "lysed", "peiP": "not_lysed",
        "note": "the only differential row: PeiW lyses it, PeiP does not. "
                "No wall chemistry published. A free hypothesis",
    },
    {
        "label": "Methanosphaera stadtmanae", "dsm": "DSM 3091",
        "species": "s__Methanosphaera stadtmanae",
        "species_confidence": "inferred_from_dsm", "verified": False,
        "wall": "pseudomurein", "peiW": "lysed", "peiP": "lysed",
        "note": "the paper attributes the lysis of a Ser-walled organism to the "
                "EgammaS-pNA activity. DSM 3091 is M. stadtmanae. Confirm "
                "against Table 4 before relying on it",
    },
    {
        "label": "strain 31A", "dsm": None, "species": None,
        "species_confidence": "explicit_strain_only", "verified": False,
        "wall": "pseudomurein", "peiW": "lysed", "peiP": "lysed", "note": "",
    },
    {
        "label": "strain BRM9", "dsm": None, "species": None,
        "species_confidence": "explicit_strain_only", "verified": False,
        "wall": "pseudomurein", "peiW": "lysed", "peiP": "lysed", "note": "",
    },
    {
        "label": "no-pseudomurein control (glycoprotein)", "dsm": "DSM 864",
        "species": None, "species_confidence": "unresolved", "verified": True,
        "wall": "glycoprotein", "peiW": "not_lysed", "peiP": "not_lysed",
        "note": "explicit negative control; the paper states it has no pseudomurein",
    },
    {
        "label": "Methanosarcina sp. CM1", "dsm": None, "species": None,
        "species_confidence": "explicit_strain_only", "verified": True,
        "wall": "methanochondroitin", "peiW": "not_lysed", "peiP": "not_lysed",
        "note": "explicit negative control; methanochondroitin, no pseudomurein",
    },
]


# ---------------------------------------------------------------------------
def predict_lysis(p1_residue, wall):
    """Predict whole-cell susceptibility to PeiW/PeiP from host wall chemistry.

    The rule is read straight off the pNA series and off M1's resistance. It is
    deliberately the *chromogenic* rule, not the isopeptide rule: what is being
    predicted is lysis of a sacculus, and the sacculus data say Thr resists.

    Returns (call, reason). `call` is one of lysed / not_lysed / unknown.
    """
    if wall in ("glycoprotein", "methanochondroitin", "no_pseudomurein",
                "protein_sheath"):
        return "not_lysed", "no pseudomurein sacculus to cut"
    if p1_residue in P1_ACCEPTED_PNA:
        return "lysed", f"{p1_residue} at P1 is cleaved (Glu-gamma-{p1_residue}-pNA)"
    if p1_residue in P1_REJECTED_PNA:
        return "not_lysed", (f"{p1_residue} at P1 is not cleaved in the pNA series, "
                             f"and the Thr-walled M1 resists both enzymes")
    return "unknown", "no P1 assignment for this host"


def check_against_reference(verified_only=True):
    """Score `predict_lysis` on the panel. Returns a DataFrame.

    This is the pipeline's only external validation with a measured phenotype.
    It has n=6 verified rows and it is still worth more than any silhouette
    score, because it is the only number here that could have come out wrong.
    """
    import pandas as pd
    from cellwall_reference import chemistry_for

    rows = []
    for e in LYSIS_PANEL:
        if verified_only and not e["verified"]:
            continue
        if e["species"]:
            p1, _p1p, src, _n = chemistry_for(e["species"])
        else:
            p1, src = "unknown", "no_species"
        # a host whose wall is not pseudomurein overrides any P1 lookup
        wall = e["wall"] if e["wall"] != "pseudomurein" else "pseudomurein"
        pred, why = predict_lysis(p1, wall)

        for enzyme in ("peiW", "peiP"):
            obs = e[enzyme]
            if pred == "unknown":
                agree = pd.NA
            else:
                agree = (pred == obs)
            rows.append({
                "strain": e["label"], "enzyme": "PeiW" if enzyme == "peiW" else "PeiP",
                "wall": e["wall"], "p1_residue": p1, "p1_source": src,
                "observed": obs, "predicted": pred, "agrees": agree,
                # A prediction for a host with no pseudomurein is not a test of
                # anything: no sacculus, no lysis, and no P1 was consulted. Only
                # the pseudomurein rows can falsify the rule.
                "nontrivial": e["wall"] == "pseudomurein",
                "reason": why, "note": e["note"],
            })
    return pd.DataFrame(rows)


def summarise(df):
    """(n_testable, n_agree, n_nontrivial, n_nontrivial_agree, disagreements)."""
    testable = df[df["agrees"].notna()]
    ok = testable["agrees"].astype(bool)
    nt = testable[testable["nontrivial"]]
    nt_ok = nt["agrees"].astype(bool)
    return (len(testable), int(ok.sum()), len(nt), int(nt_ok.sum()),
            testable[~ok])


if __name__ == "__main__":
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import pandas as pd

    pd.set_option("display.width", 200)
    print(f"Reference: {CITATION}\n")

    print("Substrate series (agarose plate, ordering only):")
    for k, v in SUBSTRATES.items():
        print(f"  {k:14s} P2={str(v['p2']):5s} P1={v['p1']:4s} "
              f"PeiW={v['peiW']:9s} PeiP={v['peiP']:9s}  {v['note']}")
    print(f"\n  P1 accepted (pNA): {sorted(P1_ACCEPTED_PNA)}")
    print(f"  P1 rejected (pNA): {sorted(P1_REJECTED_PNA)}  <- but see Wang 2025")
    print(f"  P2 accepted      : {sorted(P2_ACCEPTED)}   rejected: {sorted(P2_REJECTED)}")

    print(f"\nMetal: PeiW rescued by {'>'.join(METAL['peiW_order'])}; "
          f"PeiP by {'/'.join(METAL['peiP_order'])} only.")

    print("\nLysis panel check (verified rows only):")
    df = check_against_reference()
    print(df[["strain", "enzyme", "wall", "p1_residue", "observed",
              "predicted", "agrees", "nontrivial"]].to_string(index=False))
    n, k, nt, nt_k, bad = summarise(df)
    print(f"\n{k}/{n} testable predictions agree; {nt_k}/{nt} of the "
          f"falsifiable (pseudomurein-walled) ones.")
    if len(bad):
        print("\nDisagreements (each one is a result, not a bug):")
        print(bad[["strain", "enzyme", "observed", "predicted", "reason"]]
              .to_string(index=False))

    na = df[df["agrees"].isna()]
    if len(na):
        print("\nNo prediction attempted (this is the honest outcome, not a gap "
              "to be filled by guessing):")
        print(na[["strain", "enzyme", "observed", "p1_source"]].to_string(index=False))

    unver = [e["label"] for e in LYSIS_PANEL if not e["verified"]]
    if unver:
        print(f"\nExcluded pending confirmation against Table 4: {unver}")
