#!/usr/bin/env python3
"""What the pseudomurein-binding module actually does, per measured construct.

Single source: Visweswaran GRR, Dijkstra BW, Kok J (2011) "A Minimum of Three
Motifs Is Essential for Optimal Binding of Pseudomurein Cell Wall-Binding Domain
of Methanothermobacter thermautotrophicus", PLoS ONE 6(6):e21582.
doi:10.1371/journal.pone.0021582

They fused one, two or three PMB motifs (PF09373) from the S-layer protein
MTH719 to GFP and asked what each construct sticks to.

The abstract and the Results disagree, and the code follows the Results
---------------------------------------------------------------------
The abstract says "at least two of the three motifs ... are necessary for
binding." The Results section, the section heading, the title and the concluding
paragraph all say something stronger and more specific:

    3 motifs -> binds pseudomurein of methanogen cells   AND bacterial spheroplasts
    2 motifs -> binds bacterial spheroplasts ONLY; no binding to pseudomurein
    1 motif  -> binds nothing
    any      -> no binding to intact bacterial cells

"At least two ... necessary for binding" is true of spheroplast binding and false
of pseudomurein binding. An earlier version of `domain_arch.py` recorded the
threshold as needed for *optimal* binding, which is the abstract's wording. It is
not what was measured. Three motifs are required for the sacculus to be bound at
all.

Why this matters more than it looks
-----------------------------------
It converts a covariate into a phenotype. A C71 catalytic domain carried on fewer
than three PMB motifs is predicted to be unable to engage an intact pseudomurein
sacculus, whatever its active site looks like. The pipeline can therefore say
something falsifiable about every hit it finds, rather than merely counting
repeats.

It also makes the repeat count load-bearing. A PMB motif is 30-35 residues. At a
strict domain E-value a weak fourth or third repeat is silently dropped, the
count falls from 3 to 2, and the functional call inverts. `domain_arch.py` scans
once permissively, filters at two thresholds, and flags every protein whose class
depends on which threshold you picked.

The module reads the sugar, not the peptide
-------------------------------------------
Three-motif and two-motif constructs both bind lysozyme-treated *Lactococcus
lactis* and *Escherichia coli* spheroplasts, and the signal survives 150 mM NaCl.
Neither binds intact bacterial cells. Lysozyme exposes NAG, and NAG is the only
sugar common to murein and pseudomurein, so the authors infer the PMB domain
recognises NAG, as LysM does.

Two consequences.

  1. PMBR presence is not evidence that the substrate is pseudomurein. It is
     evidence of an NAG-binding module. A bacterial C71+PMBR protein is not
     automatically a binning artefact: it could bind exposed murein.

  2. The binding module reads the glycan backbone; the catalytic groove reads the
     peptide cross-link (the Ala-epsilon-Lys isopeptide, and P1). They are a
     priori independent axes. That is exactly why `pmbr_architecture` is a valid
     SDP partition against a barcode drawn from the catalytic site: it cannot be
     a restatement of it.

     Wang et al. 2025 are in tension with this. They report that the PB repeats
     improve recognition of Glu-gamma-Thr/Ser and Asp-beta-Ala, which are peptide
     substitutions. Either the module contacts the peptide as well as the sugar,
     or avidity alone changes the apparent preference. Unresolved, and recorded
     here as unresolved.

pH
--
The three-motif domain binds pseudomurein completely at pH 9.0 (near its pI of
9.2), partially at pH 6.5, and not at all at pH 4.0. At pH 7.0 it aggregates into
high-molecular-mass complexes; at pH 9.0 it is a 17.4 kDa monomer.

Every published Pei whole-cell lysis assay runs at pH 7.0-7.85. The binding module
is therefore operating in its aggregation-prone, partially-binding regime in all
of them. PMB domains have pIs spanning 3-10, so this is a per-protein question and
`domain_arch.py` computes the pI of each protein's own PMB region rather than
assuming MTH719's.
"""
from __future__ import annotations

CITATION = ("Visweswaran, Dijkstra & Kok 2011, PLoS ONE 6(6):e21582 "
            "(doi:10.1371/journal.pone.0021582)")

# construct -> what it bound. `pseudomurein` is the intact methanogen sacculus;
# `spheroplast` is lysozyme-treated bacterial cells with murein fragments exposed.
CONSTRUCTS = {
    "1P-GFP": {"n_motifs": 1, "pseudomurein": False, "spheroplast": False,
               "intact_bacteria": False},
    "2P-GFP": {"n_motifs": 2, "pseudomurein": False, "spheroplast": True,
               "intact_bacteria": False},
    "3P-GFP": {"n_motifs": 3, "pseudomurein": True, "spheroplast": True,
               "intact_bacteria": False},
    # PeiW's own PMB domain (four motifs), fused to GFP, also bound spheroplasts
    "PeiW-PMB-GFP": {"n_motifs": 4, "pseudomurein": True, "spheroplast": True,
                     "intact_bacteria": False},
}

MOTIFS_FOR_PSEUDOMUREIN = 3     # fewer than this: no binding to the sacculus
MOTIFS_FOR_SPHEROPLAST = 2      # fewer than this: no binding at all

# The MTH719 three-motif domain. PMB pIs range 3-10, so do not reuse these
# numbers for another protein: compute its own.
MTH719 = {
    "uniprot": "O26815", "length": 574, "role": "putative S-layer protein",
    "pmb_span": (432, 574), "n_motifs": 3, "pi": 9.2,
    "mass_kda_monomer": 17.4,
    "note": "carries PMB motifs and NO C71 catalytic domain: PF09373 is not "
            "Pei-specific",
}

PH_BINDING = {4.0: "none", 6.5: "partial", 9.0: "complete"}
PH_AGGREGATION = {7.0: "high-molecular-mass aggregate", 9.0: "17.4 kDa monomer"}
PUBLISHED_ASSAY_PH = {
    "Schofield 2015 synthetic peptide": 7.85,
    "Schofield 2015 cell suspension / plate": 7.0,
    "Morii & Koga cell wall": "6.8-7.4",
}

# Denaturation midpoint, tryptophan fluorescence: the motifs fold into a domain
# rather than acting as an unstructured tether.
GDHCL_MIDPOINT_M = 3.8

# Proteases cut once between motifs 1 and 2 and nowhere inside the motifs, even
# though sites are predicted there. The motif region is structured; the spacer is
# not. Relevant if anyone tries to define module boundaries from the alignment.
PROTEOLYSIS = ("chymotrypsin and trypsin each cut once, in the spacer between "
               "motifs 1 and 2; predicted sites inside the motifs are protected")

# Organisms known in 2011 to carry PMB motifs. Useful as sentinels: the screen
# should recover all of them, and the three bacteria are the test of whether a
# bacterial PMBR hit is real. All three also carry LysM.
KNOWN_PMBR_ARCHAEA = ["s__Methanothermobacter thermautotrophicus",
                      "s__Methanosphaera stadtmanae"]
KNOWN_PMBR_VIRUSES = ["Methanobacterium prophage", "Methanothermobacter prophage"]
KNOWN_PMBR_BACTERIA = ["s__Xanthomonas campestris",
                       "s__Granulibacter bethesdensis",
                       "s__Novosphingobium aromaticivorans"]


def rule_has_jurisdiction(n_motifs):
    """The 3-motif rule is about PF09373 domains. It has nothing to say about a
    protein that carries none.

    PeiR (UniProt D3DZZ6) has zero PMB motifs and lyses Methanobrevibacter
    ruminantium M1 -- the very host PeiW and PeiP cannot touch. Whatever PeiR uses
    to reach the sacculus, it is not a PF09373 array. Reading
    `pmbr_binding_competent == 0` as "cannot lyse" therefore contradicts the one
    protein in the family with a solved structure and a rumen phenotype.
    """
    try:
        return int(n_motifs) > 0
    except (TypeError, ValueError):
        return False


def predict_binding(n_motifs):
    """What a protein with `n_motifs` PMB repeats can bind THROUGH THAT MODULE.

    Returns (call, can_bind_sacculus_via_pmbr, reason).

    Read the middle value narrowly. It is not "this protein cannot lyse cells".
    It is "this protein's PMB array, if it has one, cannot dock on a sacculus".
    A protein with no array is outside the rule entirely: see PeiR.
    """
    try:
        n = int(n_motifs)
    except (TypeError, ValueError):
        return "unknown", False, "no repeat count"
    if n >= MOTIFS_FOR_PSEUDOMUREIN:
        return ("pseudomurein_sacculus", True,
                f"{n} motifs: >= {MOTIFS_FOR_PSEUDOMUREIN}, the minimum that bound "
                f"methanogen cells")
    if n >= MOTIFS_FOR_SPHEROPLAST:
        return ("murein_fragments_only", False,
                f"{n} motifs: the two-motif construct bound lysozyme-treated "
                f"bacterial spheroplasts but NOT pseudomurein")
    if n == 1:
        return ("none", False, "the one-motif construct bound nothing")
    return ("no_pmbr_module", False,
            "no PMB repeat detected, so the 3-motif rule does not apply. PeiR has "
            "no PMB repeat and lyses M. ruminantium M1; a zero here is not a "
            "prediction that the protein cannot lyse")


# Bjellqvist pKa values. Enough to rank PMB domains against each other and to
# say whether the published assay pH is anywhere near a given domain's pI.
_PKA_POS = {"K": 10.0, "R": 12.0, "H": 6.0}
_PKA_NEG = {"D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
_PKA_NTERM, _PKA_CTERM = 9.69, 2.34


def isoelectric_point(seq, lo=0.0, hi=14.0, tol=1e-4):
    """Bisect the net-charge curve. Returns None for an empty sequence."""
    s = "".join(c for c in str(seq).upper() if c.isalpha())
    if not s:
        return None

    def charge(ph):
        q = 1.0 / (1.0 + 10 ** (ph - _PKA_NTERM))
        q -= 1.0 / (1.0 + 10 ** (_PKA_CTERM - ph))
        for aa, pk in _PKA_POS.items():
            q += s.count(aa) / (1.0 + 10 ** (ph - pk))
        for aa, pk in _PKA_NEG.items():
            q -= s.count(aa) / (1.0 + 10 ** (pk - ph))
        return q

    while hi - lo > tol:
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def assay_ph_advice(pi):
    """The binding module works near its pI. Say so, per protein."""
    if pi is None:
        return "no PMB region; the binding-pH question does not arise"
    published = 7.0
    if abs(pi - published) < 1.0:
        return (f"PMB pI {pi:.1f}: the published lysis assay pH ({published}) is "
                f"near it, but MTH719's domain aggregates at pH 7.0. Check the "
                f"oligomeric state before reading a negative")
    return (f"PMB pI {pi:.1f}: bind at pH ~{pi:.1f}. Every published Pei lysis "
            f"assay runs at pH 7.0-7.85, where the MTH719 domain aggregates and "
            f"binds pseudomurein only partially. A negative at pH 7 may be a "
            f"binding failure, not a catalytic one")


if __name__ == "__main__":
    print(f"Reference: {CITATION}\n")
    print("Constructs:")
    for k, v in CONSTRUCTS.items():
        print(f"  {k:14s} n={v['n_motifs']}  pseudomurein={str(v['pseudomurein']):5s} "
              f"spheroplast={str(v['spheroplast']):5s} "
              f"intact_bacteria={v['intact_bacteria']}")
    print(f"\nThresholds: sacculus >= {MOTIFS_FOR_PSEUDOMUREIN} motifs, "
          f"murein fragments >= {MOTIFS_FOR_SPHEROPLAST}\n")
    for n in range(0, 6):
        call, ok, why = predict_binding(n)
        print(f"  {n} motif(s) -> {call:22s} sacculus={str(ok):5s}  {why}")

    print("\nConsistency with the measured constructs:")
    for k, v in CONSTRUCTS.items():
        _c, ok, _ = predict_binding(v["n_motifs"])
        flag = "ok " if ok == v["pseudomurein"] else "BAD"
        print(f"  {flag} {k}: predicted sacculus={ok}, observed={v['pseudomurein']}")

    print(f"\npH: {PH_BINDING}")
    print(f"Published assay pH: {PUBLISHED_ASSAY_PH}")
    print(f"\nMTH719 PMB pI (reported 9.2), recomputed from nothing: "
          f"n/a -- sequence not stored here; isoelectric_point() is applied to "
          f"each protein's own PMB span in domain_arch.py")
