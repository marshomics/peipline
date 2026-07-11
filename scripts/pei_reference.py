#!/usr/bin/env python3
"""Every experimentally characterised pseudomurein endoisopeptidase, and what it is.

This module exists because the screen was built around PeiW and PeiP, and PeiW and
PeiP are not representative of the family.

The finding that reorganised the pipeline
-----------------------------------------
PeiR (UniProt D3DZZ6, Methanobrevibacter ruminantium M1, prophage Phi-mru) is
annotated with **PF03412 Peptidase_C39 and nothing else**. No PF12386. No PF09373.
It is 228 residues, it has no pseudomurein-binding repeat at all, and it lyses
M. ruminantium cells -- the Thr-walled organism that PeiW and PeiP cannot touch.

Three consequences, each of which breaks something that used to be in this repo:

  1. A PF12386-only screen never sees PeiR. SSF54001 might: PeiR's Gene3D
     assignment is 3.90.70.10 "Cysteine proteinases", which is what SCOP 54001
     covers. So PeiR is recoverable, but only through the low-specificity model.

  2. The catalytic-triad spacing prior does not transfer. In C71 the Cys->His gap
     is 35 (PeiW C198/H233, PeiP C213/H248). PeiR's catalytic cysteine is C90
     (established by the C90A mutant, PDB 8Z4N) and residue 125 -- C90+35 -- is a
     proline. NO histidine sits at the C71 distance: PeiR's downstream histidines
     are at 104, 162, 181 and 192 (gaps 14, 72, 91, 102 from C90), and none falls
     within 35+/-8. PeiR's catalytic histidine has NOT been assigned experimentally,
     so there is no established Cys->His gap; the point is only that the C71 prior
     of 35 matches none of PeiR's histidines and so would reject it.
     `peir_his_candidates()` recomputes these gaps from the stored sequence and
     `pei_check` prints them, so this claim is checkable, not asserted.

  3. PeiR has zero PMB motifs and lyses cells. `pmbr_binding_competent == 0` is
     therefore NOT a prediction that a protein cannot lyse. The Visweswaran 2011
     three-motif rule is a statement about PF09373-containing binding domains. It
     says nothing about a protein that docks by other means, and PeiR is the
     counterexample.

There is no PeiS
----------------
Searched: UniProtKB, PubMed, the Pei literature. The characterised set is PeiW,
PeiP, PeiR, and the three CRISPRTarget-discovered viral proteins PeiG2, PeiTh2 and
PeiF3. If you have a sequence for something called PeiS, add it here; the module
will not invent one.

Substrate caution
-----------------
PeiG2*, PeiTh2* and PeiF3* were reported unable to cleave L-Ala-pNA yet able to
lyse Methanobrevibacter sp. AbM4. That is not evidence they lack the activity.
Schofield et al. 2015 showed Ala-pNA is a poor substrate for PeiW and PeiP as well:
the P2 glutamate supplies most of the binding energy, and the alanine dipeptide is
too small to engage S2. A negative on Ala-pNA is uninformative about any Pei.
"""
from __future__ import annotations

CITATIONS = {
    "PeiW/PeiP": "Luo et al. 2002, FEMS Microbiol Lett 208:47-51 (PeiW/PeiP naming); Schofield et al. 2015 (biochemistry)",
    "PeiR": ("Leahy et al. 2010, PLoS ONE 5:e8926 (genome); PDB 8Z4H, 8Z4N "
             "(Guo, Wang & Bai 2024)"),
    "structures": ("Wang et al. 2025, Int J Biol Macromol "
                   "(doi:10.1016/j.ijbiomac.2025.141813); PeiW-CD 8JX4, PeiP 8Z4F"),
    "new_viral": ("Identification and functional characterisation of novel viral "
                  "pseudomurein endoisopeptidase proteins from methanogens "
                  "(CRISPRTarget study; 21 candidates, 3 active)"),
    "binding": "Visweswaran, Dijkstra & Kok 2011, PLoS ONE 6:e21582",
    "substrates": "Schofield et al. 2015, Archaea 2015:828693",
}

# PeiR, verbatim from UniProt D3DZZ6. Kept here because every claim below about
# its residue numbering is checkable against it, and because the pipeline's
# `pei_check` rule re-derives those claims rather than trusting this comment.
PEIR_SEQUENCE = (
    "MVRFSRDMLQDGAKRMFKWLRKGEGLPNYLIMYDMDRNKEYKLVPKEYAG"
    "LYESRNIFWIKNGREPNYVTLTSVARNPLVMDYQNTNYTCCPTSLSLASQ"
    "MLYHYKSESECAKALGTSKGSGTSPAQLIANAPKLGFKIIPIKRDSKEVK"
    "KYLKKGFPVICHWQVNQSRNCKGDYTGNFGHYGLIWDMTSTHYVVADPAK"
    "GVNRKYKFSCLDNANKGYRQNYYVVCPA"
)

# family -> the two things that must never be shared across families
FAMILIES = {
    "c71": {
        "pfam": "PF12386", "merops": "C71",
        "structures": {"8JX4": "PeiW catalytic domain", "8Z4F": "PeiP full length"},
        "catalytic_cys_gap_to_his": 35,     # PeiW C198->H233; PeiP C213->H248
        "his_gap_to_asp": (17, 24),         # PeiW 17, PeiP 24. Do not narrow.
        "pei_class_applies": True,          # Wang et al.'s V252/C265 partition
    },
    "c39": {
        "pfam": "PF03412", "merops": "C39",
        "structures": {"8Z4H": "PeiR", "8Z4N": "PeiR C90A"},
        "catalytic_cys_gap_to_his": None,   # 72 in PeiR; ONE sequence is not a prior
        "his_gap_to_asp": None,
        "pei_class_applies": False,         # V252/C265 are C71 numbering
        "warning": ("MEROPS C39 proper is the bacteriocin-processing peptidase "
                    "domain of ABC transporters (double-glycine leader cleavage). "
                    "Most PF03412 hits across 350,000 proteomes will be exporters, "
                    "not Pei. The triad filter is the only thing standing between "
                    "this arm and a very large pile of transporters, so its FDR "
                    "must be reported separately per family."),
    },
}

# Characterised proteins. `pmbr_motifs` is the count of PF09373 in UniProt, which
# is a lower bound: a 30-35 residue repeat is easy for a domain E-value to drop.
CHARACTERISED = {
    "PeiW": {
        "uniprot": "Q7LYX0", "length": 284, "family": "c71",
        "pfam": ["PF12386", "PF09373"], "pmbr_motifs": 4,
        "source": "Methanothermobacter phage psiM100",
        "host": "s__Methanothermobacter wolfeii",
        "catalytic": {"cys": 198, "his": 233, "asp": 250},
        "structure": "8JX4 (catalytic domain)",
        "lyses": True, "cleaves_EgammaA_pNA": True,
    },
    "PeiP": {
        "uniprot": "Q77WJ4", "length": 305, "family": "c71",
        "pfam": ["PF12386", "PF09373"], "pmbr_motifs": 4,
        "source": "Methanobacterium phage psiM2",
        "host": "s__Methanothermobacter marburgensis",
        "catalytic": {"cys": 213, "his": 248, "asp": 272},
        "structure": "8Z4F (full length)",
        "lyses": True, "cleaves_EgammaA_pNA": True,
    },
    "PeiR": {
        "uniprot": "D3DZZ6", "length": 228, "family": "c39",
        "pfam": ["PF03412"], "pmbr_motifs": 0,
        "source": "prophage Phi-mru",
        "host": "s__Methanobrevibacter ruminantium",
        # C90 only. The His and Asp have NOT been assigned experimentally; the
        # C90A mutant is what fixes the cysteine. Do not fill these in from the
        # C39 consensus and then use them as evidence.
        "catalytic": {"cys": 90, "his": None, "asp": None},
        "pfam_domain_span": (83, 204),
        "disulfide": (171, 210),
        "structure": "8Z4H (native), 8Z4N (C90A)",
        "lyses": True, "cleaves_EgammaA_pNA": None,   # not tested
        "note": ("lyses its own Thr-walled host, which PeiW and PeiP cannot. "
                 "Carries NO pseudomurein-binding repeat."),
    },
    # The three CRISPRTarget proteins. Accessions were not given in the abstract
    # and I have not resolved them; recorded as unknown rather than guessed.
    "PeiG2": {
        "uniprot": None, "length": None, "family": "c39",
        "pfam": ["PF03412"], "pmbr_motifs": None,
        "source": "methanogen virus (CRISPRTarget)", "host": None,
        "catalytic": {"cys": None, "his": None, "asp": None},
        "lyses": True, "cleaves_L_Ala_pNA": False,
        "note": "lyses Methanobrevibacter sp. AbM4; chosen for homology to PeiR",
    },
    "PeiTh2": {
        "uniprot": None, "length": None, "family": "c39",
        "pfam": ["PF03412"], "pmbr_motifs": None,
        "source": "methanogen virus (CRISPRTarget)", "host": None,
        "catalytic": {"cys": None, "his": None, "asp": None},
        "lyses": True, "cleaves_L_Ala_pNA": False,
        "note": "more effective than PeiR on Methanobrevibacter sp. AbM4",
    },
    "PeiF3": {
        "uniprot": None, "length": None, "family": "c39",
        "pfam": ["PF03412", "PF09373?"], "pmbr_motifs": None,
        "source": "methanogen virus (CRISPRTarget)", "host": None,
        "catalytic": {"cys": None, "his": None, "asp": None},
        "lyses": True, "cleaves_L_Ala_pNA": False,
        "note": ("mutating conserved PMBR prolines killed binding; mutating the "
                 "C39 active-site residues killed lysis. So this one has both a "
                 "binding module and a C39 catalytic domain"),
    },
}

# Named in the request, absent from the literature.
NOT_FOUND = {
    "PeiS": ("no protein of this name in UniProtKB or the Pei literature. The "
             "characterised set is PeiW, PeiP, PeiR, PeiG2, PeiTh2, PeiF3. "
             "Supply a sequence if you have one; nothing here will guess."),
}

# UniProtKB proteins carrying PF12386, as of this writing. Two are bacterial and
# one is fungal, which is the whole reason for screening outside the archaea.
PF12386_CROSS_DOMAIN = {
    "A0A090AJG3": ("Thioploca ingrica", "Bacteria; Gammaproteobacteria", 121),
    "A0ABY4WHD8": ("Brevibacillus ruminantium", "Bacteria; Bacillota", 135),
    "A0A2T3YUS1": ("Trichoderma asperellum", "Eukaryota; Fungi", 989),
}


def catalytic_cys_his_gap(name):
    """Observed Cys->His gap, or None if the His was never assigned."""
    c = CHARACTERISED[name]["catalytic"]
    if c["cys"] is None or c["his"] is None:
        return None
    return c["his"] - c["cys"]


def peir_his_candidates():
    """Histidines downstream of PeiR's catalytic C90, with their gaps.

    Recomputed from the sequence every call. The point is not to assign the His;
    it is to show that no histidine sits at the C71 distance of 35.
    """
    cys = CHARACTERISED["PeiR"]["catalytic"]["cys"]
    return {i + 1: (i + 1) - cys
            for i, aa in enumerate(PEIR_SEQUENCE)
            if aa == "H" and i + 1 > cys}


def spacing_prior_admits(name, gap_1_2, tolerance):
    """Would a Cys->His spacing prior of `gap_1_2 +/- tolerance` find this protein?

    Returns (admitted, reason). For PeiR this must be False under the C71 prior,
    and `pei_check.py` turns that into a hard failure if the config ever tries.
    """
    e = CHARACTERISED[name]
    cys = e["catalytic"]["cys"]
    if cys is None:
        return None, "catalytic cysteine not assigned"
    if e["catalytic"]["his"] is not None:
        gap = e["catalytic"]["his"] - cys
        ok = abs(gap - gap_1_2) <= tolerance
        return ok, f"observed Cys->His gap {gap}; prior {gap_1_2}+/-{tolerance}"
    if name == "PeiR":
        cands = peir_his_candidates()
        hits = [h for h, g in cands.items() if abs(g - gap_1_2) <= tolerance]
        if hits:
            return True, f"His at {hits} would satisfy the prior"
        near = min(cands.values())
        return False, (f"no histidine downstream of C{cys} lies within "
                       f"{tolerance} of gap {gap_1_2}; the nearest is +{near} "
                       f"(residue {cys + near}), and C{cys}+{gap_1_2} is "
                       f"'{PEIR_SEQUENCE[cys + gap_1_2 - 1]}'")
    return None, "no His assigned and no sequence stored"


def pmbr_rule_applies(name):
    """The 3-motif binding rule is about PF09373 domains, not about lysis.

    PeiR has no PMB motif and lyses cells. Reading `pmbr_binding_competent == 0`
    as "cannot lyse" is therefore wrong for any protein with zero motifs: the rule
    has no jurisdiction there.
    """
    n = CHARACTERISED[name]["pmbr_motifs"]
    if n is None:
        return None
    return n > 0


if __name__ == "__main__":
    print("Characterised pseudomurein endoisopeptidases\n")
    hdr = f"{'name':8s} {'uniprot':10s} {'fam':5s} {'pfam':22s} {'pmbr':5s} {'lyses':6s}"
    print(hdr)
    print("-" * len(hdr))
    for n, e in CHARACTERISED.items():
        print(f"{n:8s} {str(e['uniprot'] or '-'):10s} {e['family']:5s} "
              f"{','.join(e['pfam']):22s} {str(e['pmbr_motifs']):5s} "
              f"{str(e['lyses']):6s}")

    print("\nPeiR sequence checks (recomputed, not asserted):")
    print(f"  length {len(PEIR_SEQUENCE)}  residue 90 = {PEIR_SEQUENCE[89]} "
          f"(catalytic, from the C90A mutant)")
    print(f"  downstream His and their gaps from C90: {peir_his_candidates()}")

    print("\nWould the C71 spacing prior retain each characterised protein?")
    for n in ("PeiW", "PeiP", "PeiR"):
        ok, why = spacing_prior_admits(n, FAMILIES["c71"]["catalytic_cys_gap_to_his"], 8)
        mark = {True: "keep", False: "DISCARD", None: "n/a "}[ok]
        print(f"  {n:6s} {mark:8s} {why}")

    print("\nDoes the PMBR 3-motif rule have jurisdiction?")
    for n in ("PeiW", "PeiP", "PeiR"):
        print(f"  {n:6s} {pmbr_rule_applies(n)}")

    print("\nPF12386 outside the archaea:")
    for acc, (org, lineage, ln) in PF12386_CROSS_DOMAIN.items():
        print(f"  {acc:12s} {org:26s} {lineage:32s} {ln} aa")

    print()
    for n, why in NOT_FOUND.items():
        print(f"{n}: {why}")
