#!/usr/bin/env python3
"""What the pseudomurein cross-link is actually made of, per characterised strain.

Single source: Kandler O & Koenig H (1978) "Chemical Composition of the
Peptidoglycan-Free Cell Walls of Methanogenic Bacteria", Arch Microbiol
118:141-152. Seven strains, wet chemistry, amino-acid analysis after total acid
hydrolysis, configuration checked with D-amino-acid oxidase and L-amino-acid
decarboxylases. Nothing in this module comes from a secondary summary.

Why P1 is the alanine
---------------------
Pei cleaves the epsilon-isopeptide bond between alanine and lysine, so the
alanine donates the carbonyl and is P1. Kandler & Koenig report the molar ratio
Lys : Ala : Glu as approximately 1 : 1.2 : 2. One alanine per subunit. Therefore
the alanine that M. ruminantium M1 replaces with threonine is not "an" alanine,
it is *the* alanine, and it is P1. That was the open question in an earlier
draft of this pipeline; it is now closed.

Two axes, not one
-----------------
The paper reports two independent substitutions, and they sit on opposite sides
of the scissile bond:

    P1  (acyl donor)     Ala -> Thr, completely, in M. ruminantium M1
    P1' (acyl acceptor)  Lys -> Orn, about one quarter, in M. smithii PS

An ornithine acceptor shortens the amine arm by one methylene, so the
isopeptide becomes delta(Ala)-Orn rather than epsilon(Ala)-Lys. That is an S1'
question, not an S1 question, and the published PeiW/PeiP chromogenic assay
(H-Glu-gamma-Ala-pNA) contains no acceptor residue at all, so it cannot test it.

Why genus-level assignment is not merely uncertain but wrong
------------------------------------------------------------
Both strains the 1978 paper calls "M. ruminantium" now sit in Methanobrevibacter,
in different species: M1 is the type strain of M. ruminantium, PS is the type
strain of M. smithii. A third species of the same genus, M. arboriphilus
(then "M. arbophilicum"), is Ala-type. So Methanobrevibacter contains Ala-type,
Thr-type and Orn-modified walls simultaneously. Assigning P1 from a genus label
would put every gut Methanobrevibacter into one bucket and the bucket would be
wrong. `chemistry_for()` refuses to do it.
"""
from __future__ import annotations

CITATION = ("Kandler & Koenig 1978, Arch Microbiol 118:141-152 "
            "(doi:10.1007/BF00415722)")

# modern species -> what the paper measured, on which strain
#
# `p1`  : residue donating the carbonyl of the scissile isopeptide bond
# `p1p` : residue donating the amine (epsilon-Lys, or delta-Orn)
REFERENCE = {
    "s__Methanobrevibacter ruminantium": {
        "strain": "M1", "p1": "Thr", "p1_prime": "Lys",
        "note": "alanine COMPLETELY replaced by threonine",
        "level": "type_strain",
    },
    "s__Methanobrevibacter smithii": {
        "strain": "PS", "p1": "Ala", "p1_prime": "Lys/Orn",
        "note": "~1/4 of lysine replaced by ornithine; P1 unchanged",
        "level": "type_strain",
    },
    "s__Methanobrevibacter arboriphilus": {
        "strain": "(as M. arbophilicum)", "p1": "Ala", "p1_prime": "Lys",
        "note": "", "level": "species",
    },
    "s__Methanobacterium formicicum": {
        "strain": "", "p1": "Ala", "p1_prime": "Lys", "note": "", "level": "species",
    },
    "s__Methanobacterium bryantii": {
        "strain": "M.o.H.", "p1": "Ala", "p1_prime": "Lys",
        "note": "reported as 'Methanobacterium spec. M.o.H.'", "level": "strain",
    },
    "s__Methanothermobacter thermautotrophicus": {
        "strain": "", "p1": "Ala", "p1_prime": "Lys",
        "note": "reported as M. thermoautotrophicum", "level": "species",
    },
}

# No pseudomurein sacculus at all: a protein sheath. A useful negative control,
# and a genome that must never be scored as Ala-type.
NO_PSEUDOMUREIN = {"s__Methanospirillum hungatei", "s__Methanospirillum hungatii"}

# Claims that circulated in secondary summaries and are NOT in Kandler & Koenig
# 1978. They may well be true; they need their own citation before this pipeline
# will act on them. Encountering one raises a warning rather than a silent value.
DISPUTED = {
    "s__Methanosphaera stadtmanae":
        "serine-for-alanine is not in Kandler & Koenig 1978 (the species was "
        "described in 1985). Corroborating but still secondary: Schofield et al. "
        "2015 cleave Glu-gamma-Ser-pNA with both PeiW and PeiP, and lyse "
        "DSM 3091, attributing that lysis to a serine wall. Consistent; not "
        "primary. Supply the chemistry reference before using it.",
    "s__Methanopyrus kandleri":
        "an ornithine modification is not in Kandler & Koenig 1978 (the species "
        "was described in 1991). Supply the primary reference before using it.",
}

# Residues this enzyme family has been shown to accept or reject at P1, from the
# chromogenic series of Schofield et al. 2015 (see lysis_reference.py). Kept here so
# a chemistry call of "Ser" is not treated as exotic: it is a residue PeiW and
# PeiP demonstrably cut.
#
# Thr is the interesting one. It is rejected in the pNA series, and the Thr-walled
# M. ruminantium M1 resists whole-cell lysis by both enzymes -- yet Wang et al.
# 2025 report PeiW cleaving the Glu-gamma-Thr-epsilon-Lys isopeptide. The two
# substrates differ at S1', not at S1. Nothing in this module resolves that;
# lysis_reference.py states it and assay_panel.py acts on it.
P1_CLEAVED_BY_PEIW_PEIP = {"Ala", "Ser"}
P1_NOT_CLEAVED_PNA = {"Thr"}


def genus_of(species):
    if not isinstance(species, str) or not species.startswith("s__"):
        return None
    return "g__" + species[3:].split()[0]


def _genus_index():
    idx = {}
    for sp, d in REFERENCE.items():
        idx.setdefault(genus_of(sp), set()).add((d["p1"], d["p1_prime"]))
    return idx


GENUS_CHEMISTRIES = _genus_index()


def genus_is_homogeneous(genus):
    """True only if every characterised species of the genus agrees on BOTH axes."""
    return len(GENUS_CHEMISTRIES.get(genus, set())) == 1


def chemistry_for(species):
    """Return (p1, p1_prime, source, note).

    Species-level only. A genus fallback is offered exactly when every
    characterised species of that genus agrees; for Methanobrevibacter it never
    does, which is the whole point.
    """
    if not isinstance(species, str):
        return "unknown", "unknown", "none", ""

    if species in NO_PSEUDOMUREIN:
        return "no_pseudomurein", "no_pseudomurein", "literature_species", \
            "protein sheath, no sacculus"

    if species in DISPUTED:
        return "unsupported", "unsupported", "disputed", DISPUTED[species]

    if species in REFERENCE:
        d = REFERENCE[species]
        return d["p1"], d["p1_prime"], f"literature_{d['level']}", d["note"]

    g = genus_of(species)
    if g and genus_is_homogeneous(g):
        p1, p1p = next(iter(GENUS_CHEMISTRIES[g]))
        return p1, p1p, "literature_genus_homogeneous", \
            f"no characterised strain of {species}; every characterised " \
            f"species of {g} agrees"

    if g in GENUS_CHEMISTRIES:
        return "unknown", "unknown", "genus_heterogeneous", (
            f"{g} contains more than one characterised wall chemistry "
            f"({sorted(GENUS_CHEMISTRIES[g])}); a genus-level call would be wrong")

    return "unknown", "unknown", "none", ""


def annotate(species_series):
    """Vectorised: returns a DataFrame with p1_residue, p1_prime_residue,
    p1_source, p1_note."""
    import pandas as pd
    rows = [chemistry_for(s) for s in species_series]
    return pd.DataFrame(rows, index=species_series.index,
                        columns=["p1_residue", "p1_prime_residue",
                                 "p1_source", "p1_note"])


if __name__ == "__main__":
    print(f"Reference: {CITATION}\n")
    for sp in sorted(REFERENCE) + sorted(DISPUTED) + sorted(NO_PSEUDOMUREIN):
        p1, p1p, src, note = chemistry_for(sp)
        print(f"{sp:48s} P1={p1:16s} P1'={p1p:14s} [{src}] {note}")
    print()
    for g, chems in sorted(GENUS_CHEMISTRIES.items()):
        print(f"{g:28s} homogeneous={genus_is_homogeneous(g)!s:5s} {sorted(chems)}")
