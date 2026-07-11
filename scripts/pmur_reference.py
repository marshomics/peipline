#!/usr/bin/env python3
"""Which genes actually diagnose pseudomurein, and which only look like they do.

Single source: Lupo V, Roomans C, Royen E, Ongena L, Jacquemin O, Mullender C,
Kerff F, Baurain D (2025) "Identification and characterization of archaeal
pseudomurein biosynthesis genes through pangenomics", mSystems 10(3):e...
(doi:10.1128/msystems.01201-24). Pangenomic OGs from FIVE PM-containing archaea
(four Methanobacteriales + one Methanopyrales).

Why a naive "screen for Mur ligases" is wrong
---------------------------------------------
Of the 49 orthologous groups shared by all five PM genomes:

    15  widespread across Bacteria AND Archaea   -> screening on these lights up
                                                    half the prokaryotic tree
     9  homologs in bacteria + PM-archaea only
    25  archaea-exclusive, of which only
    15  exclusive to PM-containing archaea        -> the only diagnostic tier

The Mur-domain family is the specific trap. MurT, CapB, CphA and FPGS all carry a
Mur domain and are scattered through both domains, so a bare "Mur ligase" hit
tells you close to nothing. Lupo et al. show the archaeal muramyl ligases are of
BACTERIAL origin (HGT + duplication), which is exactly why the domain alone is
not PM-specific.

The signal that IS PM-exclusive
-------------------------------
The muramyl ligases Murα/β/γ/δ, the MraY-like glycosyltransferase (OG0001163),
and OG0001014 (a small CPS) have no homologs outside Methanopyrales and
Methanobacteriales. In the ATP-grasp and MraY-like families the archaeal members
did NOT come from Bacteria. That handful is what a screen should use.

Two problems that remain even with the right subset
---------------------------------------------------
1. Circularity. The OGs were defined from five genomes already known to make PM,
   so HMMs built from n=5 encode the diversity of exactly those five organisms.
   They recover close relatives and miss anything divergent -- precisely the
   lineage a discovery screen is meant to catch. A negative outside the known
   orders is therefore weak evidence, not absence.

2. PM tracks taxonomy almost perfectly. It is restricted to two orders, so "does
   this genome make PM?" is currently near-equivalent to "is this a
   Methanobacteriales or Methanopyrales?" GTDB placement answers that faster and
   more reliably than any OG screen. The OG set only earns its keep when you are
   asking whether PM turns up somewhere UNEXPECTED, outside those two orders --
   and for that, co-occurrence in synteny is a far stronger claim than scattered
   individual hits.

None of this confirms PM. Genomic presence is a hypothesis. Confirmation is
biochemical: the wall composition itself -- N-acetyl-talosaminuronic acid (TalNAc
/ NAT) and the beta-(1->3) glycan linkage.
"""
from __future__ import annotations

CITATION = ("Lupo et al. 2025, mSystems 10:e01201-24 "
            "(doi:10.1128/msystems.01201-24)")

# The two orders that contain pseudomurein. Taxonomy is the primary PM signal
# INSIDE these; the OG screen is for testing whether PM appears OUTSIDE them.
PM_ORDERS = {"o__Methanobacteriales", "o__Methanopyrales"}

# The PM-EXCLUSIVE orthologous groups: no homologs outside the two orders. These
# are the muramyl ligases plus the MraY-like GT plus the small CPS. Roles as
# assigned by Lupo et al.; the four Mur names are their arbitrary labels for the
# Mur-domain OGs and imply no functional order.
PM_EXCLUSIVE = {
    "OG0001148": {"name": "Muralpha", "role": "muramyl ligase (Mur domain)"},
    "OG0001149": {"name": "Murbeta",  "role": "muramyl ligase (Mur domain)"},
    "OG0001150": {"name": "Murgamma", "role": "muramyl ligase (Mur domain)", "cluster": "B"},
    "OG0001473": {"name": "Murdelta", "role": "muramyl ligase (Mur domain)", "cluster": "A"},
    "OG0001163": {"name": "MraY_like", "role": "type-4 glycosyltransferase, MraY homolog",
                  "cluster": "A"},
    "OG0001014": {"name": "CPS_OG1014", "role": "smallest CPS, experimentally characterized",
                  "cluster": "A"},
}

# The muramyl ligases specifically (Muralpha..delta). "At least one of these"
# is the ligase requirement of the syntenic block.
MURAMYL_LIGASES = ["OG0001148", "OG0001149", "OG0001150", "OG0001473"]

# The two anchors that make the block PM-specific rather than "some Mur ligase".
BLOCK_ANCHORS = ["OG0001163", "OG0001014"]

# Syntenic clusters A and B (Lupo et al., Fig. 1 / synteny analysis). A genome
# that makes PM carries these arranged as clusters; scattered hits do not.
CLUSTER_A = ["OG0001014", "OG0001163", "OG0001473", "OG0001162", "OG0001472"]
CLUSTER_B = ["OG0001150", "OG0001147", "OG0001146"]
# Muralpha/beta float between A and B or sit outside, depending on species.
FLOATING = ["OG0001148", "OG0001149"]

# The five genes Lupo et al. flag as having NO known function -- new candidates
# worth characterizing, but do not read presence as proof of anything.
UNCHARACTERIZED = ["OG0001162", "OG0001472", "OG0001147", "OG0000796", "OG0000169"]

# The trap. These carry a Mur domain and are scattered across BOTH domains of
# life; a hit means "has a Mur domain", not "makes pseudomurein". Never put these
# in the marker set, and be suspicious of any user HMM that is really one of them.
MUR_DOMAIN_TRAP = {
    "MurT": "cell-wall amidotransferase, widespread in Firmicutes",
    "CapB": "capsule biosynthesis, widespread",
    "CphA": "cyanophycin synthetase, widespread",
    "FPGS": "folylpolyglutamate synthase, universal",
}


def parse_order(classification):
    """Pull the GTDB order (o__...) out of a full lineage string, or None."""
    if not isinstance(classification, str):
        return None
    for field in classification.split(";"):
        field = field.strip()
        if field.startswith("o__"):
            return field or None
    return None


def pm_expected_by_taxonomy(classification):
    """True iff the genome is in an order known to contain pseudomurein.

    This is the primary, high-confidence PM signal. It is genomic placement, not
    wall chemistry, so it still is not confirmation -- but it is more reliable
    than an n=5 OG screen and it is free.
    """
    return parse_order(classification) in PM_ORDERS


def block_status(present_ogs):
    """Given the set of PM-exclusive OGs a genome carries, describe the block.

    Returns (has_block, n_ligases, has_mray, has_cps). `has_block` is the
    CO-PRESENCE requirement: at least one muramyl ligase AND the MraY-like AND
    the CPS. It is weaker than positional synteny (which needs coordinates) but
    far stronger than "any three markers from a bag", because the two anchors are
    what make it PM-specific rather than "carries a Mur domain".
    """
    present = set(present_ogs)
    n_lig = sum(og in present for og in MURAMYL_LIGASES)
    has_mray = "OG0001163" in present
    has_cps = "OG0001014" in present
    has_block = n_lig >= 1 and has_mray and has_cps
    return has_block, n_lig, has_mray, has_cps


if __name__ == "__main__":
    print(f"Reference: {CITATION}\n")
    print("PM-exclusive OGs (the only diagnostic tier):")
    for og, d in PM_EXCLUSIVE.items():
        print(f"  {og}  {d['name']:12s} {d['role']}"
              + (f"  [cluster {d['cluster']}]" if "cluster" in d else ""))
    print(f"\nSyntenic block = >=1 of {MURAMYL_LIGASES} (muramyl ligase)")
    print(f"                 AND OG0001163 (MraY-like) AND OG0001014 (CPS)")
    print(f"\nCluster A: {CLUSTER_A}")
    print(f"Cluster B: {CLUSTER_B}")
    print(f"\nDO NOT screen on the Mur-domain trap: {list(MUR_DOMAIN_TRAP)}")
    print(f"\nPM orders (taxonomy is primary here): {sorted(PM_ORDERS)}")
    print("\nSanity:")
    for lin, want in [
        ("d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__x;g__y;s__z", True),
        ("d__Archaea;p__Methanobacteriota_B;c__x;o__Methanopyrales;f__;g__;s__", True),
        ("d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__;g__;s__", False),
    ]:
        got = pm_expected_by_taxonomy(lin)
        print(f"  {'ok ' if got == want else 'BAD'} {parse_order(lin)} -> pm_expected={got}")
    hb, nl, hm, hc = block_status({"OG0001148", "OG0001163", "OG0001014"})
    print(f"\n  block (ligase+MraY+CPS present): has_block={hb} "
          f"(ligases={nl}, mray={hm}, cps={hc})")
    hb2, *_ = block_status({"OG0001148", "OG0001149"})   # ligases only, no anchors
    print(f"  block (ligases only, no anchors): has_block={hb2}  <- must be False")
