#!/usr/bin/env python3
"""Test the taxonomy-first / syntenic-block PM logic (Lupo et al. 2025).

The advice this defends: for genomes inside the two PM orders, taxonomy beats an
OG screen; outside them, only the co-present PM-exclusive block is a signal, and
even that is a hypothesis, not confirmed pseudomurein.

Checks:

  1. an in-order genome is called pseudomurein_expected_by_taxonomy WITHOUT any
     markers -- taxonomy is primary;
  2. an out-of-order genome carrying the whole block is flagged
     pseudomurein_candidate_out_of_order;
  3. an out-of-order genome with scattered non-block hits is NOT flagged (it is
     markers_without_block) -- the anchors are what make the block PM-specific;
  4. the block requires BOTH anchors: ligases alone do not qualify;
  5. the Mur-domain trap families are recognised and refused;
  6. `classify_markers` maps OG ids and Mur names to the right roles.

Run:  python test/test_pmur_reference.py
"""
from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

FAILURES = []


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


ARCH = "d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter smithii"
PYR = "d__Archaea;p__Methanobacteriota_B;c__Methanopyri;o__Methanopyrales;f__Methanopyraceae;g__Methanopyrus;s__Methanopyrus kandleri"
BACT = "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__;g__;s__"
HALO = "d__Archaea;p__Halobacteriota;c__Halobacteria;o__Halobacteriales;f__;g__;s__"


def main() -> None:
    import pandas as pd
    import pmur_reference as PR

    # --- reference facts ----------------------------------------------------
    check(PR.PM_ORDERS == {"o__Methanobacteriales", "o__Methanopyrales"},
          "the two PM orders are wrong")
    check(PR.pm_expected_by_taxonomy(ARCH) and PR.pm_expected_by_taxonomy(PYR),
          "both PM orders must be recognised")
    check(not PR.pm_expected_by_taxonomy(BACT) and not PR.pm_expected_by_taxonomy(HALO),
          "a bacterium and a non-PM archaeon must NOT be PM-expected")
    check(set(PR.MURAMYL_LIGASES) == {"OG0001148", "OG0001149", "OG0001150", "OG0001473"},
          "the four muramyl-ligase OGs are wrong")
    check(PR.block_status({"OG0001148", "OG0001163", "OG0001014"})[0],
          "ligase + MraY + CPS must satisfy the block")
    check(not PR.block_status({"OG0001148", "OG0001149"})[0],
          "ligases without the anchors must NOT satisfy the block")
    check(not PR.block_status({"OG0001163", "OG0001014"})[0],
          "anchors without a ligase must NOT satisfy the block")
    print("1. reference:        two PM orders, four ligases, block needs a ligase "
          "AND both anchors")

    # --- classify_markers ---------------------------------------------------
    from cellwall_genotype import classify_markers
    roles = classify_markers(
        ["OG0001148", "Murgamma", "OG0001163", "OG0001014", "MurT", "some_hyp"])
    check(roles["OG0001148"] == "muramyl_ligase", f"OG0001148 -> {roles['OG0001148']}")
    check(roles["Murgamma"] == "muramyl_ligase", f"Murgamma -> {roles['Murgamma']}")
    check(roles["OG0001163"] == "mray_like", f"OG0001163 -> {roles['OG0001163']}")
    check(roles["OG0001014"] == "cps", f"OG0001014 -> {roles['OG0001014']}")
    check(roles["MurT"] == "trap", f"MurT must be flagged as a trap, got {roles['MurT']}")
    check(roles["some_hyp"] == "other", f"unknown marker -> {roles['some_hyp']}")
    print("2. classify_markers: OG ids and Mur names map to roles; MurT is a trap")

    # --- end-to-end call logic: drive the REAL cellwall_genotype code --------
    from cellwall_genotype import annotate_pathway_calls
    # marker set = the six PM-exclusive OGs; ligase/anchor roles as classified
    mk = ["OG0001148", "OG0001163", "OG0001014", "MurT_stray"]
    lig, mray, cps = ["OG0001148"], ["OG0001163"], ["OG0001014"]

    # five genomes, one per scenario, with explicit pmur_<marker> columns
    rows = [
        # in-order, NO markers at all -> taxonomy still calls it
        {"classification": ARCH, "completeness": 100,
         "pmur_OG0001148": 0, "pmur_OG0001163": 0, "pmur_OG0001014": 0, "pmur_MurT_stray": 0},
        # out-of-order (bacterium) with the FULL block -> candidate
        {"classification": BACT, "completeness": 100,
         "pmur_OG0001148": 1, "pmur_OG0001163": 1, "pmur_OG0001014": 1, "pmur_MurT_stray": 0},
        # out-of-order with scattered non-block hits (a ligase + a stray) -> weak
        {"classification": BACT, "completeness": 100,
         "pmur_OG0001148": 1, "pmur_OG0001163": 0, "pmur_OG0001014": 0, "pmur_MurT_stray": 1},
        # complete non-PM archaeon, nothing -> clean negative
        {"classification": HALO, "completeness": 100,
         "pmur_OG0001148": 0, "pmur_OG0001163": 0, "pmur_OG0001014": 0, "pmur_MurT_stray": 0},
        # incomplete, nothing -> indeterminate, not negative
        {"classification": HALO, "completeness": 70,
         "pmur_OG0001148": 0, "pmur_OG0001163": 0, "pmur_OG0001014": 0, "pmur_MurT_stray": 0},
    ]
    g = pd.DataFrame(rows)
    g = annotate_pathway_calls(g, lig, mray, cps, block_ok=True,
                               min_markers=3, cls_col="classification")
    got = list(g["pathway_call"])
    want = ["pseudomurein_expected_by_taxonomy",
            "pseudomurein_candidate_out_of_order",
            "markers_without_block",
            "no_pathway_detected",
            "indeterminate_low_completeness"]
    check(got == want, f"real annotate_pathway_calls gave:\n  {got}\n  want {want}")
    # and the block flag is only set for the one out-of-order block carrier
    check(list(g["pm_block_present"]) == [0, 1, 0, 0, 0],
          f"pm_block_present wrong: {list(g['pm_block_present'])}")
    check(int(g["pm_expected_by_taxonomy"].iloc[0]) == 1 and
          int(g["pm_expected_by_taxonomy"].iloc[1]) == 0,
          "taxonomy flag wrong")
    print("3. real code:        annotate_pathway_calls reproduces all five "
          "scenarios (taxonomy-first, block gates out-of-order)")

    # --- the discovery case is a HYPOTHESIS, never confirmed ----------------
    c = g["pathway_call"].iloc[1]
    check("candidate" in c and "expected" not in c,
          "an out-of-order genome must be a 'candidate', never confirmed/expected")
    # and taxonomy alone never says 'confirmed' either
    check("expected_by_taxonomy" in g["pathway_call"].iloc[0],
          "an in-order genome is 'expected_by_taxonomy', flagging it as genomic")
    print("4. humility:         out-of-order PM is 'candidate'; in-order is "
          "'expected_by_taxonomy'; neither claims confirmed PM")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
