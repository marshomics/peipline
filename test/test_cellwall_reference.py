#!/usr/bin/env python3
"""The wall-chemistry table is the only place in this pipeline where a claim from
the literature turns into a value a model will condition on. It gets a test.

Everything asserted here is checked against Kandler & Koenig 1978, Arch Microbiol
118:141-152, and nothing else.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "scripts"))

from cellwall_reference import (chemistry_for, genus_is_homogeneous,  # noqa: E402
                               genus_of)

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print(f"!! {msg}")
    else:
        print(f"ok  {msg}")


def main():
    # --- what the paper actually says ---------------------------------------
    p1, p1p, src, _ = chemistry_for("s__Methanobrevibacter ruminantium")
    check((p1, p1p) == ("Thr", "Lys"),
          "M. ruminantium (strain M1): Ala completely replaced by Thr at P1")
    check(src == "literature_type_strain",
          "M. ruminantium call is flagged as type-strain evidence")

    p1, p1p, src, _ = chemistry_for("s__Methanobrevibacter smithii")
    check((p1, p1p) == ("Ala", "Lys/Orn"),
          "M. smithii (strain PS): P1 stays Ala; ~1/4 of Lys -> Orn at P1'")

    p1, p1p, _, _ = chemistry_for("s__Methanobrevibacter arboriphilus")
    check((p1, p1p) == ("Ala", "Lys"), "M. arboriphilus: Ala / Lys")

    p1, p1p, _, _ = chemistry_for("s__Methanothermobacter thermautotrophicus")
    check((p1, p1p) == ("Ala", "Lys"),
          "M. thermautotrophicus: Ala / Lys (the PeiW/PeiP host genus)")

    # --- the genus trap ------------------------------------------------------
    check(not genus_is_homogeneous("g__Methanobrevibacter"),
          "Methanobrevibacter is NOT homogeneous: it holds Ala/Lys, Ala/Orn and Thr/Lys")
    p1, p1p, src, note = chemistry_for("s__Methanobrevibacter woesei")
    check((p1, p1p, src) == ("unknown", "unknown", "genus_heterogeneous"),
          "an uncharacterised Methanobrevibacter gets 'unknown', never a genus guess")
    check("would be wrong" in note, "the refusal explains itself")

    check(genus_is_homogeneous("g__Methanobacterium"),
          "Methanobacterium is homogeneous over its characterised species")
    p1, p1p, src, _ = chemistry_for("s__Methanobacterium congolense")
    check((p1, p1p, src) == ("Ala", "Lys", "literature_genus_homogeneous"),
          "an uncharacterised Methanobacterium may inherit its genus, and is labelled so")

    # --- claims this pipeline will not act on --------------------------------
    p1, _, src, note = chemistry_for("s__Methanosphaera stadtmanae")
    check((p1, src) == ("unsupported", "disputed"),
          "Ser-for-Ala in M. stadtmanae is NOT in Kandler & Koenig 1978")
    check("1985" in note, "the refusal says why (species described after the paper)")

    p1, _, src, _ = chemistry_for("s__Methanopyrus kandleri")
    check((p1, src) == ("unsupported", "disputed"),
          "an Orn modification in M. kandleri is NOT in Kandler & Koenig 1978")

    # --- negative control ----------------------------------------------------
    p1, _, _, _ = chemistry_for("s__Methanospirillum hungatei")
    check(p1 == "no_pseudomurein",
          "M. hungatei has a protein sheath, not a sacculus; never Ala-type")

    # --- plumbing ------------------------------------------------------------
    check(genus_of("s__Methanobrevibacter smithii") == "g__Methanobrevibacter",
          "genus parsed from a GTDB species string")
    check(chemistry_for(None)[0] == "unknown", "a null species is 'unknown', not a crash")
    check(chemistry_for("s__Escherichia coli")[0] == "unknown",
          "a bacterium has no pseudomurein call at all")

    print("\nRESULT:", "PASS" if not FAIL else f"FAIL ({len(FAIL)})")
    sys.exit(0 if not FAIL else 1)


if __name__ == "__main__":
    main()
