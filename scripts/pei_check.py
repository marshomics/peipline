#!/usr/bin/env python3
"""Refuse to start a run whose rules would discard the proteins it is looking for.

A screen for Pei enzymes that throws away PeiR is not a screen for Pei enzymes.
This runs in milliseconds, reads only `config.yaml` and `pei_reference.py`, and
exits non-zero before 350,000 proteomes are searched.

What it asserts, in order of how badly it used to fail:

  1. Every characterised Pei is assigned to a family, and that family's profile is
     configured. PeiR is PF03412; if the C39 arm is missing, PeiR cannot be found.

  2. No family's spacing prior discards its own seeds. The C71 prior (Cys->His
     gap 35) rejects PeiR, whose gap is 72. If someone points the C39 arm at the
     C71 prior -- the obvious, wrong shortcut -- this fails loudly.

  3. `pei_class` is not applied outside C71. The four-class partition is defined
     by V252 and C265 in PeiW numbering. Those positions do not exist in a C39
     protein, and assigning classes to one would be reading a coordinate system
     onto a protein that has none.

  4. The C39 arm declares no spacing prior. One characterised sequence is not a
     prior; the spacing must be learned from the hits and reported.

  5. Each family's alignment scaffold is its own specific model, never the
     SSF54001 fold model, and never another family's model.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pei_reference import (  # noqa: E402
    CHARACTERISED, NOT_FOUND, PF12386_CROSS_DOMAIN, peir_his_candidates,
    pmbr_rule_applies, spacing_prior_admits)
from utils import load_config  # noqa: E402

PROBLEMS = []


def fail(msg):
    PROBLEMS.append(msg)
    print(f"  FAIL  {msg}", file=sys.stderr)


def ok(msg):
    print(f"  ok    {msg}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", action="store_true")
    a = ap.parse_args()

    cfg = load_config(a.config)
    fams = cfg.get("families") or {}
    profiles = cfg.get("profiles") or {}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    rows = []
    unfindable = []

    # --- 1. every seed has a configured home ---------------------------------
    print("\n[pei_check] family arms configured:", file=sys.stderr)
    for name, e in CHARACTERISED.items():
        fam = e["family"]
        if fam not in fams:
            fail(f"{name} belongs to family '{fam}', which is not in config.families. "
                 f"It cannot be found by this screen.")
            rows.append((name, fam, "no_family_arm", ""))
            continue
        want = fams[fam]["specific_profile"]
        if want not in profiles:
            fail(f"{name}: family '{fam}' names specific_profile '{want}', which is "
                 f"not in config.profiles.")
            rows.append((name, fam, "profile_missing", want))
            continue
        # PeiR carries PF03412 and nothing else; PeiW/PeiP carry PF12386.
        if e["pfam"] and want not in [p.rstrip("?") for p in e["pfam"]]:
            fail(f"{name} is annotated {e['pfam']} but family '{fam}' searches with "
                 f"'{want}'. The model will not find the protein.")
            rows.append((name, fam, "profile_mismatch", want))
            continue
        # A declared-but-disabled model is not a screen. Say so, once per protein,
        # rather than certifying an arm that will never run.
        if not profiles[want].get("enabled", True):
            unfindable.append((name, want))
            print(f"  WARN  {name} ({e['uniprot'] or 'no accession'}) needs '{want}', "
                  f"which is declared but DISABLED. It will not be found.",
                  file=sys.stderr)
            rows.append((name, fam, "profile_disabled", want))
            continue
        ok(f"{name} ({e['uniprot'] or 'no accession'}) -> family {fam} -> {want}")
        rows.append((name, fam, "found_by", want))

    # --- 2. no family's spacing prior discards its own seeds ------------------
    print("\n[pei_check] spacing priors vs their own seeds:", file=sys.stderr)
    for fam, fc in fams.items():
        g12, tol = fc.get("expected_gap_1_2"), fc.get("gap_tolerance")
        for name in fc.get("seeds", []):
            if name not in CHARACTERISED:
                fail(f"family '{fam}' seeds an unknown protein '{name}'")
                continue
            if g12 is None or tol is None:
                ok(f"{fam}/{name}: no spacing prior declared; it will be learned")
                rows.append((name, fam, "prior_learned", ""))
                continue
            admitted, why = spacing_prior_admits(name, int(g12), int(tol))
            if admitted is False:
                fail(f"family '{fam}' prior (Cys->His {g12}+/-{tol}) DISCARDS its own "
                     f"seed {name}: {why}")
            elif admitted is None:
                ok(f"{fam}/{name}: {why}")
            else:
                ok(f"{fam}/{name}: retained ({why})")
            rows.append((name, fam, f"prior_{admitted}", why))

    # The specific trap: the C71 prior applied to the C39 arm.
    c71 = fams.get("c71", {})
    if c71.get("expected_gap_1_2") is not None:
        adm, why = spacing_prior_admits("PeiR", int(c71["expected_gap_1_2"]),
                                        int(c71.get("gap_tolerance") or 8))
        if adm:
            fail("the C71 spacing prior admits PeiR. That should be impossible "
                 "(no PeiR histidine lies within 35+/-8 of C90); the reference "
                 "data are wrong.")
        else:
            ok(f"the C71 prior correctly rejects PeiR -- which is exactly why the "
               f"C39 arm must not borrow it. Gaps from C90: {peir_his_candidates()}")

    # --- 3. pei_class does not escape C71 ------------------------------------
    print("\n[pei_check] pei_class jurisdiction:", file=sys.stderr)
    for fam, fc in fams.items():
        applies = bool(fc.get("pei_class_applies"))
        if fam != "c71" and applies:
            fail(f"family '{fam}' sets pei_class_applies. The four-class partition "
                 f"is defined by V252 and C265 in PeiW (C71) numbering; those "
                 f"positions do not exist in a {fam.upper()} protein.")
        elif fam == "c71" and not applies:
            fail("family 'c71' has pei_class_applies false; the partition was "
                 "derived on C71 and should be used there.")
        else:
            ok(f"{fam}: pei_class_applies={applies}")

    # --- 4. C39 declares no borrowed prior -----------------------------------
    c39 = fams.get("c39", {})
    if c39:
        if c39.get("expected_gap_1_2") is not None:
            fail("family 'c39' declares expected_gap_1_2. Only one C39 Pei has an "
                 "assigned catalytic residue (PeiR, C90, and its His was never "
                 "assigned). One sequence is not a prior. Set it null and learn it.")
        else:
            ok("c39 declares no spacing prior; it will be learned and reported")

    # --- 5. scaffolds are specific models ------------------------------------
    print("\n[pei_check] alignment scaffolds:", file=sys.stderr)
    for fam, fc in fams.items():
        scaf = fc.get("align_profile")
        role = (profiles.get(scaf) or {}).get("role")
        if role != "specific":
            fail(f"family '{fam}' aligns to '{scaf}', whose role is '{role}'. The "
                 f"scaffold must be a specific model: match-state coordinates from "
                 f"a fold model are not comparable across families.")
        elif (profiles[scaf].get("family") or fam) != fam:
            fail(f"family '{fam}' aligns to '{scaf}', which belongs to family "
                 f"'{profiles[scaf].get('family')}'.")
        else:
            ok(f"{fam} -> {scaf} (role={role})")

    # --- reporting -----------------------------------------------------------
    print("\n[pei_check] PMBR rule jurisdiction:", file=sys.stderr)
    for name in CHARACTERISED:
        j = pmbr_rule_applies(name)
        n = CHARACTERISED[name]["pmbr_motifs"]
        if j is False:
            print(f"  note  {name} has {n} PMB motifs and lyses cells. "
                  f"`pmbr_binding_competent == 0` is NOT a prediction that it "
                  f"cannot lyse.", file=sys.stderr)
        rows.append((name, CHARACTERISED[name]["family"], f"pmbr_rule_applies_{j}",
                     str(n)))

    print("\n[pei_check] PF12386 outside the archaea (why the screen is "
          "cross-domain):", file=sys.stderr)
    for acc, (org, lineage, ln) in PF12386_CROSS_DOMAIN.items():
        print(f"  note  {acc} {org} ({lineage}) {ln} aa", file=sys.stderr)

    for n, why in NOT_FOUND.items():
        print(f"\n[pei_check] note  {n}: {why}", file=sys.stderr)

    with open(a.out, "w") as fh:
        fh.write("protein\tfamily\tcheck\tdetail\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    if unfindable:
        names = ", ".join(f"{n} (needs {p})" for n, p in unfindable)
        print(f"\n[pei_check] ===================================================\n"
              f"[pei_check] THIS RUN CANNOT FIND: {names}\n"
              f"[pei_check] The model is declared in config.profiles but disabled, "
              f"because searching it before the per-family alignment path exists "
              f"would pool C39 hits with C71 hits and score them against C71 triad "
              f"columns. No PeiR histidine lies at the C71 gap of 35 (its "
              f"downstream His sit at gaps 14/72/91/102 from C90).\n"
              f"[pei_check] Found and mangled is worse than not found. Enable it "
              f"only once triad_detect_filter.py aligns each family to its own "
              f"scaffold.\n"
              f"[pei_check] ===================================================",
              file=sys.stderr)

    print(f"\n[pei_check] {len(PROBLEMS)} problem(s), {len(unfindable)} unfindable "
          f"protein(s); wrote {a.out}", file=sys.stderr)
    if PROBLEMS and a.strict:
        sys.exit("[pei_check] the configured rules would discard proteins this "
                 "screen exists to find. Fix config.yaml before running.")


if __name__ == "__main__":
    main()
