#!/usr/bin/env python3
"""Test that the screen cannot be configured to discard the proteins it seeks.

The premise this file defends: PeiR (UniProt D3DZZ6) is a pseudomurein
endoisopeptidase that lyses Methanobrevibacter ruminantium M1, carries PF03412
Peptidase_C39 and nothing else, has no PF09373 binding repeat, and has a catalytic
cysteine at C90 whose nearest downstream histidine is 72 residues away. Every
assumption the C71 arm makes is false for it.

Checks:

  1. the PeiR sequence really does say what pei_reference claims: C at 90, no His
     at 90+35, cysteines at the annotated disulfide;
  2. the C71 spacing prior DISCARDS PeiR -- if it ever admits it, the reference
     data have been corrupted and every "PeiR is different" claim is hollow;
  3. `pei_check --strict` passes on the shipped config;
  4. it FAILS if the C39 arm borrows the C71 spacing prior;
  5. it FAILS if the C39 arm is deleted (PeiR then has no home);
  6. it FAILS if pei_class is switched on for C39;
  7. it FAILS if a family aligns to the SSF54001 fold model;
  8. the PMBR rule declares no jurisdiction over PeiR, and `predict_binding(0)`
     does not claim the protein cannot lyse;
  9. nothing in the repo asserts PeiS exists.

Run:  python test/test_pei_reference.py
"""
from __future__ import annotations

import copy
import os
import re
import subprocess
import sys
import tempfile

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

FAILURES = []


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def run_check(cfg_obj, td, name):
    p = os.path.join(td, f"{name}.yaml")
    yaml.safe_dump(cfg_obj, open(p, "w"))
    return subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "pei_check.py"),
         "--config", p, "--out", os.path.join(td, f"{name}.tsv"), "--strict"],
        capture_output=True, text=True)


def main() -> None:
    import pei_reference as PR
    import pmbr_reference as MR

    # --- 1. the sequence says what the module claims -------------------------
    seq = PR.PEIR_SEQUENCE
    e = PR.CHARACTERISED["PeiR"]
    check(len(seq) == e["length"] == 228, f"PeiR length {len(seq)} != 228")
    check(seq[e["catalytic"]["cys"] - 1] == "C",
          "residue 90 of PeiR is not a cysteine; the C90A mutant says it must be")
    d1, d2 = e["disulfide"]
    check(seq[d1 - 1] == "C" and seq[d2 - 1] == "C",
          f"the annotated disulfide {d1}-{d2} does not connect two cysteines")
    check(e["catalytic"]["his"] is None and e["catalytic"]["asp"] is None,
          "PeiR's His and Asp were never assigned experimentally; recording a "
          "guess as data is how a consensus becomes evidence")
    gaps = PR.peir_his_candidates()
    check(min(gaps.values()) == 14 and 72 in gaps.values(),
          f"unexpected His gaps from C90: {gaps}")
    check(seq[90 + 35 - 1] == "P",
          "C90+35, where a C71 histidine would sit, must be the proline at 125")
    print(f"1. PeiR sequence:    C90 confirmed, disulfide C{d1}-C{d2}, "
          f"His gaps {gaps}, residue 125 = {seq[124]}")

    # --- 2. the C71 prior must reject PeiR -----------------------------------
    admitted, why = PR.spacing_prior_admits("PeiR", 35, 8)
    check(admitted is False,
          f"the C71 prior must DISCARD PeiR, got admitted={admitted}: {why}")
    for n in ("PeiW", "PeiP"):
        a, _ = PR.spacing_prior_admits(n, 35, 8)
        check(a is True, f"the C71 prior must retain {n}")
    # widen the tolerance far enough and it should admit PeiR: the check is a
    # statement about the prior, not an artefact of the function
    wide, _ = PR.spacing_prior_admits("PeiR", 72, 2)
    check(wide is True,
          "a prior centred on 72 must admit PeiR, or spacing_prior_admits is "
          "simply returning False for everything")
    print(f"2. C71 prior:        rejects PeiR ({why[:58]}...), keeps PeiW/PeiP")

    # --- 3-7. pei_check on the shipped config and on four sabotaged ones ------
    base = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    with tempfile.TemporaryDirectory() as td:
        r = run_check(base, td, "shipped")
        check(r.returncode == 0,
              f"pei_check --strict must pass on the shipped config:\n{r.stderr[-1200:]}")
        check("PeiR" in r.stderr and "PF03412" in r.stderr,
              "the check must show PeiR routed to PF03412")

        # PF03412 ships disabled. The check must SAY SO rather than certify a
        # screen that cannot find PeiR. Silence here would be the worst outcome:
        # a green run that never looked for half the family.
        c39_enabled = base["profiles"]["PF03412"].get("enabled", True)
        if not c39_enabled:
            check("CANNOT FIND" in r.stderr and "PeiR" in r.stderr,
                  "with PF03412 disabled, pei_check must announce that PeiR "
                  "cannot be found by this run")
            print("3. shipped config:   pei_check passes AND warns that the C39 "
                  "arm is disabled, so PeiR is unfindable")
        else:
            # If someone enables it, the execution path must be split first.
            tdf = open(os.path.join(SCRIPTS, "triad_detect_filter.py")).read()
            check("families" in tdf,
                  "PF03412 is enabled but triad_detect_filter.py still reads a "
                  "single align_profile. C39 hits would be aligned to PF12386 and "
                  "scored against C71 triad columns. Split the path first.")
            print("3. shipped config:   pei_check --strict passes, C39 arm enabled")

        # a disabled profile must never reach the searcher
        sb = open(os.path.join(SCRIPTS, "search_batch.sh")).read()
        check('d.get("enabled", True)' in sb,
              "search_batch.sh must skip disabled profiles")
        snake = open(os.path.join(ROOT, "Snakefile")).read()
        check('d.get("enabled", True)' in snake,
              "the Snakefile must exclude disabled profiles from PROFILE_NAMES")

        # 4. C39 borrows the C71 prior -- the obvious, wrong shortcut
        bad = copy.deepcopy(base)
        bad["families"]["c39"]["expected_gap_1_2"] = 35
        bad["families"]["c39"]["gap_tolerance"] = 8
        r = run_check(bad, td, "borrowed_prior")
        check(r.returncode != 0, "pei_check must fail when C39 borrows the C71 prior")
        check("DISCARDS its own seed PeiR" in r.stderr,
              f"the failure must name PeiR:\n{r.stderr[-800:]}")
        print("4. borrowed prior:   rejected, and it names PeiR")

        # 5. the C39 arm is deleted
        bad = copy.deepcopy(base)
        del bad["families"]["c39"]
        r = run_check(bad, td, "no_c39")
        check(r.returncode != 0, "pei_check must fail when the C39 arm is missing")
        check("cannot be found by this screen" in r.stderr,
              f"the failure must say PeiR is unfindable:\n{r.stderr[-800:]}")
        print("5. missing C39 arm:  rejected; PeiR would be unfindable")

        # 6. pei_class escapes C71
        bad = copy.deepcopy(base)
        bad["families"]["c39"]["pei_class_applies"] = True
        r = run_check(bad, td, "class_escape")
        check(r.returncode != 0, "pei_check must fail when pei_class is applied to C39")
        check("V252" in r.stderr, "the failure must name the C71-numbered positions")
        print("6. pei_class escape: rejected (V252/C265 are PeiW numbering)")

        # 7. a family aligns to the fold model
        bad = copy.deepcopy(base)
        bad["families"]["c71"]["align_profile"] = "SSF54001"
        r = run_check(bad, td, "fold_scaffold")
        check(r.returncode != 0,
              "pei_check must fail when a family aligns to the SSF54001 fold model")
        check("role" in r.stderr and "sensitivity" in r.stderr,
              "the failure must name the role of the offending model")
        print("7. fold scaffold:    rejected (match states must come from a "
              "specific model)")

    # --- 8. the PMBR rule knows its own limits -------------------------------
    check(MR.rule_has_jurisdiction(0) is False, "no PMB array -> no jurisdiction")
    check(MR.rule_has_jurisdiction(2) is True, "2 motifs -> the rule applies")
    call, can_bind, why = MR.predict_binding(0)
    check(call == "no_pmbr_module", f"got {call}")
    check("PeiR" in why and "not a prediction" in why,
          "predict_binding(0) must say a zero is not a prediction that the "
          "protein cannot lyse, and cite PeiR")
    check(PR.pmbr_rule_applies("PeiR") is False,
          "the PMBR rule must declare no jurisdiction over PeiR")
    check(PR.CHARACTERISED["PeiR"]["lyses"] is True,
          "PeiR lyses cells; that is the whole point of the exclusion")
    # and the panel must not collapse the three states into two
    panel_src = open(os.path.join(SCRIPTS, "assay_panel.py")).read()
    check("np.select" in panel_src and "pmbr_rule_applies" in panel_src,
          "assay_panel must branch three ways: no array / cannot dock / docks")
    print("8. PMBR jurisdiction: zero motifs != cannot lyse (PeiR is the "
          "counterexample)")

    # --- 9. PeiS is not invented ---------------------------------------------
    check("PeiS" in PR.NOT_FOUND, "PeiS must be recorded as not found")
    check("PeiS" not in PR.CHARACTERISED, "PeiS must not appear as characterised")
    offenders = []
    for dirpath, _, files in os.walk(ROOT):
        if any(p in dirpath for p in (".git", "__pycache__")):
            continue
        for f in files:
            if not f.endswith((".py", ".yaml", ".md")) and f != "Snakefile":
                continue
            p = os.path.join(dirpath, f)
            if os.path.basename(p) in ("pei_reference.py", "pei_check.py",
                                       os.path.basename(__file__)):
                continue
            txt = open(p, errors="ignore").read()
            if re.search(r"\bPeiS\b", txt):
                offenders.append(os.path.relpath(p, ROOT))
    check(not offenders,
          f"PeiS does not exist in UniProtKB or the literature; these files "
          f"reference it as if it does: {offenders}")
    print("9. PeiS:             recorded as not found, referenced nowhere else")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
