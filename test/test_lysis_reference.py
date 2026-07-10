#!/usr/bin/env python3
"""Test the one part of this pipeline that could be proven wrong by an experiment.

Subedi et al. 2015 measured which methanogens PeiW and PeiP lyse. The pipeline
predicts lysis from host wall chemistry (Kandler & Koenig 1978) via a rule
derived from the pNA substrate series. Neither table saw the other.

What is checked:

  1. the rule reproduces every falsifiable row of the panel;
  2. the rule is FALSIFIABLE -- corrupt the chemistry table and it must disagree.
     A rule that agrees with everything has not been tested, it has been decorated;
  3. Ser is accepted at P1 and Thr is rejected in the pNA series;
  4. the rule refuses to predict where it has no chemistry, rather than guessing;
  5. `lysis_check.py --strict` exits non-zero when the rule breaks;
  6. no file in the repo still claims that Glu-gamma-Thr-pNA was never made. It
     was. It was assayed. Neither enzyme cleaves it. An earlier version of
     assay_panel.py asserted the opposite, and this test exists so that the claim
     cannot come back.

Run:  python test/test_lysis_reference.py
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

FAILURES = []


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def main() -> None:
    import lysis_reference as LR
    import cellwall_reference as CW

    # --- 1. the rule reproduces the panel -----------------------------------
    df = LR.check_against_reference()
    n, k, nt, nt_k, bad = LR.summarise(df)
    check(nt >= 4, f"only {nt} falsifiable rows; the panel lost its content")
    check(nt_k == nt,
          f"the host-chemistry rule contradicts the measured panel on "
          f"{nt - nt_k} of {nt} falsifiable rows:\n{bad.to_string(index=False)}")
    print(f"1. panel:            {nt_k}/{nt} falsifiable predictions agree "
          f"({k}/{n} including the trivial no-pseudomurein negatives)")

    # The row that carries the hypothesis. Match on the full label: "M1" is a
    # substring of "CM1", the Methanosarcina negative control, and matching it
    # would quietly fold a no-pseudomurein control into the Thr test.
    m1 = df[df["strain"] == "Methanobrevibacter ruminantium M1"]
    check(len(m1) == 2 and set(m1["p1_residue"]) == {"Thr"},
          "M. ruminantium M1 must carry Thr at P1 (Kandler & Koenig, strain M1)")
    check(set(m1["observed"]) == {"not_lysed"},
          "M1 must be recorded as resistant to both enzymes")

    # --- 2. falsifiability ---------------------------------------------------
    # If M1's wall were Ala, the rule would predict lysis, and the panel says it
    # is not lysed. The test must therefore FAIL. If it passes, the rule is not
    # reading the chemistry at all.
    saved = dict(CW.REFERENCE["s__Methanobrevibacter ruminantium"])
    try:
        CW.REFERENCE["s__Methanobrevibacter ruminantium"]["p1"] = "Ala"
        df_bad = LR.check_against_reference()
        _, _, _, nt_k_bad, bad_bad = LR.summarise(df_bad)
        check(nt_k_bad < nt,
              "corrupting M1's P1 residue to Ala did not produce a disagreement. "
              "The rule is not actually consulting the chemistry table, so the "
              "4/4 above means nothing")
        check(any("M1" in s for s in bad_bad["strain"]),
              "the disagreement should be on the M1 rows")
        print(f"2. falsifiable:      corrupting M1 Thr->Ala breaks "
              f"{nt - nt_k_bad} prediction(s), as it must")
    finally:
        CW.REFERENCE["s__Methanobrevibacter ruminantium"] = saved

    # --- 3. the substrate series --------------------------------------------
    check(LR.SUBSTRATES["EgammaT-pNA"]["peiW"] == "inactive",
          "Glu-gamma-Thr-pNA is NOT cleaved by PeiW")
    check(LR.SUBSTRATES["EgammaT-pNA"]["peiP"] == "inactive",
          "Glu-gamma-Thr-pNA is NOT cleaved by PeiP")
    check(LR.SUBSTRATES["EgammaS-pNA"]["peiW"] == "active",
          "Glu-gamma-Ser-pNA IS cleaved: Ser substitutes for Ala at P1")
    check(LR.SUBSTRATES["DbetaA-pNA"]["peiW"] == "poor",
          "Asp does not substitute for Glu at P2")
    check(LR.P1_ACCEPTED_PNA == {"Ala", "Ser"}, "P1 accepted set is wrong")
    check(LR.P1_REJECTED_PNA == {"Thr"}, "P1 rejected set is wrong")
    check(CW.P1_CLEAVED_BY_PEIW_PEIP == LR.P1_ACCEPTED_PNA,
          "cellwall_reference and lysis_reference disagree about which P1 "
          "residues are cleaved")
    # every substrate has pNA where the acceptor belongs -- that is why the
    # series cannot probe P1'
    check(all(v["acceptor"] == "pNA" for v in LR.SUBSTRATES.values()),
          "a substrate in the pNA series claims a real acyl acceptor")
    print("3. substrates:       Ala/Ser cleaved, Thr not, Asp not at P2, "
          "S1' empty throughout")

    # --- 4. metal ------------------------------------------------------------
    check(LR.METAL["peiP_order"] == ["Ca"],
          "PeiP is rescued by Ca alone; the rest give <15%")
    check(len(LR.METAL["peiW_order"]) == 5,
          "PeiW is rescued by Ca, Mn, Mg, Ba and Ni")
    print(f"4. metal:            PeiW {'>'.join(LR.METAL['peiW_order'])}; "
          f"PeiP {LR.METAL['peiP_order'][0]} only")

    # --- 5. refusal to guess -------------------------------------------------
    # SM9 is lysed by PeiW and not by PeiP. No wall chemistry is published. The
    # rule must say "unknown", not split the difference.
    sm9 = df[df["strain"].str.contains("SM9")]
    check(len(sm9) == 2 and sm9["agrees"].isna().all(),
          "SM9 has no published chemistry; the rule must decline to predict it")
    check(set(sm9["observed"]) == {"lysed", "not_lysed"},
          "SM9 is the only differential row in the panel and must stay so")
    pred, _ = LR.predict_lysis("unknown", "pseudomurein")
    check(pred == "unknown", "predict_lysis must return 'unknown' without a P1")
    print("5. refusal:          SM9 (PeiW lyses, PeiP does not) left unpredicted")

    # --- 6. lysis_check.py --strict ------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        r = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "lysis_check.py"),
             "--out", os.path.join(td, "check.tsv"),
             "--out-substrates", os.path.join(td, "sub.tsv"), "--strict"],
            capture_output=True, text=True)
        check(r.returncode == 0,
              f"lysis_check --strict should pass on the real tables:\n{r.stderr}")
        check(os.path.exists(os.path.join(td, "check.tsv")), "no check table written")

        # and it must fail loudly when the rule is broken
        broken = os.path.join(td, "broken.py")
        src = open(os.path.join(SCRIPTS, "cellwall_reference.py")).read()
        src = src.replace('"strain": "M1", "p1": "Thr"', '"strain": "M1", "p1": "Ala"')
        os.makedirs(os.path.join(td, "scripts"), exist_ok=True)
        for f in ("lysis_reference.py", "lysis_check.py"):
            open(os.path.join(td, "scripts", f), "w").write(
                open(os.path.join(SCRIPTS, f)).read())
        open(os.path.join(td, "scripts", "cellwall_reference.py"), "w").write(src)
        r2 = subprocess.run(
            [sys.executable, os.path.join(td, "scripts", "lysis_check.py"),
             "--out", os.path.join(td, "c2.tsv"),
             "--out-substrates", os.path.join(td, "s2.tsv"), "--strict"],
            capture_output=True, text=True)
        check(r2.returncode != 0,
              "lysis_check --strict must exit non-zero when the chemistry table "
              "contradicts the measured lysis panel")
        check("DISAGREEMENT" in r2.stderr.upper(),
              "the failure must name the disagreement")
        print(f"6. strict mode:      exit {r2.returncode} on a corrupted "
              f"chemistry table")

    # --- 7. the retracted claim must not return ------------------------------
    bad_patterns = [
        r"no\s+Thr-pNA\s+exists",
        r"no\s+Thr-pNA\s+was\s+(ever\s+)?made",
        r"Thr-pNA\s+was\s+never\s+made",
    ]
    offenders = []
    for dirpath, _, files in os.walk(ROOT):
        if any(p in dirpath for p in (".git", "__pycache__")):
            continue
        for f in files:
            if not f.endswith((".py", ".yaml", ".md", ".smk")) and f != "Snakefile":
                continue
            p = os.path.join(dirpath, f)
            if os.path.abspath(p) == os.path.abspath(__file__):
                continue
            text = open(p, errors="ignore").read()
            for pat in bad_patterns:
                if re.search(pat, text, re.I):
                    offenders.append(f"{os.path.relpath(p, ROOT)}: /{pat}/")
    check(not offenders,
          "a file still claims Glu-gamma-Thr-pNA was never synthesised. It was "
          "(JPT Peptide Technologies), it was assayed, and neither enzyme "
          "cleaves it:\n  " + "\n  ".join(offenders))
    print("7. retracted claim:  no file asserts Thr-pNA was never made")

    # --- 8. assay_panel routes Thr to an isopeptide --------------------------
    import assay_panel as AP
    check(AP.PNA_VALID_P1 == {"Ala", "Ser"},
          "the pNA assay is valid only for Ala and Ser at P1")
    check("isopeptide" in AP.P1_PREDICTION["Thr"],
          "a Thr panel member must be sent to the isopeptide assay")
    check("Do NOT use" in AP.P1_PREDICTION["Thr"],
          "the Thr prediction must warn against EgammaT-pNA explicitly")
    check("pNA" in AP.P1_PREDICTION["Ser"] and "isopeptide" not in AP.P1_PREDICTION["Ser"],
          "Ser is cleaved in the pNA series; the cheap assay is valid")
    print("8. assay routing:    Thr -> isopeptide, Ala/Ser -> pNA")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
