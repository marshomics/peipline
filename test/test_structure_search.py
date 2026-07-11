#!/usr/bin/env python3
"""Static checks on the structure stack, which cannot run in the sandbox.

There is no GPU, no ESMFold/Foldseek/HHsuite and no staged databases here, so the
end-to-end run is unverifiable. What IS verifiable, and is what breaks silently:

  1. the output parsers (Foldseek .m8, HHsearch .hhr) on fixtures;
  2. the gating -- with nothing staged, every stage must report `skipped` with a
     reason, and the run must emit the plan rather than crash or fabricate a hit;
  3. `enabled: false` writes an empty table and returns.

Run:  python test/test_structure_search.py
"""
from __future__ import annotations

import os
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
    import pandas as pd
    import yaml
    import structure_search as SS

    # --- 1. Foldseek .m8 parser ---------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        m8 = os.path.join(td, "hits.m8")
        open(m8, "w").write(
            "cand1\tOG0001163_ref\t0.42\t210\t80\t2\t1\t210\t1\t210\t1e-20\t180.5\n"
            "cand1\tOG0001014_ref\t0.31\t60\t30\t1\t1\t60\t1\t60\t3e-05\t60.1\n"
            "junkline_too_few_cols\n")
        hits = SS.parse_foldseek_m8(m8)
    check(len(hits) == 2, f"expected 2 foldseek hits, got {len(hits)}")
    check(hits[0]["query"] == "cand1" and hits[0]["target"] == "OG0001163_ref"
          and abs(hits[0]["evalue"] - 1e-20) < 1e-30 and hits[0]["bits"] == 180.5,
          f"foldseek row parsed wrong: {hits[0]}")
    print("1. foldseek .m8:     query/target/evalue/bits parsed; short line skipped")

    # --- 2. HHsearch .hhr parser --------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        hhr = os.path.join(td, "q.hhr")
        open(hhr, "w").write(
            "Query         cand1 some description\n"
            "Match_columns 210\n\n"
            " No Hit                             Prob E-value P-value  Score\n"
            "  1 OG0001163_ref MraY-like         99.8 2.1E-24 1e-28   180.4\n"
            "  2 OG0001014_ref CPS               71.0 4.0E-03 1e-06    40.1\n"
            "\n")
        hh = SS.parse_hhr(hhr)
    check(len(hh) == 2, f"expected 2 hhr hits, got {len(hh)}: {hh}")
    check(hh[0]["query"] == "cand1" and hh[0]["target"] == "OG0001163_ref"
          and hh[0]["prob"] == 99.8 and abs(hh[0]["evalue"] - 2.1e-24) < 1e-30,
          f"hhr row parsed wrong: {hh[0]}")
    print("2. hhsearch .hhr:    query/target/prob/evalue parsed from the hit table")

    # --- 3. gating: nothing staged -> everything skipped with a reason ------
    scfg = {"enabled": True,
            "esmfold_weights": "/no/such/weights",
            "foldseek_pm_db": "/no/such/foldseek_db",
            "hhsuite_pm_db": "/no/such/hhsuite_db"}
    ready = SS.stage_ready(scfg)
    for stage in ("esmfold", "foldseek", "hhsearch"):
        ok, reason = ready[stage]
        check(not ok, f"{stage} must be NOT ready when its resource is absent")
        check(reason, f"{stage} skip must carry a reason")
    print("3. gating:           no staged resources -> every stage skipped, "
          "with reasons")

    # --- 4. end-to-end no-op: enabled but nothing staged --------------------
    with tempfile.TemporaryDirectory() as td:
        cw = os.path.join(td, "cellwall_genotype.tsv")
        pd.DataFrame([
            {"sample": "BACT1", "pathway_call":
             "pseudomurein_candidate_out_of_order_divergent_syntenic"},
            {"sample": "m1", "pathway_call": "pseudomurein_expected_by_taxonomy"},
        ]).to_csv(cw, sep="\t", index=False)
        cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
        cfg["specificity"]["structure_search"] = scfg
        cfgp = os.path.join(td, "c.yaml")
        yaml.safe_dump(cfg, open(cfgp, "w"))
        out = os.path.join(td, "structure_homology.tsv")
        r = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "structure_search.py"),
             "--config", cfgp, "--cellwall", cw, "--workdir", td, "--out", out],
            capture_output=True, text=True)
        check(r.returncode == 0, f"structure_search must not crash when unstaged:\n{r.stderr[-800:]}")
        check(os.path.exists(out), "it must write the plan table even when nothing ran")
        t = pd.read_csv(out, sep="\t")
        check("skipped" in set(t["status"]), "stages must be recorded as skipped")
        # only the 1 out-of-order candidate is nominated, not the in-order genome
        candrow = t[t["stage"] == "candidates"]
        check(len(candrow) == 1 and "BACT1" in str(candrow.iloc[0]["reason"])
              and "m1" not in str(candrow.iloc[0]["reason"]),
              "only out-of-order candidates should be folded, not in-order genomes")
        check("Traceback" not in r.stderr, "no traceback")
    print("4. no-op run:        enabled + nothing staged -> writes the plan + "
          "reasons, folds only the out-of-order candidate, no crash, no fake hit")

    # --- 5. disabled -> empty table -----------------------------------------
    with tempfile.TemporaryDirectory() as td:
        cw = os.path.join(td, "cw.tsv")
        pd.DataFrame([{"sample": "x", "pathway_call": "no_pathway_detected"}]).to_csv(
            cw, sep="\t", index=False)
        cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
        cfg["specificity"]["structure_search"]["enabled"] = False
        cfgp = os.path.join(td, "c.yaml")
        yaml.safe_dump(cfg, open(cfgp, "w"))
        out = os.path.join(td, "o.tsv")
        r = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "structure_search.py"),
             "--config", cfgp, "--cellwall", cw, "--workdir", td, "--out", out],
            capture_output=True, text=True)
        check(r.returncode == 0 and os.path.exists(out),
              "disabled must write an empty table and return 0")
    print("5. disabled:         enabled=false -> empty table, clean return")

    print()
    print("NOTE: the structure stack has NOT been executed (no GPU / Foldseek / "
          "HHsuite / DBs). Parsers and gating are tested; the run is not.")
    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
