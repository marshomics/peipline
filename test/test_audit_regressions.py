#!/usr/bin/env python3
"""Regressions for defects the July audit found. Each locks a fix that was, before
the audit, either a guaranteed crash on real data or a silently wrong number.

1. triad low-coverage gating: a fragment that carries C/H/D at the three triad
   columns but covers < min_match_coverage of the profile must be `low_coverage`,
   NOT `triad_positive`, and must not trip the "coordinate inconsistent" assert
   that used to abort the whole run.
2. module_trees tanglegram title: the `p` name was undefined and crashed every
   real run at the plotting step. The source must reference p_rf, not a bare p.
3. cellwall_genotype bacteria check: referenced a column `pmur_pathway` that is
   never created (real name `pmur_count_pathway`) -> KeyError at the end of every
   run. The wrong name must be gone.
4. convergence statistic names: the old columns were mislabels (a parsimony change
   count reported as "origins"; a home-rolled index reported as Fritz-Purvis D).
   The output must use the corrected names parsimony_changes / clustering_index.

Run:  python test/test_audit_regressions.py
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")

FAILURES = []


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def build_sto(path, L=60, ci=5, hi=25, di=45):
    """One full-coverage specific seq defines the columns; one FRAGMENT carries
    C/H/D at the same three columns but is otherwise all gaps (low coverage)."""
    def row(fill, triad, gaps=()):
        s = [fill] * L
        for p, aa in triad.items():
            s[p] = aa
        for g in gaps:
            s[g] = "-"
        return "".join(s)
    seqs = {}
    for n in range(4):
        seqs[f"spec{n}"] = row("V", {ci: "C", hi: "H", di: "D"})
    # fragment: C/H/D at the triad columns, everything else a gap -> coverage 3/L
    frag = ["-"] * L
    frag[ci], frag[hi], frag[di] = "C", "H", "D"
    seqs["frag_lowcov"] = "".join(frag)
    with open(path, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n")
        w = max(len(k) for k in seqs) + 2
        for k, v in seqs.items():
            fh.write(f"{k:<{w}}{v}\n")
        fh.write(f"{'#=GC RF':<{w}}{'x' * L}\n//\n")
    return list(seqs), (ci, hi, di)


def main() -> None:
    import pandas as pd
    import yaml

    # --- 1. triad low-coverage gating ---------------------------------------
    with tempfile.TemporaryDirectory() as td:
        sto = os.path.join(td, "hits.sto")
        names, (ci, hi, di) = build_sto(sto)
        pd.DataFrame({"seq_id": names,
                      "evidence": ["specific"] * (len(names) - 1) + ["ssf_only"],
                      "sample": "S", "protein_id": names, "faa": "x"}).to_csv(
            os.path.join(td, "idmap.tsv"), sep="\t", index=False)
        pd.DataFrame({"seq_id": names, "weight": 1.0, "cluster": names,
                      "cluster_size": 1}).to_csv(
            os.path.join(td, "w.tsv"), sep="\t", index=False)
        open(os.path.join(td, "c.txt"), "w").write("x\n")
        cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
        cfg["triad"].update(expected_gap_1_2=hi - ci, expected_gap_2_3=di - hi,
                            spacing_tolerance=5, min_occupancy=0.3,
                            min_residue_freq=0.3, min_match_coverage=0.5,
                            alt_third_residue=None)
        cfgp = os.path.join(td, "c.yaml")
        yaml.safe_dump(cfg, open(cfgp, "w"))
        out = {k: os.path.join(td, k) for k in
               ("cands", "chosen", "keep", "afa", "colstats", "tiers", "outcomes")}
        r = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "triad_detect_filter.py"),
             "--sto", sto, "--config", cfgp, "--combined", os.path.join(td, "c.txt"),
             "--idmap", os.path.join(td, "idmap.tsv"), "--weights", os.path.join(td, "w.tsv"),
             "--out-candidates", out["cands"], "--out-chosen", out["chosen"],
             "--out-keep", out["keep"], "--out-afa", out["afa"],
             "--out-colstats", out["colstats"], "--out-tiers", out["tiers"],
             "--out-outcomes", out["outcomes"]],
            capture_output=True, text=True)
        check(r.returncode == 0,
              f"triad must not crash on a low-coverage triad-bearing fragment:\n{r.stderr[-800:]}")
        check("coordinate system is inconsistent" not in r.stderr,
              "the low-coverage fragment tripped the run-ending assertion")
        if os.path.exists(out["outcomes"]):
            o = pd.read_csv(out["outcomes"], sep="\t").set_index("seq_id")["outcome"]
            check(o.get("frag_lowcov") == "low_coverage",
                  f"fragment must be low_coverage, got {o.get('frag_lowcov')}")
            keep = open(out["keep"]).read().split()
            check("frag_lowcov" not in keep,
                  "a low-coverage fragment must NOT enter the triad-positive set / c71.faa")
    print("1. triad low-cov:    fragment with C/H/D at the columns -> low_coverage, "
          "kept out of c71.faa, no crash")

    # --- 2. module_trees uses p_rf, not undefined p -------------------------
    mt = open(os.path.join(SCRIPTS, "module_trees.py")).read()
    check(not re.search(r"\(p\s*=\s*\{p:\.4g\}\)", mt),
          "module_trees still formats an undefined `p` in the tanglegram title")
    check("tip-shuffle p = {p_rf" in mt or "{p_rf:" in mt,
          "module_trees tanglegram must report p_rf (the tip-shuffle p)")
    print("2. module_trees:     tanglegram title references p_rf, not undefined p")

    # --- 3. cellwall_genotype no phantom column -----------------------------
    cw = open(os.path.join(SCRIPTS, "cellwall_genotype.py")).read()
    check('g["pmur_pathway"]' not in cw,
          "cellwall_genotype still references the non-existent column pmur_pathway")
    check('g["pmur_count_pathway"]' in cw,
          "cellwall_genotype must use pmur_count_pathway")
    print("3. cellwall:         bacteria check uses pmur_count_pathway (no KeyError)")

    # --- 4. convergence renamed statistics ----------------------------------
    cv = open(os.path.join(SCRIPTS, "convergence.py")).read()
    old_origins = "observed_" + "origins"      # split so this test's own sed can't rewrite it
    old_d = "fritz_purvis_" + "D"
    check("parsimony_changes" in cv and old_origins not in cv,
          "convergence must rename the origins column to parsimony_changes")
    check("clustering_index" in cv and old_d not in cv,
          "convergence must rename the D column to clustering_index (not caper's D)")
    print("4. convergence:      parsimony_changes / clustering_index (honest names)")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
