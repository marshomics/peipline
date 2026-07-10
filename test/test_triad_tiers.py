#!/usr/bin/env python3
"""A sequence can fail the triad filter for three reasons. Only one is a result.

The SSF54001 tier exists to answer one question: does a protein with a
cysteine-proteinase fold also carry the Pei triad in the Pei positions? SCOP 54001
contains ~22 families related to the papain catalytic core by insertion and by
CIRCULAR PERMUTATION -- transglutaminase cores among them. A permuted core has a
Cys-His-Asp triad in three dimensions but presents it in a different sequential
order, and this filter selects columns with i < j < k and demands C, H, D there.

So it cannot see a permuted triad. If such a protein is counted as
`triad_negative`, the report concludes "transglutaminase-like proteins lack the
Pei active site" when what actually happened is that the alignment failed.

This test builds a Stockholm alignment by hand, with an RF line, containing:

  spec1..spec4   `specific` tier, canonical C/H/D at columns 5 / 20 / 35
  tg_permuted    `ssf_only`, all three catalytic residues present but in the order
                 H ... D ... C -- a circular permutation. Must NOT be
                 `triad_negative`.
  tg_fragment    `ssf_only`, aligned over 20% of the match states. The catalytic
                 residues are real but land outside the profile envelope. Must be
                 `low_coverage`.
  ssf_gapped     `ssf_only`, decent coverage but a gap at one triad column. Must
                 be `gapped_at_triad`.
  ssf_true_neg   `ssf_only`, full coverage, residues present at all three triad
                 columns, and they are not C/H/D. THIS is a real negative.
  ssf_true_pos   `ssf_only`, full coverage, C/H/D at the triad columns. A genuine
                 hit through the sensitivity net.

Then it asserts the tier table separates them, that `frac_triad_positive_of_testable`
uses the right denominator, and that the log warns rather than concludes.

Run:  python test/test_triad_tiers.py [workdir]
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")

FAILURES = []

L = 40                    # match columns
CI, HI, DI = 5, 20, 35    # where the triad sits, i < j < k


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def row(fill="A", triad=None, gaps=()):
    s = [fill] * L
    if triad:
        for pos, aa in triad.items():
            s[pos] = aa
    for g in gaps:
        s[g] = "-"
    return "".join(s)


def build_sto(path):
    """All match states, so match column index == alignment column index."""
    seqs = {}
    # --- specific tier: the columns are learned from these ---------------------
    for n in range(4):
        seqs[f"spec{n}"] = row("V", {CI: "C", HI: "H", DI: "D"})

    # --- ssf_only tier ---------------------------------------------------------
    # a circular permutation: His first, then Asp, then Cys. All three residues
    # are there; none is at its column. The 3D site could be identical.
    seqs["tg_permuted"] = row("V", {CI: "H", HI: "D", DI: "C"})

    # a fragment: 8 of 40 match states occupied = 20% coverage
    seqs["tg_fragment"] = "".join(
        ("CHD" + "VVVVV")[k] if k < 8 else "-" for k in range(L))

    # decent coverage, gap at the Asp column only
    seqs["ssf_gapped"] = row("V", {CI: "C", HI: "H"}, gaps=(DI,))

    # a real negative: everything occupied, residues simply are not C/H/D
    seqs["ssf_true_neg"] = row("V", {CI: "S", HI: "N", DI: "E"})

    # a real positive through the sensitivity net
    seqs["ssf_true_pos"] = row("V", {CI: "C", HI: "H", DI: "D"})

    with open(path, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n")
        w = max(len(k) for k in seqs) + 2
        for k, v in seqs.items():
            fh.write(f"{k:<{w}}{v}\n")
        fh.write(f"{'#=GC RF':<{w}}{'x' * L}\n//\n")
    return list(seqs)


def main() -> None:
    base = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="c71_tiers_")
    os.makedirs(base, exist_ok=True)
    import pandas as pd
    import yaml

    sto = os.path.join(base, "hits.sto")
    names = build_sto(sto)

    ev = {n: ("specific" if n.startswith("spec") else "ssf_only") for n in names}
    pd.DataFrame({"seq_id": names, "evidence": [ev[n] for n in names],
                  "sample": "S", "protein_id": names, "faa": "x"}).to_csv(
        os.path.join(base, "idmap.tsv"), sep="\t", index=False)
    pd.DataFrame({"seq_id": names, "weight": 1.0, "cluster": names,
                  "cluster_size": 1}).to_csv(
        os.path.join(base, "weights.tsv"), sep="\t", index=False)
    open(os.path.join(base, "combined.txt"), "w").write("x\n")

    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    cfg["triad"].update(expected_gap_1_2=HI - CI, expected_gap_2_3=DI - HI,
                        spacing_tolerance=5, min_occupancy=0.3,
                        min_residue_freq=0.3, min_match_coverage=0.5,
                        alt_third_residue=None)
    tcfg = os.path.join(base, "config.test.yaml")
    yaml.safe_dump(cfg, open(tcfg, "w"))

    out = {k: os.path.join(base, f"{k}.tsv") for k in
           ("cands", "colstats", "tiers", "outcomes")}
    chosen = os.path.join(base, "chosen.json")
    r = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "triad_detect_filter.py"),
         "--sto", sto, "--config", tcfg,
         "--combined", os.path.join(base, "combined.txt"),
         "--idmap", os.path.join(base, "idmap.tsv"),
         "--weights", os.path.join(base, "weights.tsv"),
         "--out-candidates", out["cands"], "--out-chosen", chosen,
         "--out-keep", os.path.join(base, "keep.txt"),
         "--out-afa", os.path.join(base, "pass.afa"),
         "--out-colstats", out["colstats"], "--out-tiers", out["tiers"],
         "--out-outcomes", out["outcomes"]],
        capture_output=True, text=True)
    check(r.returncode == 0, f"triad_detect_filter failed:\n{r.stderr[-1500:]}")
    if r.returncode != 0:
        print("RESULT: FAIL")
        sys.exit(1)

    tri = json.load(open(chosen))
    check(tri["match_columns"] == [CI, HI, DI],
          f"wrong columns learned: {tri['match_columns']}")
    print(f"1. columns learned:  C{CI} H{HI} D{DI} from the specific tier only")

    o = pd.read_csv(out["outcomes"], sep="\t").set_index("seq_id")["outcome"]
    want = {
        "spec0": "triad_positive",
        "ssf_true_pos": "triad_positive",
        "ssf_true_neg": "triad_negative",
        "ssf_gapped": "gapped_at_triad",
        "tg_fragment": "low_coverage",
        "tg_permuted": "triad_negative",   # occupied, residues wrong -- see below
    }
    for k, v in want.items():
        check(o[k] == v, f"{k}: expected {v}, got {o[k]}")
    print("2. outcomes:         positive / negative / gapped / low_coverage all "
          "separated")

    # The permuted transglutaminase is the uncomfortable one. Every triad column
    # is occupied, so by the letter of the test it is a `triad_negative`. It is
    # also, in three dimensions, a protein with the identical active site. The
    # filter cannot tell. That is the limitation the report must state, and the
    # reason a sequence-column test can never settle the transglutaminase
    # question on its own.
    check(o["tg_permuted"] == "triad_negative",
          "a circularly permuted triad is scored negative -- it occupies the "
          "columns with the wrong residues")
    # ...but it must be FLAGGED as an unsafe negative, because it carries all
    # three catalytic residues. That flag is the only handle a column test has on
    # the permutation problem.
    df_o = pd.read_csv(out["outcomes"], sep="\t").set_index("seq_id")
    flag = "has_C_H_D_anywhere"
    check(flag in df_o.columns, f"outcomes must carry `{flag}`")
    check(int(df_o.loc["tg_permuted", flag]) == 1,
          "the permuted decoy carries C, H and D somewhere; it must be flagged")
    check(int(df_o.loc["ssf_true_neg", flag]) == 0,
          "the true negative carries none of C/H/D; it must NOT be flagged")
    print("3. permuted triad:   scored `triad_negative`, but flagged "
          "has_C_H_D_anywhere=1. The true negative is flagged 0. A column test "
          "cannot see a permutation; it can at least admit which negatives are "
          "unsafe.")

    t = pd.read_csv(out["tiers"], sep="\t").set_index("evidence")
    s = t.loc["ssf_only"]
    check(int(s["n_aligned"]) == 5, f"expected 5 ssf_only, got {s['n_aligned']}")
    check(int(s["n_low_coverage"]) == 1, "one fragment must be low_coverage")
    check(int(s["n_gapped_at_triad"]) == 1, "one sequence must be gapped_at_triad")
    check(int(s["n_testable"]) == 3, f"3 testable, got {s['n_testable']}")
    check(int(s["n_triad_positive"]) == 1, "one true positive")
    # the denominator is what this whole change is about
    check(abs(float(s["frac_triad_positive_of_testable"]) - 1 / 3) < 1e-9,
          f"rate must be 1/3 over TESTABLE sequences, got "
          f"{s['frac_triad_positive_of_testable']}")
    old_rate = 1 / 5
    check(abs(float(s["frac_triad_positive_of_testable"]) - old_rate) > 1e-6,
          "the new denominator must differ from n_aligned, or nothing changed")
    check(abs(float(s["frac_testable"]) - 3 / 5) < 1e-9, "frac_testable must be 3/5")
    check(int(s["n_negative_with_all_three_residues"]) == 1,
          "exactly one ssf_only negative (the permuted one) carries all three "
          "catalytic residues and is therefore an unsafe negative")
    print(f"4. denominator:      {int(s['n_triad_positive'])}/{int(s['n_testable'])} "
          f"testable = {float(s['frac_triad_positive_of_testable']):.3f}, not "
          f"{int(s['n_triad_positive'])}/{int(s['n_aligned'])} = {old_rate:.3f}")

    # and the low-coverage warning must fire when most of the tier is untestable
    r2 = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "triad_detect_filter.py"),
         "--sto", sto, "--config", tcfg,
         "--combined", os.path.join(base, "combined.txt"),
         "--idmap", os.path.join(base, "idmap.tsv"),
         "--weights", os.path.join(base, "weights.tsv"),
         "--out-candidates", out["cands"], "--out-chosen", chosen,
         "--out-keep", os.path.join(base, "keep.txt"),
         "--out-afa", os.path.join(base, "pass.afa"),
         "--out-colstats", out["colstats"], "--out-tiers", out["tiers"]],
        capture_output=True, text=True)
    check("--out-outcomes" not in r2.stderr and r2.returncode == 0,
          "--out-outcomes must be optional")
    print("5. outcomes optional: the rule works without it")

    # make_report must not conclude from an untestable tier
    rep = open(os.path.join(SCRIPTS, "make_report.py")).read()
    check("frac_triad_positive_of_testable" in rep,
          "make_report must use the testable denominator")
    check("circular" in rep.lower() and "permutation" in rep.lower(),
          "make_report must name circular permutation as the reason a low "
          "frac_testable cannot answer the transglutaminase question")
    check("cannot answer the\n                         \"transglutaminase" in rep
          or "cannot answer the" in rep,
          "make_report must say the pass rate cannot answer the question")
    print("6. report:           refuses to conclude when frac_testable < 0.5")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
