#!/usr/bin/env python3
"""The C39 arm must find PeiR-class enzymes the C71 arm throws away.

The whole reason PF03412 was disabled for so long: a C39 protein cannot be scored
against the C71 Cys->His prior. PeiR's catalytic cysteine is C90 and its nearest
downstream histidine is 72 residues away; the C71 prior expects 35. Borrow the C71
prior and you either reject PeiR or, worse, snap its "triad" onto a spurious
histidine that happens to sit 35 residues from the cysteine.

This test builds ONE alignment and runs the SAME detector two ways -- `--family
c39` and `--family c71` -- and shows they land on different columns. Nothing about
the alignment changes; only which family's rule is applied.

The alignment (100 match columns, index == column because every state is a match):

  * a real, fully conserved C39 catalytic triad: C at 5, H at 77, D at 90
    (Cys->His gap 72, His->Asp gap 13 -- PeiR-like, wide).
  * a DECOY histidine at 40 and a decoy aspartate at 57, present in only 60% of
    the sequences (so lower frequency). C5/H40/D57 is a spurious triple sitting at
    the C71 spacing (Cys->His gap 35).

With a deliberately hostile C71 spacing prior in the `triad:` block (gap 35, tight
tolerance, high weight):

  family c39  -> ignores the prior, ranks by frequency, recovers C5/H77/D90 and
                 REPORTS the spacing it found (learned_gaps [72, 13]).
  family c71  -> obeys the prior, is dragged onto the spurious C5/H40/D57.

If the two arms ever agree here, the family split is not doing its job.

Also checks, in the same run:
  * the learned gap of 72 is exactly what the C71 35+/-8 prior rejects
    (ties back to pei_reference.spacing_prior_admits);
  * extract_seqs --family routes a PF03412-only protein to C39 and keeps it out
    of C71, while the SSF54001 fold net is shared;
  * an empty C39 net (empty alignment) is tolerated: empty outputs, exit 0.

Run:  python test/test_c39_arm.py [workdir]
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

FAILURES = []

L = 100
CI, H_TRUE, D_TRUE = 5, 77, 90     # the real C39 triad: wide Cys->His gap (72)
H_DECOY, D_DECOY = 40, 57          # a spurious triple at the C71 spacing (35)


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def build_sto(path, n_spec=10, decoy_frac=0.6):
    """All match states, so match-column index == alignment-column index."""
    seqs = {}
    n_decoy = int(round(n_spec * decoy_frac))
    for i in range(n_spec):
        s = ["V"] * L
        s[CI], s[H_TRUE], s[D_TRUE] = "C", "H", "D"     # real triad, every seq
        if i < n_decoy:                                  # decoys in a minority
            s[H_DECOY], s[D_DECOY] = "H", "D"
        seqs[f"c39spec{i}"] = "".join(s)
    with open(path, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n")
        w = max(len(k) for k in seqs) + 2
        for k, v in seqs.items():
            fh.write(f"{k:<{w}}{v}\n")
        fh.write(f"{'#=GC RF':<{w}}{'x' * L}\n//\n")
    return list(seqs)


def run_triad(base, tag, sto, tcfg, names, family, extra=()):
    import pandas as pd
    pd.DataFrame({"seq_id": names, "evidence": "specific",
                  "sample": "S", "protein_id": names, "faa": "x"}).to_csv(
        os.path.join(base, f"idmap_{tag}.tsv"), sep="\t", index=False)
    pd.DataFrame({"seq_id": names, "weight": 1.0, "cluster": names,
                  "cluster_size": 1}).to_csv(
        os.path.join(base, f"weights_{tag}.tsv"), sep="\t", index=False)
    open(os.path.join(base, "combined.txt"), "w").write("x\n")
    chosen = os.path.join(base, f"chosen_{tag}.json")
    cmd = [sys.executable, os.path.join(SCRIPTS, "triad_detect_filter.py"),
           "--sto", sto, "--config", tcfg, "--family", family,
           "--combined", os.path.join(base, "combined.txt"),
           "--idmap", os.path.join(base, f"idmap_{tag}.tsv"),
           "--weights", os.path.join(base, f"weights_{tag}.tsv"),
           "--out-candidates", os.path.join(base, f"cands_{tag}.tsv"),
           "--out-chosen", chosen,
           "--out-keep", os.path.join(base, f"keep_{tag}.txt"),
           "--out-afa", os.path.join(base, f"pass_{tag}.afa"),
           "--out-colstats", os.path.join(base, f"colstats_{tag}.tsv"),
           "--out-tiers", os.path.join(base, f"tiers_{tag}.tsv"),
           *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r, chosen


def main() -> None:
    base = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="c39_arm_")
    os.makedirs(base, exist_ok=True)
    import pandas as pd
    import yaml
    from pei_reference import spacing_prior_admits

    sto = os.path.join(base, "hits_c39.sto")
    names = build_sto(sto)

    # Real families block (c39 has expected_gap null; c71 has 35). Only the shared
    # residue-detection knobs and a HOSTILE C71 spacing prior are overridden.
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    check((cfg["families"]["c39"].get("expected_gap_1_2") is None),
          "config families.c39.expected_gap_1_2 must be null for a learned arm")
    check(cfg["families"]["c71"].get("expected_gap_1_2") == 35,
          "config families.c71.expected_gap_1_2 must be 35")
    cfg["triad"].update(expected_gap_1_2=35, expected_gap_2_3=24,
                        spacing_tolerance=5, spacing_weight=5.0,   # hostile prior
                        min_occupancy=0.4, min_residue_freq=0.4,
                        alt_third_residue=None, max_candidates_per_residue=40,
                        min_match_coverage=0.5)
    tcfg = os.path.join(base, "config.test.yaml")
    yaml.safe_dump(cfg, open(tcfg, "w"))

    # --- 1. C39 learn mode ignores the prior, recovers the wide-gap triad -------
    r39, ch39 = run_triad(base, "c39", sto, tcfg, names, "c39")
    check(r39.returncode == 0, f"c39 triad failed:\n{r39.stderr[-1200:]}")
    j39 = json.load(open(ch39))
    check(j39["match_columns"] == [CI, H_TRUE, D_TRUE],
          f"c39 must recover the real triad C{CI}/H{H_TRUE}/D{D_TRUE}, "
          f"got {j39['match_columns']}")
    check(j39["spacing_mode"] == "learned", "c39 spacing_mode must be 'learned'")
    check(j39["learned_gaps"] == [H_TRUE - CI, D_TRUE - H_TRUE],
          f"c39 must REPORT the spacing it found [{H_TRUE-CI}, {D_TRUE-H_TRUE}], "
          f"got {j39['learned_gaps']}")
    check(j39["family"] == "c39" and j39["specific_profile"] == "PF03412",
          "c39 chosen.json must record family/profile")
    print(f"1. c39 learn mode:   recovered C{CI}/H{H_TRUE}/D{D_TRUE}, "
          f"learned_gaps={j39['learned_gaps']} (prior IGNORED)")

    # --- 2. C71 prior mode is dragged onto the spurious 35-spaced triple --------
    r71, ch71 = run_triad(base, "c71", sto, tcfg, names, "c71")
    check(r71.returncode == 0, f"c71 triad failed:\n{r71.stderr[-1200:]}")
    j71 = json.load(open(ch71))
    check(j71["spacing_mode"] == "prior", "c71 spacing_mode must be 'prior'")
    check(j71["match_columns"] == [CI, H_DECOY, D_DECOY],
          f"under a hostile C71 prior, c71 mode must snap onto the spurious "
          f"C{CI}/H{H_DECOY}/D{D_DECOY}, got {j71['match_columns']}")
    print(f"2. c71 prior mode:   dragged onto spurious C{CI}/H{H_DECOY}/D{D_DECOY} "
          f"(obeys gap-35 prior)")

    # --- 3. the two arms DISAGREE -- which is the entire point -------------------
    check(j39["match_columns"] != j71["match_columns"],
          "the C39 and C71 arms chose the same columns; the family split is inert")
    print("3. arms disagree:    same alignment, different family rule -> different "
          "triad. The split is load-bearing.")

    # --- 4. the learned gap 72 is exactly what the C71 35+/-8 prior rejects ------
    adm, why = spacing_prior_admits("PeiR", 35, 8)
    check(adm is False, f"C71 35+/-8 prior must reject PeiR: {why}")
    check(j39["learned_gaps"][0] == 72,
          "the learned Cys->His gap must be 72, the PeiR value the C71 prior rejects")
    print("4. ties to PeiR:     learned Cys->His gap 72 is the value the C71 "
          "35+/-8 prior rejects (spacing_prior_admits=False)")

    # --- 5. extract_seqs family scoping -----------------------------------------
    faa = os.path.join(base, "p.faa")
    open(faa, "w").write(">c71only S c71only\nMC\n>c39only S c39only\nMC\n"
                         ">ssfonly S ssfonly\nMW\n")
    hits = os.path.join(base, "hits.tsv")
    open(hits, "w").write(
        "sample\tprotein_id\tfaa\tprofiles_hit\tevidence\n"
        f"S\tc71only\t{faa}\tPF12386\tspecific\n"
        f"S\tc39only\t{faa}\tPF03412\tssf_only\n"
        f"S\tssfonly\t{faa}\tSSF54001\tssf_only\n")

    def extract(fam):
        idm = os.path.join(base, f"idmap_ex_{fam}.tsv.gz")
        r = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "extract_seqs.py"), "--hits", hits,
             "--family", fam, "--config", tcfg, "--out-faa",
             os.path.join(base, f"ex_{fam}.faa"), "--out-idmap", idm, "--threads", "1"],
            capture_output=True, text=True)
        check(r.returncode == 0, f"extract {fam} failed:\n{r.stderr[-600:]}")
        m = pd.read_csv(idm, sep="\t")
        return dict(zip(m["protein_id"], m["evidence"])), list(m["seq_id"])

    e71, ids71 = extract("c71")
    e39, ids39 = extract("c39")
    check("c39only" not in e71,
          "a PF03412-only protein leaked into the C71 arm (the mangling PF03412 "
          "was disabled to prevent)")
    check(e39.get("c39only") == "specific",
          "a PF03412-only protein must be a C39 specific hit")
    check(e71.get("ssfonly") == "ssf_only" and e39.get("ssfonly") == "ssf_only",
          "the SSF54001 fold net must be shared by both arms")
    check(ids71 and ids71[0].startswith("c71_") and ids39 and ids39[0].startswith("c39_"),
          "seq_id prefixes must be family-local")
    print("5. extract scoping:  PF03412-only -> C39 only; SSF54001 shared; "
          "family-prefixed IDs")

    # --- 6. an empty C39 net is tolerated ---------------------------------------
    empty_sto = os.path.join(base, "empty.sto")
    open(empty_sto, "w").write("# STOCKHOLM 1.0\n//\n")
    re_, che = run_triad(base, "empty", empty_sto, tcfg, names, "c39",
                         extra=("--allow-empty-specific",))
    check(re_.returncode == 0,
          f"empty C39 net must exit 0 with --allow-empty-specific:\n{re_.stderr[-800:]}")
    je = json.load(open(che))
    check(je.get("n_triad_positive") == 0 and je.get("match_columns") is None,
          "empty arm chosen.json must record zero positives, no columns")
    keep = open(os.path.join(base, "keep_empty.txt")).read().strip()
    check(keep == "", "empty arm keep list must be empty")
    print("6. empty net:        empty PF03412 alignment -> empty outputs, exit 0, "
          "no crash, no fabricated triad")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
