#!/usr/bin/env python3
"""`c71.faa` pools two evidence tiers. Make sure the pipeline says so, and tests it.

A sequence enters `c71.faa` if it carries the triad, whether it cleared PF12386's
curated gathering threshold (`specific`) or only the SSF54001 fold model
(`ssf_only`). After that point nothing distinguishes them: the tree, the
subgroups, the SSN, the coupling, dN/dS and the prevalence model all treat them
identically.

That is defensible if the ssf_only sequences really are C71 proteins that
PF12386's threshold was too strict to catch. It is not defensible if they are a
different family that happens to carry Cys/His/Asp at the same columns -- and
SSF54001 spans ~22 families.

The gene tree can tell the difference, and `convergence.py` now asks it:

    ssf_only tips CLUMPED     -> their own lineage -> probably a different family
    ssf_only tips INTERSPERSED -> PF12386's threshold is just conservative

This test plants both worlds on the same tree topology and asserts the answer
flips. It also checks that `extract_seqs.py` writes the tier beside the FASTA and
warns when ssf_only sequences are present.

Run:  python test/test_evidence_pooling.py [workdir]
"""
from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")

FAILURES = []
N = 32                     # tips, a balanced binary tree


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def balanced(ids, depth=0):
    if len(ids) == 1:
        return f"{ids[0]}:0.1"
    h = len(ids) // 2
    return f"({balanced(ids[:h], depth + 1)},{balanced(ids[h:], depth + 1)}):0.1"


def run_convergence(base, tag, evidence_of, tips):
    import pandas as pd
    import yaml

    tree = os.path.join(base, f"{tag}.treefile")
    open(tree, "w").write(balanced(tips) + ";\n")

    assign = os.path.join(base, f"assign_{tag}.tsv")
    pd.DataFrame({
        "seq_id": tips,
        # one subgroup for everyone: the subgroup traits are not what is on trial
        "subgroup": 0,
        "barcode": "X", "weight": 1.0, "cluster": tips,
        "evidence": [evidence_of[t] for t in tips],
    }).to_csv(assign, sep="\t", index=False)

    reps = os.path.join(base, f"reps_{tag}.tsv")
    pd.DataFrame({"seq_id": tips}).to_csv(reps, sep="\t", index=False)

    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    cfg["convergence"] = {"n_permutations": 300, "n_brownian": 200, "seed": 7}
    tcfg = os.path.join(base, f"cfg_{tag}.yaml")
    yaml.safe_dump(cfg, open(tcfg, "w"))

    out = os.path.join(base, f"convergence_{tag}.tsv")
    r = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "convergence.py"),
         "--tree", tree, "--assign", assign, "--reps", reps,
         "--config", tcfg, "--figdir", os.path.join(base, "fig"), "--out", out],
        capture_output=True, text=True,
        env=dict(os.environ, MPLBACKEND="Agg"))
    if r.returncode != 0:
        print(r.stderr[-1500:])
    return r, (pd.read_csv(out, sep="\t") if os.path.exists(out) else None)


def main() -> None:
    base = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="c71_ev_")
    os.makedirs(base, exist_ok=True)
    import pandas as pd

    tips = [f"c71_{i:07d}" for i in range(N)]

    # --- world A: the ssf_only tips are one contiguous clade -----------------
    # In a balanced binary tree, the first 8 tips are a clade. If SSF54001 let in a
    # different family, this is what it would look like.
    clumped = {t: ("ssf_only" if i < 8 else "specific") for i, t in enumerate(tips)}
    rA, dA = run_convergence(base, "clumped", clumped, tips)
    check(rA.returncode == 0, f"convergence failed on the clumped world:\n{rA.stderr[-800:]}")
    check("ssf_only" in rA.stderr, "convergence must announce the ssf_only trait")

    rowA = dA[dA["subgroup"] == "evidence:ssf_only"]
    check(len(rowA) == 1, "the evidence trait must appear in convergence.tsv")
    dA_val = float(rowA.iloc[0]["fritz_purvis_D"])
    origins_A = int(rowA.iloc[0]["observed_origins"])
    check(origins_A == 1, f"a contiguous clade must have 1 origin, got {origins_A}")
    check(dA_val < 0.25, f"a clade must give D < 0.25, got {dA_val}")
    check("own lineage" in rowA.iloc[0]["interpretation"],
          f"the interpretation must warn that c71.faa pools two families, got "
          f"{rowA.iloc[0]['interpretation']!r}")
    print(f"1. clumped world:    {origins_A} origin, D = {dA_val:.3f} -> "
          f"'own lineage' warning fires")

    # --- world B: same tree, ssf_only scattered across it --------------------
    # Every fourth tip. If PF12386's gathering threshold is merely conservative,
    # this is what it would look like.
    inter = {t: ("ssf_only" if i % 4 == 0 else "specific") for i, t in enumerate(tips)}
    rB, dB = run_convergence(base, "interspersed", inter, tips)
    check(rB.returncode == 0, f"convergence failed on the interspersed world:\n"
                              f"{rB.stderr[-800:]}")
    rowB = dB[dB["subgroup"] == "evidence:ssf_only"]
    check(len(rowB) == 1, "the evidence trait must appear in convergence.tsv")
    dB_val = float(rowB.iloc[0]["fritz_purvis_D"])
    origins_B = int(rowB.iloc[0]["observed_origins"])
    check(origins_B == 8, f"8 scattered tips must give 8 origins, got {origins_B}")
    check(dB_val > dA_val,
          f"scattered tips must give a higher D than a clade: {dB_val} vs {dA_val}")
    check("interspersed" in rowB.iloc[0]["interpretation"],
          f"got {rowB.iloc[0]['interpretation']!r}")
    print(f"2. interspersed:     {origins_B} origins, D = {dB_val:.3f} -> "
          f"'threshold is conservative' reading")

    check(dA_val < 0.25 < dB_val,
          "the test must SEPARATE the two worlds. If both land on the same side "
          "of the threshold the statistic is not doing anything.")
    print(f"3. discrimination:   D flips across 0.25 ({dA_val:.3f} -> {dB_val:.3f}) "
          f"on the same topology, so the trait is read, not the tree")

    # --- extract_seqs must write the tier and warn ---------------------------
    idmap = os.path.join(base, "idmap.tsv.gz")
    faa_src = os.path.join(base, "src.faa")
    with open(faa_src, "w") as fh:
        for t in tips:
            fh.write(f">{t}_prot\nMCVKHAAD\n")
    pd.DataFrame({"seq_id": tips, "sample": "S", "protein_id": [f"{t}_prot" for t in tips],
                  "faa": faa_src, "profiles_hit": "x",
                  "evidence": [inter[t] for t in tips]}).to_csv(
        idmap, sep="\t", index=False, compression="gzip")
    keep = os.path.join(base, "keep.txt")
    open(keep, "w").write("\n".join(tips) + "\n")

    ev_out = os.path.join(base, "c71_evidence.tsv")
    r = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "extract_seqs.py"),
         "--keep-ids", keep, "--idmap", idmap,
         "--out-faa", os.path.join(base, "c71.faa"),
         "--out-evidence", ev_out, "--threads", "1"],
        capture_output=True, text=True)
    check(r.returncode == 0, f"extract_seqs failed:\n{r.stderr[-800:]}")
    check(os.path.exists(ev_out), "c71_evidence.tsv must be written")
    e = pd.read_csv(ev_out, sep="\t")
    check(set(e["evidence"]) == {"specific", "ssf_only"}, "both tiers must be recorded")
    check(int((e["evidence"] == "ssf_only").sum()) == 8, "8 ssf_only sequences")
    check("composition" in r.stderr, "the composition must be printed")
    check("NO downstream analysis conditions on this" in r.stderr,
          "extract_seqs must say plainly that nothing downstream stratifies on the "
          "tier; that is the whole risk")
    check("convergence.tsv" in r.stderr,
          "the warning must point at the check that can settle it")
    print("4. c71.faa tier:     written to c71_evidence.tsv; the warning names the "
          "check that settles it")

    # and it must stay quiet when there is nothing to warn about
    pd.DataFrame({"seq_id": tips, "sample": "S",
                  "protein_id": [f"{t}_prot" for t in tips], "faa": faa_src,
                  "profiles_hit": "x", "evidence": "specific"}).to_csv(
        idmap, sep="\t", index=False, compression="gzip")
    r = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "extract_seqs.py"),
         "--keep-ids", keep, "--idmap", idmap,
         "--out-faa", os.path.join(base, "c71b.faa"),
         "--out-evidence", os.path.join(base, "ev2.tsv"), "--threads", "1"],
        capture_output=True, text=True)
    check("NO downstream analysis" not in r.stderr,
          "with a pure-specific c71.faa there is nothing to warn about")
    print("5. no false alarm:   a pure PF12386 c71.faa produces no warning")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
