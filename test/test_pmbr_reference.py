#!/usr/bin/env python3
"""Test the binding-module rule and the machinery that keeps it honest.

What is checked:

  1. `predict_binding` reproduces every GFP-fusion construct Visweswaran et al.
     2011 actually measured (1P, 2P, 3P, and PeiW's own four-motif domain);
  2. the rule is FALSIFIABLE -- move the threshold and the constructs disagree;
  3. the two-motif case is the interesting one: it binds spheroplasts and NOT
     pseudomurein, so `pmbr_binding_competent` must be 0 while the protein is
     still a plausible murein binder;
  4. the repeat count is fragile, and `domain_arch.py` says so. A synthetic
     protein with one strong and two weak PMB repeats must be called
     `pmbr_count_ambiguous`, not `reduced_pmbr`;
  5. `sdp.py` refuses `pmbr_partition_mode: both`, because the count and the
     binarized partition are nested and would double-count one observation;
  6. the pI is computed per protein, not borrowed from MTH719;
  7. no file still records the threshold as needed for "optimal" binding. Three
     motifs are needed for ANY binding to pseudomurein. The abstract of that
     paper says otherwise; the Results, heading, title and conclusion do not.

Run:  python test/test_pmbr_reference.py [workdir]
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


# A domtblout line as hmmscan writes it: target model first, query protein after.
def domtbl_line(model, acc, tlen, prot, qlen, ievalue, score, hf, ht, ef, et):
    return (f"{model} {acc} {tlen} {prot} - {qlen} 1e-30 100.0 0.0 1 1 "
            f"{ievalue} {ievalue} {score} 0.0 {hf} {ht} {ef} {et} {ef} {et} 0.95 -")


def main() -> None:
    import pmbr_reference as PR

    # --- 1. the measured constructs -----------------------------------------
    for name, c in PR.CONSTRUCTS.items():
        call, can_bind, _ = PR.predict_binding(c["n_motifs"])
        check(can_bind == c["pseudomurein"],
              f"{name}: predicted sacculus binding {can_bind}, observed "
              f"{c['pseudomurein']}")
        if c["n_motifs"] == 2:
            check(call == "murein_fragments_only",
                  "the two-motif construct bound spheroplasts and not "
                  "pseudomurein; the call must record both halves")
        if c["n_motifs"] == 1:
            check(call == "none", "the one-motif construct bound nothing")
    print(f"1. constructs:       {len(PR.CONSTRUCTS)}/{len(PR.CONSTRUCTS)} "
          f"reproduced (1P none, 2P spheroplast-only, 3P and PeiW-PMB sacculus)")

    # --- 2. falsifiability ---------------------------------------------------
    saved = PR.MOTIFS_FOR_PSEUDOMUREIN
    try:
        PR.MOTIFS_FOR_PSEUDOMUREIN = 2      # the abstract's (wrong) reading
        bad = [n for n, c in PR.CONSTRUCTS.items()
               if PR.predict_binding(c["n_motifs"])[1] != c["pseudomurein"]]
        check(bad,
              "moving the threshold to 2 motifs did not contradict any construct. "
              "The rule is not reading the data and the agreement above is empty")
        check("2P-GFP" in bad,
              "the 2-motif construct is what discriminates the two readings; it "
              "must be the one that breaks")
        print(f"2. falsifiable:      threshold=2 contradicts {bad}, as it must")
    finally:
        PR.MOTIFS_FOR_PSEUDOMUREIN = saved

    # --- 3. the pI is per protein -------------------------------------------
    acidic = PR.isoelectric_point("DDDDEEEEDDDD")
    basic = PR.isoelectric_point("KKKKRRRRKKKK")
    check(acidic is not None and acidic < 4.5, f"poly-D/E pI came out {acidic}")
    check(basic is not None and basic > 10.0, f"poly-K/R pI came out {basic}")
    check(PR.isoelectric_point("") is None, "empty sequence must give None")
    check("9.2" not in PR.assay_ph_advice(basic),
          "assay_ph_advice must use the protein's own pI, not MTH719's 9.2")
    check("aggregat" in PR.assay_ph_advice(7.2),
          "a domain with pI near 7 must be warned about aggregation at pH 7.0")
    print(f"3. pI:               poly-D/E {acidic}, poly-K/R {basic}, per protein")

    # --- 4. domain_arch: fragile counts -------------------------------------
    base = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="c71_pmbr_")
    os.makedirs(base, exist_ok=True)
    import yaml

    faa = os.path.join(base, "prot.faa")
    with open(faa, "w") as fh:
        # solid: 3 strong PMB repeats + catalytic domain
        fh.write(">solid\n" + "A" * 40 + "K" * 200 + "\n")
        # fragile: 1 strong + 2 weak repeats -> 1 at strict E, 3 at permissive E
        fh.write(">fragile\n" + "D" * 40 + "E" * 200 + "\n")
        # duo: exactly 2 strong repeats -> murein fragments only, not ambiguous
        fh.write(">duo\n" + "R" * 40 + "R" * 200 + "\n")

    dt = os.path.join(base, "scan.domtbl")
    L = ["# hmmscan domtblout stub"]
    for prot, strong, weak in (("solid", 3, 0), ("fragile", 1, 2), ("duo", 2, 0)):
        pos = 1
        for _ in range(strong):
            L.append(domtbl_line("PMBR", "PF09373.1", 32, prot, 240, "1e-12",
                                 60.0, 1, 32, pos, pos + 31))
            pos += 35
        for _ in range(weak):
            # clears 1e-2 but not 1e-5: exactly the repeat that flips the class
            L.append(domtbl_line("PMBR", "PF09373.1", 32, prot, 240, "1e-3",
                                 12.0, 1, 32, pos, pos + 31))
            pos += 35
        L.append(domtbl_line("Peptidase_C71", "PF12386.1", 150, prot, 240,
                             "1e-40", 200.0, 1, 150, 150, 240))
    open(dt, "w").write("\n".join(L) + "\n")

    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    cfg["specificity"]["pmbr"] = {"count_evalue_strict": 1e-5,
                                  "count_evalue_permissive": 1e-2,
                                  "report_pi": True}
    tcfg = os.path.join(base, "config.test.yaml")
    yaml.safe_dump(cfg, open(tcfg, "w"))

    out_arch = os.path.join(base, "arch.tsv")
    r = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "domain_arch.py"),
         "--faa", faa, "--config", tcfg, "--domtbl", dt, "--skip-scan",
         "--out-arch", out_arch, "--out-domains", os.path.join(base, "dom.tsv")],
        capture_output=True, text=True)
    check(r.returncode == 0, f"domain_arch.py failed:\n{r.stderr[-1500:]}")

    import pandas as pd
    A = pd.read_csv(out_arch, sep="\t").set_index("seq_id")

    check(A.loc["solid", "n_pmbr"] == 3 and A.loc["solid", "pmbr_binding_competent"] == 1,
          f"3 strong repeats must bind the sacculus, got "
          f"n={A.loc['solid', 'n_pmbr']} competent={A.loc['solid', 'pmbr_binding_competent']}")
    check(A.loc["solid", "architecture_class"] == "canonical_pei",
          f"got {A.loc['solid', 'architecture_class']}")

    check(A.loc["duo", "n_pmbr"] == 2 and A.loc["duo", "pmbr_binding_competent"] == 0,
          "2 motifs must NOT be called sacculus-competent")
    check(A.loc["duo", "predicted_binding"] == "murein_fragments_only",
          f"2 motifs bind lysozyme-exposed murein only, got "
          f"{A.loc['duo', 'predicted_binding']}")
    check(A.loc["duo", "pmbr_count_fragile"] == 0,
          "a clean 2-motif protein is not ambiguous, it is a negative")
    check(A.loc["duo", "architecture_class"] == "reduced_pmbr",
          f"got {A.loc['duo', 'architecture_class']}")

    check(A.loc["fragile", "n_pmbr"] == 1 and A.loc["fragile", "n_pmbr_permissive"] == 3,
          f"strict count should be 1 and permissive 3, got "
          f"{A.loc['fragile', 'n_pmbr']} / {A.loc['fragile', 'n_pmbr_permissive']}")
    check(A.loc["fragile", "pmbr_count_fragile"] == 1,
          "a protein whose count straddles the 3-motif cliff must be flagged")
    check(A.loc["fragile", "architecture_class"] == "pmbr_count_ambiguous",
          f"a fragile protein must not be assigned to either side, got "
          f"{A.loc['fragile', 'architecture_class']}")
    print("4. count fragility:  1 strong + 2 weak repeats -> pmbr_count_ambiguous, "
          "not reduced_pmbr")

    # the pI must differ between the acidic and basic synthetic PMB spans
    pi_f = float(A.loc["fragile", "pmbr_pi"])
    pi_d = float(A.loc["duo", "pmbr_pi"])
    check(pi_f < 5.0 < pi_d,
          f"per-protein pI is not being computed from the PMB span: "
          f"fragile(poly-D)={pi_f}, duo(poly-R)={pi_d}")
    print(f"5. per-protein pI:   poly-D span {pi_f}, poly-R span {pi_d}")

    # --- 6. sdp.py rejects a nested partition mode ---------------------------
    src = open(os.path.join(SCRIPTS, "sdp.py")).read()
    check("pmbr_partition_mode" in src,
          "sdp.py must read pmbr_partition_mode")
    check(re.search(r"nested", src, re.I),
          "sdp.py must explain why count and binding_competent cannot both be "
          "used: they are nested and would double-count one observation")
    check("'count', 'binding_competent'" in src or
          '("count", "binding_competent")' in src,
          "sdp.py must validate pmbr_partition_mode against exactly two values")
    print("6. partition nesting: sdp.py validates the mode and refuses both")

    # --- 7. the weaker wording must not return -------------------------------
    # Grep for "optimal binding" anywhere, having first removed the paper's own
    # title -- which contains the phrase and is the only legitimate use of it.
    # A narrower regex missed the inline comment `# >=3 motifs: optimal binding`,
    # which is exactly the form the mistake took the first time.
    TITLE = ("A Minimum of Three Motifs Is Essential for Optimal Binding of "
             "Pseudomurein Cell Wall-Binding Domain")
    offenders = []
    for dirpath, _, files in os.walk(ROOT):
        if any(p in dirpath for p in (".git", "__pycache__")):
            continue
        for f in files:
            if not f.endswith((".py", ".yaml", ".md")) and f != "Snakefile":
                continue
            p = os.path.join(dirpath, f)
            if os.path.abspath(p) == os.path.abspath(__file__):
                continue
            # Collapse whitespace FIRST: the title is line-wrapped in the
            # docstring, so stripping it before unwrapping matches nothing and
            # the legitimate citation gets reported as an offender.
            text = re.sub(r"\s+", " ", open(p, errors="ignore").read())
            text = re.sub(re.escape(TITLE), "", text, flags=re.I)
            for m in re.finditer(r".{0,60}optimal\s+binding.{0,30}", text, re.I):
                offenders.append(f"{os.path.relpath(p, ROOT)}: ...{m.group(0).strip()}...")
    check(not offenders,
          "a file still frames the 3-motif threshold as needed for *optimal* "
          "binding. That is the abstract's wording. The Results, heading, title "
          "and conclusion all say three motifs are needed for ANY binding to "
          "pseudomurein; two bind only lysozyme-treated bacterial spheroplasts:\n  "
          + "\n  ".join(offenders))
    print("7. wording:          no file records the threshold as merely 'optimal'")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
