#!/usr/bin/env python3
"""Static checks on phyloglm.R, which no test in this repo can execute.

There is no R in the development sandbox, so the phylogenetic regression has
never been run here. That is exactly why it accumulated a silent, total failure:

    lapply(setNames(numcov, numcov), function(x) mean(x, na.rm = TRUE))

`setNames(numcov, numcov)` is a character vector of column NAMES. lapply passed
each name as a string, so R evaluated mean("completeness") -> NA with a warning.
Every covariate became NA for every species, complete.cases() dropped every row,
and the headline inference never ran. Every sg_* indicator became 0 the same way,
so the subgroup models were all skipped for "0 positives". The script exited 0.

These checks cannot prove the R is correct. They can stop the specific mistakes
from coming back, and they cost nothing.

Run:  python test/test_phyloglm_static.py
"""
from __future__ import annotations

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "scripts", "phyloglm.R")

FAILURES = []


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def strip_comments(src):
    """Drop R comments. The fix documents the old idiom verbatim, and a naive
    grep would flag the explanation as if it were the bug."""
    out = []
    for line in src.split("\n"):
        in_str, esc, cut = None, False, len(line)
        for i, ch in enumerate(line):
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif in_str:
                if ch == in_str:
                    in_str = None
            elif ch in "\"'":
                in_str = ch
            elif ch == "#":
                cut = i
                break
        out.append(line[:cut])
    return "\n".join(out)


def main() -> None:
    src = strip_comments(open(R).read())

    # 1. the aggregation bug must not return
    check(not re.search(r"lapply\s*\(\s*setNames\s*\(\s*(numcov|sgcols)",
                        src),
          "phyloglm.R aggregates with lapply(setNames(<colnames>, ...), f). "
          "lapply passes the NAME, not the column: mean('completeness') is NA. "
          "Use lapply(.SD, f) with .SDcols=")
    check(".SDcols" in src,
          "the aggregation must use .SD/.SDcols so the columns themselves are "
          "passed to mean()/any()")
    print("1. aggregation:      uses .SD/.SDcols, not lapply(setNames(names))")

    # 2. the all-NA guard turns a silent failure into a loud one
    check(re.search(r"all\(is\.na\(agg\[\[cc\]\]\)\)", src),
          "phyloglm.R must stop when a covariate is all-NA after aggregation; "
          "that is an aggregation bug, not missing data")
    print("2. all-NA guard:     stops rather than silently dropping every species")

    # 3. phylolm aligns by row name, so row names must be tip labels
    check(re.search(r"rownames\(ddf\)\s*<-\s*ddf\$tree_tip", src),
          "the data frame handed to phyloglm() must carry tip labels as row "
          "names. as.data.frame() on a data.table gives integer row names, and "
          "phylolm aligns by name -- silently mismatching covariates to tips, or "
          "erroring into the try() and leaving only the non-phylogenetic glm")
    print("3. tip alignment:    row names set to tree_tip before phyloglm()")

    # 4. non-convergent / boundary fits must be discarded, not pooled
    check("at_bound" in src and "n_bad" in src,
          "phyloglm returns a completed object when alpha hits the btol bound; "
          "its standard errors are then untrustworthy and nothing throws. Such "
          "replicates must be discarded, not pooled")
    check("phyloglm_replicates_discarded" in src,
          "the number of discarded replicates must reach the output table")
    print("4. convergence:      boundary/non-convergent replicates discarded and "
          "counted")

    # 5. Rubin's rules need a t reference with the Rubin df, not a normal
    check("rubin_df" in src and re.search(r"pt\(-abs\(tstat\)", src),
          "Rubin's rules pool with a t reference on the Rubin df. Using pnorm "
          "treats the variance as known when it was estimated from m replicates")
    print("5. pooling:          t reference with the Rubin df")

    # 6. separation in the taxon-adjusted prevalence model must be handled
    check("brglm2" in src or "logistf" in src,
          "has_c71 is rare and taxon is a many-level factor, so taxa with zero "
          "positives are perfectly separated: glm returns a divergent coefficient "
          "and plogis() pins the adjusted prevalence at 0 or 1 with a [0,1] CI. "
          "Penalise the likelihood (Firth)")
    check("separated" in src,
          "separated taxa must be flagged in the output, whether or not a "
          "penalised fit was available")
    print("6. separation:       Firth when available, separated taxa flagged")

    # 7. rough balance check: R has no compiler here, so at least catch a
    #    grossly unbalanced edit
    for opener, closer in (("(", ")"), ("{", "}"), ("[", "]")):
        n_o, n_c = src.count(opener), src.count(closer)
        check(n_o == n_c,
              f"unbalanced {opener}{closer} in phyloglm.R: {n_o} vs {n_c}. "
              f"There is no R here to parse it, so this is the only structural "
              f"check available")
    print("7. balance:          parentheses, braces and brackets balanced")

    print()
    print("NOTE: phyloglm.R has NOT been executed. These are static checks only. "
          "Run it once on the cluster before trusting any coefficient.")
    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
