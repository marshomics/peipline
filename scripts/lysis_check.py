#!/usr/bin/env python3
"""Score the pipeline's host-chemistry rule against the only measured phenotype.

Subedi et al. 2015 put purified PeiW and PeiP on plates of eleven methanogens and
recorded which lysed. That table is the one place in this project where a
prediction can be wrong in a way that shows.

The rule under test is stated in lysis_reference.predict_lysis and derives from
two independent observations: the pNA substrate series (Ala and Ser at P1 are
cleaved, Thr is not) and M. ruminantium M1 (Thr wall, resistant). It sees only
the host's wall chemistry from Kandler & Koenig 1978. It never sees the panel.

Rows whose host has no pseudomurein are reported but flagged `nontrivial=False`:
predicting that an enzyme cannot cut a wall that does not exist is not a test.

This rule runs before the screen and does not depend on it. If it fails, the
chemistry table or the substrate table is wrong, and every specificity claim
downstream inherits the error. Better to find out in ten milliseconds.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lysis_reference import (CITATION, LYSIS_PANEL, METAL,  # noqa: E402
                             SUBSTRATES, check_against_reference, summarise)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-substrates", required=True)
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any falsifiable prediction disagrees")
    a = ap.parse_args()

    import pandas as pd
    for p in (a.out, a.out_substrates):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    df = check_against_reference()
    df.to_csv(a.out, sep="\t", index=False)

    sub = pd.DataFrame(SUBSTRATES).T.rename_axis("substrate").reset_index()
    sub.to_csv(a.out_substrates, sep="\t", index=False)

    n, k, nt, nt_k, bad = summarise(df)
    print(f"[lysis] reference: {CITATION}", file=sys.stderr)
    print(f"[lysis] {k}/{n} predictions agree; {nt_k}/{nt} of the falsifiable "
          f"(pseudomurein-walled) ones", file=sys.stderr)
    print(f"[lysis] PeiW metal rescue {'>'.join(METAL['peiW_order'])}; "
          f"PeiP {'/'.join(METAL['peiP_order'])} only", file=sys.stderr)

    na = df[df["agrees"].isna()]
    for _, r in na.iterrows():
        print(f"[lysis] no prediction for {r['strain']} / {r['enzyme']} "
              f"(p1_source={r['p1_source']}); observed {r['observed']}",
              file=sys.stderr)

    unver = [e["label"] for e in LYSIS_PANEL if not e["verified"]]
    if unver:
        print(f"[lysis] {len(unver)} panel rows are excluded because the species "
              f"name could not be read from the machine-readable text: {unver}. "
              f"Confirm against Table 4 and set verified=True.", file=sys.stderr)

    if len(bad):
        print("\n[lysis] DISAGREEMENTS:", file=sys.stderr)
        print(bad[["strain", "enzyme", "observed", "predicted", "reason"]]
              .to_string(index=False), file=sys.stderr)
        if a.strict:
            sys.exit("[lysis] the host-chemistry rule contradicts the measured "
                     "lysis panel. Fix the rule before trusting any specificity "
                     "call downstream.")
    print(f"[lysis] {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
