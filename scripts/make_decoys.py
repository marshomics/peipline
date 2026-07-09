#!/usr/bin/env python3
"""Reverse every sequence to build a decoy database.

Reversal preserves length and amino-acid composition (and therefore the
low-complexity and biased-composition regions that inflate HMMER bit scores)
while destroying homology. The fraction of decoy hits at a given score is a
direct empirical estimate of the false-discovery rate at that score, which is
the honest way to defend a bit-score cutoff of 25 against ~10^9 target
sequences rather than trusting the E-value's asymptotics.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import read_fasta, write_fasta  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-faa", required=True)
    ap.add_argument("--out-faa", required=True)
    a = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(a.out_faa)), exist_ok=True)
    n = 0
    with open(a.out_faa, "w") as fh:
        for header, seq in read_fasta(a.in_faa):
            write_fasta(fh, header, seq[::-1])
            n += 1
    print(f"[decoys] {n} reversed sequences -> {a.out_faa}", file=sys.stderr)


if __name__ == "__main__":
    main()
