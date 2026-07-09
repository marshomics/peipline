#!/usr/bin/env python3
"""Left-merge the filtered HMM hit table onto the sample metadata table.

Left = the hit table, so every hit survives and the metadata (taxonomy etc.)
is annotation. The join key is `sample` on both sides -- the hit table carries
it because batch_faa.py wrote it into the FASTA description field.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--sample-col", default="sample")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    hits = pd.read_csv(a.combined, sep="\t", dtype={"sample": str, "protein_id": str})
    meta = pd.read_csv(a.table, sep="\t", dtype={a.sample_col: str}, low_memory=False)

    if meta[a.sample_col].duplicated().any():
        sys.exit(f"[merge] duplicated '{a.sample_col}' in {a.table}: a left merge "
                 f"would multiply hit rows. Deduplicate first.")

    meta = meta.rename(columns={a.sample_col: "sample"})
    # `faa` already lives in the hit table; keep the metadata copy under a
    # distinct name rather than letting pandas append _x/_y suffixes.
    overlap = (set(meta.columns) & set(hits.columns)) - {"sample"}
    meta = meta.rename(columns={c: f"{c}_meta" for c in overlap})

    merged = hits.merge(meta, on="sample", how="left",
                        validate="many_to_one", indicator="_merge")

    unmatched = int((merged["_merge"] == "left_only").sum())
    merged = merged.drop(columns="_merge")
    print(f"[merge] {len(hits)} hit rows -> {len(merged)} merged rows; "
          f"{unmatched} with no metadata match", file=sys.stderr)
    if len(merged) != len(hits):
        sys.exit("[merge] row count changed: the join was not many-to-one.")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    merged.to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
