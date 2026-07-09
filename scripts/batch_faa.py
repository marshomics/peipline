#!/usr/bin/env python3
"""Concatenate a batch's faa files into one, renaming sequences to short IDs.

Design note. At ~350k proteomes this step handles on the order of 1e9 protein
sequences, so we must not persist a per-sequence mapping table -- that would be
tens of gigabytes of bookkeeping to recover a few thousand hits.

Instead each sequence is written as

    >b00042_s00000123 <protein_id> <sample>

The synthetic name is collision-free by construction (no delimiter can clash
with a sample or protein ID) and short, which keeps the domtblout small. HMMER
copies everything after the name into the `description of target` field of the
domtblout, so the provenance of a *hit* comes back for free and only the hits
cost us anything. The per-batch map file therefore stores one row per sample
(proteome size, status), not one row per sequence.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import read_fasta, seq_id, write_fasta  # noqa: E402

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUO")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", required=True)
    ap.add_argument("--batch-id", required=True)
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--allow-missing", action="store_true")
    a = ap.parse_args()

    df = pd.read_csv(a.batch, sep="\t", dtype=str)
    os.makedirs(os.path.dirname(a.out_faa), exist_ok=True)
    os.makedirs(os.path.dirname(a.out_map), exist_ok=True)

    n_seq = 0
    n_bad = 0

    with open(a.out_faa, "w") as fo, gzip.open(a.out_map, "wt") as fm:
        fm.write("sample\tfaa\tstatus\tn_proteins\tn_residues\n")
        for sample, faa in zip(df["sample"], df["faa"]):
            if not isinstance(faa, str) or not isinstance(sample, str):
                fm.write(f"{sample}\t{faa}\tnull_path\t0\t0\n")
                n_bad += 1
                continue
            if any(c.isspace() for c in sample):
                sys.exit(f"[batch_faa] sample ID contains whitespace: {sample!r}. "
                         f"Provenance is carried in the FASTA description field, "
                         f"which is whitespace-delimited. Rename the sample.")
            if not os.path.exists(faa):
                fm.write(f"{sample}\t{faa}\tmissing\t0\t0\n")
                n_bad += 1
                continue
            if os.path.getsize(faa) == 0:
                fm.write(f"{sample}\t{faa}\tempty\t0\t0\n")
                n_bad += 1
                continue

            try:
                got = res = 0
                for header, seq in read_fasta(faa):
                    if not seq:
                        continue
                    # HMMER rejects '*' (stop codons from Prodigal) and other
                    # stray characters in an amino-acid alphabet.
                    clean = "".join(c for c in seq.upper() if c in VALID_AA)
                    if not clean:
                        continue
                    sid = f"b{a.batch_id}_s{n_seq:08d}"
                    write_fasta(fo, f"{sid} {seq_id(header)} {sample}", clean)
                    n_seq += 1
                    got += 1
                    res += len(clean)
            except (OSError, EOFError, gzip.BadGzipFile) as e:
                fm.write(f"{sample}\t{faa}\tunreadable:{type(e).__name__}\t0\t0\n")
                n_bad += 1
                continue

            fm.write(f"{sample}\t{faa}\t{'ok' if got else 'no_seqs'}\t{got}\t{res}\n")
            if not got:
                n_bad += 1

    if n_bad and not a.allow_missing:
        sys.exit(f"[batch_faa] batch {a.batch_id}: {n_bad}/{len(df)} samples unusable. "
                 f"See {a.out_map}. Set allow_missing_faa: true to tolerate.")

    print(f"[batch_faa] batch {a.batch_id}: {n_seq} sequences from "
          f"{len(df) - n_bad}/{len(df)} samples ({n_bad} skipped)", file=sys.stderr)


if __name__ == "__main__":
    main()
