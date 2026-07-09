#!/usr/bin/env python3
"""Call genes on a chunk of archaeal assemblies with Prodigal.

`-p single` trains a gene model on the genome itself, which is what you want for
a complete or near-complete assembly. It needs ~20 kb of sequence to train on;
below that Prodigal refuses. Rather than let a short or heavily fragmented MAG
kill the array job, this falls back to `-p meta` (pre-trained models) and records
which mode each genome got, because that choice affects gene-call sensitivity
and therefore belongs in the detection-bias bookkeeping downstream.
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config  # noqa: E402


def genome_size(path: str) -> int:
    op = gzip.open if path.endswith(".gz") else open
    n = 0
    with op(path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n


def count_faa(path: str) -> int:
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def run_prodigal(fna, faa, ffn, mode, table):
    cmd = ["prodigal", "-i", fna, "-a", faa, "-d", ffn,
           "-p", mode, "-g", str(table), "-m", "-q", "-o", "/dev/null"]
    return subprocess.run(cmd, capture_output=True, text=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fna", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)["prodigal"]
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.manifest)), exist_ok=True)

    rows = []
    for fna in a.fna:
        stem = os.path.basename(fna)
        for ext in (".fna.gz", ".fna", ".fa", ".fasta"):
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        faa = os.path.join(a.outdir, f"{stem}.faa")
        ffn = os.path.join(a.outdir, f"{stem}.ffn")

        if os.path.exists(faa) and os.path.getsize(faa) > 0:
            rows.append((stem, fna, faa, "cached", count_faa(faa), genome_size(fna)))
            continue

        size = genome_size(fna)
        mode = cfg["mode"] if size >= int(cfg["min_len_for_single"]) else "meta"
        r = run_prodigal(fna, faa, ffn, mode, cfg["translation_table"])

        if r.returncode != 0 and mode == "single":
            print(f"[prodigal] {stem}: single-mode failed, retrying with meta\n"
                  f"{r.stderr.strip()[:300]}", file=sys.stderr)
            mode = "meta"
            r = run_prodigal(fna, faa, ffn, mode, cfg["translation_table"])

        if r.returncode != 0:
            print(f"[prodigal] {stem}: FAILED\n{r.stderr.strip()[:500]}", file=sys.stderr)
            rows.append((stem, fna, "", "failed", 0, size))
            continue

        rows.append((stem, fna, faa, mode, count_faa(faa), size))

    with open(a.manifest, "w") as fh:
        fh.write("stem\tfna\tfaa\tprodigal_mode\tn_proteins\tgenome_size\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    ok = sum(1 for r in rows if r[3] != "failed")
    print(f"[prodigal] {ok}/{len(rows)} genomes called -> {a.manifest}", file=sys.stderr)
    if ok == 0 and rows:
        sys.exit("[prodigal] every genome in this chunk failed")


if __name__ == "__main__":
    main()
