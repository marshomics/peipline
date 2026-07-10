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
    """Write to .tmp, rename on success, then drop a .done marker.

    Prodigal used to write straight to the final path. The per-genome .faa/.ffn
    are not Snakemake outputs, so a job killed by walltime or OOM left a
    truncated .faa behind; on the rerun the cache check saw a non-empty file,
    called it "cached", and a partial proteome went into hmmsearch. One killed
    chunk out of ~220 would be invisible in a 10-genome test and silent in
    production.

    os.replace is atomic within a filesystem, so anything at the final path is
    complete. `.done` is written last, so a half-renamed pair is detectable too.
    """
    tmp_faa, tmp_ffn, done = faa + ".tmp", ffn + ".tmp", faa + ".done"
    for p in (tmp_faa, tmp_ffn, done):
        if os.path.exists(p):
            os.remove(p)
    cmd = ["prodigal", "-i", fna, "-a", tmp_faa, "-d", tmp_ffn,
           "-p", mode, "-g", str(table), "-m", "-q", "-o", "/dev/null"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode == 0 and os.path.exists(tmp_faa) and os.path.getsize(tmp_faa) > 0:
        os.replace(tmp_faa, faa)
        if os.path.exists(tmp_ffn):
            os.replace(tmp_ffn, ffn)
        with open(done, "w") as fh:
            fh.write("ok\n")
    else:
        for p in (tmp_faa, tmp_ffn):
            if os.path.exists(p):
                os.remove(p)
    return r


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

        # The per-genome .faa/.ffn are not Snakemake outputs, so a job killed by
        # walltime or OOM leaves a TRUNCATED .faa behind. `size > 0` then accepts
        # it as "cached" on the rerun and a partial proteome enters hmmsearch.
        # run_prodigal() now writes to .tmp and os.replace()s on success, so any
        # file at the final path is complete by construction. A stale truncated
        # file from an older run is still possible: `.done` is the marker.
        if os.path.exists(faa) and os.path.getsize(faa) > 0 and \
                os.path.exists(faa + ".done"):
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
