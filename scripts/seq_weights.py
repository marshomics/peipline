#!/usr/bin/env python3
"""Redundancy weights: w_i = 1 / |cluster(i)| at 90% identity.

Sequence databases are not random samples of protein space. A genus with 8,000
sequenced isolates contributes 8,000 near-identical active sites; a genus with
one contributes one. Every unweighted residue frequency, sequence logo, entropy,
mutual information and consensus in this pipeline would then describe the
sequencing effort rather than the enzyme family.

Weighting each sequence by the reciprocal of its cluster size makes each cluster
contribute exactly 1 to any frequency, which is the same idea as Henikoff
position-based weighting but computed on whole sequences and tied to a threshold
you can state in a methods section. The weights sum to the number of clusters,
so weighted counts read as "effective number of independent sequences".
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, read_fasta, seq_id  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--tmpdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)["redundancy"]
    os.makedirs(a.tmpdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    ids = [seq_id(h) for h, _ in read_fasta(a.faa)]
    if not ids:
        sys.exit("[weights] hits.faa is empty")

    if not cfg.get("enabled", True):
        with open(a.out, "w") as fh:
            fh.write("seq_id\tcluster\tcluster_size\tweight\n")
            for i in ids:
                fh.write(f"{i}\t{i}\t1\t1.0\n")
        print("[weights] weighting disabled; all weights = 1", file=sys.stderr)
        return

    if not shutil.which("mmseqs"):
        sys.exit("[weights] mmseqs not on PATH")

    pref = os.path.join(a.tmpdir, "clu")
    with tempfile.TemporaryDirectory(dir=a.tmpdir) as mm_tmp:
        subprocess.run(
            ["mmseqs", "easy-cluster", a.faa, pref, mm_tmp,
             "--min-seq-id", str(cfg["cluster_min_seq_id"]),
             "-c", str(cfg["cluster_coverage"]), "--cov-mode", "0",
             "--threads", str(a.threads), "-v", "1"], check=True)

    member_of = {}
    with open(f"{pref}_cluster.tsv") as fh:
        for line in fh:
            r, m = line.rstrip("\n").split("\t")[:2]
            member_of[m] = r

    missing = [i for i in ids if i not in member_of]
    if missing:
        sys.exit(f"[weights] {len(missing)} sequences absent from the MMseqs2 clustering, "
                 f"e.g. {missing[:5]}")

    size = Counter(member_of[i] for i in ids)
    with open(a.out, "w") as fh:
        fh.write("seq_id\tcluster\tcluster_size\tweight\n")
        for i in ids:
            c = member_of[i]
            fh.write(f"{i}\t{c}\t{size[c]}\t{1.0 / size[c]:.8f}\n")

    print(f"[weights] {len(ids)} sequences -> {len(size)} clusters at "
          f"{cfg['cluster_min_seq_id']} identity; effective n = {len(size)}",
          file=sys.stderr)
    top = size.most_common(3)
    print(f"[weights] largest clusters: {top}", file=sys.stderr)


if __name__ == "__main__":
    main()
