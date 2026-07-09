#!/usr/bin/env python3
"""Prepare the alignment that goes into IQ-TREE.

Three reductions, each recorded in tree_representatives.tsv so nothing about
the final tree's tip set is a mystery six months from now:

  1. Dereplication of identical alignment rows. Identical rows carry zero
     phylogenetic signal and only cost runtime and zero-length branches.
  2. If still above `tree.max_seqs`, MMseqs2 clustering of the *ungapped*
     sequences at `cluster_min_seq_id`, keeping one representative per cluster.
     ML inference on >5k tips is where IQ-TREE stops being a weekend job.
  3. trimAl column trimming (-gappyout by default).

If you need every tip, set tree.max_seqs very high and tree.mode: fast.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, read_fasta, write_fasta  # noqa: E402


def sh(cmd):
    print("+ " + " ".join(map(str, cmd)), file=sys.stderr)
    subprocess.run(list(map(str, cmd)), check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--afa", required=True)
    ap.add_argument("--out-aln", required=True)
    ap.add_argument("--out-reps", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--tmpdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    cfg = load_config(a.config)["tree"]
    os.makedirs(a.tmpdir, exist_ok=True)
    for p in (a.out_aln, a.out_reps):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    aln = list(read_fasta(a.afa))
    if len(aln) < 4:
        sys.exit(f"[tree_input] only {len(aln)} triad-positive sequences; "
                 f"a phylogeny is not meaningful. Nothing to do.")

    # --- 1. dereplicate identical rows ---------------------------------------
    rep_of: dict[str, str] = {}      # member -> representative
    seen: dict[str, str] = {}        # hash -> representative
    order = []
    if cfg["dereplicate"]:
        for name, seq in aln:
            h = hashlib.blake2b(seq.encode(), digest_size=16).hexdigest()
            if h in seen:
                rep_of[name] = seen[h]
            else:
                seen[h] = name
                rep_of[name] = name
                order.append(name)
        print(f"[tree_input] dereplication: {len(aln)} -> {len(order)} unique rows",
              file=sys.stderr)
    else:
        order = [n for n, _ in aln]
        rep_of = {n: n for n in order}

    seqs = dict(aln)
    derep = os.path.join(a.tmpdir, "derep.afa")
    with open(derep, "w") as fh:
        for n in order:
            write_fasta(fh, n, seqs[n])

    # --- 2. cluster if still too many ----------------------------------------
    keep = order
    clustered = False
    if len(order) > int(cfg["max_seqs"]):
        if not shutil.which("mmseqs"):
            sys.exit("[tree_input] need MMseqs2 to reduce the tip count; not on PATH")
        clustered = True
        ung = os.path.join(a.tmpdir, "derep.ungapped.faa")
        with open(ung, "w") as fh:
            for n in order:
                write_fasta(fh, n, seqs[n].replace("-", ""))

        pref = os.path.join(a.tmpdir, "clu")
        with tempfile.TemporaryDirectory(dir=a.tmpdir) as mm_tmp:
            sh(["mmseqs", "easy-cluster", ung, pref, mm_tmp,
                "--min-seq-id", cfg["cluster_min_seq_id"],
                "-c", cfg["cluster_coverage"], "--cov-mode", 0,
                "--threads", a.threads, "-v", 1])

        # <pref>_cluster.tsv: representative <TAB> member
        member_of = {}
        with open(f"{pref}_cluster.tsv") as fh:
            for line in fh:
                r, m = line.split("\t")[:2]
                member_of[m.strip()] = r.strip()
        for name in list(rep_of):
            rep_of[name] = member_of.get(rep_of[name], rep_of[name])
        keep = sorted(set(member_of.values()))
        print(f"[tree_input] clustered at {cfg['cluster_min_seq_id']} id: "
              f"{len(order)} -> {len(keep)} representatives", file=sys.stderr)

        if len(keep) > int(cfg["max_seqs"]):
            print(f"[tree_input] WARNING: {len(keep)} reps still exceeds "
                  f"max_seqs={cfg['max_seqs']}. Lower cluster_min_seq_id, or accept "
                  f"the runtime.", file=sys.stderr)

    reps_afa = os.path.join(a.tmpdir, "reps.afa")
    keep_set = set(keep)
    with open(reps_afa, "w") as fh:
        for n in order:
            if n in keep_set:
                write_fasta(fh, n, seqs[n])

    # --- 3. trim -------------------------------------------------------------
    if not shutil.which("trimal"):
        sys.exit("[tree_input] trimal not on PATH")
    mode = cfg["trimal_mode"].split()
    sh(["trimal", "-in", reps_afa, "-out", a.out_aln, "-fasta", *mode])

    n_out = sum(1 for _ in read_fasta(a.out_aln))
    if n_out < 4:
        sys.exit(f"[tree_input] trimAl left {n_out} sequences; refusing to build a tree")
    ncol = len(next(iter(read_fasta(a.out_aln)))[1])
    print(f"[tree_input] tree input: {n_out} tips x {ncol} columns", file=sys.stderr)

    # --- provenance ----------------------------------------------------------
    with open(a.out_reps, "w") as fh:
        fh.write("sequence\trepresentative\tis_tip\tdereplicated\tclustered\n")
        for name, _ in aln:
            r = rep_of[name]
            fh.write(f"{name}\t{r}\t{int(name in keep_set)}\t"
                     f"{int(cfg['dereplicate'])}\t{int(clustered)}\n")


if __name__ == "__main__":
    main()
