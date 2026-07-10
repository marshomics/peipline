#!/usr/bin/env python3
"""Stand-ins for cluster binaries the test sandbox does not have.

These are NOT part of the pipeline. They exist so the pipeline's plumbing can be
exercised end to end without HMMER, MMseqs2, DIAMOND, seqkit or trimAl on PATH.
Each one implements only the exact invocation the pipeline makes.
"""
from __future__ import annotations

import os
import stat
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)


def _write(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


# ---------------------------------------------------------------------------
def seqkit_main(argv):
    """seqkit rmdup -s -o OUT IN"""
    from utils import read_fasta, write_fasta
    assert argv[0] == "rmdup"
    out = argv[argv.index("-o") + 1]
    inp = argv[-1]
    seen, n = set(), 0
    with open(out, "w") as fh:
        for h, s in read_fasta(inp):
            if s in seen:
                continue
            seen.add(s)
            write_fasta(fh, h, s)
            n += 1
    print(f"[stub seqkit] {n} unique", file=sys.stderr)


def mmseqs_main(argv):
    """mmseqs easy-cluster IN PREFIX TMP --min-seq-id X -c Y ..."""
    from utils import read_fasta, seq_id, write_fasta
    if argv[0] == "version":
        print("stub"); return
    assert argv[0] == "easy-cluster"
    inp, pref = argv[1], argv[2]
    minid = float(argv[argv.index("--min-seq-id") + 1])

    recs = [(seq_id(h), s) for h, s in read_fasta(inp)]
    recs.sort(key=lambda r: -len(r[1]))

    # Greedy clustering on 5-mer Jaccard, as a fast stand-in for MMseqs2's
    # identity. Good enough to reproduce cluster structure on the test data.
    def km(s):
        return set(s[i:i + 5] for i in range(len(s) - 4))

    reps, assign = [], {}
    rep_km = []
    for name, seq in recs:
        k = km(seq)
        placed = False
        for (rname, _), rk in zip(reps, rep_km):
            j = len(k & rk) / max(len(k | rk), 1)
            if j >= minid ** 5:            # Jaccard of k-mers ~ identity^k
                assign[name] = rname
                placed = True
                break
        if not placed:
            reps.append((name, seq))
            rep_km.append(k)
            assign[name] = name

    with open(f"{pref}_cluster.tsv", "w") as fh:
        for m, r in assign.items():
            fh.write(f"{r}\t{m}\n")
    with open(f"{pref}_rep_seq.fasta", "w") as fh:
        for name, seq in reps:
            write_fasta(fh, name, seq)
    print(f"[stub mmseqs] {len(recs)} -> {len(reps)} clusters", file=sys.stderr)


def diamond_main(argv):
    """diamond makedb ... | diamond blastp --query Q --db D --out M8 ...

    Scores are a shared-4-mer proxy, not Smith-Waterman: an exact all-vs-all with
    Biopython's aligner takes minutes on even 76 sequences, and the SSN code only
    needs a monotone similarity to turn into an E-value. Real runs use DIAMOND.
    """
    import math
    from utils import read_fasta, seq_id
    if argv[0] == "version":
        print("stub"); return
    if argv[0] == "makedb":
        return
    assert argv[0] == "blastp"
    q = argv[argv.index("--query") + 1]
    out = argv[argv.index("--out") + 1]
    emax = float(argv[argv.index("--evalue") + 1])

    recs = [(seq_id(h), s) for h, s in read_fasta(q)]
    kmers = [set(s[i:i + 4] for i in range(len(s) - 3)) for _, s in recs]
    N = sum(len(s) for _, s in recs)
    LAMB, K = 0.267, 0.041                     # BLOSUM62, gap 11/1

    rows = []
    for i, (ni, si) in enumerate(recs):
        for j in range(i, len(recs)):
            nj, sj = recs[j]
            shared = len(si) if i == j else len(kmers[i] & kmers[j])
            score = 2.0 * shared
            e = K * len(si) * N * math.exp(-LAMB * score)
            if e > emax:
                continue
            pid = 100.0 * shared / max(len(kmers[i] | kmers[j]), 1)
            ln = min(len(si), len(sj))
            rows.append((ni, nj, pid, ln, 0, 0, 1, len(si), 1, len(sj), e, score,
                         len(si), len(sj)))
            if i != j:
                rows.append((nj, ni, pid, ln, 0, 0, 1, len(sj), 1, len(si), e, score,
                             len(sj), len(si)))
    with open(out, "w") as fh:
        for r in rows:
            fh.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x)
                               for x in r) + "\n")
    print(f"[stub diamond] {len(rows)} alignments", file=sys.stderr)


def hmmsearch_main(argv):
    """hmmsearch --cpu N --noali {--cut_ga | -T x --domT x --incT x --incdomT x}
                 --domtblout D --tblout T -o /dev/null HMM SEQDB

    Backed by pyhmmer, which is the HMMER3 code, so the domtblout this writes is
    the real format. Only the argument surface that search_batch.sh uses.
    """
    import pyhmmer
    from pyhmmer.easel import Alphabet, SequenceFile
    from pyhmmer.plan7 import HMMFile

    dom = tbl = None
    T = None
    use_ga = False
    pos = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--domtblout":
            dom = argv[i + 1]; i += 2
        elif a == "--tblout":
            tbl = argv[i + 1]; i += 2
        elif a in ("-o", "--cpu"):
            i += 2
        elif a == "--cut_ga":
            use_ga = True; i += 1
        elif a in ("-T", "--domT", "--incT", "--incdomT"):
            T = float(argv[i + 1]); i += 2
        elif a == "--noali":
            i += 1
        else:
            pos.append(a); i += 1

    hmm_path, db = pos[0], pos[1]
    with HMMFile(hmm_path) as hf:
        hmms = list(hf)
    if use_ga:
        ga = hmms[0].cutoffs.gathering
        if ga is None:
            sys.stderr.write(f"hmmsearch: GA bit thresholds unavailable on {hmm_path}\n")
            sys.exit(1)
        T = float(ga[0])
    if T is None:
        sys.exit("stub hmmsearch: no threshold given")

    abc = Alphabet.amino()
    with SequenceFile(db, digital=True, alphabet=abc) as sf:
        seqs = sf.read_block()

    with open(dom, "wb") as fd, open(tbl, "wb") as ft:
        first = True
        for hits in pyhmmer.hmmsearch(hmms, seqs, T=T, domT=T, incT=T, incdomT=T, cpus=1):
            hits.write(fd, format="domains", header=first)
            hits.write(ft, format="targets", header=first)
            first = False


def hmmalign_main(argv):
    """hmmalign --trim --amino --outformat Stockholm -o OUT <hmm> <faa>

    pyhmmer *is* the HMMER3 code, so this is the real algorithm behind a CLI
    shim rather than an approximation of it. The `#=GC RF` line is what matters:
    groove_map.py and triad_detect_filter.py use the profile match states as
    their coordinate system, and an alignment written without RF silently has no
    coordinates at all.
    """
    import pyhmmer
    from pyhmmer.easel import Alphabet, SequenceFile
    from pyhmmer.plan7 import HMMFile

    out, trim, pos = None, False, []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "-o":
            out, i = argv[i + 1], i + 2
        elif a == "--trim":
            trim, i = True, i + 1
        elif a == "--outformat":
            i += 2
        elif a.startswith("-"):
            i += 1
        else:
            pos.append(a)
            i += 1
    if len(pos) != 2 or out is None:
        sys.exit(f"hmmalign stub: expected '<hmm> <faa>' plus -o, got {argv}")

    with HMMFile(pos[0]) as hf:
        hmm = next(iter(hf))
    with SequenceFile(pos[1], digital=True, alphabet=Alphabet.amino()) as sf:
        seqs = sf.read_block()
    msa = pyhmmer.hmmalign(hmm, seqs, trim=trim)
    msa.name = b"aln"
    with open(out, "wb") as fh:
        msa.write(fh, "stockholm")


DISPATCH = {"seqkit": seqkit_main, "mmseqs": mmseqs_main, "diamond": diamond_main,
            "hmmsearch": hmmsearch_main, "hmmalign": hmmalign_main}


def install(bindir):
    os.makedirs(bindir, exist_ok=True)
    for tool in DISPATCH:
        _write(os.path.join(bindir, tool), textwrap.dedent(f"""\
            #!/usr/bin/env bash
            exec {sys.executable} {os.path.abspath(__file__)} {tool} "$@"
            """))
    # trimAl: our alignments are already trimmed to match columns
    _write(os.path.join(bindir, "trimal"), textwrap.dedent("""\
        #!/usr/bin/env bash
        while [ $# -gt 0 ]; do case "$1" in
          -in) IN=$2; shift 2;; -out) OUT=$2; shift 2;; *) shift;; esac; done
        cp "$IN" "$OUT"
        """))
    # prodigal: pyrodigal's CLI is argument-compatible apart from -q
    _write(os.path.join(bindir, "prodigal"), textwrap.dedent(f"""\
        #!/usr/bin/env bash
        if [ "$1" = "-v" ]; then echo "prodigal stub (pyrodigal)"; exit 0; fi
        ARGS=()
        for a in "$@"; do [ "$a" = "-q" ] || ARGS+=("$a"); done
        exec {sys.executable} -m pyrodigal "${{ARGS[@]}}" 2>/dev/null
        """))
    return bindir


if __name__ == "__main__":
    DISPATCH[sys.argv[1]](sys.argv[2:])
