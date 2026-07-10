#!/usr/bin/env python3
"""Genotype each genome's pseudomurein biosynthesis pathway.

Pei cleaves the epsilon-isopeptide bond between alanine and lysine. The alanine
is P1, and P1 is not fixed across hosts: it is threonine in Methanobrevibacter
ruminantium, serine in M. stadtmanae, and M. smithii / M. kandleri carry an
ornithine modification of the cross-link. If Pei subgroups mean anything
functional, they should track the substrate the host actually builds.

This screens each genome against a supplied set of Pmur marker HMMs (PmurB,
PmurC, PmurE, the MurE/F-homologous peptide ligases, and whatever else you put
in the directory) and emits a presence/absence genotype vector plus a coarse
call.

Three things this deliberately does not do.

It does not call a genome "Ala-type" because a marker is missing. A marker can
be absent because the pathway is absent, because the MAG is 70% complete, or
because the HMM is bad. Genomes below `pmur_min_markers` are called
`no_pathway_detected`, and the completeness of every genome is carried forward
so the regression can condition on it.

It does not infer P1 chemistry from the marker set, because no marker in the
published set is known to determine P1. The P1 call comes from
`cellwall_reference.py`, which encodes Kandler & Koenig 1978 at the level the
paper actually supports: species, and in two cases a single type strain.

It does not assign P1 from a genus. Methanobrevibacter contains M. ruminantium
(Thr at P1), M. smithii (Ala at P1, ~1/4 Orn at P1') and M. arboriphilus (Ala,
Lys). A genus-level call would put every gut Methanobrevibacter in one bucket
and the bucket would be wrong.

It does not silently use a bit-score cutoff when a model carries a GA line.
"""
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cellwall_reference import CITATION, annotate, genus_of  # noqa: E402
from utils import load_config, parse_domtblout, read_fasta  # noqa: E402

# The same alphabet batch_faa.py enforces. HMMER rejects "*", which Prodigal
# appends to every protein it calls -- and Prodigal calls the archaea, the only
# genomes that can make pseudomurein.
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUO")


def hmm_has_ga(path):
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("HMM "):
                return False
            if line.startswith("GA "):
                return True
    return False


def scan(hmm, faa, out, threshold, threads):
    flags = ["--cut_ga"] if hmm_has_ga(hmm) else \
            ["-T", str(threshold), "--domT", str(threshold),
             "--incT", str(threshold), "--incdomT", str(threshold)]
    subprocess.run(["hmmsearch", "--cpu", str(threads), "--noali", *flags,
                    "--domtblout", out, "-o", os.devnull, hmm, faa], check=True)
    return "cut_ga" if "--cut_ga" in flags else str(threshold)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="unified sample table")
    ap.add_argument("--genomes", required=True, help="genome_level_table.tsv")
    ap.add_argument("--config", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-markers", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--only-with-c71", action="store_true",
                    help="screen only genomes carrying a C71 (much cheaper)")
    a = ap.parse_args()

    cfg = load_config(a.config)["specificity"]
    hmms = sorted(glob.glob(os.path.join(cfg["pmur_hmm_dir"], "*.hmm")))
    if not hmms:
        sys.exit(f"[cellwall] no .hmm in {cfg['pmur_hmm_dir']}")
    markers = [os.path.splitext(os.path.basename(h))[0] for h in hmms]
    print(f"[cellwall] {len(markers)} markers: {markers}", file=sys.stderr)

    os.makedirs(a.workdir, exist_ok=True)
    for p in (a.out, a.out_markers):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    tab = pd.read_csv(a.table, sep="\t", dtype={"sample": str}, low_memory=False)
    gen = pd.read_csv(a.genomes, sep="\t", dtype={"sample": str}, low_memory=False)
    gen = gen[["sample", "has_c71", "completeness", "contamination", "domain"]]
    tab = tab.merge(gen, on="sample", how="inner", suffixes=("", "_g"))

    if a.only_with_c71:
        tab = tab[tab["has_c71"] == 1]
    # Only archaea can make pseudomurein. Screening 340k bacteria for Pmur is a
    # negative control worth running once, not every time.
    print(f"[cellwall] screening {len(tab):,} genomes "
          f"({int((tab['domain'] == 'Archaea').sum()):,} archaea)", file=sys.stderr)
    if tab.empty:
        sys.exit("[cellwall] no genomes to screen")

    # One concatenated faa, sequences prefixed by sample, so we scan once.
    #
    # This used to be a plain `open(faa, errors="replace")` line loop. Two silent
    # failures, both on inputs this pipeline actually produces:
    #
    #   * a gzipped proteome was opened as text, decode errors were SWALLOWED by
    #     errors="replace", no '>' was ever seen, and the gzip bytes were written
    #     out as "sequence". The genotype for that genome was quietly wrong.
    #   * Prodigal's `-a` output appends '*' to every protein. batch_faa.py strips
    #     it because HMMER rejects it; this did not. Prodigal calls the ARCHAEA,
    #     which are the only genomes that can make pseudomurein.
    #
    # Use the same gzip-aware reader and the same alphabet filter as the screen.
    cat = os.path.join(a.workdir, "genomes.faa")
    n_seq, n_missing = 0, 0
    with open(cat, "w") as out:
        for sample, faa in zip(tab["sample"], tab["faa"]):
            if not isinstance(faa, str) or not os.path.exists(faa):
                n_missing += 1
                continue
            if "|" in str(sample):
                sys.exit(f"[cellwall] sample id {sample!r} contains '|', which is "
                         f"the provenance delimiter. Hits would be misattributed.")
            for _hdr, seq in read_fasta(faa):
                clean = "".join(c for c in seq.upper() if c in VALID_AA)
                if not clean:
                    continue
                out.write(f">{sample}|{n_seq}\n{clean}\n")
                n_seq += 1
    if n_missing:
        print(f"[cellwall] {n_missing} genomes have no readable faa and were skipped",
              file=sys.stderr)
    n_screened = len(tab) - n_missing
    print(f"[cellwall] {n_seq:,} proteins from {n_screened:,} genomes", file=sys.stderr)

    hits = defaultdict(set)
    thresholds = {}
    for hmm, m in zip(hmms, markers):
        with tempfile.NamedTemporaryFile(suffix=".domtbl", dir=a.workdir,
                                         delete=False) as tf:
            dom = tf.name
        thresholds[m] = scan(hmm, cat, dom, cfg["pmur_score_threshold"], a.threads)
        n = 0
        for r in parse_domtblout(dom):
            hits[m].add(r["target_name"].split("|", 1)[0])
            n += 1
        os.remove(dom)
        print(f"[cellwall] {m}: {n:,} domain hits in {len(hits[m]):,} genomes "
              f"(threshold {thresholds[m]})", file=sys.stderr)

    pd.DataFrame({"marker": markers,
                  "threshold": [thresholds[m] for m in markers],
                  "n_genomes_hit": [len(hits[m]) for m in markers],
                  "frac_genomes_hit": [len(hits[m]) / len(tab) for m in markers],
                  "has_ga_line": [int(hmm_has_ga(h)) for h in hmms]}
                 ).to_csv(a.out_markers, sep="\t", index=False)

    g = tab[["sample", "domain", "species", "completeness", "contamination",
             "has_c71"]].copy()
    for m in markers:
        g[f"pmur_{m}"] = g["sample"].isin(hits[m]).astype(int)
    g["n_pmur_markers"] = g[[f"pmur_{m}" for m in markers]].sum(axis=1)
    g["pmur_pathway"] = (g["n_pmur_markers"] >= int(cfg["pmur_min_markers"])).astype(int)

    # An absent marker in a 65%-complete MAG is not an absent gene.
    g["pathway_call"] = "no_pathway_detected"
    g.loc[g["pmur_pathway"] == 1, "pathway_call"] = "pseudomurein"
    low = (g["pmur_pathway"] == 0) & (g["completeness"] < 90)
    g.loc[low, "pathway_call"] = "indeterminate_low_completeness"
    print(f"[cellwall] {int(low.sum()):,} genomes are indeterminate because they "
          f"are <90% complete; they are not called negative", file=sys.stderr)

    # Literature wall chemistry, species-level, never inferred from markers and
    # never inherited from a heterogeneous genus. See cellwall_reference.py.
    g = pd.concat([g, annotate(g["species"])], axis=1)
    g["genus"] = g["species"].map(genus_of)
    g["p1_citation"] = CITATION

    g.to_csv(a.out, sep="\t", index=False)
    print("\n[cellwall] pathway calls:", file=sys.stderr)
    print(g["pathway_call"].value_counts().to_string(), file=sys.stderr)
    print(f"\n[cellwall] wall chemistry, from {CITATION}", file=sys.stderr)
    print("\n  P1 (acyl donor; the residue Pei's S1 pocket reads):", file=sys.stderr)
    print(g["p1_residue"].value_counts().to_string(), file=sys.stderr)
    print("\n  P1' (acyl acceptor; Lys epsilon-amine, or Orn delta-amine):",
          file=sys.stderr)
    print(g["p1_prime_residue"].value_counts().to_string(), file=sys.stderr)
    print("\n  provenance:", file=sys.stderr)
    print(g["p1_source"].value_counts().to_string(), file=sys.stderr)

    n_het = int((g["p1_source"] == "genus_heterogeneous").sum())
    if n_het:
        print(f"\n[cellwall] {n_het:,} genomes are in a genus whose characterised "
              f"species disagree about the wall (Methanobrevibacter has Ala/Lys, "
              f"Ala/Orn and Thr/Lys members). They are 'unknown', not guessed.",
              file=sys.stderr)
    disp = g.loc[g["p1_source"] == "disputed", "species"].value_counts()
    if len(disp):
        print(f"\n[cellwall] {int(disp.sum()):,} genomes belong to species whose "
              f"wall chemistry is asserted in secondary sources but not in the "
              f"primary reference this pipeline uses:", file=sys.stderr)
        print(disp.to_string(), file=sys.stderr)
        print("  Supply the primary citation and add them to "
              "cellwall_reference.REFERENCE before relying on them.", file=sys.stderr)
    n_nops = int((g["p1_residue"] == "no_pseudomurein").sum())
    if n_nops:
        print(f"\n[cellwall] {n_nops:,} genomes have no pseudomurein sacculus at "
              f"all (protein sheath). They are a negative control, not Ala-type.",
              file=sys.stderr)

    bac = g[(g["domain"] == "Bacteria") & (g["pmur_pathway"] == 1)]
    if len(bac):
        print(f"\n[cellwall] WARNING: {len(bac)} BACTERIA appear to carry the "
              f"pseudomurein pathway. Bacteria do not make pseudomurein. Either "
              f"the markers cross-react with MurC/MurE (they are homologous), or "
              f"those genomes are contaminated. Check before using this column.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
