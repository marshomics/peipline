#!/usr/bin/env python3
"""Build a miniature dataset with a known answer, structured like the real one.

Contents
  60 bacterial proteomes (.faa), as provided files
  20 archaeal assemblies (.fna), for Prodigal to call genes on
  a PF12386-like model, WITH a GA line, built from the Pei family only
  an SSF54001-like model, WITHOUT a GA line, built from a broad fold sample

Two protein families share one structural template:
  Pei family   catalytic triad planted at match columns 112 / 147 / 171,
               three active-site subgroups with distinct flanking motifs
  fold family  same fold, 55% diverged, triad at 120 / 160 / 200 instead

The fold family is the point of the test. It clears SSF54001 at score 25 (it is
a cysteine-protease fold) but not PF12386's gathering threshold, so it must land
in the `ssf_only` tier, must not contribute to the triad column call, and must
almost entirely fail the triad filter. Sequences with a deliberately broken
catalytic residue (C->S, H->Y, D->E) must fail too.
"""
from __future__ import annotations

import os
import random
import sys

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
DOMAIN_LEN = 300
TRIAD = {"C": 112, "H": 147, "D": 171}
FOLD_TRIAD = {"C": 120, "H": 160, "D": 200}

SUBGROUPS = {
    0: {("C", -3): "G", ("C", -2): "S", ("C", 2): "W", ("H", -1): "N", ("D", 1): "Y"},
    1: {("C", -3): "A", ("C", -2): "T", ("C", 2): "F", ("H", -1): "Q", ("D", 1): "L"},
    2: {("C", -3): "P", ("C", -2): "V", ("C", 2): "Y", ("H", -1): "K", ("D", 1): "M"},
}

PHYLA = ["Bacillota", "Pseudomonadota", "Actinomycetota", "Bacteroidota"]
AR_PHYLA = ["Methanobacteriota", "Halobacteriota", "Thermoproteota"]

def _codon_table():
    """All synonymous codons per residue. A one-codon-per-residue table makes the
    DNA so low-complexity that Prodigal trains on nonsense and calls genes on the
    reverse strand."""
    import collections
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    t = collections.defaultdict(list)
    for c, a in zip(codons, aas):
        if a != "*":
            t[a].append(c)
    return t


CODON = _codon_table()


def rand_prot(n, rng):
    return "".join(rng.choice(list(AA), int(n)))


def mutate(template, rng, mut):
    s = list(template)
    for i in range(len(s)):
        if rng.random() < mut:
            s[i] = rng.choice(list(AA))
    return s


def make_pei(template, rng, subgroup, mut=0.30, break_residue=None):
    s = mutate(template, rng, mut)
    for r, pos in TRIAD.items():
        s[pos] = r
    for (r, off), aa in SUBGROUPS[subgroup].items():
        s[TRIAD[r] + off] = aa
    if break_residue:
        s[TRIAD[break_residue]] = {"C": "S", "H": "Y", "D": "E"}[break_residue]
    return "".join(s)


def make_fold(template, rng, mut=0.55):
    """Same fold, different family: triad in different places."""
    s = mutate(template, rng, mut)
    for r, pos in FOLD_TRIAD.items():
        s[pos] = r
    for pos in TRIAD.values():                 # make sure it is NOT a Pei triad
        if s[pos] in "CHD":
            s[pos] = rng.choice(list("AGILVST"))
    return "".join(s)


def to_dna(prot, rng):
    return "ATG" + "".join(rng.choice(CODON[c]) for c in prot[1:]) + "TAA"


def main() -> None:
    out = sys.argv[1] if len(sys.argv) > 1 else "/tmp/c71_test"
    rng = np.random.default_rng(7)
    random.seed(7)

    faa_dir = os.path.join(out, "faa")
    fna_dir = os.path.join(out, "gtdb_archaea", "seqs")
    hmm_dir = os.path.join(out, "hmms")
    for d in (faa_dir, fna_dir, hmm_dir):
        os.makedirs(d, exist_ok=True)

    template = rand_prot(DOMAIN_LEN, rng)
    truth, bac_rows = [], []

    # ---- bacteria: provided .faa -------------------------------------------
    n_bac = 60
    for si in range(n_bac):
        sample = f"SAMPLE_{si:03d}"
        path = os.path.join(faa_dir, f"{sample}.faa")
        bac_rows.append((sample, path, PHYLA[si % len(PHYLA)]))
        with open(path, "w") as fh:
            for pi in range(25):
                fh.write(f">{sample}_prot{pi:04d} hypothetical protein\n"
                         f"{rand_prot(rng.integers(120, 400), rng)}\n")

            for ri in range(int(rng.choice([0, 1, 1, 2], p=[0.15, 0.55, 0.20, 0.10]))):
                sg = int(rng.integers(0, 3))
                u = rng.random()
                brk = "D" if u < 0.08 else "C" if u < 0.12 else "H" if u < 0.15 else None
                dom = make_pei(template, rng, sg, break_residue=brk)
                seq = rand_prot(rng.integers(10, 40), rng) + dom + rand_prot(rng.integers(10, 40), rng)
                pid = f"{sample}_pei_{ri}"
                fh.write(f">{pid} peptidase\n{seq}\n")
                truth.append((sample, pid, "pei", sg, brk or "none"))

            # the fold family: SSF54001-positive, PF12386-negative
            for fi in range(int(rng.choice([0, 1, 2], p=[0.45, 0.40, 0.15]))):
                dom = make_fold(template, rng)
                seq = rand_prot(rng.integers(10, 40), rng) + dom + rand_prot(rng.integers(10, 40), rng)
                pid = f"{sample}_fold_{fi}"
                fh.write(f">{pid} cysteine protease\n{seq}\n")
                truth.append((sample, pid, "fold", -1, "none"))

    # a missing and an empty file, to exercise allow_missing_faa
    bac_rows.append(("SAMPLE_MISSING", os.path.join(faa_dir, "nope.faa"), "Bacillota"))
    empty = os.path.join(faa_dir, "SAMPLE_EMPTY.faa")
    open(empty, "w").close()
    bac_rows.append(("SAMPLE_EMPTY", empty, "Bacillota"))

    with open(os.path.join(out, "faa_sample_table.tsv"), "w") as fh:
        fh.write("sample\tfaa\n")
        for s, p, _ in bac_rows:
            fh.write(f"{s}\t{p}\n")

    # Mirrors the real table: gtdb_* rank columns (no lineage string), a bare
    # `phylum` host-organism column that is always empty, an N50 column misnamed
    # ctg_L50, and 10% archaea hiding among the "bacterial" proteomes.
    bac_species = {}
    with open(os.path.join(out, "bacteria_qc.tsv"), "w") as fh:
        fh.write("genome_id\tgtdb_domain\tgtdb_phylum\tgtdb_species\tphylum\t"
                 "Completeness\tContamination\tctg_L50\tn_contigs\n")
        for i, (s, _, ph) in enumerate(bac_rows):
            arch = (i % 10 == 3)
            dom = "d__Archaea" if arch else "d__Bacteria"
            phy = "p__Methanobacteriota" if arch else f"p__{ph}"
            sp = f"s__Genus{i % 12} species{i % 12}"
            bac_species[s] = (dom, sp)
            fh.write(f"{s}\t{dom}\t{phy}\t{sp}\t\t"
                     f"{rng.uniform(70, 100):.2f}\t{rng.uniform(0, 5):.2f}\t"
                     f"{int(rng.integers(5_000, 250_000))}\t{int(rng.integers(1, 400))}\n")

    # ---- archaea: .fna for prodigal ----------------------------------------
    n_arc = 20
    ar_rows = []
    for ai in range(n_arc):
        acc = f"GCA_{900000000 + ai}.1"
        path = os.path.join(fna_dir, f"{acc}.fna")
        ph = AR_PHYLA[ai % len(AR_PHYLA)]
        ar_rows.append((acc, ph))

        genes = ["M" + rand_prot(rng.integers(150, 300), rng) for _ in range(40)]
        n_real = int(rng.choice([0, 1, 2], p=[0.2, 0.5, 0.3]))
        for ri in range(n_real):
            sg = int(rng.integers(0, 3))
            brk = "D" if rng.random() < 0.10 else None
            dom = make_pei(template, rng, sg, break_residue=brk)
            # embed the domain inside a longer ORF: Prodigal is free to pick a
            # different start codon, as it does on real assemblies
            genes.append("M" + rand_prot(25, rng) + dom + rand_prot(25, rng))
            truth.append((acc, f"__prodigal_pei_{ri}", "pei", sg, brk or "none"))
        if rng.random() < 0.5:
            genes.append("M" + rand_prot(20, rng) + make_fold(template, rng)
                         + rand_prot(20, rng))

        rng.shuffle(genes)
        dna = []
        for g in genes:
            # AT-rich intergenic spacer, then a Shine-Dalgarno-like motif
            dna.append("".join(rng.choice(list("ATATGCAT"), int(rng.integers(60, 150)))))
            dna.append("AGGAGGT" + "".join(rng.choice(list("ACGT"), 6)))
            dna.append(to_dna(g, rng))
        seq = "".join(dna)
        with open(path, "w") as fh:
            fh.write(f">{acc}_contig1\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    # Archaeal metadata: species reps flagged, so non-reps inherit a rep's tip.
    ar_species, ar_reps = {}, {}
    with open(os.path.join(out, "ar53_metadata.tsv"), "w") as fh:
        fh.write("accession\tcheckm2_completeness\tcheckm2_contamination\tn50_contigs\t"
                 "contig_count\tgtdb_taxonomy\tgtdb_representative\n")
        for i, (acc, ph) in enumerate(ar_rows):
            sp = f"s__Archaeon{i % 6} sp{i % 6}"
            isrep = sp not in ar_reps
            if isrep:
                ar_reps[sp] = f"GB_{acc}"
            ar_species[acc] = sp
            fh.write(f"GB_{acc}\t{rng.uniform(80, 100):.2f}\t{rng.uniform(0, 3):.2f}\t"
                     f"{int(rng.integers(20_000, 900_000))}\t{int(rng.integers(1, 60))}\t"
                     f"d__Archaea;p__{ph};c__X;o__Y;f__Z;g__A{i % 5};{sp}\t"
                     f"{'t' if isrep else 'f'}\n")
        # the stray archaea inside the bacterial table get species too
        for i, (s_, (d_, sp_)) in enumerate(bac_species.items()):
            if d_ != "d__Archaea":
                continue
            acc = f"GCA_{800000000 + i}.1"
            if sp_ not in ar_reps:
                ar_reps[sp_] = f"GB_{acc}"
                fh.write(f"GB_{acc}\t95.0\t1.0\t100000\t20\t"
                         f"d__Archaea;p__Methanobacteriota;c__X;o__Y;f__Z;g__M;{sp_}\tt\n")

    # ---- trees --------------------------------------------------------------
    def random_newick(tips, rng):
        nodes = [f"{t}:{rng.uniform(0.01, 0.2):.4f}" for t in tips]
        while len(nodes) > 1:
            i, j = sorted(rng.choice(len(nodes), 2, replace=False))
            b = nodes.pop(j); a = nodes.pop(i)
            nodes.append(f"({a},{b}):{rng.uniform(0.01, 0.2):.4f}")
        return nodes[0].rsplit(":", 1)[0] + ";"

    bac_tips = sorted({sp for (d, sp) in bac_species.values() if d == "d__Bacteria"})
    bac_tips.append("s__")          # GTDB-Tk really does emit this; must be ignored
    with open(os.path.join(out, "bacteria.tree"), "w") as fh:
        fh.write(random_newick(bac_tips, rng) + "\n")
    with open(os.path.join(out, "ar53.tree"), "w") as fh:
        fh.write(random_newick(sorted(set(ar_reps.values())), rng) + "\n")

    with open(os.path.join(out, "truth.tsv"), "w") as fh:
        fh.write("sample\tprotein_id\tfamily\tsubgroup\tbroken_residue\n")
        for t in truth:
            fh.write("\t".join(map(str, t)) + "\n")

    # ---- profiles -----------------------------------------------------------
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextMSA, TextSequence

    abc = Alphabet.amino()

    def build(name, seqs):
        msa = TextMSA(name=name.encode(),
                      sequences=[TextSequence(name=f"t{i}".encode(), sequence=s)
                                 for i, s in enumerate(seqs)]).digitize(abc)
        hmm, _, _ = pyhmmer.plan7.Builder(abc).build_msa(
            msa, pyhmmer.plan7.Background(abc))
        hmm.name = name.encode()
        hmm.accession = name.encode()
        return hmm

    pei_train = [make_pei(template, rng, i % 3, mut=0.30) for i in range(60)]
    fold_train = [make_fold(template, rng) for i in range(30)]

    pf = build("Peptidase_C71", pei_train)
    ssf = build("SSF54001", pei_train[:30] + fold_train)

    # Choose PF12386's gathering threshold empirically: above the best score any
    # fold-family sequence achieves, below the worst score a real Pei achieves.
    def scores(hmm, seqs):
        blk = pyhmmer.easel.DigitalSequenceBlock(
            abc, [TextSequence(name=f"q{i}".encode(), sequence=s).digitize(abc)
                  for i, s in enumerate(seqs)])
        best = {}
        for hits in pyhmmer.hmmsearch([hmm], blk, T=0.0, domT=0.0, incT=0.0, incdomT=0.0):
            for h in hits:
                nm = h.name.decode() if isinstance(h.name, bytes) else str(h.name)
                best[nm] = h.score
        return np.array([best.get(f"q{i}", 0.0) for i in range(len(seqs))])

    pei_probe = [make_pei(template, rng, i % 3) for i in range(40)]
    fold_probe = [make_fold(template, rng) for i in range(40)]
    s_pei, s_fold = scores(pf, pei_probe), scores(pf, fold_probe)
    ga = float(np.floor((s_fold.max() + s_pei.min()) / 2))
    if not (s_fold.max() < ga < s_pei.min()):
        ga = float(np.floor(s_fold.max() + 1))
    print(f"[testdata] PF12386 scores: Pei {s_pei.min():.0f}-{s_pei.max():.0f}, "
          f"fold {s_fold.min():.0f}-{s_fold.max():.0f}  ->  GA = {ga:.0f}")
    pf.cutoffs.gathering = (ga, ga)

    with open(os.path.join(hmm_dir, "PF12386.hmm"), "wb") as fh:
        pf.write(fh)
    with open(os.path.join(hmm_dir, "SSF54001.hmm"), "wb") as fh:
        ssf.write(fh)      # no GA line, by design

    n_pei = sum(1 for t in truth if t[2] == "pei")
    print(f"[testdata] {n_bac} bacteria, {n_arc} archaea")
    print(f"[testdata] {n_pei} Pei proteins ({sum(1 for t in truth if t[2]=='pei' and t[4]=='none')} intact), "
          f"{sum(1 for t in truth if t[2] == 'fold')} fold-family decoys")
    print(f"[testdata] ground-truth triad columns: {TRIAD}")
    print(f"[testdata] -> {out}")


if __name__ == "__main__":
    main()
