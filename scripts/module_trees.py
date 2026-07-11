#!/usr/bin/env python3
"""Do the binding module and the catalytic module share a history?

If Pei retargets by swapping its pseudomurein-binding repeat array rather than
by retuning its catalytic groove, then the PMBR tree and the C71 tree will
disagree, and the disagreement will not be noise. That is a testable statement
and this script tests it three ways:

  Mantel   correlation of the two patristic distance matrices over shared tips,
           with a permutation null. Answers "are the two modules evolving at
           proportional rates over the same tips".
  RF       normalised Robinson-Foulds between the two topologies, against a null
           built from random trees on the same tip set. Answers "do they induce
           the same splits".
  cophylo  a tanglegram, because a number that nobody can look at is a number
           nobody will believe.

Both distances are computed on the SAME tips, from the SAME proteins, so a
difference cannot be a sampling artefact. Tips are the redundancy-dereplicated
representatives, so an over-sequenced genus cannot manufacture congruence.

A caveat that has to travel with the result: two modules of one protein are
physically linked, so they share a gene tree unless recombination has separated
them. Congruence is the null. Incongruence is the finding. The RF null therefore
uses random trees, which is conservative in the right direction.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from Bio import Phylo

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (PALETTE, load_config, read_fasta, savefig, seq_id,  # noqa: E402
                   set_style, write_fasta)

import matplotlib.pyplot as plt  # noqa: E402


def sh(cmd, **kw):
    print("+ " + " ".join(map(str, cmd)), file=sys.stderr)
    subprocess.run(list(map(str, cmd)), check=True, **kw)


def extract_module(faa, arch, column, out_faa, min_len=40):
    """Slice each protein down to one module, using the envelope coordinates
    hmmscan reported. Repeats are concatenated in N->C order."""
    spans = dict(zip(arch["seq_id"], arch[column].fillna("")))
    n = 0
    with open(out_faa, "w") as fh:
        for h, s in read_fasta(faa):
            sid = seq_id(h)
            sp = spans.get(sid, "")
            if not sp:
                continue
            frag = "".join(s[int(a) - 1:int(b)] for a, b in
                           (x.split("-") for x in sp.split(";")))
            if len(frag) >= min_len:
                write_fasta(fh, sid, frag)
                n += 1
    return n


def build_tree(faa, prefix, mode, threads, seed):
    aln = prefix + ".aln"
    sh(["mafft", "--auto", "--anysymbol", "--thread", threads, faa],
       stdout=open(aln, "w"))
    trimmed = prefix + ".trim.aln"
    sh(["trimal", "-in", aln, "-out", trimmed, "-fasta", "-gappyout"])
    if mode == "full":
        sh(["iqtree2", "-s", trimmed, "--prefix", prefix, "--seqtype", "AA",
            "-m", "MFP", "-B", "1000", "-T", threads, "-seed", seed, "-redo"])
    else:
        sh(["iqtree2", "-s", trimmed, "--prefix", prefix, "--seqtype", "AA",
            "-m", "LG+F+G4", "--fast", "-T", threads, "-seed", seed, "-redo"])
    return prefix + ".treefile"


def patristic(tree, tips):
    """Pairwise patristic distances. Bio.Phylo's distance() is O(n) per pair, so
    do it once via depth-to-root plus MRCA depth."""
    idx = {t: i for i, t in enumerate(tips)}
    depths = tree.depths(unit_branch_lengths=False)
    term = {t.name: t for t in tree.get_terminals() if t.name in idx}
    D = np.zeros((len(tips), len(tips)))
    for i, a in enumerate(tips):
        for j in range(i + 1, len(tips)):
            b = tips[j]
            mrca = tree.common_ancestor(term[a], term[b])
            d = depths[term[a]] + depths[term[b]] - 2 * depths[mrca]
            D[i, j] = D[j, i] = d
    return D


def mantel(A, B, n_perm, rng):
    iu = np.triu_indices_from(A, k=1)
    a, b = A[iu], B[iu]
    a = (a - a.mean()) / (a.std() or 1)
    b = (b - b.mean()) / (b.std() or 1)
    r = float((a * b).mean())
    n = A.shape[0]
    ge = 0
    for _ in range(n_perm):
        p = rng.permutation(n)
        bp = B[np.ix_(p, p)][iu]
        bp = (bp - bp.mean()) / (bp.std() or 1)
        if float((a * bp).mean()) >= r:
            ge += 1
    return r, (ge + 1) / (n_perm + 1)


def splits(tree, tips):
    """Set of non-trivial bipartitions, each as a frozenset of the smaller side."""
    tipset = set(tips)
    out = set()
    for cl in tree.get_nonterminals():
        s = frozenset(t.name for t in cl.get_terminals()) & tipset
        if 1 < len(s) < len(tipset) - 1:
            out.add(min(s, frozenset(tipset - s), key=lambda x: (len(x), sorted(x))))
    return out


def rf_distance(t1, t2, tips):
    s1, s2 = splits(t1, tips), splits(t2, tips)
    return len(s1 ^ s2), len(s1), len(s2)


def rf_null(t1, t2, tips, n_perm, rng):
    """Distribution of RF between the two trees when one has its tip labels shuffled.

    This is the null the docstring always promised and the code never computed.
    It is the honest question: is the observed RF smaller than what you would get
    if the two topologies were unrelated but had the same shape and the same tip
    set?

    It replaces the Mantel test, which was doing the significance work. A Mantel
    test on two patristic distance matrices over the same taxa has badly inflated
    type I error: the matrices share divergence-time structure whatever the
    topologies do, and permuting rows and columns does not restore exchangeability
    under that autocorrelation. Congruence was the easy call, and congruence is
    the boring answer -- so the inflated test biased the pipeline away from
    detecting the module swapping it exists to find.
    """
    tips = list(tips)
    obs, _, _ = rf_distance(t1, t2, set(tips))
    # Relabelling the tips of a fixed topology samples uniformly from the
    # topologies with that shape, which is exactly "unrelated but same shape".
    original = {}
    for term in t2.get_terminals():
        original[id(term)] = term.name
    terms = list(t2.get_terminals())
    null = np.empty(n_perm, dtype=int)
    names = [original[id(x)] for x in terms]
    for b in range(n_perm):
        perm = rng.permutation(len(names))
        for term, j in zip(terms, perm):
            term.name = names[j]
        null[b], _, _ = rf_distance(t1, t2, set(tips))
    for term in terms:                       # restore, this tree is used later
        term.name = original[id(term)]
    # fewer splits differing = more congruent
    p = (1 + int((null <= obs).sum())) / (n_perm + 1)
    return obs, null, p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True)
    ap.add_argument("--arch", required=True)
    ap.add_argument("--reps", required=True, help="tree_representatives.tsv")
    ap.add_argument("--config", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    cfg = load_config(a.config)
    scfg = cfg["specificity"]
    fmts, dpi = tuple(cfg["plots"]["formats"]), cfg["plots"]["dpi"]
    set_style()
    os.makedirs(a.workdir, exist_ok=True)
    os.makedirs(a.figdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    arch = pd.read_csv(a.arch, sep="\t", dtype={"seq_id": str})
    reps = pd.read_csv(a.reps, sep="\t")
    keep = set(reps.loc[reps["is_tip"] == 1, "sequence"])

    both = arch[(arch["pmbr_span"].notna()) & (arch["pmbr_span"] != "") &
                (arch["catalytic_span"].notna()) & (arch["catalytic_span"] != "")]
    both = both[both["seq_id"].isin(keep)]
    n = len(both)
    print(f"[module_trees] {n} dereplicated proteins carry BOTH modules", file=sys.stderr)

    result = {"n_proteins_both_modules": n}
    if n < 20:
        print("[module_trees] too few proteins with both modules for a congruence "
              "test; writing an empty result rather than a meaningless one",
              file=sys.stderr)
        pd.Series(result).rename("value").rename_axis("metric").to_csv(
            a.out, sep="\t")
        return

    if n > scfg["module_max_tips"]:
        both = both.sample(scfg["module_max_tips"], random_state=1)
        print(f"[module_trees] subsampled to {len(both)} tips", file=sys.stderr)

    sub = os.path.join(a.workdir, "both.faa")
    ids = set(both["seq_id"])
    with open(sub, "w") as fh:
        for h, s in read_fasta(a.faa):
            if seq_id(h) in ids:
                write_fasta(fh, seq_id(h), s)

    trees = {}
    for mod, col in (("catalytic", "catalytic_span"), ("pmbr", "pmbr_span")):
        f = os.path.join(a.workdir, f"{mod}.faa")
        k = extract_module(sub, both, col, f)
        print(f"[module_trees] {mod}: {k} sequences", file=sys.stderr)
        trees[mod] = build_tree(f, os.path.join(a.workdir, mod),
                                scfg["module_tree_mode"], a.threads,
                                cfg["tree"]["seed"])

    t1 = Phylo.read(trees["catalytic"], "newick")
    t2 = Phylo.read(trees["pmbr"], "newick")
    for t in (t1, t2):
        try:
            t.root_at_midpoint()
        except Exception:  # noqa: BLE001
            pass

    tips = sorted({x.name for x in t1.get_terminals()} &
                  {x.name for x in t2.get_terminals()})
    print(f"[module_trees] {len(tips)} shared tips", file=sys.stderr)
    if len(tips) < 20:
        sys.exit("[module_trees] fewer than 20 shared tips")

    rng = np.random.default_rng(cfg["tree"]["seed"])

    # Mantel r is kept as a DESCRIPTIVE effect size. Its p-value is not used and
    # is not reported: a Mantel test on two patristic matrices over the same taxa
    # has inflated type I error, because the matrices share divergence-time
    # structure regardless of topology and the row/column permutation does not
    # restore exchangeability under that autocorrelation.
    D1, D2 = patristic(t1, tips), patristic(t2, tips)
    r, _p_mantel_unused = mantel(D1, D2, 99, rng)

    # The significance test is on the topology, against the tip-shuffle null the
    # docstring always described.
    n_perm = int(scfg.get("rf_permutations", scfg["mantel_permutations"]))
    n_perm = min(n_perm, 999)              # each replicate re-derives all splits
    rf, n1, n2 = rf_distance(t1, t2, tips)
    rf_max = n1 + n2
    rf_norm = rf / rf_max if rf_max else np.nan
    _obs, rf_nulldist, p_rf = rf_null(t1, t2, tips, n_perm, rng)
    rf_null_mean = float(rf_nulldist.mean())
    rf_z = ((rf - rf_null_mean) / rf_nulldist.std()
            if rf_nulldist.std() > 0 else np.nan)

    # Congruence is now a claim about the null, not about a hardcoded 0.3.
    if p_rf < 0.05:
        interp = (f"modules share a history: RF {rf} is smaller than "
                  f"{100 * (1 - p_rf):.0f}% of tip-shuffled topologies "
                  f"(null mean {rf_null_mean:.0f})")
    elif rf >= rf_null_mean:
        interp = ("modules are INCONGRUENT: the observed RF is no smaller than a "
                  "tip-shuffled null, which is what binding-module swapping looks "
                  "like")
    else:
        interp = ("ambiguous: RF is below the null mean but not significantly so. "
                  "More tips, or a bootstrap-based test, would settle it")

    result.update({
        "n_shared_tips": len(tips),
        "mantel_r_descriptive_only": r,
        "robinson_foulds": rf,
        "rf_max": rf_max,
        "rf_normalised": rf_norm,
        "rf_null_mean": rf_null_mean,
        "rf_z": rf_z,
        "rf_p_tip_shuffle": p_rf,
        "rf_permutations": n_perm,
        "catalytic_splits": n1,
        "pmbr_splits": n2,
        "interpretation": interp,
    })
    pd.Series(result).rename("value").rename_axis("metric").to_csv(a.out, sep="\t")
    for k, v in result.items():
        print(f"  {k}: {v}", file=sys.stderr)

    # --- tanglegram ----------------------------------------------------------
    def order(tree):
        tree.ladderize()
        return [t.name for t in tree.get_terminals() if t.name in set(tips)]

    o1, o2 = order(t1), order(t2)
    y1 = {t: i for i, t in enumerate(o1)}
    y2 = {t: i for i, t in enumerate(o2)}

    fig, ax = plt.subplots(figsize=(5.6, max(2.4, 0.05 * len(tips))))
    arch_cls = dict(zip(both["seq_id"], both["architecture_class"]))
    classes = sorted(set(arch_cls.values()))
    cmap = dict(zip(classes, PALETTE))
    for t in tips:
        ax.plot([0, 1], [y1[t], y2[t]], lw=0.3,
                color=cmap.get(arch_cls.get(t), "0.7"), alpha=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["C71 catalytic tree", "PMBR tree"])
    ax.set_yticks([])
    # Mantel r is descriptive only (its permutation p is inflated on patristic
    # matrices and is deliberately NOT reported); significance is the tip-shuffle
    # RF null. `p` was never in scope here -- this line used to crash.
    ax.set_title(f"Mantel r = {r:.3f} (descriptive); norm. RF = {rf_norm:.3f}, "
                 f"tip-shuffle p = {p_rf:.4g}", loc="left")
    for c in classes:
        ax.plot([], [], color=cmap[c], lw=1.2, label=c)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    savefig(fig, a.figdir, "23_module_tanglegram", fmts, dpi)


if __name__ == "__main__":
    main()
