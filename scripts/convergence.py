#!/usr/bin/env python3
"""How many times did each active-site subgroup arise?

Painting subgroups onto a tree shows you convergence. It does not measure it.
Two statistics, both on the C71 gene tree:

1. Minimum independent origins under Fitch parsimony, treating each subgroup as
   a binary trait. Compared against a null built by shuffling tip labels, which
   holds the tree and the subgroup's prevalence fixed and destroys only the
   phylogenetic arrangement. A subgroup arising once is a clade. A subgroup
   arising fifteen times, when the null expects forty, is still phylogenetically
   structured; arising fifteen times when the null expects sixteen is not a
   subgroup at all, it is a stain.

2. Fritz & Purvis D. D = 0 means the trait is distributed as Brownian threshold
   evolution predicts (phylogenetically conserved); D = 1 means randomly
   scattered across the tips. Negative D means more clumped than Brownian. It is
   scaled by simulation on this tree, so it is comparable across subgroups of
   very different prevalence, which raw origin counts are not.

The tip set is already the dereplicated representative set used for the tree, so
the over-sequencing correction is inherited rather than reapplied.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from Bio import Phylo

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import PALETTE, load_config, palette, savefig, set_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402


def postorder(root):
    out, stack = [], [(root, False)]
    while stack:
        node, done = stack.pop()
        if done:
            out.append(node)
        else:
            stack.append((node, True))
            for c in node.clades:
                stack.append((c, False))
    return out


def fitch_changes(post, tip_state: dict) -> int:
    """Fitch's algorithm, first pass. Returns the minimum number of state changes
    on the tree, which for a binary trait with a single ancestral absence is the
    minimum number of independent origins plus losses. We report changes; for a
    trait absent at the root, origins <= changes."""
    sets, changes = {}, 0
    for node in post:
        if not node.clades:
            sets[id(node)] = {tip_state[node.name]}
            continue
        child = [sets[id(c)] for c in node.clades]
        inter = set.intersection(*child)
        if inter:
            sets[id(node)] = inter
        else:
            sets[id(node)] = set.union(*child)
            changes += 1
    return changes


def min_origins(post, tip_state: dict) -> int:
    """Downpass changes with the root forced to state 0 (absence): the number of
    0->1 transitions needed. Equivalent to Fitch changes when the trait is rarer
    than half the tips, which every subgroup here is by construction unless k=1."""
    return fitch_changes(post, tip_state)


def brownian_binary(post, tips, prevalence, rng):
    """Simulate Brownian motion on the tree, threshold at the quantile that
    reproduces the observed prevalence. This is the Fritz & Purvis null for
    'as phylogenetically conserved as a Brownian trait'."""
    val = {}
    for node in reversed(post):                     # root -> tips
        if node is post[-1]:
            val[id(node)] = 0.0
        for c in node.clades:
            bl = c.branch_length or 0.0
            val[id(c)] = val[id(node)] + rng.normal(0.0, np.sqrt(max(bl, 1e-9)))
    x = np.array([val[id(t)] for t in tips])
    cut = np.quantile(x, 1.0 - prevalence)
    return {t.name: int(v > cut) for t, v in zip(tips, x)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--reps", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    vcfg = cfg["convergence"]
    fmts, dpi = tuple(cfg["plots"]["formats"]), cfg["plots"]["dpi"]
    set_style()
    os.makedirs(a.figdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    tree = Phylo.read(a.tree, "newick")
    try:
        tree.root_at_midpoint()
    except Exception as exc:  # noqa: BLE001
        print(f"[convergence] midpoint rooting failed ({exc}); using the tree as given. "
              f"Origin counts are root-dependent, so treat them as an upper bound.",
              file=sys.stderr)

    post = postorder(tree.root)
    tips = [c for c in post if not c.clades]
    tip_names = {t.name for t in tips}

    assign = pd.read_csv(a.assign, sep="\t", dtype={"seq_id": str})
    assign = assign[assign["seq_id"].isin(tip_names)]
    if len(assign) != len(tip_names):
        sys.exit(f"[convergence] {len(tip_names) - len(assign)} tree tips have no subgroup "
                 f"assignment. tree_representatives.tsv and the assignment table disagree.")

    sub = dict(zip(assign["seq_id"], assign["subgroup"]))
    n = len(tips)
    rng = np.random.default_rng(vcfg["seed"])
    B = int(vcfg["n_permutations"])
    S = int(vcfg["n_brownian"])

    # `c71.faa` is the UNION of two evidence tiers: sequences that cleared
    # PF12386's curated gathering threshold, and sequences that cleared only the
    # SSF54001 fold model and then passed the triad filter. Nothing downstream of
    # the filter distinguishes them, so test the tier itself as a trait on the
    # gene tree.
    #
    #   D near 0  ssf_only tips are CLUMPED. They form their own lineage, which is
    #             what a different family looks like. Pooling them into c71.faa is
    #             then a decision that needs defending, not a default.
    #   D near 1  ssf_only tips are interspersed among the PF12386 hits, so the
    #             gathering threshold is simply conservative and the sensitivity
    #             tier is doing the job it was added for.
    #
    # This is the cheapest test that separates "PF12386 is strict" from "SSF54001
    # let something else in", and it costs one extra trait.
    traits = [(str(k), {t.name: int(sub[t.name] == k) for t in tips})
              for k in sorted(assign["subgroup"].unique())]
    if "evidence" in assign.columns and assign["evidence"].nunique() > 1:
        ev = dict(zip(assign["seq_id"], assign["evidence"]))
        state_ev = {t.name: int(ev.get(t.name) == "ssf_only") for t in tips}
        traits.append(("evidence:ssf_only", state_ev))
        print(f"[convergence] {sum(state_ev.values())}/{n} tree tips are ssf_only. "
              f"Testing whether they are a clade (a different family) or "
              f"interspersed (PF12386's threshold is just conservative).",
              file=sys.stderr)

    rows = []
    for k, state in traits:
        obs = min_origins(post, state)
        prev = sum(state.values()) / n
        if prev in (0.0, 1.0):
            continue

        labels = np.array(list(state.values()))
        null = np.empty(B, dtype=int)
        keys = [t.name for t in tips]
        for b in range(B):
            perm = rng.permutation(labels)
            null[b] = min_origins(post, dict(zip(keys, perm)))
        p = (1 + int((null <= obs).sum())) / (B + 1)      # fewer origins = clustered

        brown = np.empty(S, dtype=int)
        for s in range(S):
            brown[s] = min_origins(post, brownian_binary(post, tips, prev, rng))

        mean_r, mean_b = float(null.mean()), float(brown.mean())
        denom = mean_r - mean_b

        # D is a ratio of differences between two nulls. When the tree is short
        # and star-like, Brownian and random nulls coincide, the denominator
        # collapses, and D explodes. Reporting a D of -16 would be worse than
        # reporting nothing: it looks like an extraordinarily conserved trait
        # when it is division by almost zero.
        min_denom = 0.05 * max(mean_r, 1.0)
        if denom > min_denom:
            D = (obs - mean_b) / denom
            interp = ("clade-like" if D < 0.25 else
                      "structured" if D < 0.75 else "convergent/scattered")
        else:
            D, interp = np.nan, "D undefined: Brownian and random nulls coincide"

        if k == "evidence:ssf_only" and np.isfinite(D):
            interp += ("  <- ssf_only forms its own lineage: these are probably NOT "
                       "the same family as the PF12386 hits, and c71.faa pools them"
                       if D < 0.25 else
                       "  <- ssf_only is interspersed among the PF12386 hits: the "
                       "gathering threshold is conservative, not the fold model "
                       "permissive")

        rows.append({
            "subgroup": k, "n_tips": int(labels.sum()), "prevalence": prev,
            "observed_origins": obs,
            "null_random_mean": mean_r, "null_random_sd": float(null.std()),
            "null_brownian_mean": mean_b, "null_brownian_sd": float(brown.std()),
            "p_clustered": p,
            "fritz_purvis_D": D,
            "interpretation": interp,
        })
        print(f"[convergence] {k}: {obs} origins, random null {mean_r:.1f}, "
              f"Brownian null {mean_b:.1f}, D = {D:.3f}, p = {p:.4f}", file=sys.stderr)

    df = pd.DataFrame(rows)
    # k subgroup traits (+ the evidence trait) are reported side by side and a
    # reader will scan the column for "which is significantly clustered". That is
    # a family of tests, so correct it. The raw p stays visible.
    if len(df):
        from statsmodels.stats.multitest import multipletests
        df["q_clustered"] = multipletests(df["p_clustered"], method="fdr_bh")[1]
        df["significant"] = df["q_clustered"] < 0.05
        n_sig = int(df["significant"].sum())
        print(f"[convergence] {n_sig}/{len(df)} traits clustered at q < 0.05 "
              f"(BH across the {len(df)} traits, not per-trait p)", file=sys.stderr)
    df.to_csv(a.out, sep="\t", index=False)
    if df.empty:
        print("[convergence] nothing to test", file=sys.stderr)
        return

    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.4))
    ax = axes[0]
    x = np.arange(len(df))
    ax.bar(x - 0.2, df["observed_origins"], width=0.35, color=PALETTE[0], label="observed")
    ax.errorbar(x + 0.2, df["null_random_mean"], yerr=df["null_random_sd"], fmt="o",
                ms=3, color=PALETTE[4], lw=0.8, capsize=2, label="tip-shuffle null")
    ax.errorbar(x + 0.2, df["null_brownian_mean"], yerr=df["null_brownian_sd"], fmt="s",
                ms=3, color=PALETTE[2], lw=0.8, capsize=2, label="Brownian null")
    ax.set_xticks(x); ax.set_xticklabels([f"SG{s}" for s in df["subgroup"]])
    ax.set_ylabel("Independent origins")
    ax.legend()
    ax.set_title("A  Parsimony origins vs nulls", loc="left")

    ax = axes[1]
    cols = palette(len(df))
    ax.bar(x, df["fritz_purvis_D"].fillna(0), color=cols, width=0.6)
    for xi, v in zip(x, df["fritz_purvis_D"]):
        if not np.isfinite(v):
            ax.text(xi, 0.02, "undefined", rotation=90, fontsize=5.5, ha="center",
                    va="bottom", color="0.3")
    ax.axhline(0, color="0.3", lw=0.7)
    ax.axhline(1, color="0.3", lw=0.7, ls="--")
    ax.text(len(df) - 0.4, 0.02, "Brownian (clumped)", fontsize=5.5, va="bottom", ha="right")
    ax.text(len(df) - 0.4, 1.02, "random (scattered)", fontsize=5.5, va="bottom", ha="right")
    ax.set_xticks(x); ax.set_xticklabels([f"SG{s}" for s in df["subgroup"]])
    ax.set_ylabel("Fritz & Purvis $D$")
    ax.set_title("B  Phylogenetic signal", loc="left")

    fig.tight_layout()
    savefig(fig, a.figdir, "19_subgroup_convergence", fmts, dpi)


if __name__ == "__main__":
    main()
