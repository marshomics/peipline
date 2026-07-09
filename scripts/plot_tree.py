#!/usr/bin/env python3
"""Rectangular and circular phylogenies of the C71 sequences.

Drawn from scratch on a matplotlib Axes rather than via ete3/ggtree so that
every element stays a vector path and every label stays live text in the SVG.
Tips are coloured by active-site subgroup and, in a second pair of panels, by
taxonomy. Nodes with UFBoot >= 95 (or the configured cutoff) get a support dot.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from Bio import Phylo

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, palette, savefig, set_style, top_n_with_other  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

SUPPORT_CUTOFF = 95.0


def parse_support(clade):
    """IQ-TREE writes 'SH-aLRT/UFboot' as the internal node label. Take UFboot."""
    for src in (clade.confidence, clade.name):
        if src is None:
            continue
        s = str(src)
        try:
            return float(s.split("/")[-1])
        except ValueError:
            continue
    return None


def _postorder(root):
    """Iterative post-order. Bio.Phylo's own traversals recurse, which blows the
    stack on the caterpillar-shaped trees you get from ladderizing thousands of
    near-identical sequences."""
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


def layout(tree):
    """Return (xs, ys, terminals, internals) keyed by clade id()."""
    tree.ladderize(reverse=True)
    post = _postorder(tree.root)

    terms = [c for c in post if not c.clades]
    internals = [c for c in post if c.clades]

    ys = {}
    i = 0
    for c in post:
        if not c.clades:
            ys[id(c)] = float(i)
            i += 1
    # terminals in the post-order sweep appear left-to-right, so this is the
    # same tip order a recursive draw would give.
    for c in post:                       # children precede parents in post-order
        if c.clades:
            vals = [ys[id(k)] for k in c.clades]
            ys[id(c)] = 0.5 * (min(vals) + max(vals))

    xs = {id(tree.root): 0.0}
    stack = [tree.root]
    while stack:
        node = stack.pop()
        for c in node.clades:
            bl = c.branch_length if c.branch_length is not None else 0.0
            xs[id(c)] = xs[id(node)] + bl
            stack.append(c)

    terms.sort(key=lambda t: ys[id(t)])
    return xs, ys, terms, internals


def _segments(tree, xs, ys, internals):
    """Yield ('h', x0, x1, y) and ('v', x, y0, y1) for a rectangular cladogram."""
    for clade in internals:
        childs = clade.clades
        y0 = min(ys[id(c)] for c in childs)
        y1 = max(ys[id(c)] for c in childs)
        yield ("v", xs[id(clade)], y0, y1)
        for c in childs:
            yield ("h", xs[id(clade)], xs[id(c)], ys[id(c)])
    yield ("h", 0.0, xs[id(tree.root)], ys[id(tree.root)])   # root stub


def draw_rectangular(ax, tree, xs, ys, terms, internals, tip_colors, lw=0.25):
    from matplotlib.collections import LineCollection
    lines = []
    for kind, *v in _segments(tree, xs, ys, internals):
        if kind == "h":
            x0, x1, y = v
            lines.append([(x0, y), (x1, y)])
        else:
            x, y0, y1 = v
            lines.append([(x, y0), (x, y1)])
    ax.add_collection(LineCollection(lines, colors="0.25", linewidths=lw,
                                     rasterized=len(lines) > 4000))

    tx = [xs[id(t)] for t in terms]
    ty = [ys[id(t)] for t in terms]
    tc = [tip_colors.get(t.name, "0.7") for t in terms]
    ax.scatter(tx, ty, s=2.0, c=tc, linewidths=0, zorder=3,
               rasterized=len(terms) > 4000)

    sup = [(xs[id(c)], ys[id(c)]) for c in internals
           if (s := parse_support(c)) is not None and s >= SUPPORT_CUTOFF]
    if sup:
        ax.scatter(*zip(*sup), s=1.2, c="k", linewidths=0, zorder=4,
                   rasterized=len(sup) > 4000)

    ax.set_xlim(-0.02 * max(tx or [1]), 1.06 * max(tx or [1]))
    ax.set_ylim(-1, len(terms))
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Substitutions per site")


def draw_circular(ax, tree, xs, ys, terms, internals, tip_colors, lw=0.25,
                  span=0.97, ring=True):
    from matplotlib.collections import LineCollection
    n = len(terms)
    rmax = max(xs.values()) or 1.0

    def theta(y):
        return 2 * np.pi * span * (y / max(n, 1)) + np.pi / 2

    segs = []
    for kind, *v in _segments(tree, xs, ys, internals):
        if kind == "h":                       # radial
            x0, x1, y = v
            t = theta(y)
            segs.append([(x0 * np.cos(t), x0 * np.sin(t)),
                         (x1 * np.cos(t), x1 * np.sin(t))])
        else:                                  # arc
            x, y0, y1 = v
            t0, t1 = theta(y0), theta(y1)
            k = max(4, int(abs(t1 - t0) * 60))
            tt = np.linspace(t0, t1, k)
            segs.append(list(zip(x * np.cos(tt), x * np.sin(tt))))
    ax.add_collection(LineCollection(segs, colors="0.25", linewidths=lw,
                                     rasterized=len(segs) > 4000))

    tt = np.array([theta(ys[id(t)]) for t in terms])
    rr = np.array([xs[id(t)] for t in terms])
    tc = [tip_colors.get(t.name, "0.7") for t in terms]
    ax.scatter(rr * np.cos(tt), rr * np.sin(tt), s=1.6, c=tc, linewidths=0,
               zorder=3, rasterized=n > 4000)

    if ring:
        r0, r1 = rmax * 1.04, rmax * 1.10
        for t, c in zip(tt, tc):
            dt = 2 * np.pi * span / max(n, 1)
            a = np.linspace(t - dt / 2, t + dt / 2, 3)
            ax.fill(np.concatenate([r0 * np.cos(a), (r1 * np.cos(a))[::-1]]),
                    np.concatenate([r0 * np.sin(a), (r1 * np.sin(a))[::-1]]),
                    color=c, lw=0, zorder=2, rasterized=n > 4000)

    sup = [(xs[id(c)], theta(ys[id(c)])) for c in internals
           if (s := parse_support(c)) is not None and s >= SUPPORT_CUTOFF]
    if sup:
        r, t = np.array([s[0] for s in sup]), np.array([s[1] for s in sup])
        ax.scatter(r * np.cos(t), r * np.sin(t), s=1.0, c="k", linewidths=0, zorder=4,
                   rasterized=len(sup) > 4000)

    # scale bar
    bar = 10 ** np.floor(np.log10(rmax / 3))
    ax.plot([-rmax * 1.15, -rmax * 1.15 + bar], [-rmax * 1.15] * 2, color="k", lw=0.8)
    ax.text(-rmax * 1.15 + bar / 2, -rmax * 1.19, f"{bar:g}", ha="center", va="top",
            fontsize=6)

    lim = rmax * 1.22
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")


def colour_map(labels):
    uniq = [u for u in pd.unique(labels.dropna()) if u is not None]
    uniq = sorted(uniq, key=lambda x: str(x))
    return dict(zip(uniq, palette(len(uniq))))


def one_figure(tree_path, tip_labels, title, out, figdir, fmts, dpi, legend_title):
    tree = Phylo.read(tree_path, "newick")
    try:
        tree.root_at_midpoint()
    except Exception as e:  # noqa: BLE001 - midpoint rooting can fail on odd trees
        print(f"[tree] midpoint rooting failed ({e}); leaving tree as-is", file=sys.stderr)

    xs, ys, terms, internals = layout(tree)
    cmap = colour_map(tip_labels)
    tip_colors = {k: cmap.get(v, "0.7") for k, v in tip_labels.items()}

    handles = [Line2D([0], [0], marker="o", ls="", ms=3, color=c, label=str(k))
               for k, c in cmap.items()]
    handles.append(Line2D([0], [0], marker="o", ls="", ms=2, color="k",
                          label=f"UFBoot $\\geq$ {SUPPORT_CUTOFF:g}"))

    fig, ax = plt.subplots(figsize=(4.0, 5.6))
    draw_rectangular(ax, tree, xs, ys, terms, internals, tip_colors)
    ax.set_title(title, loc="left")
    ax.legend(handles=handles, title=legend_title, loc="upper left",
              bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    fig.tight_layout()
    savefig(fig, figdir, f"{out}_rectangular", fmts, dpi)

    fig, ax = plt.subplots(figsize=(5.6, 5.6))
    draw_circular(ax, tree, xs, ys, terms, internals, tip_colors)
    ax.set_title(title, loc="left")
    ax.legend(handles=handles, title=legend_title, loc="center left",
              bbox_to_anchor=(1.0, 0.5), borderaxespad=0)
    fig.tight_layout()
    savefig(fig, figdir, f"{out}_circular", fmts, dpi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    fmts = tuple(cfg["plots"]["formats"])
    dpi = cfg["plots"]["dpi"]
    set_style()

    assign = pd.read_csv(a.assign, sep="\t", dtype=str)
    assign = assign.set_index("seq_id")

    sub = assign["subgroup"].map(lambda s: f"Subgroup {s}")
    one_figure(a.tree, sub, "C71 phylogeny", "08_tree_by_subgroup", a.figdir, fmts, dpi,
               "Active-site subgroup")

    if "taxon" in assign.columns and assign["taxon"].notna().any():
        tax = top_n_with_other(assign["taxon"], cfg["plots"]["max_taxa_in_legend"])
        one_figure(a.tree, tax, "C71 phylogeny", "09_tree_by_taxonomy", a.figdir, fmts, dpi,
                   cfg["plots"]["taxonomy_rank"].capitalize())
    else:
        print("[tree] no taxon column in the subgroup table; skipping taxonomy tree",
              file=sys.stderr)

    print(f"[tree] figures -> {a.figdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
