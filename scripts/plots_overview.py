#!/usr/bin/env python3
"""Overview figures: score distributions, profile agreement, coverage,
prevalence across samples, taxonomy composition, and the alignment
conservation profile with the catalytic triad marked."""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (PALETTE, load_config, palette, resolve_taxonomy,  # noqa: E402
                   savefig, set_style, top_n_with_other)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def fig_score_distributions(hits, figdir, fmts, dpi):
    profiles = sorted(hits["profile"].unique())
    cols = palette(len(profiles))
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.2))

    ax = axes[0]
    for p, c in zip(profiles, cols):
        s = hits.loc[hits["profile"] == p, "full_score"]
        ax.hist(s, bins=60, histtype="step", lw=1.0, color=c, label=f"{p} (n={len(s):,})")
    ax.set_xlabel("Full-sequence bit score")
    ax.set_ylabel("Proteins")
    ax.set_yscale("log")
    ax.legend(loc="upper right")
    ax.set_title("A  Score distribution", loc="left")

    ax = axes[1]
    for p, c in zip(profiles, cols):
        s = np.sort(hits.loc[hits["profile"] == p, "full_score"].to_numpy())
        ax.plot(s, 1 - np.arange(len(s)) / len(s), color=c, lw=1.0, label=p)
    ax.set_xlabel("Full-sequence bit score")
    ax.set_ylabel("Fraction of hits $\\geq$ score")
    ax.set_yscale("log")
    ax.set_title("B  Survival", loc="left")

    ax = axes[2]
    e = hits["full_evalue"].replace(0, np.nextafter(0, 1))
    for p, c in zip(profiles, cols):
        m = hits["profile"] == p
        ax.scatter(-np.log10(e[m]), hits.loc[m, "full_score"], s=1.5, alpha=0.25,
                   color=c, rasterized=True, linewidths=0)
    ax.set_xlabel("$-\\log_{10}$ E-value")
    ax.set_ylabel("Bit score")
    ax.set_title("C  Score vs E-value", loc="left")

    fig.tight_layout()
    savefig(fig, figdir, "01_score_distributions", fmts, dpi)


def fig_decoy_fdr(fdr_path, figdir, fmts, dpi):
    """Reversed-decoy false-discovery rate. This is the defence of the cutoffs."""
    fdr = pd.read_csv(fdr_path, sep="\t")
    if fdr.empty:
        print("[plots] decoy FDR disabled; skipping", file=sys.stderr)
        return
    profiles = sorted(fdr["profile"].unique())
    cols = palette(len(profiles))
    fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.2))

    ax = axes[0]
    for p, c in zip(profiles, cols):
        d = fdr[fdr["profile"] == p]
        ax.plot(d["bit_score"], d["fdr"].clip(lower=1e-6), color=c, lw=1.0,
                label=f"{p} ({d['applied_threshold'].iloc[0]})")
    ax.set_yscale("log")
    ax.set_xlabel("Full-sequence bit score")
    ax.set_ylabel("Decoy FDR")
    ax.legend()
    ax.set_title("A  Empirical FDR from reversed decoys", loc="left")

    ax = axes[1]
    for p, c in zip(profiles, cols):
        d = fdr[fdr["profile"] == p]
        ax.plot(d["bit_score"], d["n_target"].clip(lower=0.5), color=c, lw=1.0, label=f"{p} target")
        ax.plot(d["bit_score"], d["n_decoy"].clip(lower=0.5), color=c, lw=1.0, ls="--",
                label=f"{p} decoy")
    ax.set_yscale("log")
    ax.set_xlabel("Full-sequence bit score")
    ax.set_ylabel("Hits $\\geq$ score")
    ax.legend(fontsize=5.5)
    ax.set_title("B  Target vs decoy", loc="left")

    fig.tight_layout()
    savefig(fig, figdir, "20_decoy_fdr", fmts, dpi)


def fig_evidence_tiers(tiers_path, figdir, fmts, dpi):
    """Does an SSF54001-only hit also carry the Pei triad in the Pei positions?"""
    t = pd.read_csv(tiers_path, sep="\t")
    if t.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.2))
    x = np.arange(len(t))
    cols = palette(len(t))

    ax = axes[0]
    ax.bar(x - 0.2, t["n_aligned"], width=0.35, color="0.75", label="aligned")
    ax.bar(x + 0.2, t["n_triad_positive"], width=0.35, color=PALETTE[0],
           label="full triad")
    ax.set_xticks(x); ax.set_xticklabels(t["evidence"])
    ax.set_yscale("log"); ax.set_ylabel("Sequences")
    ax.legend()
    ax.set_title("A  Triad filter by evidence tier", loc="left")

    ax = axes[1]
    ax.bar(x, 100 * t["frac_triad_positive"], color=cols, width=0.55)
    for i, v in enumerate(100 * t["frac_triad_positive"]):
        ax.text(i, v, f"{v:.1f}%", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x); ax.set_xticklabels(t["evidence"])
    ax.set_ylabel("% carrying the full triad")
    ax.set_title("B  Pass rate", loc="left")

    fig.tight_layout()
    savefig(fig, figdir, "21_evidence_tiers", fmts, dpi)


def fig_adjusted_prevalence(prev_path, figdir, fmts, dpi):
    """Raw vs detection-bias-adjusted prevalence, per taxon."""
    p = pd.read_csv(prev_path, sep="\t")
    if p.empty or "adjusted_prevalence" not in p.columns:
        print("[plots] no adjusted prevalence table; skipping", file=sys.stderr)
        return
    p = p.sort_values("adjusted_prevalence", ascending=True)
    fig, ax = plt.subplots(figsize=(4.4, max(2.0, 0.16 * len(p) + 0.8)))
    y = np.arange(len(p))
    ax.hlines(y, p["lo"], p["hi"], color="0.6", lw=0.8)
    ax.scatter(p["adjusted_prevalence"], y, s=10, color=PALETTE[0], zorder=3,
               label="adjusted (complete genome, median sequencing effort)")
    ax.scatter(p["raw_prevalence"], y, s=10, marker="x", color=PALETTE[4], zorder=3,
               label="raw")
    labels = p["taxon"].astype(str)
    if "n_species" in p.columns:
        labels = labels + " (" + p["n_species"].astype(int).astype(str) + " spp.)"
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel("P(species carries C71)")
    ax.legend(loc="lower right")
    ax.set_title("Detection-bias-adjusted prevalence, per species", loc="left")
    fig.tight_layout()
    savefig(fig, figdir, "22_adjusted_prevalence", fmts, dpi)


def fig_profile_agreement(hits, figdir, fmts, dpi):
    wide = hits.pivot_table(index=["sample", "protein_id"], columns="profile",
                            values="full_score", aggfunc="max")
    profiles = list(wide.columns)
    fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.4))

    ax = axes[0]
    counts = [
        int((wide.notna().sum(axis=1) == len(profiles)).sum()),
        *[int((wide[p].notna() & wide.drop(columns=p).isna().all(axis=1)).sum())
          for p in profiles],
    ]
    labels = ["Both"] + [f"{p} only" for p in profiles]
    ax.bar(labels, counts, color=palette(len(counts)), width=0.6)
    for i, v in enumerate(counts):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=6)
    ax.set_ylabel("Proteins")
    ax.set_title("A  Profile overlap", loc="left")
    ax.tick_params(axis="x", rotation=20)

    ax = axes[1]
    if len(profiles) == 2 and wide.dropna().shape[0] > 1:
        sub = wide.dropna()
        ax.scatter(sub[profiles[0]], sub[profiles[1]], s=1.5, alpha=0.25,
                   color=PALETTE[0], rasterized=True, linewidths=0)
        lims = [0, float(np.nanmax(sub.to_numpy())) * 1.05]
        ax.plot(lims, lims, color="0.4", lw=0.6, ls="--")
        r = float(np.corrcoef(sub[profiles[0]], sub[profiles[1]])[0, 1])
        ax.set_xlabel(f"{profiles[0]} bit score")
        ax.set_ylabel(f"{profiles[1]} bit score")
        ax.text(0.03, 0.95, f"Pearson $r$ = {r:.3f}\nn = {len(sub):,}",
                transform=ax.transAxes, va="top", fontsize=6)
    ax.set_title("B  Score concordance", loc="left")

    fig.tight_layout()
    savefig(fig, figdir, "02_profile_agreement", fmts, dpi)


def fig_coverage_length(hits, figdir, fmts, dpi):
    profiles = sorted(hits["profile"].unique())
    cols = palette(len(profiles))
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.2))

    ax = axes[0]
    data = [hits.loc[hits["profile"] == p, "profile_coverage"].dropna() for p in profiles]
    parts = ax.violinplot(data, showextrema=False, widths=0.8)
    for b, c in zip(parts["bodies"], cols):
        b.set_facecolor(c); b.set_alpha(0.7); b.set_edgecolor("none")
    for i, d in enumerate(data, 1):
        ax.boxplot(d, positions=[i], widths=0.12, showfliers=False,
                   medianprops=dict(color="k", lw=0.8),
                   boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    ax.set_xticks(range(1, len(profiles) + 1)); ax.set_xticklabels(profiles)
    ax.set_ylabel("Profile coverage of best domain")
    ax.set_title("A  Domain coverage", loc="left")

    ax = axes[1]
    lens = hits.drop_duplicates(["sample", "protein_id"])["protein_len"]
    ax.hist(lens, bins=80, color=PALETTE[0], alpha=0.85)
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Proteins")
    ax.set_title("B  Hit protein length", loc="left")

    ax = axes[2]
    d = hits.drop_duplicates(["sample", "protein_id"])
    ax.scatter(d["protein_len"], d["full_score"], s=1.5, alpha=0.2, color=PALETTE[2],
               rasterized=True, linewidths=0)
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Best bit score")
    ax.set_title("C  Length vs score", loc="left")

    fig.tight_layout()
    savefig(fig, figdir, "03_coverage_and_length", fmts, dpi)


def fig_prevalence(hits, stats, figdir, tabdir, fmts, dpi):
    per_sample = hits.drop_duplicates(["sample", "protein_id"]).groupby("sample").size()
    n_searched = int(stats.loc["samples_searched", "value"])
    n_with = int(per_sample.shape[0])

    fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.2))

    ax = axes[0]
    ax.bar(["With $\\geq$1 hit", "No hit"], [n_with, n_searched - n_with],
           color=[PALETTE[0], "0.8"], width=0.55)
    for i, v in enumerate([n_with, n_searched - n_with]):
        ax.text(i, v, f"{v:,}\n({100 * v / n_searched:.1f}%)", ha="center",
                va="bottom", fontsize=6)
    ax.set_ylabel("Samples")
    ax.set_title(f"A  Prevalence (n={n_searched:,})", loc="left")

    ax = axes[1]
    vc = per_sample.value_counts().sort_index()
    ax.bar(vc.index, vc.values, color=PALETTE[1], width=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Hit proteins per sample")
    ax.set_ylabel("Samples")
    ax.set_title("B  Copy number", loc="left")

    fig.tight_layout()
    savefig(fig, figdir, "04_prevalence", fmts, dpi)
    per_sample.rename("n_hit_proteins").to_csv(os.path.join(tabdir, "hits_per_sample.tsv"),
                                              sep="\t")


def fig_taxonomy(merged, cfg, figdir, tabdir, fmts, dpi):
    rank = cfg["plots"]["taxonomy_rank"]
    tax = resolve_taxonomy(merged, cfg["plots"]["taxonomy_col"], rank)
    if tax is None:
        print(f"[plots] no taxonomy column found for rank '{rank}'; skipping", file=sys.stderr)
        return
    merged = merged.assign(_tax=tax)
    per_prot = merged.drop_duplicates(["sample", "protein_id"])
    lab = top_n_with_other(per_prot["_tax"], cfg["plots"]["max_taxa_in_legend"])
    counts = lab.value_counts()

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8),
                             gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    y = np.arange(len(counts))
    ax.barh(y, counts.values, color=palette(len(counts)), height=0.7)
    ax.set_yticks(y); ax.set_yticklabels(counts.index)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Hit proteins")
    ax.set_title(f"A  Hits by {rank}", loc="left")

    ax = axes[1]
    grp = merged.assign(_lab=top_n_with_other(merged["_tax"],
                                              cfg["plots"]["max_taxa_in_legend"]))
    order = counts.index.tolist()
    data = [grp.loc[grp["_lab"] == t, "full_score"].dropna().to_numpy() for t in order]
    bp = ax.boxplot(data, vert=False, showfliers=False, widths=0.6, patch_artist=True,
                    medianprops=dict(color="k", lw=0.8))
    for b, c in zip(bp["boxes"], palette(len(order))):
        b.set_facecolor(c); b.set_alpha(0.75); b.set_linewidth(0.5)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Bit score")
    ax.set_title("B  Score by taxon", loc="left")

    fig.tight_layout()
    savefig(fig, figdir, "05_taxonomy", fmts, dpi)
    counts.rename("n_hit_proteins").rename_axis(rank).to_csv(
        os.path.join(tabdir, f"hits_by_{rank}.tsv"), sep="\t")


def fig_conservation(colstats, chosen, figdir, fmts, dpi):
    cs = pd.read_csv(colstats, sep="\t")
    with open(chosen) as fh:
        tri = json.load(fh)
    cols = tri["match_columns"]
    res = tri["residues"]

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 3.2), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1]})
    ax = axes[0]
    ax.fill_between(cs["match_col"], cs["occupancy"], color=PALETTE[0], alpha=0.5, lw=0)
    ax.set_ylabel("Occupancy")
    ax.set_ylim(0, 1.02)
    ax.set_title("A  Column occupancy across the profile", loc="left")

    ax = axes[1]
    ax.plot(cs["match_col"], cs["entropy_bits"], color="0.3", lw=0.6)
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_xlabel("Profile match-state column")
    ax.set_title("B  Column entropy; catalytic triad marked", loc="left")
    ax.invert_yaxis()

    for c, r, col in zip(cols, res, [PALETTE[3], PALETTE[4], PALETTE[2]]):
        for A in axes:
            A.axvline(c, color=col, lw=0.9, alpha=0.9)
        ax.annotate(f"{r}{c}", xy=(c, ax.get_ylim()[1]), xytext=(0, -2),
                    textcoords="offset points", ha="center", va="top",
                    fontsize=6, color=col, fontweight="bold")

    handles = [Line2D([0], [0], color=col, lw=1.2, label=f"{r} @ col {c} "
                      f"(freq {f:.2f})")
               for c, r, f, col in zip(cols, res, tri["residue_frequencies"],
                                       [PALETTE[3], PALETTE[4], PALETTE[2]])]
    axes[0].legend(handles=handles, loc="lower right", ncol=3)

    fig.tight_layout()
    savefig(fig, figdir, "06_alignment_conservation", fmts, dpi)


def fig_triad_funnel(stats, chosen, figdir, fmts, dpi):
    with open(chosen) as fh:
        tri = json.load(fh)
    steps = [
        ("Samples searched", int(stats.loc["samples_searched", "value"])),
        ("Proteins searched", int(stats.loc["proteins_searched", "value"])),
        ("Proteins with HMM hit", int(stats.loc["unique_proteins_with_hit", "value"])),
        ("Aligned", tri["n_input_sequences"]),
        ("Full catalytic triad", tri["n_triad_positive"]),
    ]
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    y = np.arange(len(steps))
    vals = [v for _, v in steps]
    ax.barh(y, vals, color=palette(len(steps)), height=0.65)
    ax.set_yticks(y); ax.set_yticklabels([s for s, _ in steps])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Count (log scale)")
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v:,}", va="center", fontsize=6)
    ax.set_title("Screening funnel", loc="left")
    fig.tight_layout()
    savefig(fig, figdir, "07_funnel", fmts, dpi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--search-stats", required=True)
    ap.add_argument("--fdr", required=True)
    ap.add_argument("--colstats", required=True)
    ap.add_argument("--chosen", required=True)
    ap.add_argument("--tiers", required=True)
    ap.add_argument("--prevalence", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tabdir", required=True)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    fmts = tuple(cfg["plots"]["formats"])
    dpi = cfg["plots"]["dpi"]
    set_style()
    os.makedirs(a.figdir, exist_ok=True)
    os.makedirs(a.tabdir, exist_ok=True)

    hits = pd.read_csv(a.combined, sep="\t", dtype={"sample": str, "protein_id": str})
    stats = pd.read_csv(a.search_stats, sep="\t", index_col="metric")

    if hits.empty:
        sys.exit("[plots] the combined hit table is empty; nothing to plot")

    fig_score_distributions(hits, a.figdir, fmts, dpi)
    fig_profile_agreement(hits, a.figdir, fmts, dpi)
    fig_coverage_length(hits, a.figdir, fmts, dpi)
    fig_prevalence(hits, stats, a.figdir, a.tabdir, fmts, dpi)
    fig_conservation(a.colstats, a.chosen, a.figdir, fmts, dpi)
    fig_triad_funnel(stats, a.chosen, a.figdir, fmts, dpi)
    fig_decoy_fdr(a.fdr, a.figdir, fmts, dpi)
    fig_evidence_tiers(a.tiers, a.figdir, fmts, dpi)
    fig_adjusted_prevalence(a.prevalence, a.figdir, fmts, dpi)

    merged = pd.read_csv(a.merged, sep="\t", dtype={"sample": str, "protein_id": str},
                         low_memory=False)
    fig_taxonomy(merged, cfg, a.figdir, a.tabdir, fmts, dpi)

    print(f"[plots] overview figures -> {a.figdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
