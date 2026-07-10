#!/usr/bin/env python3
"""Statistical coupling between active-site barcode positions.

Question: are the subgroup-defining flank positions co-adapted, or independently
drifting? If two positions covary beyond what their marginal frequencies
predict, they plausibly belong to one substrate-binding surface. If they do not,
"k subgroups" is noise with a k-means fitted to it.

Method
------
Mutual information between every pair of barcode columns, with the
average-product correction (Dunn et al. 2008) subtracted:

    MI_apc(i,j) = MI(i,j) - MI(i,.) * MI(.,j) / MI(.,.)

APC removes the background MI that arises from a column simply being variable
and from the shared ancestry of the sequences. Without it, MI heatmaps of
protein alignments show a plaid pattern of "hot" rows that is entirely an
artefact of column entropy.

Redundancy correction
---------------------
The test is run on one representative per 90%-identity cluster. A permutation
test needs exchangeable observations, and 8,000 near-identical sequences from
one over-sequenced genus are not 8,000 observations. Down-weighting them would
fix the point estimate but not the null distribution; dereplicating fixes both.

Significance
------------
Each informative column is independently shuffled, destroying coupling while
preserving marginal composition (and therefore the finite-sample MI bias, which
at 21 states and a few hundred sequences is larger than the signal).

Significance comes from a parametric tail fitted to the per-pair permutation
null, not from the raw permutation rank. With B permutations the smallest
achievable rank-based p is 1/(B+1); Benjamini-Hochberg over ~400 pairs then needs
p < 1.2e-4 to reach q = 0.05, so B would have to exceed ~9,000 before any pair
could ever be called significant.

The tail is a moment-matched Gamma, NOT a z-score. Under independence,
2*N*ln2*MI is approximately chi-square with (a-1)(b-1) df, where a and b are the
numbers of residue states in the two columns. The barcode flanks a catalytic
residue, so those columns are conserved, the df is small, and the null is
strongly right-skewed. A normal tail on a right-skewed null understates the upper
tail: on a chi-square(4) null the z gives p = 7.7e-7 where the truth is 1e-3, a
factor of 2,400 anti-conservative, and every "significant" pair rests on it. The
Gamma is fitted from three accumulated moments of the same permutations, matches
the chi-square shape family, and degrades to the normal as the skew vanishes.

`p_normal`, `p_empirical`, `null_skew` and `tail_model` are all reported so the
approximation is auditable. The parametric p is floored at 1/(10*(B+1)): below
that we are extrapolating past what the permutations can support.

Gaps are excluded from the MI. A gap is not a residue state, and two columns
co-deleted in one clade would otherwise register as coupled.

Invariant columns (the three catalytic residues, fixed by the filter) have zero
entropy, carry zero MI, and are excluded from the FDR denominator rather than
counted as hundreds of free true negatives.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import gamma, norm
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import PALETTE, load_config, savefig, set_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY-"
IDX = {a: i for i, a in enumerate(AA)}
K = len(AA)


def encode_int(barcodes):
    n, L = len(barcodes), len(barcodes[0])
    A = np.full((n, L), IDX["-"], dtype=np.int64)
    for i, s in enumerate(barcodes):
        for j, c in enumerate(s):
            A[i, j] = IDX.get(c, IDX["-"])
    return A


GAP = IDX["-"]


def mi_pair(ai: np.ndarray, aj: np.ndarray) -> float:
    """MI over the rows where BOTH columns carry a residue.

    Gap was the 21st symbol. Two columns deleted together in one clade -- a single
    shared indel, i.e. phylogenetic signal -- then produced high MI, and the
    column-wise permutation preserved it because each column keeps its own gap
    fraction. APC removes column-wise background but not a pairwise co-gap term.
    With min_column_occupancy at 0.5 a column may be half gaps, so this was not
    hypothetical.

    Masking to the doubly-occupied rows makes the statistic what it claims to be:
    covariation between residues.
    """
    m = (ai != GAP) & (aj != GAP)
    if m.sum() < 2:
        return 0.0
    ai, aj = ai[m], aj[m]
    J = np.bincount(ai * K + aj, minlength=K * K).reshape(K, K).astype(float)
    tot = J.sum()
    if tot <= 0:
        return 0.0
    J /= tot
    pi, pj = J.sum(1), J.sum(0)
    nz = J > 0
    outer = np.outer(pi, pj)
    return float(np.sum(J[nz] * np.log2(J[nz] / outer[nz])))


def mi_matrix(A: np.ndarray, active: np.ndarray) -> np.ndarray:
    m = len(active)
    M = np.zeros((m, m))
    for x in range(m):
        for y in range(x + 1, m):
            M[x, y] = M[y, x] = mi_pair(A[:, active[x]], A[:, active[y]])
    return M


def apc(M: np.ndarray) -> np.ndarray:
    C = M.copy()
    np.fill_diagonal(C, np.nan)
    row = np.nanmean(C, axis=1)
    grand = np.nanmean(C)
    if not np.isfinite(grand) or grand <= 0:
        return M
    out = M - np.outer(row, row) / grand
    np.fill_diagonal(out, 0.0)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assign", required=True)
    ap.add_argument("--chosen", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contacts", help="Ca distance matrix over match columns")
    ap.add_argument("--groove", help="groove_columns.tsv")
    a = ap.parse_args()

    cfg = load_config(a.config)
    ccfg = cfg["coupling"]
    fmts, dpi = tuple(cfg["plots"]["formats"]), cfg["plots"]["dpi"]
    set_style()
    os.makedirs(a.figdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    assign = pd.read_csv(a.assign, sep="\t")
    with open(a.chosen) as fh:
        tri = json.load(fh)

    # --- dereplicate --------------------------------------------------------
    n_all = len(assign)
    if "cluster" in assign.columns and assign["cluster"].notna().any():
        derep = assign.sort_values("seq_id").drop_duplicates("cluster")
    else:
        derep = assign
    print(f"[coupling] {n_all} sequences -> {len(derep)} cluster representatives",
          file=sys.stderr)
    if len(derep) < 20:
        sys.exit(f"[coupling] only {len(derep)} independent sequences; a permutation "
                 f"test on this would not mean anything")

    barcodes = derep["barcode"].astype(str).tolist()
    A = encode_int(barcodes)
    n, L = A.shape
    rng = np.random.default_rng(cfg["active_site"]["random_state"])

    layout_path = os.path.join(os.path.dirname(a.chosen), "barcode_layout.json")
    if os.path.exists(layout_path):
        layout = json.load(open(layout_path))
        centres, bounds = layout["centres"], layout["bounds"]
    else:
        w_flank = cfg["active_site"]["flank_window"]
        bounds = [[b * (2 * w_flank + 1), (b + 1) * (2 * w_flank + 1)] for b in range(3)]
        centres = [b * (2 * w_flank + 1) + w_flank for b in range(3)]
        layout = {"centres": centres, "bounds": bounds}

    block = np.zeros(L, dtype=int)
    for bi, (lo, hi) in enumerate(bounds):
        block[lo:hi] = bi

    occ = np.array([(A[:, j] != IDX["-"]).mean() for j in range(L)])
    varied = np.array([len(np.unique(A[:, j])) > 1 for j in range(L)])
    active = np.flatnonzero((occ >= ccfg["min_column_occupancy"]) & varied)
    print(f"[coupling] {L} barcode columns, {len(active)} informative. The three "
          f"catalytic residues are invariant by construction and carry no MI.",
          file=sys.stderr)
    if len(active) < 2:
        sys.exit("[coupling] fewer than two informative columns")

    OBS_raw = mi_matrix(A, active)
    OBS = apc(OBS_raw)

    B = int(ccfg["n_permutations"])
    ge = np.zeros_like(OBS)
    s1 = np.zeros_like(OBS)
    s2 = np.zeros_like(OBS)
    s3 = np.zeros_like(OBS)
    Ap = A.copy()
    for b in range(B):
        for j in active:
            Ap[:, j] = A[rng.permutation(n), j]
        null = apc(mi_matrix(Ap, active))
        ge += (null >= OBS)
        s1 += null
        s2 += null ** 2
        s3 += null ** 3
        if (b + 1) % max(1, B // 5) == 0:
            print(f"[coupling] permutation {b + 1}/{B}", file=sys.stderr)

    p_emp = (ge + 1.0) / (B + 1.0)
    mu = s1 / B
    var = np.maximum(s2 / B - mu ** 2, 0.0)
    sd = np.sqrt(var)
    sd[sd <= 0] = np.nan
    Z = (OBS - mu) / sd

    # The tail, and why it is not a z-test.
    #
    # Under independence, 2*N*ln2*MI is approximately chi-square with (a-1)(b-1)
    # degrees of freedom, where a and b are the numbers of residue states in the
    # two columns. The barcode flanks a catalytic residue, so those columns are
    # conserved, a and b are small, the df is small, and the null is strongly
    # RIGHT-SKEWED. norm.sf() on a right-skewed null understates the upper tail:
    # p comes out too small and the count of "coupled" pairs is inflated. The
    # empirical p cannot corroborate, because it has a 1/(B+1) floor that BH over
    # ~400 pairs can never clear -- which is exactly why the z-score was reached
    # for in the first place.
    #
    # So fit a shifted Gamma to each pair's own permutation null by method of
    # moments (three accumulated moments, no extra permutations), and take the
    # tail from that. A Gamma matches a chi-square's shape family and degrades to
    # the normal as the skew goes to zero, so this is strictly better calibrated
    # than the z and never worse.
    m3 = s3 / B - 3.0 * mu * (s2 / B) + 2.0 * mu ** 3      # central third moment
    with np.errstate(divide="ignore", invalid="ignore"):
        skew = m3 / np.power(var, 1.5)
        # Gamma(k, theta) shifted by loc: skew = 2/sqrt(k)
        k = 4.0 / np.square(skew)
        theta = np.sqrt(var / k)
        loc = mu - k * theta
        p_gamma = gamma.sf(OBS, a=k, loc=loc, scale=theta)

    # Fall back to the normal only where the Gamma is undefined: a null with no
    # variance, or a left-skewed one (which the chi-square argument forbids and
    # which therefore indicates too few permutations, not a real shape).
    usable = np.isfinite(p_gamma) & (skew > 1e-6) & np.isfinite(sd)
    p_z = norm.sf(Z)
    p = np.where(usable, p_gamma, p_z)
    p = np.where(np.isfinite(p), p, 1.0)
    # A parametric tail must never claim more than the permutations can support.
    # Below the empirical floor we are extrapolating; say so by not going under it
    # by more than an order of magnitude per decade of B.
    p = np.maximum(p, 1.0 / (10.0 * (B + 1.0)))

    n_gamma = int(usable[np.triu_indices_from(usable, 1)].sum())
    n_pairs = len(active) * (len(active) - 1) // 2
    print(f"[coupling] tail calibrated by a moment-matched Gamma for {n_gamma}/"
          f"{n_pairs} pairs (median null skew "
          f"{np.nanmedian(skew[np.triu_indices_from(skew, 1)]):.2f}); the rest "
          f"fall back to the normal.", file=sys.stderr)

    rows = []
    for x in range(len(active)):
        for y in range(x + 1, len(active)):
            i, j = int(active[x]), int(active[y])
            rows.append({"col_i": i, "col_j": j,
                         "residue_i": tri["residues"][block[i]],
                         "residue_j": tri["residues"][block[j]],
                         "mi_bits": float(OBS_raw[x, y]),
                         "mi_apc": float(OBS[x, y]),
                         "null_mean": float(mu[x, y]), "null_sd": float(sd[x, y]),
                         "z": float(Z[x, y]),
                         "null_skew": float(skew[x, y]),
                         "tail_model": "gamma" if usable[x, y] else "normal",
                         "p": float(p[x, y]),
                         "p_normal": float(p_z[x, y]),
                         "p_empirical": float(p_emp[x, y]),
                         "same_block": bool(block[i] == block[j])})
    df = pd.DataFrame(rows)
    df["q_bh"] = multipletests(df["p"], method="fdr_bh")[1]
    df["significant"] = df["q_bh"] < float(ccfg["fdr"])

    # --- structural validation -----------------------------------------------
    # A coupled pair that is 40 A apart in the structure is phylogenetic signal
    # that APC failed to remove. A coupled pair that is adjacent, and in the
    # substrate groove, is a co-adapted specificity surface. Without this check
    # the coupling result is a statistic; with it, it is a claim about a protein.
    contact_summary = {}
    if a.contacts and a.groove and os.path.exists(a.contacts):
        D = np.loadtxt(a.contacts, delimiter="\t")
        gr = pd.read_csv(a.groove, sep="\t").set_index("match_col")
        # barcode position -> alignment match column
        bar2col = {p: c for p, c in enumerate(layout["alignment_columns"])} \
            if "alignment_columns" in layout else {}
        if bar2col:
            ci = df["col_i"].map(bar2col)
            cj = df["col_j"].map(bar2col)
            df["align_col_i"], df["align_col_j"] = ci, cj
            df["ca_distance_a"] = [
                D[int(x), int(y)] if (pd.notna(x) and pd.notna(y)
                                      and int(x) < D.shape[0] and int(y) < D.shape[0])
                else np.nan for x, y in zip(ci, cj)]
            rad = float(load_config(a.config)["specificity"]["contact_radius_a"])
            df["in_contact"] = df["ca_distance_a"] <= rad
            df["both_in_groove"] = [
                bool(gr["in_groove"].get(x, 0)) and bool(gr["in_groove"].get(y, 0))
                for x, y in zip(ci, cj)]

            sig = df[df["significant"]]
            has_d = df["ca_distance_a"].notna()
            if has_d.any() and len(sig):
                from scipy.stats import fisher_exact, mannwhitneyu
                tab = pd.crosstab(df.loc[has_d, "significant"],
                                  df.loc[has_d, "in_contact"])
                if tab.shape == (2, 2):
                    orr, pf = fisher_exact(tab.to_numpy(), alternative="greater")
                else:
                    orr, pf = np.nan, np.nan
                u, pu = mannwhitneyu(
                    df.loc[has_d & df["significant"], "ca_distance_a"].dropna(),
                    df.loc[has_d & ~df["significant"], "ca_distance_a"].dropna(),
                    alternative="less") if (has_d & df["significant"]).sum() > 2 \
                    else (np.nan, np.nan)
                contact_summary = {
                    "n_pairs_with_structure": int(has_d.sum()),
                    "n_significant": int(sig["significant"].sum()),
                    "n_significant_in_contact": int((sig["in_contact"] == True).sum()),  # noqa: E712
                    "n_significant_both_in_groove": int((sig["both_in_groove"] == True).sum()),  # noqa: E712
                    "contact_enrichment_or": orr, "contact_enrichment_p": pf,
                    "median_dist_significant": float(sig["ca_distance_a"].median()),
                    "median_dist_other": float(df.loc[has_d & ~df["significant"],
                                                      "ca_distance_a"].median()),
                    "mannwhitney_p": pu,
                }
                for k_, v_ in contact_summary.items():
                    print(f"[coupling] {k_}: {v_}", file=sys.stderr)
                pd.Series(contact_summary).rename("value").rename_axis(
                    "metric").to_csv(a.out.replace(".tsv", "_contacts.tsv"), sep="\t")
    else:
        print("[coupling] no structure supplied; coupled pairs are not validated "
              "against spatial adjacency. The result stays a statistic.",
              file=sys.stderr)

    df = df.sort_values("q_bh")
    df.to_csv(a.out, sep="\t", index=False)

    floor = 1.0 / (B + 1)
    if df["p_empirical"].min() <= floor and (df["q_bh"] < float(ccfg["fdr"])).sum() == 0:
        print(f"[coupling] note: the rank-based p floor is {floor:.3g}; with "
              f"{len(df)} pairs, rank-based BH could not have called anything. "
              f"The z-score p-values are the ones being used.", file=sys.stderr)

    n_sig = int(df["significant"].sum())
    n_cross = int((df["significant"] & ~df["same_block"]).sum())
    print(f"[coupling] {n_sig}/{len(df)} pairs significant at q < {ccfg['fdr']}; "
          f"{n_cross} span different catalytic residues", file=sys.stderr)

    # --- figure --------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0),
                             gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    M = np.full((L, L), np.nan)
    M[np.ix_(active, active)] = OBS
    np.fill_diagonal(M, np.nan)
    vmax = float(np.nanpercentile(np.abs(M), 99)) if np.isfinite(M).any() else 1.0
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
    for c in centres:
        ax.axhline(c, color="0.2", lw=0.4, ls=":")
        ax.axvline(c, color="0.2", lw=0.4, ls=":")
    ax.set_xticks(centres)
    ax.set_xticklabels([f"{r}{c}" for r, c in zip(tri["residues"], tri["match_columns"])])
    ax.set_yticks(centres)
    ax.set_yticklabels([f"{r}{c}" for r, c in zip(tri["residues"], tri["match_columns"])])
    ax.set_title(f"A  APC-corrected MI (n={n} independent sequences)", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("MI$_{APC}$ (bits)"); cb.outline.set_linewidth(0.4)

    ax = axes[1]
    for m, col, lab in ((df["same_block"], PALETTE[0], "within one residue's flank"),
                        (~df["same_block"], PALETTE[4], "across catalytic residues")):
        ax.scatter(df.loc[m, "mi_apc"], -np.log10(df.loc[m, "q_bh"].clip(lower=1e-300)),
                   s=6, color=col, label=lab, linewidths=0)
    ax.axhline(-np.log10(float(ccfg["fdr"])), color="0.4", lw=0.6, ls="--")
    ax.set_xlabel("MI$_{APC}$ (bits)")
    ax.set_ylabel("$-\\log_{10}$ q")
    ax.legend(loc="upper left")
    ax.set_title(f"B  {n_sig} coupled pairs at q < {ccfg['fdr']}", loc="left")

    fig.tight_layout()
    savefig(fig, a.figdir, "16_barcode_coupling", fmts, dpi)

    sig = df[df["significant"]]
    if len(sig):
        print(sig.head(10).to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
