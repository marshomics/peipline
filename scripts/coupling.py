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

Significance is a z-score against the per-pair permutation null, not the raw
permutation rank. With B permutations the smallest achievable rank-based p is
1/(B+1); Benjamini-Hochberg over ~400 pairs then needs p < 1.2e-4 to reach
q = 0.05, so B would have to exceed ~9,000 before any pair could ever be called
significant. The z-score uses the null's mean and standard deviation rather than
its tail rank, which gets resolution beyond 1/B out of the same B permutations.
This is the standard treatment of MI_apc in the coevolution literature. The raw
empirical p is reported alongside so the approximation is visible.

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
from scipy.stats import norm
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


def mi_pair(ai: np.ndarray, aj: np.ndarray) -> float:
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
    Ap = A.copy()
    for b in range(B):
        for j in active:
            Ap[:, j] = A[rng.permutation(n), j]
        null = apc(mi_matrix(Ap, active))
        ge += (null >= OBS)
        s1 += null
        s2 += null ** 2
        if (b + 1) % max(1, B // 5) == 0:
            print(f"[coupling] permutation {b + 1}/{B}", file=sys.stderr)

    p_emp = (ge + 1.0) / (B + 1.0)
    mu = s1 / B
    sd = np.sqrt(np.maximum(s2 / B - mu ** 2, 0.0))
    sd[sd <= 0] = np.nan
    Z = (OBS - mu) / sd
    p_z = norm.sf(Z)                       # one-sided: coupling, not repulsion
    p_z = np.where(np.isfinite(p_z), p_z, 1.0)

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
                         "p": float(p_z[x, y]),
                         "p_empirical": float(p_emp[x, y]),
                         "same_block": bool(block[i] == block[j])})
    df = pd.DataFrame(rows)
    df["q_bh"] = multipletests(df["p"], method="fdr_bh")[1]
    df["significant"] = df["q_bh"] < float(ccfg["fdr"])
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
