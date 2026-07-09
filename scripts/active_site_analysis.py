#!/usr/bin/env python3
"""Compare active sites across the triad-positive sequences.

What is being asked: how many defensible sub-groups do the catalytic-site
neighbourhoods fall into, and what distinguishes them?

Approach
--------
The three triad columns are, by construction, invariant across every retained
sequence. All the signal therefore lives in the flanks. We take a window of
`flank_window` match columns either side of each triad residue, concatenate the
three windows into one active-site "barcode", and analyse that:

  * BLOSUM62 embedding rather than one-hot. One-hot treats D->E and D->W as
    equally distant; BLOSUM62 does not, and substitution-aware distances are
    what you want when the question is functional rather than phylogenetic.
  * Columns are weighted by conservation, w_j = 1 - H_j/log2(20). Without this
    the barcode is dominated by fast-evolving positions: a 33-column barcode
    with five functionally-constrained sites and twenty-eight free ones will
    cluster on the noise. Weighting is the difference between recovering the
    real subgroups and recovering nothing.
  * PCA before clustering, keeping enough components for `pca_variance` of the
    variance (capped at `pca_components`), so k-means operates on a compact
    basis where Euclidean distance is meaningful. Retaining every component of
    a wide, small-n matrix re-admits the noise that weighting removed.
  * k is chosen by a sweep over silhouette, Calinski-Harabasz and
    Davies-Bouldin. A density-based HDBSCAN run is reported alongside as an
    independent check that does not presuppose k, or convex clusters.
  * Subgroup identity is then tested against taxonomy: Fisher exact per
    subgroup x taxon cell, Benjamini-Hochberg corrected. Subgroups are also
    painted onto the phylogeny (plot_tree.py) so that convergence -- the same
    active-site barcode arising in distant clades -- is visible rather than
    assumed away.

Outputs: subgroup assignments, per-subgroup sequence logos, entropy profile,
pairwise-identity heatmap, PCA scatter, model-selection curves, a consensus
motif table, and the taxonomy enrichment table.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import fisher_exact
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import (calinski_harabasz_score, davies_bouldin_score,
                             silhouette_score)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (PALETTE, load_config, palette, read_fasta,  # noqa: E402
                   resolve_taxonomy, savefig, set_style, top_n_with_other)

import matplotlib.pyplot as plt  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}


# ---------------------------------------------------------------------------
def blosum_embedding():
    """20x20 BLOSUM62 rows, z-scored, plus a zero row for gaps."""
    from Bio.Align import substitution_matrices
    m = substitution_matrices.load("BLOSUM62")
    M = np.array([[float(m[a, b]) for b in AA] for a in AA], dtype=np.float32)
    M = (M - M.mean()) / M.std()
    return M


def encode(barcodes, M, weights=None):
    n, L = len(barcodes), len(barcodes[0])
    X = np.zeros((n, L * 20), dtype=np.float32)
    for i, s in enumerate(barcodes):
        for j, c in enumerate(s):
            k = AA_IDX.get(c)
            if k is not None:
                X[i, j * 20:(j + 1) * 20] = M[k] if weights is None else M[k] * weights[j]
    return X


def conservation_weights(barcodes, seq_w=None):
    """w_j = 1 - H_j / log2(20), in [0, 1]. Invariant columns get 1, columns at
    background composition get ~0. H is redundancy-weighted."""
    arr = np.array([list(b) for b in barcodes])
    h = np.array([shannon(arr[:, j], seq_w) for j in range(arr.shape[1])])
    return np.clip(1.0 - h / np.log2(20), 0.0, 1.0).astype(np.float32)


def barcode_blocks(cols, w, L):  # noqa: D401
    """Concatenate a +/-w window around each triad column.

    Returns (idx, centres, bounds): the alignment columns making up the
    barcode, the barcode-relative position of each catalytic residue, and the
    (start, stop) of each block. Blocks are kept separate rather than merged
    into a sorted set, so a barcode position always maps to exactly one
    catalytic residue even if two windows overlap.
    """
    idx, centres, bounds = [], [], []
    for c in cols:
        lo, hi = max(0, c - w), min(L, c + w + 1)
        bounds.append((len(idx), len(idx) + hi - lo))
        centres.append(len(idx) + (c - lo))
        idx.extend(range(lo, hi))
    return idx, centres, bounds


def shannon(col: np.ndarray, w: np.ndarray | None = None) -> float:
    """Redundancy-weighted Shannon entropy over non-gap residues, in bits."""
    if w is None:
        w = np.ones(len(col))
    c: dict = {}
    tot = 0.0
    for x, wi in zip(col, w):
        if x == "-":
            continue
        c[x] = c.get(x, 0.0) + wi
        tot += wi
    if tot <= 0 or not c:
        return 0.0
    p = np.array(list(c.values())) / tot
    return float(-(p * np.log2(p)).sum())


def weighted_counts(barcodes, weights, sel=None):
    """L x 20 weighted residue-count matrix."""
    L = len(barcodes[0])
    M = np.zeros((L, 20), dtype=float)
    it = range(len(barcodes)) if sel is None else sel
    for i in it:
        wi = weights[i]
        for j, ch in enumerate(barcodes[i]):
            k = AA_IDX.get(ch)
            if k is not None:
                M[j, k] += wi
    return M


# ---------------------------------------------------------------------------
def n_components_parallel_analysis(X, cap, seed, n_perm=25, q=95, max_rows=3000):
    """Horn's parallel analysis: keep the leading components whose eigenvalue
    beats the 95th percentile of the eigenvalue you would get from the same
    matrix with every feature independently permuted. Data-driven, and it does
    not require guessing a variance cutoff -- 90% of the variance of a wide,
    nearly-isotropic matrix is mostly noise."""
    rng = np.random.default_rng(seed)
    Xs = X[rng.choice(len(X), max_rows, replace=False)] if len(X) > max_rows else X
    cap = int(min(cap, Xs.shape[0] - 1, Xs.shape[1]))
    obs = PCA(n_components=cap, random_state=seed).fit(Xs).explained_variance_

    null = np.zeros((n_perm, cap))
    Xp = np.array(Xs, copy=True)
    for p in range(n_perm):
        for j in range(Xp.shape[1]):
            rng.shuffle(Xp[:, j])
        null[p] = PCA(n_components=cap, random_state=seed).fit(Xp).explained_variance_
    thr = np.percentile(null, q, axis=0)

    # keep a prefix: stop at the first component that fails
    failed = np.flatnonzero(obs <= thr)
    n = int(failed[0]) if failed.size else cap
    return max(2, min(n, cap))


def gap_statistic(Z, kmin, kmax, seed, B=20, max_rows=3000):
    """Tibshirani's gap statistic against a uniform reference in the bounding
    box. Returns the smallest k with Gap(k) >= Gap(k+1) - s(k+1)."""
    rng = np.random.default_rng(seed)
    Zs = Z[rng.choice(len(Z), max_rows, replace=False)] if len(Z) > max_rows else Z
    lo, hi = Zs.min(0), Zs.max(0)
    gap, sd = {}, {}
    for k in range(max(1, kmin - 1), kmax + 1):
        w = np.log(KMeans(k, n_init=10, random_state=seed).fit(Zs).inertia_ + 1e-12)
        refs = np.array([
            np.log(KMeans(k, n_init=5, random_state=seed)
                   .fit(rng.uniform(lo, hi, size=Zs.shape)).inertia_ + 1e-12)
            for _ in range(B)])
        gap[k], sd[k] = refs.mean() - w, refs.std() * np.sqrt(1 + 1 / B)
    for k in sorted(gap)[:-1]:
        if gap[k] >= gap[k + 1] - sd[k + 1]:
            return k, gap
    return kmax, gap


def choose_k(Z, kmin, kmax, seed, criterion="silhouette_cosine", sil_sample=5000,
             sample_weight=None):
    """Sweep k and score each solution four ways.

    The primary criterion is the *cosine* silhouette. On a BLOSUM-embedded
    barcode the informative signal is the direction of the residue vector, not
    its magnitude; Euclidean silhouette systematically prefers k=2 because it
    is dominated by vector length, which tracks how many non-gap positions a
    sequence has rather than which subgroup it belongs to. The other three
    metrics are reported so a disagreement is visible rather than hidden.
    """
    rng = np.random.default_rng(seed)
    sub = (rng.choice(len(Z), sil_sample, replace=False)
           if len(Z) > sil_sample else np.arange(len(Z)))
    rows, labels = [], {}
    for k in range(kmin, min(kmax, len(Z) - 1) + 1):
        km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(
            Z, sample_weight=sample_weight)
        lab = km.labels_
        labels[k] = lab
        if len(set(lab[sub])) < 2:
            continue
        rows.append({
            "k": k,
            "silhouette_cosine": silhouette_score(Z[sub], lab[sub], metric="cosine"),
            "silhouette_euclidean": silhouette_score(Z[sub], lab[sub]),
            "calinski_harabasz": calinski_harabasz_score(Z, lab),
            "davies_bouldin": davies_bouldin_score(Z, lab),
            "inertia": km.inertia_,
        })
    m = pd.DataFrame(rows)
    if m.empty:
        raise SystemExit("[active_site] k sweep produced no valid solution")

    k_gap, _ = gap_statistic(Z, kmin, min(kmax, len(Z) - 2), seed)
    m["k_gap_statistic"] = k_gap

    if criterion == "davies_bouldin":
        best = int(m.loc[m["davies_bouldin"].idxmin(), "k"])
    elif criterion == "gap":
        best = int(k_gap)
    else:
        best = int(m.loc[m[criterion].idxmax(), "k"])
    return best, m, labels, k_gap


def run_hdbscan(Z, min_cluster_size):
    try:
        from sklearn.cluster import HDBSCAN
    except ImportError:
        print("[active_site] sklearn too old for HDBSCAN; skipping", file=sys.stderr)
        return None
    h = HDBSCAN(min_cluster_size=min_cluster_size).fit(Z)
    lab = h.labels_
    n = len(set(lab) - {-1})
    print(f"[active_site] HDBSCAN: {n} clusters, "
          f"{(lab == -1).sum()} noise points", file=sys.stderr)
    return lab


# ---------------------------------------------------------------------------
def fig_model_selection(metrics, k_best, k_gap, hdb_k, figdir, fmts, dpi):
    cols = ["silhouette_cosine", "silhouette_euclidean", "calinski_harabasz",
            "davies_bouldin", "inertia"]
    labs = ["Silhouette, cosine (higher)", "Silhouette, Euclidean (higher)",
            "Calinski-Harabasz (higher)", "Davies-Bouldin (lower)", "Within-cluster SS"]
    fig, axes = plt.subplots(1, 5, figsize=(9.0, 1.9))
    for ax, col, lab in zip(axes, cols, labs):
        ax.plot(metrics["k"], metrics[col], "o-", ms=2.5, color=PALETTE[0])
        ax.axvline(k_best, color=PALETTE[4], lw=0.8, ls="--")
        ax.axvline(k_gap, color=PALETTE[2], lw=0.8, ls=":")
        ax.set_xlabel("k"); ax.set_ylabel(lab)
    axes[0].set_title(f"k = {k_best} (orange); gap statistic k = {k_gap} (green)"
                      + (f"; HDBSCAN {hdb_k}" if hdb_k is not None else ""),
                      loc="left")
    fig.tight_layout()
    savefig(fig, figdir, "10_subgroup_model_selection", fmts, dpi)


def fig_pca(Z, labels, hdb, figdir, fmts, dpi, var):
    ncol = 2 if hdb is None else 3
    fig, axes = plt.subplots(1, ncol, figsize=(2.5 * ncol, 2.5))
    uniq = sorted(set(labels))
    cols = palette(len(uniq))
    for ax, (x, y) in zip(axes[:2], [(0, 1), (0, 2)]):
        if Z.shape[1] <= max(x, y):
            ax.axis("off"); continue
        for u, c in zip(uniq, cols):
            m = labels == u
            ax.scatter(Z[m, x], Z[m, y], s=2, color=c, alpha=0.5, linewidths=0,
                       rasterized=True, label=f"SG{u} (n={int(m.sum()):,})")
        ax.set_xlabel(f"PC{x + 1} ({100 * var[x]:.1f}%)")
        ax.set_ylabel(f"PC{y + 1} ({100 * var[y]:.1f}%)")
    axes[0].legend(loc="best", markerscale=2)
    axes[0].set_title("A  Active-site PCA, k-means", loc="left")
    if hdb is not None:
        ax = axes[2]
        hu = sorted(set(hdb))
        for u, c in zip(hu, palette(len(hu))):
            m = hdb == u
            ax.scatter(Z[m, 0], Z[m, 1], s=2,
                       color=("0.8" if u == -1 else c), alpha=0.5, linewidths=0,
                       rasterized=True)
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        ax.set_title("B  HDBSCAN (grey = noise)", loc="left")
    fig.tight_layout()
    savefig(fig, figdir, "11_active_site_pca", fmts, dpi)


def fig_logos(barcodes, labels, seq_w, cols, res, centres, bounds, figdir, fmts, dpi):
    try:
        import logomaker
    except ImportError:
        print("[active_site] logomaker missing; skipping logos", file=sys.stderr)
        return
    uniq = sorted(set(labels))
    L = len(barcodes[0])
    fig, axes = plt.subplots(len(uniq), 1, figsize=(7.2, 1.05 * len(uniq) + 0.4),
                             sharex=True, squeeze=False)
    for ax, u in zip(axes[:, 0], uniq):
        sel = np.flatnonzero(np.asarray(labels) == u)
        eff = float(seq_w[sel].sum())
        counts = pd.DataFrame(weighted_counts(barcodes, seq_w, sel),
                              index=range(L), columns=list(AA))
        # All-gap columns would divide by zero; make them uniform (= 0 bits).
        counts.loc[counts.sum(axis=1) == 0] = 1.0
        # pseudocount=0. logomaker's default of 1-per-residue adds 20 phantom
        # observations, which flattens an invariant column from 4.3 bits to
        # under 1 whenever a subgroup has only a few dozen members.
        info = logomaker.transform_matrix(counts, from_type="counts",
                                          to_type="information", pseudocount=0)
        logomaker.Logo(info, ax=ax, color_scheme="chemistry", show_spines=False)
        ax.set_ylabel(f"SG{u}\n(n={len(sel):,}, eff {eff:.0f})", fontsize=6)
        ax.set_ylim(0, 4.32)
        for c in centres:
            ax.axvline(c, color=PALETTE[4], lw=0.6, alpha=0.6)
        for lo, _ in bounds[1:]:
            ax.axvline(lo - 0.5, color="0.6", lw=0.6, ls=":")
    axes[-1, 0].set_xticks(centres)
    axes[-1, 0].set_xticklabels([f"{r} (col {c})" for r, c in zip(res, cols)])
    axes[-1, 0].set_xlabel("Active-site barcode: flanking columns around each "
                           "catalytic residue")
    axes[0, 0].set_title("Redundancy-weighted sequence logo per active-site subgroup",
                         loc="left")
    fig.tight_layout()
    savefig(fig, figdir, "12_subgroup_logos", fmts, dpi)


def fig_entropy(barcodes, seq_w, cols, res, centres, figdir, tabdir, fmts, dpi):
    arr = np.array([list(b) for b in barcodes])
    ent = np.array([shannon(arr[:, j], seq_w) for j in range(arr.shape[1])])
    fig, ax = plt.subplots(figsize=(7.2, 1.9))
    x = np.arange(len(ent))
    colors = [PALETTE[4] if j in set(centres) else PALETTE[0] for j in x]
    ax.bar(x, ent, color=colors, width=0.8)
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_xticks(centres)
    ax.set_xticklabels([f"{r}{c}" for r, c in zip(res, cols)])
    ax.set_xlabel("Position in active-site barcode")
    ax.set_title("Per-position variability; catalytic residues in orange", loc="left")
    fig.tight_layout()
    savefig(fig, figdir, "13_active_site_entropy", fmts, dpi)
    pd.DataFrame({"barcode_position": x, "entropy_bits": ent}).to_csv(
        os.path.join(tabdir, "active_site_entropy.tsv"), sep="\t", index=False)


def fig_heatmap_dendrogram(barcodes, labels, seed, figdir, fmts, dpi, max_n=400):
    rng = np.random.default_rng(seed)
    idx = (rng.choice(len(barcodes), max_n, replace=False)
           if len(barcodes) > max_n else np.arange(len(barcodes)))
    B = np.array([list(barcodes[i]) for i in idx])
    lab = np.asarray(labels)[idx]

    D = pdist(B.view(np.uint32).reshape(len(B), -1), metric="hamming")
    Zl = linkage(D, method="average")
    order = dendrogram(Zl, no_plot=True)["leaves"]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0),
                             gridspec_kw={"width_ratios": [1, 1.4]})
    ax = axes[0]
    dendrogram(Zl, ax=ax, no_labels=True, color_threshold=0,
               link_color_func=lambda _: "0.35")
    ax.set_ylabel("Average-linkage Hamming distance")
    ax.set_title(f"A  Barcode dendrogram (n={len(B)} sampled)", loc="left")

    ax = axes[1]
    from scipy.spatial.distance import squareform
    S = 1 - squareform(D)
    im = ax.imshow(S[np.ix_(order, order)], cmap="magma", vmin=0, vmax=1,
                   interpolation="nearest", rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("B  Pairwise barcode identity", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("Identity")
    cb.outline.set_linewidth(0.4)

    uniq = sorted(set(lab))
    cmap = dict(zip(uniq, palette(len(uniq))))
    for i, o in enumerate(order):
        ax.add_patch(plt.Rectangle((-len(order) * 0.03, i - 0.5), len(order) * 0.025, 1,
                                   color=cmap[lab[o]], clip_on=False, lw=0))
    fig.tight_layout()
    savefig(fig, figdir, "14_barcode_similarity", fmts, dpi)


def taxonomy_enrichment(df, tabdir):
    """Fisher-exact subgroup x taxon, run on the 90%-identity cluster
    representatives rather than on every sequence.

    Running it on all sequences treats a genus with 8,000 near-identical
    proteins as 8,000 independent observations, which makes the p-values
    meaningless. Dereplicating is the cheap fix; the phylogenetic logistic
    regression (phyloglm.R) is the expensive, better one, and the two are meant
    to be read together.
    """
    if "taxon" not in df.columns or df["taxon"].isna().all():
        return None
    from statsmodels.stats.multitest import multipletests

    derep = df.sort_values("weight", ascending=False).drop_duplicates("cluster") \
        if "cluster" in df.columns else df
    tab = pd.crosstab(derep["subgroup"], derep["taxon"])
    rows = []
    total = tab.to_numpy().sum()
    for sg in tab.index:
        for tx in tab.columns:
            a = tab.loc[sg, tx]
            b = tab.loc[sg].sum() - a
            c = tab[tx].sum() - a
            d = total - a - b - c
            odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({"subgroup": sg, "taxon": tx, "n_derep": int(a),
                         "odds_ratio": odds, "p": p})
    out = pd.DataFrame(rows)
    out["q_bh"] = multipletests(out["p"], method="fdr_bh")[1]
    out = out.sort_values("q_bh")
    out.attrs["n_derep"] = len(derep)
    out.to_csv(os.path.join(tabdir, "subgroup_taxonomy_enrichment.tsv"),
               sep="\t", index=False)
    tab.to_csv(os.path.join(tabdir, "subgroup_by_taxon_counts.tsv"), sep="\t")
    print(f"[active_site] taxonomy enrichment on {len(derep)} dereplicated sequences "
          f"(of {len(df)})", file=sys.stderr)
    return tab


def fig_taxonomy_composition(tab, figdir, fmts, dpi):
    if tab is None or tab.empty:
        return
    frac = tab.div(tab.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(5.2, 2.6))
    bottom = np.zeros(len(frac))
    cols = palette(frac.shape[1])
    for c, col in zip(frac.columns, cols):
        ax.bar(frac.index.astype(str), frac[c], bottom=bottom, color=col,
               width=0.7, label=str(c))
        bottom += frac[c].to_numpy()
    ax.set_ylabel("Fraction of sequences")
    ax.set_xlabel("Active-site subgroup")
    ax.set_ylim(0, 1)
    ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", borderaxespad=0)
    ax.set_title("Taxonomic composition of each subgroup", loc="left")
    fig.tight_layout()
    savefig(fig, figdir, "15_subgroup_taxonomy", fmts, dpi)


def consensus_table(barcodes, labels, seq_w, cols, res, bounds, tabdir):
    rows = []
    for u in sorted(set(labels)):
        sel = np.flatnonzero(np.asarray(labels) == u)
        M = weighted_counts(barcodes, seq_w, sel)      # L x 20, weighted
        eff = float(seq_w[sel].sum())
        cons, conf = [], []
        for j in range(M.shape[0]):
            tot = M[j].sum()
            if tot <= 0:
                cons.append("-"); conf.append(0.0); continue
            k = int(np.argmax(M[j]))
            frac = M[j, k] / tot
            # upper = majority residue carries >=50% of the weighted mass
            cons.append(AA[k] if frac >= 0.5 else AA[k].lower())
            conf.append(frac)
        motifs = ["".join(cons[lo:hi]) for lo, hi in bounds]
        rows.append({
            "subgroup": u, "n": len(sel), "effective_n": round(eff, 2),
            **{f"motif_{r}{c}": m for r, c, m in zip(res, cols, motifs)},
            "consensus_barcode": "".join(cons),
            "mean_consensus_support": float(np.mean(conf)),
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(tabdir, "subgroup_consensus_motifs.tsv"), sep="\t", index=False)
    return df


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--afa", required=True)
    ap.add_argument("--chosen", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--idmap", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tabdir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-assign", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    acfg = cfg["active_site"]
    fmts = tuple(cfg["plots"]["formats"])
    dpi = cfg["plots"]["dpi"]
    seed = acfg["random_state"]
    set_style()
    os.makedirs(a.figdir, exist_ok=True)
    os.makedirs(a.tabdir, exist_ok=True)

    with open(a.chosen) as fh:
        tri = json.load(fh)
    cols, res = tri["match_columns"], tri["residues"]
    w = acfg["flank_window"]

    names, rows = zip(*read_fasta(a.afa))
    L = len(rows[0])
    idx, centres, bounds = barcode_blocks(cols, w, L)
    if len(idx) != 3 * (2 * w + 1):
        print(f"[active_site] a window was clipped at the profile edge: {len(idx)} "
              f"columns instead of {3 * (2 * w + 1)}", file=sys.stderr)
    barcodes = ["".join(r[j] for j in idx) for r in rows]

    # Downstream scripts must not re-derive the barcode geometry from the flank
    # window: a window clipped at a profile edge breaks that assumption.
    with open(os.path.join(os.path.dirname(a.chosen), "barcode_layout.json"), "w") as fh:
        json.dump({"alignment_columns": idx, "centres": centres, "bounds": bounds,
                   "flank_window": w, "triad_columns": cols, "residues": res}, fh, indent=2)

    wdf = pd.read_csv(a.weights, sep="\t").set_index("seq_id")
    seq_cluster = wdf["cluster"].reindex(list(names))
    seq_w = wdf["weight"].reindex(list(names)).to_numpy(dtype=float)
    if np.isnan(seq_w).any():
        sys.exit(f"[active_site] {int(np.isnan(seq_w).sum())} sequences lack a weight")

    print(f"[active_site] {len(barcodes)} sequences (effective {seq_w.sum():.1f} after "
          f"redundancy weighting), barcode length {len(idx)}", file=sys.stderr)

    # --- embed + PCA ---------------------------------------------------------
    M = blosum_embedding()
    weights = (conservation_weights(barcodes, seq_w)
               if acfg.get("entropy_weighting", True) else None)
    if weights is not None:
        print("[active_site] conservation weights: "
              + " ".join(f"{w:.2f}" for w in weights), file=sys.stderr)

    rng = np.random.default_rng(seed)
    n = len(barcodes)
    fit_idx = (rng.choice(n, acfg["max_cluster_seqs"], replace=False)
               if n > acfg["max_cluster_seqs"] else np.arange(n))

    X_fit = encode([barcodes[i] for i in fit_idx], M, weights)
    cap = int(min(acfg["pca_components"], X_fit.shape[0] - 1, X_fit.shape[1]))
    sel = acfg.get("pca_selection", "parallel_analysis")
    if sel == "parallel_analysis":
        ncomp = n_components_parallel_analysis(X_fit, cap, seed)
    elif sel == "variance":
        cum = np.cumsum(PCA(n_components=cap, random_state=seed)
                        .fit(X_fit).explained_variance_ratio_)
        ncomp = max(2, min(int(np.searchsorted(cum, acfg.get("pca_variance", 0.90)) + 1), cap))
    else:
        ncomp = cap
    pca = PCA(n_components=ncomp, random_state=seed).fit(X_fit)
    Z_fit = pca.transform(X_fit)
    print(f"[active_site] PCA ({sel}): {ncomp} components, "
          f"{100 * pca.explained_variance_ratio_.sum():.1f}% variance", file=sys.stderr)

    w_fit = seq_w[fit_idx]
    k_best, metrics, _, k_gap = choose_k(Z_fit, acfg["k_min"], acfg["k_max"], seed,
                                         acfg.get("k_criterion", "silhouette_cosine"),
                                         sample_weight=w_fit)
    km = KMeans(n_clusters=k_best, n_init=10, random_state=seed).fit(
        Z_fit, sample_weight=w_fit)

    hdb = run_hdbscan(Z_fit, acfg["hdbscan_min_cluster_size"]) if acfg["run_hdbscan"] else None
    hdb_k = None if hdb is None else len(set(hdb) - {-1})

    # Assign every sequence. Sequences held out of the PCA/k-means fit are
    # projected into the same basis and given their nearest centroid, which is
    # exactly what KMeans.predict does.
    Z_all = np.empty((n, ncomp), dtype=np.float32)
    step = 20000
    for s in range(0, n, step):
        Z_all[s:s + step] = pca.transform(encode(barcodes[s:s + step], M, weights))
    labels = km.predict(Z_all)

    metrics.to_csv(os.path.join(a.tabdir, "subgroup_model_selection.tsv"),
                   sep="\t", index=False)
    print(f"[active_site] k = {k_best} subgroups "
          f"({acfg.get('k_criterion', 'silhouette_cosine')}); "
          f"gap statistic says {k_gap}; HDBSCAN says {hdb_k}", file=sys.stderr)

    # --- annotate ------------------------------------------------------------
    idmap = pd.read_csv(a.idmap, sep="\t", dtype=str).set_index("seq_id")
    merged = pd.read_csv(a.merged, sep="\t", dtype={"sample": str, "protein_id": str},
                         low_memory=False)
    tax = resolve_taxonomy(merged, cfg["plots"]["taxonomy_col"], cfg["plots"]["taxonomy_rank"])
    if tax is not None:
        merged = merged.assign(_tax=tax)
        smp2tax = merged.drop_duplicates("sample").set_index("sample")["_tax"]
    else:
        smp2tax = None

    assign = pd.DataFrame({
        "seq_id": list(names),
        "subgroup": labels,
        "barcode": barcodes,
        "weight": seq_w,
        "cluster": seq_cluster.to_numpy(),
    })
    # -2 marks sequences that were not part of the HDBSCAN subsample; -1 is
    # HDBSCAN's own noise label.
    hdb_full = np.full(n, -2, dtype=int)
    if hdb is not None:
        hdb_full[fit_idx] = hdb
    assign["hdbscan"] = hdb_full
    assign["sample"] = assign["seq_id"].map(idmap["sample"])
    assign["protein_id"] = assign["seq_id"].map(idmap["protein_id"])
    assign["evidence"] = assign["seq_id"].map(idmap["evidence"])
    if smp2tax is not None:
        assign["taxon"] = top_n_with_other(assign["sample"].map(smp2tax),
                                           cfg["plots"]["max_taxa_in_legend"])
    assign.to_csv(a.out_assign, sep="\t", index=False)

    # --- figures -------------------------------------------------------------
    fig_model_selection(metrics, k_best, k_gap, hdb_k, a.figdir, fmts, dpi)
    fig_pca(Z_fit, km.labels_, hdb, a.figdir, fmts, dpi, pca.explained_variance_ratio_)
    fig_logos(barcodes, labels, seq_w, cols, res, centres, bounds, a.figdir, fmts, dpi)
    fig_entropy(barcodes, seq_w, cols, res, centres, a.figdir, a.tabdir, fmts, dpi)
    fig_heatmap_dendrogram(barcodes, labels, seed, a.figdir, fmts, dpi)

    tab = taxonomy_enrichment(assign, a.tabdir)
    fig_taxonomy_composition(tab, a.figdir, fmts, dpi)
    cons = consensus_table(barcodes, labels, seq_w, cols, res, bounds, a.tabdir)

    print(cons.to_string(index=False), file=sys.stderr)
    print(f"[active_site] outputs -> {a.tabdir}, {a.figdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
