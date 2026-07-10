#!/usr/bin/env python3
"""Specificity-determining positions, called so that the answer cannot be an
artefact of how the groups were defined.

The circularity problem
-----------------------
The active-site subgroups in this pipeline are k-means on a barcode made of the
columns flanking the catalytic triad. Running an SDP test on those same columns,
against those same groups, is guaranteed to rediscover them. The result would be
beautiful, reproducible, and empty.

So SDPs are called under partitions that never saw the barcode:

    tree_clade         clades cut from the C71 gene tree by subtree size
    ssn_cluster        connected components of the sequence similarity network
    pmbr_architecture  whether the binding module can dock on an intact sacculus
                       (>= 3 PMB motifs), or the repeat count, per
                       `pmbr_partition_mode`. Never both: they are nested.
    host_genus         GTDB genus of the genome the protein came from
    pei_class          the published four-class partition (Wang et al. 2025),
                       defined by two structural positions. Its own two defining
                       columns are excluded from the test, because one of them
                       (V252) lies inside the triad-flank barcode.

`pmbr_architecture` is a legitimate partition against a barcode drawn from the
catalytic site because the two modules read different things. The PMB domain binds
NAG in the glycan backbone -- it sticks to lysozyme-treated bacterial spheroplasts
and survives 150 mM NaCl (Visweswaran et al. 2011) -- while the groove cleaves the
Ala-epsilon-Lys isopeptide in the peptide cross-link. Glycan and peptide are a
priori independent, so agreement between them is evidence rather than tautology.

That independence is not settled. Wang et al. 2025 report the PB repeats improving
recognition of Glu-gamma-Thr/Ser and Asp-beta-Ala, which are peptide substitutions.
Either the module touches the peptide or avidity alone shifts the apparent
preference. If PMBR turns out to read the peptide, this partition is no longer
independent of the barcode and the result would need re-reading.

A position is called an SDP only if it clears FDR under at least
`sdp_min_partitions` of them. A position significant under one partition and no
other is a property of that partition, and the 4-way concordance matrix says so.

The statistic
-------------
Type II SDP score, in the GroupSim / SDPpred sense: a position is a candidate
when it is conserved *within* groups and different *between* them. Implemented
as the redundancy-weighted mutual information between the residue at a column
and the group label,

    I(residue ; group)

with the same average-product correction used elsewhere, and with the
finite-sample bias handled by the permutation null rather than by an analytic
correction (21 residue states x k groups is well outside the regime where the
Miller-Madow correction is trustworthy).

The null
--------
Group labels are shuffled WITHIN tree clades, not globally. A global shuffle
destroys phylogenetic structure, so every phylogenetically clustered position
looks like an SDP and the FDR is meaningless. Shuffling inside clades preserves
the tree and asks the only question worth asking: given how these sequences are
related, does this column know more about the group label than it should?

For `tree_clade` itself the within-clade null is degenerate, so that partition
uses a global shuffle and is flagged as such in the output. It is the weakest of
the four for exactly this reason, which is why replication is required.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import PALETTE, load_config, palette, read_fasta, savefig, set_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY-"
IDX = {a: i for i, a in enumerate(AA)}
K = len(AA)


# ---------------------------------------------------------------------------
def tree_clades(tree_path, tips, target_size):
    """Cut the tree into clades of at most `target_size` tips.

    Cutting at a fixed relative depth is what a first draft does, and on a
    star-like or very unbalanced tree it hands back 240 singleton clades. A
    within-clade permutation over singleton clades is the identity permutation:
    the null becomes the observation, every z is zero, and the test can never
    reject anything. It fails silently and looks like "no signal".

    Recursing on subtree size cannot produce that outcome unless the tree really
    is a star, and `null_is_degenerate` detects that case explicitly.
    """
    from Bio import Phylo
    t = Phylo.read(tree_path, "newick")
    try:
        t.root_at_midpoint()
    except Exception:  # noqa: BLE001
        pass

    labels, k = {}, 0
    stack = [t.root]
    while stack:
        node = stack.pop()
        leaves = node.get_terminals()
        if len(leaves) <= target_size or not node.clades:
            for leaf in leaves:
                labels[leaf.name] = f"clade{k}"
            k += 1
        else:
            stack.extend(node.clades)
    return pd.Series({t_: labels.get(t_) for t_ in tips}, dtype=object)


def null_is_degenerate(clade, labels, min_permutable=0.5):
    """Fraction of labelled tips sitting in a clade with >= 2 labelled tips.

    Below `min_permutable`, a within-clade shuffle barely moves anything and the
    null collapses onto the observation.
    """
    ok = pd.notna(labels)
    cl = clade[ok]
    if not len(cl):
        return True, 0.0
    sizes = cl.value_counts()
    permutable = sizes[sizes >= 2].sum() / len(cl)
    return permutable < min_permutable, float(permutable)


def encode(rows):
    n, L = len(rows), len(rows[0])
    A = np.full((n, L), IDX["-"], dtype=np.int64)
    for i, s in enumerate(rows):
        for j, c in enumerate(s):
            A[i, j] = IDX.get(c, IDX["-"])
    return A


def mi_col_group(col, grp, w, n_grp):
    """Weighted mutual information between residue and group label, in bits."""
    J = np.zeros((K, n_grp))
    np.add.at(J, (col, grp), w)
    tot = J.sum()
    if tot <= 0:
        return 0.0
    J /= tot
    pi, pj = J.sum(1), J.sum(0)
    nz = J > 0
    outer = np.outer(pi, pj)
    return float(np.sum(J[nz] * np.log2(J[nz] / outer[nz])))


def score_all(A, grp, w, n_grp, active):
    return np.array([mi_col_group(A[:, c], grp, w, n_grp) for c in active])


def apc(scores, A, grp, w, n_grp, active):
    """Average-product correction against the column's own entropy: subtract
    what a column of this variability would score against a random label."""
    # column entropy proxy
    ent = []
    for c in active:
        cnt = Counter(A[:, c].tolist())
        p = np.array(list(cnt.values()), float)
        p /= p.sum()
        ent.append(-(p * np.log2(p)).sum())
    ent = np.array(ent)
    if ent.sum() <= 0 or scores.mean() <= 0:
        return scores
    return scores - (ent * scores.mean() / (ent.mean() or 1))


def within_clade_permute(grp, clade, rng):
    out = grp.copy()
    for c in np.unique(clade):
        m = clade == c
        if m.sum() > 1:
            out[m] = grp[m][rng.permutation(m.sum())]
    return out


# ---------------------------------------------------------------------------
def call_sdps(A, labels, clade, w, active, n_perm, fdr, rng, null_mode):
    ok = pd.notna(labels)
    A_, lab, cl, w_ = A[ok.values], labels[ok].to_numpy(), clade[ok].to_numpy(), w[ok.values]
    groups = sorted(set(lab))
    gmap = {g: i for i, g in enumerate(groups)}
    grp = np.array([gmap[x] for x in lab])
    n_grp = len(groups)

    obs = apc(score_all(A_, grp, w_, n_grp, active), A_, grp, w_, n_grp, active)

    s1 = np.zeros(len(active))
    s2 = np.zeros(len(active))
    for _ in range(n_perm):
        gp = (within_clade_permute(grp, cl, rng) if null_mode == "within_clade"
              else grp[rng.permutation(len(grp))])
        nul = apc(score_all(A_, gp, w_, n_grp, active), A_, gp, w_, n_grp, active)
        s1 += nul
        s2 += nul ** 2
    mu = s1 / n_perm
    sd = np.sqrt(np.maximum(s2 / n_perm - mu ** 2, 1e-12))

    from scipy.stats import norm
    z = (obs - mu) / sd
    p = norm.sf(z)
    q = multipletests(p, method="fdr_bh")[1]
    return pd.DataFrame({"match_col": active, "mi_bits": obs, "null_mean": mu,
                         "null_sd": sd, "z": z, "p": p, "q_bh": q,
                         "significant": q < fdr, "n_groups": n_grp,
                         "n_sequences": int(ok.sum())})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--afa", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--arch", required=True)
    ap.add_argument("--ssn", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--idmap", required=True)
    ap.add_argument("--groove", required=True)
    ap.add_argument("--groove-json")
    ap.add_argument("--pei-class")
    ap.add_argument("--weights", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-concordance", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    scfg = cfg["specificity"]
    fmts, dpi = tuple(cfg["plots"]["formats"]), cfg["plots"]["dpi"]
    set_style()
    os.makedirs(a.figdir, exist_ok=True)
    for p in (a.out, a.out_concordance):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    names, rows = zip(*read_fasta(a.afa))
    A = encode(rows)
    L = A.shape[1]

    wdf = pd.read_csv(a.weights, sep="\t").set_index("seq_id")
    w = wdf["weight"].reindex(names).to_numpy(float)
    cluster = wdf["cluster"].reindex(names)

    assign = pd.read_csv(a.assign, sep="\t", dtype={"seq_id": str}).set_index("seq_id")
    arch = pd.read_csv(a.arch, sep="\t", dtype={"seq_id": str}).set_index("seq_id")
    ssn = pd.read_csv(a.ssn, sep="\t", dtype={"seq_id": str}).set_index("seq_id")
    idmap = pd.read_csv(a.idmap, sep="\t", dtype=str).set_index("seq_id")
    merged = pd.read_csv(a.merged, sep="\t", dtype={"sample": str}, low_memory=False)

    tips = list(names)
    target = scfg.get("sdp_clade_target_size") or max(5, len(tips) // 40)
    target = int(target)
    clade = tree_clades(a.tree, tips, target)
    # sequences absent from the tree (not representatives) get their cluster as a
    # stand-in clade, so the within-clade null still has structure to preserve
    clade = clade.fillna(pd.Series({t: f"clu_{cluster.get(t)}" for t in tips}))
    csz = clade.value_counts()
    print(f"[sdp] {clade.nunique()} clades, target size {target}, "
          f"median size {int(csz.median())}, {int((csz == 1).sum())} singletons",
          file=sys.stderr)

    # --- the four partitions -------------------------------------------------
    genus = None
    if "gtdb_genus" in merged.columns:
        s2g = merged.drop_duplicates("sample").set_index("sample")["gtdb_genus"]
        genus = pd.Series({t: s2g.get(idmap["sample"].get(t)) for t in tips})

    # Two ways to partition on the binding module, and they must not both be used.
    #
    #   count             one group per repeat number
    #   binding_competent two groups, split at the 3-motif cliff where the
    #                     Visweswaran 2011 constructs stop binding pseudomurein
    #
    # These are nested. A column that separates 2-motif from 3-motif proteins
    # separates `pmbr2` from `pmbr3` too, so scoring it as replicated across both
    # partitions would count one observation twice and inflate `n_partitions`.
    # `sdp_min_partitions` would then be satisfied by a single piece of evidence.
    mode = str(scfg.get("pmbr_partition_mode", "binding_competent"))
    if mode not in ("count", "binding_competent"):
        sys.exit(f"[sdp] pmbr_partition_mode must be 'count' or "
                 f"'binding_competent', got {mode!r}")

    def _accessory(t):
        v = arch["accessory_binding_domains"].get(t)
        return ("+" + v) if isinstance(v, str) and v else ""

    if mode == "count":
        pmbr_arch = pd.Series({
            t: f"pmbr{int(arch['n_pmbr'].get(t, 0))}{_accessory(t)}"
            if t in arch.index else None for t in tips})
    else:
        if "pmbr_binding_competent" not in arch.columns:
            sys.exit("[sdp] domain_architecture.tsv has no "
                     "`pmbr_binding_competent` column. Re-run domain_arch.py: the "
                     "3-motif binding threshold is where the biology is.")
        pmbr_arch = pd.Series({
            t: (("sacculus_binder" if int(arch["pmbr_binding_competent"].get(t, 0))
                 else "cannot_bind_sacculus") + _accessory(t))
            if t in arch.index else None for t in tips})
        # Proteins whose count straddles the cliff have no defensible group.
        if "pmbr_count_fragile" in arch.columns:
            frag = [t for t in tips
                    if t in arch.index and int(arch["pmbr_count_fragile"].get(t, 0))]
            for t in frag:
                pmbr_arch[t] = None
            if frag:
                print(f"[sdp] {len(frag)} proteins dropped from the "
                      f"pmbr_architecture partition: their repeat count changes "
                      f"binding class between the strict and permissive E-values, "
                      f"so assigning them to either group would be assigning them "
                      f"to a threshold, not to a phenotype.", file=sys.stderr)

    print(f"[sdp] pmbr_architecture partition mode: {mode}", file=sys.stderr)

    parts = {
        "tree_clade": clade,
        "ssn_cluster": pd.Series({t: ssn["ssn_cluster"].get(t) for t in tips}),
        "pmbr_architecture": pmbr_arch,
        "host_genus": genus if genus is not None else pd.Series(dtype=object),
    }

    # The published four-class partition (Wang et al. 2025), if it was computed.
    # It is defined by exactly two alignment columns, and one of them (V252) is
    # two residues from the catalytic Asp, so it lies inside the +/-5 triad-flank
    # barcode. Using it as an SDP partition without excluding its own defining
    # columns would rediscover them. `excluded` carries that through.
    excluded = set()
    if a.pei_class and os.path.exists(a.pei_class):
        pcl = pd.read_csv(a.pei_class, sep="\t", dtype={"seq_id": str}).set_index("seq_id")
        if "pei_class" in pcl.columns:
            lab = pd.Series({t: pcl["pei_class"].get(t) for t in tips})
            parts["pei_class"] = lab.where(lab != "unassigned")
        gdef = json.load(open(a.groove_json)) if a.groove_json else {}
        for k in ("position_1", "position_2"):
            cpos = (gdef.get("class_positions") or {}).get(k)
            if cpos:
                excluded.add(int(cpos["match_col"]))
        if excluded:
            print(f"[sdp] excluding the class-defining columns {sorted(excluded)} from "
                  f"the SDP test: pei_class is defined by them, so leaving them in "
                  f"would be circular", file=sys.stderr)

    parts = {k: v for k, v in parts.items() if k in scfg["sdp_partitions"]}

    # drop groups that are too small, then partitions with too few groups
    usable = {}
    for name, lab in parts.items():
        if lab.empty or lab.isna().all():
            print(f"[sdp] partition '{name}' unavailable; skipping", file=sys.stderr)
            continue
        vc = lab.value_counts()
        big = set(vc[vc >= scfg["sdp_min_group_size"]].index)
        if len(big) > scfg["sdp_max_groups"]:
            big = set(vc[vc >= scfg["sdp_min_group_size"]].index[:scfg["sdp_max_groups"]])
        lab = lab.where(lab.isin(big))
        if lab.notna().sum() < 30 or lab.nunique() < 2:
            print(f"[sdp] partition '{name}': {lab.nunique()} usable groups; skipping",
                  file=sys.stderr)
            continue
        usable[name] = lab
        print(f"[sdp] partition '{name}': {lab.nunique()} groups, "
              f"{int(lab.notna().sum())} sequences", file=sys.stderr)

    if len(usable) < scfg["sdp_min_partitions"]:
        sys.exit(f"[sdp] only {len(usable)} usable partitions; "
                 f"sdp_min_partitions is {scfg['sdp_min_partitions']}. "
                 f"Replication is the whole point; refusing to report single-"
                 f"partition SDPs.")

    # informative columns only
    occ = np.array([(A[:, c] != IDX["-"]).mean() for c in range(L)])
    varied = np.array([len(np.unique(A[:, c])) > 1 for c in range(L)])
    active = np.flatnonzero((occ >= 0.5) & varied)
    print(f"[sdp] {len(active)}/{L} informative columns", file=sys.stderr)

    rng = np.random.default_rng(cfg["active_site"]["random_state"])
    per = {}
    for name, lab in usable.items():
        if name == "tree_clade":
            # a within-clade shuffle of clade labels is the identity by definition
            null_mode = "global"
        else:
            null_mode = scfg["sdp_null"]
            if null_mode == "within_clade":
                degen, frac = null_is_degenerate(clade, lab)
                if degen:
                    print(f"[sdp] partition '{name}': only {100 * frac:.0f}% of tips "
                          f"sit in a clade with >=2 members. A within-clade shuffle "
                          f"would be close to the identity and the test would have "
                          f"no power. Falling back to a GLOBAL null, which does not "
                          f"control for phylogeny; treat these SDPs as upper bounds.",
                          file=sys.stderr)
                    null_mode = "global_fallback"
        # a partition never gets to test the columns that define it
        cols_here = np.array([c for c in active if not (name == "pei_class"
                                                        and c in excluded)])
        df = call_sdps(A, lab, clade, w, cols_here, int(scfg["sdp_permutations"]),
                       float(scfg["sdp_fdr"]), rng,
                       "within_clade" if null_mode == "within_clade" else "global")
        df["partition"] = name
        df["null_mode"] = null_mode
        df["columns_excluded"] = (",".join(map(str, sorted(excluded)))
                                  if name == "pei_class" else "")
        per[name] = df
        print(f"[sdp] {name}: {int(df['significant'].sum())} significant columns "
              f"(null={null_mode})", file=sys.stderr)

    allsdp = pd.concat(per.values(), ignore_index=True)

    # --- replication ---------------------------------------------------------
    hits = (allsdp[allsdp["significant"]]
            .groupby("match_col")["partition"].agg(list))
    rep = pd.DataFrame({
        "match_col": hits.index.to_numpy(),
        "n_partitions": hits.map(len).to_numpy(),
        "partitions": hits.map(lambda x: ",".join(sorted(x))).to_numpy(),
    })
    if rep.empty:
        rep = pd.DataFrame(columns=["match_col", "n_partitions", "partitions"])
    rep["replicated"] = (rep["n_partitions"] >= int(scfg["sdp_min_partitions"])
                         if len(rep) else pd.Series(dtype=bool))

    groove = pd.read_csv(a.groove, sep="\t")
    gcols = ["match_col", "in_groove", "is_triad", "structure_resnum",
             "structure_residue", "distance_to_seed_a"]
    gcols += [c for c in ("in_groove_literature", "in_groove_any",
                          "metal_coordinating", "in_metal_shell",
                          "distance_to_metal_a") if c in groove.columns]
    rep = rep.merge(groove[gcols], on="match_col", how="left")
    rep = rep.sort_values(["replicated", "n_partitions"], ascending=False)
    rep.to_csv(a.out, sep="\t", index=False)

    allsdp.to_csv(a.out.replace(".tsv", "_per_partition.tsv"), sep="\t", index=False)

    # concordance matrix
    names_ = sorted(usable)
    C = pd.DataFrame(0, index=names_, columns=names_)
    sig = {n: set(per[n].loc[per[n]["significant"], "match_col"]) for n in names_}
    for i in names_:
        for j in names_:
            u = len(sig[i] | sig[j])
            C.loc[i, j] = len(sig[i] & sig[j]) / u if u else np.nan
    C.to_csv(a.out_concordance, sep="\t")

    n_rep = int(rep["replicated"].sum())
    n_groove = int((rep["replicated"] & (rep["in_groove"] == 1)).sum())
    n_tri = int((rep["replicated"] & (rep["is_triad"] == 1)).sum())
    print(f"\n[sdp] {n_rep} columns replicate in >= {scfg['sdp_min_partitions']} "
          f"partitions; {n_groove} of them are in the substrate groove, "
          f"{n_tri} are the catalytic triad itself", file=sys.stderr)

    # enrichment of replicated SDPs in a structural region: Fisher, columns as
    # units. Run once per region. The groove is where substrate specificity
    # should live; the metal shell is where the PeiW/PeiP cation difference
    # should live. They are different predictions and get separate tests, so a
    # hit in one is not laundered into evidence for the other.
    from scipy.stats import fisher_exact

    is_rep = pd.Series(0, index=active)
    is_rep.loc[is_rep.index.isin(rep.loc[rep["replicated"], "match_col"])] = 1

    def enrich(region):
        if region not in groove.columns:
            return np.nan, np.nan, 0
        g = groove.set_index("match_col")[region].reindex(active).fillna(0)
        n_in = int((is_rep.to_numpy() & g.to_numpy().astype(int)).sum())
        tab = pd.crosstab(is_rep, g)
        if tab.shape != (2, 2):
            return np.nan, np.nan, n_in
        odds, p = fisher_exact(tab.to_numpy(), alternative="greater")
        return odds, p, n_in

    odds, p, _ = enrich("in_groove")
    if np.isfinite(p):
        print(f"[sdp] replicated SDPs are enriched in the groove: "
              f"OR = {odds:.2f}, p = {p:.3g}", file=sys.stderr)

    m_odds, m_p, n_metal = enrich("in_metal_shell")
    n_shell = int(groove["in_metal_shell"].sum()) if "in_metal_shell" in groove else 0
    summary = {"n_replicated": n_rep, "n_replicated_in_groove": n_groove,
               "n_replicated_is_triad": n_tri, "groove_enrichment_or": odds,
               "groove_enrichment_p": p, "n_partitions_used": len(usable),
               "partitions": ",".join(names_)}
    if n_shell == 0:
        print("[sdp] no metal shell in the structure, so no metal-shell "
              "enrichment test. PeiW and PeiP both need a divalent cation and "
              "differ in which ones work (Subedi et al. 2015); the site is "
              "simply not resolved in the deposited coordinates.",
              file=sys.stderr)
        summary.update({"n_metal_shell_columns": 0,
                        "metal_shell_enrichment_or": "not_tested",
                        "metal_shell_enrichment_p": "not_tested"})
    else:
        print(f"[sdp] {n_metal} replicated SDPs lie in the {n_shell}-column "
              f"metal shell: OR = {m_odds:.2f}, p = {m_p:.3g}", file=sys.stderr)
        summary.update({"n_replicated_in_metal_shell": n_metal,
                        "n_metal_shell_columns": n_shell,
                        "metal_shell_enrichment_or": m_odds,
                        "metal_shell_enrichment_p": m_p})

    pd.Series(summary).rename("value").rename_axis(
        "metric").to_csv(a.out.replace(".tsv", "_summary.tsv"), sep="\t")

    # --- figure --------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 3.4), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    ax = axes[0]
    cols = palette(len(names_))
    for (n_, c) in zip(names_, cols):
        d = per[n_]
        ax.plot(d["match_col"], -np.log10(d["q_bh"].clip(lower=1e-300)), lw=0.7,
                color=c, label=f"{n_} ({int(d['significant'].sum())})")
    ax.axhline(-np.log10(float(scfg["sdp_fdr"])), color="0.4", lw=0.6, ls="--")
    ax.set_ylabel("$-\\log_{10}$ q")
    ax.legend(ncol=2, fontsize=6)
    ax.set_title("A  Specificity-determining positions, per partition", loc="left")

    ax = axes[1]
    gcols = groove.loc[groove["in_groove"] == 1, "match_col"]
    ax.vlines(gcols, 0, 1, color=PALETTE[5], lw=0.6, alpha=0.5, label="groove")
    r = rep[rep["replicated"]]
    ax.vlines(r["match_col"], 0, 1, color=PALETTE[4], lw=1.2,
              label=f"replicated SDP ({n_rep})")
    ax.vlines(groove.loc[groove["is_triad"] == 1, "match_col"], 0, 1,
              color="k", lw=1.2, label="catalytic triad")
    ax.set_yticks([])
    ax.set_xlabel("PF12386 match-state column")
    ax.legend(ncol=3, fontsize=6, loc="upper right")
    ax.set_title(f"B  Groove overlap (Fisher OR = {odds:.2f}, p = {p:.3g})", loc="left")
    fig.tight_layout()
    savefig(fig, a.figdir, "24_specificity_determining_positions", fmts, dpi)

    fig, ax = plt.subplots(figsize=(3.0, 2.6))
    im = ax.imshow(C.to_numpy(), cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(names_))); ax.set_xticklabels(names_, rotation=45, ha="right")
    ax.set_yticks(range(len(names_))); ax.set_yticklabels(names_)
    for i in range(len(names_)):
        for j in range(len(names_)):
            ax.text(j, i, f"{C.iloc[i, j]:.2f}", ha="center", va="center",
                    fontsize=6, color="w" if C.iloc[i, j] < 0.6 else "k")
    fig.colorbar(im, ax=ax, fraction=0.046).set_label("Jaccard of SDP sets")
    ax.set_title("SDP concordance across partitions", loc="left")
    fig.tight_layout()
    savefig(fig, a.figdir, "25_sdp_concordance", fmts, dpi)


if __name__ == "__main__":
    main()
