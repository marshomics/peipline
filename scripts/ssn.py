#!/usr/bin/env python3
"""Sequence similarity network, EFI-EST convention.

Edge weight is the alignment score, -log10(E-value), from an all-vs-all DIAMOND
search. Raising the threshold prunes edges; the network fragments; the number of
connected components rises. The point of the sweep is that there is no single
"right" threshold, and a subgrouping that only exists at one threshold is not a
subgrouping.

The knee of the components-vs-threshold curve is picked automatically (maximum
distance from the chord joining the curve's endpoints, i.e. the Kneedle rule)
and can be overridden in config.

Why bother when there is already a tree and a k-means?
Because all three make different assumptions. The tree assumes a substitution
model and vertical descent. The k-means assumes convex clusters in a BLOSUM
embedding. The SSN assumes only that homologous sequences align well. Where all
three agree the subgroup is real; where they disagree, the disagreement is the
result. `ssn_concordance.tsv` is that comparison.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import PALETTE, load_config, palette, savefig, set_style, top_n_with_other  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

M8_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
FLOOR = 1e-300     # DIAMOND reports E=0 for identical sequences


def read_nodes(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def load_edges(path, nodes):
    df = pd.read_csv(path, sep="\t", names=M8_COLS, dtype={"qseqid": str, "sseqid": str})
    df = df[df["qseqid"] != df["sseqid"]]
    keep = set(nodes)
    df = df[df["qseqid"].isin(keep) & df["sseqid"].isin(keep)]
    # Undirected: keep one row per pair, the best (lowest E).
    pair = df[["qseqid", "sseqid"]]
    df = df.assign(_u=pair.min(axis=1), _v=pair.max(axis=1))
    df = df.sort_values("evalue").drop_duplicates(["_u", "_v"], keep="first")
    df["alignment_score"] = -np.log10(df["evalue"].clip(lower=FLOOR))
    return df[["_u", "_v", "alignment_score", "pident", "bitscore"]].rename(
        columns={"_u": "source", "_v": "target"})


def components(nodes, edges):
    """Union-find. igraph would do this too, but this keeps the sweep dependency
    free and it is O(E alpha(V))."""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    return {n: find(n) for n in nodes}


def write_graphml(path, nodes, clusters, edges):
    """Minimal GraphML, so the deliverable exists even without python-igraph."""
    import xml.sax.saxutils as sx
    attrs = [c for c in ("ssn_cluster", "subgroup", "evidence") if c in clusters.columns]
    with open(path, "w") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n'
                 '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
        for i, k in enumerate(attrs):
            fh.write(f'<key id="d{i}" for="node" attr.name="{k}" attr.type="string"/>\n')
        fh.write('<key id="e0" for="edge" attr.name="alignment_score" attr.type="double"/>\n')
        fh.write('<graph edgedefault="undirected">\n')
        by_id = clusters.set_index("seq_id")
        for n in nodes:
            fh.write(f'<node id="{sx.quoteattr(n)[1:-1]}">')
            for i, k in enumerate(attrs):
                fh.write(f'<data key="d{i}">{sx.escape(str(by_id.at[n, k]))}</data>')
            fh.write("</node>\n")
        for j, (u, v, s) in enumerate(zip(edges["source"], edges["target"],
                                          edges["alignment_score"])):
            fh.write(f'<edge id="e{j}" source="{sx.quoteattr(u)[1:-1]}" '
                     f'target="{sx.quoteattr(v)[1:-1]}">'
                     f'<data key="e0">{s:.4f}</data></edge>\n')
        fh.write("</graph>\n</graphml>\n")


def knee(x, y):
    """Kneedle: the point of maximum distance from the chord between endpoints."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) < 3:
        return x[0]
    xn = (x - x.min()) / max(float(np.ptp(x)), 1e-12)
    yn = (y - y.min()) / max(float(np.ptp(y)), 1e-12)
    d = np.abs(yn - (xn * (yn[-1] - yn[0]) + yn[0]))
    return float(x[int(np.argmax(d))])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--m8", required=True)
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-clusters", required=True)
    ap.add_argument("--out-sweep", required=True)
    ap.add_argument("--out-graphml", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tabdir", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    scfg = cfg["ssn"]
    fmts, dpi = tuple(cfg["plots"]["formats"]), cfg["plots"]["dpi"]
    set_style()
    for p in (a.out_clusters, a.out_sweep, a.out_graphml):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    nodes = read_nodes(a.nodes)
    edges = load_edges(a.m8, nodes)
    print(f"[ssn] {len(nodes)} nodes, {len(edges)} unique edges "
          f"(E <= {scfg['evalue_max']})", file=sys.stderr)
    if len(edges) == 0:
        sys.exit("[ssn] no edges: every sequence is its own component. Loosen evalue_max.")

    # --- threshold sweep -----------------------------------------------------
    ts = np.arange(scfg["threshold_min"], scfg["threshold_max"] + 1e-9,
                   scfg["threshold_step"])
    sw = []
    for t in ts:
        e = edges[edges["alignment_score"] >= t]
        comp = components(nodes, zip(e["source"], e["target"]))
        sizes = pd.Series(list(comp.values())).value_counts()
        sw.append({"threshold": float(t), "n_edges": len(e),
                   "n_components": int(sizes.shape[0]),
                   "largest_component": int(sizes.iloc[0]),
                   "n_singletons": int((sizes == 1).sum())})
    sweep = pd.DataFrame(sw)
    sweep.to_csv(a.out_sweep, sep="\t", index=False)

    thr = scfg["threshold"]
    if thr is None:
        # A knee on a flat curve can land where every node is isolated, which is
        # not a clustering. Restrict the search to thresholds that leave at least
        # one real component behind.
        usable = sweep[sweep["largest_component"] >= 2]
        if usable.empty:
            sys.exit("[ssn] every threshold leaves only singletons. The all-vs-all "
                     "found no significant similarity between any two sequences: "
                     "check ssn.evalue_max and that the aligner actually ran.")
        thr = knee(usable["threshold"], usable["n_components"])
        print(f"[ssn] knee at alignment score {thr:g}", file=sys.stderr)
    thr = float(thr)

    e = edges[edges["alignment_score"] >= thr]
    comp = components(nodes, zip(e["source"], e["target"]))
    sizes = pd.Series(list(comp.values())).value_counts()
    remap = {r: i for i, r in enumerate(sizes.index)}       # 0 = largest
    clusters = pd.DataFrame({"seq_id": nodes,
                             "ssn_cluster": [remap[comp[n]] for n in nodes]})
    clusters["ssn_cluster_size"] = clusters["ssn_cluster"].map(
        clusters["ssn_cluster"].value_counts())
    clusters["ssn_threshold"] = thr

    assign = pd.read_csv(a.assign, sep="\t", dtype={"seq_id": str})
    clusters = clusters.merge(
        assign[[c for c in ("seq_id", "subgroup", "evidence", "taxon") if c in assign.columns]],
        on="seq_id", how="left")
    clusters.to_csv(a.out_clusters, sep="\t", index=False)
    n_singleton = int((sizes == 1).sum())
    print(f"[ssn] {sizes.shape[0]} clusters at score >= {thr:g}; "
          f"largest holds {int(sizes.iloc[0])} nodes, {n_singleton} singletons",
          file=sys.stderr)
    if n_singleton > 0.8 * len(nodes):
        print("[ssn] WARNING: over 80% of nodes are singletons at the chosen "
              "threshold. The network is not telling you anything; lower "
              "ssn.threshold or loosen ssn.evalue_max.", file=sys.stderr)

    # --- concordance ---------------------------------------------------------
    if "subgroup" in clusters.columns and clusters["subgroup"].notna().any():
        ct = pd.crosstab(clusters["ssn_cluster"], clusters["subgroup"])
        ct.to_csv(os.path.join(a.tabdir, "ssn_concordance.tsv"), sep="\t")
        from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
        ok = clusters["subgroup"].notna()
        ari = adjusted_rand_score(clusters.loc[ok, "ssn_cluster"], clusters.loc[ok, "subgroup"])
        ami = adjusted_mutual_info_score(clusters.loc[ok, "ssn_cluster"],
                                         clusters.loc[ok, "subgroup"])
        pd.Series({"adjusted_rand": ari, "adjusted_mutual_info": ami,
                   "n_ssn_clusters": int(sizes.shape[0]),
                   "n_active_site_subgroups": int(clusters["subgroup"].nunique()),
                   "threshold": thr}).rename("value").rename_axis("metric").to_csv(
            os.path.join(a.tabdir, "ssn_vs_subgroup_agreement.tsv"), sep="\t")
        print(f"[ssn] SSN clusters vs active-site subgroups: ARI={ari:.3f} AMI={ami:.3f}",
              file=sys.stderr)

    # --- graphml + layout ----------------------------------------------------
    idx = {n: i for i, n in enumerate(nodes)}
    coords = None
    try:
        import igraph as ig
        g = ig.Graph(n=len(nodes),
                     edges=[(idx[u], idx[v]) for u, v in zip(e["source"], e["target"])])
        g.vs["name"] = nodes
        g.vs["ssn_cluster"] = clusters["ssn_cluster"].tolist()
        if "subgroup" in clusters.columns:
            g.vs["subgroup"] = clusters["subgroup"].fillna(-1).astype(int).tolist()
        if "evidence" in clusters.columns:
            g.vs["evidence"] = clusters["evidence"].fillna("NA").astype(str).tolist()
        g.es["alignment_score"] = e["alignment_score"].tolist()
        g.write_graphml(a.out_graphml)
        # Seed igraph's RNG so the force-directed layout (drl) is reproducible;
        # otherwise the published SSN figure's node coordinates change run to run.
        _seed = int((cfg.get("active_site") or {}).get("random_state", 12345))
        try:
            ig.set_random_number_generator(random.Random(_seed))
        except Exception:  # noqa: BLE001
            random.seed(_seed)
        coords = np.array(g.layout(scfg["layout"]).coords)
    except Exception as exc:  # noqa: BLE001
        print(f"[ssn] igraph unavailable or layout failed ({exc}); writing plain GraphML "
              f"and skipping the network plot. Open it in Cytoscape.", file=sys.stderr)
        write_graphml(a.out_graphml, nodes, clusters, e)

    # --- figures -------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.2))
    ax = axes[0]
    ax.plot(sweep["threshold"], sweep["n_components"], "-", color=PALETTE[0])
    ax.axvline(thr, color=PALETTE[4], ls="--", lw=0.8)
    ax.set_xlabel("Alignment score, $-\\log_{10}E$")
    ax.set_ylabel("Connected components")
    ax.set_title("A  Threshold sweep", loc="left")

    ax = axes[1]
    ax.plot(sweep["threshold"], sweep["largest_component"], "-", color=PALETTE[0],
            label="largest component")
    ax.plot(sweep["threshold"], sweep["n_singletons"], "-", color=PALETTE[1],
            label="singletons")
    ax.axvline(thr, color=PALETTE[4], ls="--", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Alignment score")
    ax.set_ylabel("Nodes")
    ax.legend()
    ax.set_title("B  Fragmentation", loc="left")

    ax = axes[2]
    cs = clusters["ssn_cluster_size"].drop_duplicates().sort_values(ascending=False)
    ax.plot(np.arange(1, len(cs) + 1), cs.to_numpy(), "o-", ms=2, color=PALETTE[2])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Cluster rank"); ax.set_ylabel("Cluster size")
    ax.set_title(f"C  {len(cs)} clusters at score {thr:g}", loc="left")
    fig.tight_layout()
    savefig(fig, a.figdir, "17_ssn_threshold_sweep", fmts, dpi)

    if coords is not None and len(nodes) > 1:
        for colour_by in ("subgroup", "ssn_cluster"):
            if colour_by not in clusters.columns:
                continue
            lab = clusters[colour_by].fillna(-1)
            if colour_by == "ssn_cluster":
                lab = top_n_with_other(lab.astype(str), cfg["plots"]["max_taxa_in_legend"])
            uniq = sorted(map(str, pd.unique(lab)))
            cmap = dict(zip(uniq, palette(len(uniq))))
            fig, ax = plt.subplots(figsize=(4.6, 4.6))
            src = np.array([[coords[idx[u]], coords[idx[v]]]
                            for u, v in zip(e["source"], e["target"])])
            from matplotlib.collections import LineCollection
            ax.add_collection(LineCollection(src, colors="0.85", linewidths=0.15,
                                             rasterized=True, zorder=1))
            for u in uniq:
                m = (lab.astype(str) == u).to_numpy()
                ax.scatter(coords[m, 0], coords[m, 1], s=4, color=cmap[u], linewidths=0,
                           label=f"{u} (n={int(m.sum()):,})", zorder=2, rasterized=True)
            ax.set_aspect("equal"); ax.axis("off")
            ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), markerscale=2,
                      title=colour_by.replace("_", " "))
            ax.set_title(f"SSN, alignment score $\\geq$ {thr:g}", loc="left")
            fig.tight_layout()
            savefig(fig, a.figdir, f"18_ssn_by_{colour_by}", fmts, dpi)


if __name__ == "__main__":
    main()
