#!/usr/bin/env python3
"""Selection on the substrate groove versus the domain core.

The hypothesis is a phage-host arms race over cell-wall chemistry: the catalytic
triad should be under strong purifying selection (it cannot change), and the
groove residues that read P1 should show episodic diversifying selection (they
must change when the host changes its cross-link).

Pipeline: recover the CDS for every triad-positive protein, build a codon
alignment by threading the nucleotides through the existing protein alignment
(pal2nal semantics, implemented here so the match-state coordinates are
preserved exactly), then hand it to HyPhy.

  FEL    site-wise dN/dS, both directions. Used for the groove-vs-core contrast.
  MEME   episodic diversifying selection at individual sites. The right test for
         an arms race, where positive selection is transient and lineage-specific.
  RELAX  NOT IMPLEMENTED. It is described in the literature as the right test for
         "are prophage-encoded copies under different pressure from host-encoded
         ones", and `specificity.selection.run_relax` used to sit in config.yaml
         as if it did something. It never ran: the HyPhy loop iterates only
         `hyphy_methods`, `--prophage` was declared and never read, and RELAX
         needs an explicit `--test` branch set that was never built. Rather than
         leave a flag that misrepresents the run, `run_relax: true` is now a hard
         error. To implement it: label tips prophage/host from geNomad, write the
         branch set, and call `hyphy relax --test <set> --reference <complement>`.
         The partition is defined by genomic context, independent of dN/dS, so the
         design is sound -- it is simply absent.

Three honest limits, all reported rather than buried.

CDS recovery is by exact translated match: a CDS is accepted only if its
translation equals the protein. Prokka and Prodigal both emit `.ffn` keyed on
the same locus tag, but a genome re-annotated since the `.faa` was written will
not match, and a silent mismatch here would corrupt every dN/dS estimate.

Sequences are dereplicated at the nucleotide level. Identical CDS contribute no
substitutions, and HyPhy's runtime is superlinear in tips.

The groove/core contrast is a comparison of two site classes within one
alignment, so it is a paired test on the same tree and the same model. It is not
a comparison across two separate HyPhy runs, which would confound the site class
with the model fit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (PALETTE, load_config, match_columns, read_fasta,  # noqa: E402
                   read_stockholm, savefig, seq_id,
                   set_style, write_fasta)

import matplotlib.pyplot as plt  # noqa: E402

CODON_TABLE_11_STOPS = {"TAA", "TAG", "TGA"}
GENCODE = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in zip((a + b + c for a in _bases for b in _bases for c in _bases), _aas):
    GENCODE[_i] = _a


def translate(nt):
    nt = nt.upper().replace("U", "T")
    return "".join(GENCODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def sh(cmd, **kw):
    print("+ " + " ".join(map(str, cmd)), file=sys.stderr)
    return subprocess.run(list(map(str, cmd)), check=True, **kw)


def ffn_for(faa, prodigal_dir):
    cand = os.path.splitext(faa)[0] + ".ffn"
    if os.path.exists(cand):
        return cand
    stem = os.path.splitext(os.path.basename(faa))[0]
    cand = os.path.join(prodigal_dir, stem + ".ffn")
    return cand if os.path.exists(cand) else None


def recover_cds(idmap, keep_ids, prodigal_dir, prot_seq):
    """protein seq_id -> CDS, accepted only on an exact translated match."""
    want = defaultdict(dict)
    for _, r in idmap[idmap["seq_id"].isin(keep_ids)].iterrows():
        want[r["faa"]][r["protein_id"]] = r["seq_id"]

    out, stats = {}, {"no_ffn": 0, "not_found": 0, "mismatch": 0, "ok": 0,
                      "internal_stop": 0}
    for faa, m in want.items():
        ffn = ffn_for(faa, prodigal_dir)
        if not ffn:
            stats["no_ffn"] += len(m)
            continue
        found = set()
        for h, s in read_fasta(ffn):
            pid = seq_id(h)
            if pid not in m:
                continue
            found.add(pid)
            sid = m[pid]
            nt = s.upper().replace("U", "T")
            if len(nt) % 3:
                nt = nt[:len(nt) - len(nt) % 3]
            aa = translate(nt).rstrip("*")
            if "*" in aa:
                stats["internal_stop"] += 1
                continue
            if aa != prot_seq.get(sid, ""):
                stats["mismatch"] += 1
                continue
            out[sid] = nt[:3 * len(aa)]
            stats["ok"] += 1
        stats["not_found"] += len(set(m) - found)
    return out, stats


def matchcol_residue_index(sto_path, full_prot):
    """Per sequence: match column -> 0-based residue index in the FULL protein.

    This is the piece that was missing, and its absence made the whole module a
    no-op. `triad_pass_matchcols.afa` holds only PF12386 match states: insert
    residues and anything `hmmalign --trim` clipped from the termini are gone.
    Two consequences, both fatal:

      1. Ungapping a match-column row does NOT give you the protein, so an exact
         translated match against the CDS can never succeed.
      2. Threading a full CDS through a match-column row by consuming one codon
         per non-gap column silently misassigns codons for any protein with an
         insertion. In frame, wrong position.

    So walk the Stockholm instead. Insert residues (lowercase, or in non-RF
    columns) consume a residue but occupy no match column. The trimmed fragment
    is located inside the full protein by substring search, exactly as
    groove_map.py does, because `--trim` clips the ends.
    """
    order, seqs, rf = read_stockholm(sto_path)
    if not rf:
        sys.exit("[selection] no '#=GC RF' line in the Stockholm file.")
    cols = match_columns(rf)
    matchset = {c: i for i, c in enumerate(cols)}
    L = len(cols)

    out, problems = {}, []
    for name in order:
        prot = full_prot.get(name)
        if prot is None:
            continue
        aligned = seqs[name]
        frag = "".join(ch for ch in aligned if ch not in "-.").upper()
        off = prot.find(frag)
        if off < 0:
            problems.append(name)
            continue
        idx = [None] * L
        k = off
        for pos, ch in enumerate(aligned):
            if ch in "-.":
                continue
            if pos in matchset:
                idx[matchset[pos]] = k
            k += 1                      # insert residues consume the protein too
        out[name] = idx
    if problems:
        print(f"[selection] {len(problems)} aligned sequences are not substrings of "
              f"their own protein (e.g. {problems[:3]}). hmmalign rewrote them, or "
              f"c71.faa and hits.sto disagree.", file=sys.stderr)
    return out, L


def codon_align(colmap, cds, L):
    """Build a codon alignment on PF12386 match columns.

    Each match column takes the codon of the residue that actually occupies it,
    looked up by residue index. Never by sequential consumption.
    """
    out, dropped = {}, 0
    for sid, idx in colmap.items():
        nt = cds.get(sid)
        if nt is None:
            continue
        row, ok = [], True
        for c in range(L):
            i = idx[c]
            if i is None:
                row.append("---")
                continue
            cod = nt[3 * i:3 * i + 3]
            if len(cod) < 3:
                ok = False
                break
            row.append(cod)
        if ok:
            out[sid] = "".join(row)
        else:
            dropped += 1
    if dropped:
        print(f"[selection] {dropped} sequences dropped: a match column indexes a "
              f"residue past the end of its CDS.", file=sys.stderr)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sto", required=True,
                    help="hits.sto -- the ONLY place match columns can be mapped "
                         "back to residue indices")
    ap.add_argument("--faa", required=True,
                    help="c71.faa -- full protein sequences, not match-column rows")
    ap.add_argument("--idmap", required=True)
    ap.add_argument("--keep", required=True)
    ap.add_argument("--groove", required=True)
    ap.add_argument("--chosen", required=True)
    ap.add_argument("--prophage", default=None,
                    help="geNomad prophage calls. Reserved for RELAX, which is "
                         "not implemented; supplying this changes nothing.")
    ap.add_argument("--config", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    cfg = load_config(a.config)
    scfg = cfg["specificity"]["selection"]
    fmts, dpi = tuple(cfg["plots"]["formats"]), cfg["plots"]["dpi"]
    set_style()
    os.makedirs(a.workdir, exist_ok=True)
    os.makedirs(a.figdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    if not scfg.get("enabled", True):
        pd.DataFrame().to_csv(a.out, sep="\t", index=False)
        return

    # A config flag that claims an analysis ran is worse than no flag.
    if scfg.get("run_relax"):
        sys.exit("[selection] specificity.selection.run_relax is true, but RELAX is "
                 "not implemented: no test/reference branch partition is ever built "
                 "and `hyphy relax` is never invoked. Set run_relax: false, or "
                 "implement it. Do not leave a flag that says an analysis happened.")
    bad = set(scfg["hyphy_methods"]) - {"FEL", "MEME"}
    if bad:
        sys.exit(f"[selection] hyphy_methods contains {sorted(bad)}, which this "
                 f"script cannot run. Only FEL and MEME are wired: the generic "
                 f"invocation has no --test/--branches and RELAX/aBSREL would "
                 f"either error or prompt interactively.")

    # The FULL protein, from c71.faa. Not an ungapped match-column row: those
    # drop insert residues and anything hmmalign --trim clipped, so they can
    # never equal the translation of a CDS.
    full_prot = {seq_id(h): s.upper().rstrip("*") for h, s in read_fasta(a.faa)}
    keep_ids = {l.strip() for l in open(a.keep) if l.strip()}
    idmap = pd.read_csv(a.idmap, sep="\t", dtype=str)
    prodigal_dir = load_config(a.config)["outputs"]["prodigal_dir"]

    colmap, L = matchcol_residue_index(a.sto, full_prot)
    print(f"[selection] mapped {len(colmap)} sequences onto {L} match columns",
          file=sys.stderr)
    if not colmap:
        sys.exit("[selection] no sequence could be mapped from match columns back "
                 "to residue indices. hits.sto and c71.faa disagree.")

    cds, stats = recover_cds(idmap, keep_ids, prodigal_dir, full_prot)
    print(f"[selection] CDS recovery: {stats}", file=sys.stderr)
    if stats["ok"] < 30:
        n_try = stats["ok"] + stats["mismatch"] + stats["internal_stop"]
        sys.exit(f"[selection] only {stats['ok']} of {n_try} CDS recovered by exact "
                 f"translated match. Not enough for dN/dS.\n"
                 f"  mismatch={stats['mismatch']} means the .ffn translation differs "
                 f"from the protein in c71.faa: different annotation runs, or a "
                 f"genetic-code mismatch (archaea on table 11).\n"
                 f"  Refusing to emit a dN/dS table from {stats['ok']} sequences. "
                 f"Previously this returned a stats file and looked like a result.")

    # dereplicate identical CDS: they carry no substitutions
    seen, uniq = {}, {}
    for sid, nt in cds.items():
        h = hashlib.blake2b(nt.encode(), digest_size=16).hexdigest()
        if h not in seen:
            seen[h] = sid
            uniq[sid] = nt
    print(f"[selection] {len(cds)} CDS -> {len(uniq)} unique", file=sys.stderr)

    if len(uniq) > scfg["max_seqs"]:
        rng = np.random.default_rng(12345)
        pick = set(rng.choice(sorted(uniq), scfg["max_seqs"], replace=False))
        uniq = {k: v for k, v in uniq.items() if k in pick}
        print(f"[selection] subsampled to {len(uniq)} tips", file=sys.stderr)

    ca = codon_align({k: colmap[k] for k in uniq if k in colmap}, uniq, L)
    if not ca:
        sys.exit("[selection] the codon alignment is empty.")
    print(f"[selection] codon alignment: {len(ca)} x "
          f"{len(next(iter(ca.values()))) // 3} codons", file=sys.stderr)

    aln = os.path.join(a.workdir, "codon.fasta")
    with open(aln, "w") as fh:
        for k, v in ca.items():
            write_fasta(fh, k, v)

    # The protein alignment the tree is built from must be the exact translation
    # of the codon alignment, column for column. Reading it from a separate file
    # would let the two drift apart, and HyPhy would then map sites onto a tree
    # inferred from different data.
    prot = os.path.join(a.workdir, "prot.fasta")
    with open(prot, "w") as fh:
        for k, v in ca.items():
            aa = "".join("-" if v[i:i + 3] == "---" else translate(v[i:i + 3])
                         for i in range(0, len(v), 3))
            write_fasta(fh, k, aa)
    tree = os.path.join(a.workdir, "codon")
    sh(["iqtree2", "-s", prot, "--prefix", tree, "--seqtype", "AA",
        "-m", "LG+F+G4", "--fast", "-T", str(a.threads), "-redo"])
    treefile = tree + ".treefile"

    # --- HyPhy ---------------------------------------------------------------
    results = {}
    for method in scfg["hyphy_methods"]:
        outj = os.path.join(a.workdir, f"{method}.json")
        cmd = ["hyphy", method.lower(), "--alignment", aln, "--tree", treefile,
               "--output", outj, "--code", "Universal"]
        try:
            sh(cmd)
        except subprocess.CalledProcessError:
            print(f"[selection] {method} failed; skipping", file=sys.stderr)
            continue
        results[method] = json.load(open(outj))

    groove = pd.read_csv(a.groove, sep="\t").set_index("match_col")
    with open(a.chosen) as fh:
        tri = json.load(fh)

    def _check_site_count(d, method):
        # HyPhy emits one MLE row per codon/match column, in order. If that count
        # ever differs from L (the match-column count the CDS was aligned to), the
        # positional match_col = i mapping below is silently off and every dN/dS is
        # attributed to the wrong column -- the coordinate-drift class this pipeline
        # guards against everywhere else. Fail instead of misattributing.
        if len(d) != L:
            sys.exit(f"[selection] {method} returned {len(d)} sites but the "
                     f"alignment has {L} match columns; the site->column mapping "
                     f"would be wrong. Aborting rather than misplacing selection.")

    rows = []
    if "FEL" in results:
        d = results["FEL"]["MLE"]["content"]["0"]
        _check_site_count(d, "FEL")
        hdr = [h[0] for h in results["FEL"]["MLE"]["headers"]]
        for i, r in enumerate(d):
            rec = dict(zip(hdr, r))
            rows.append({"match_col": i, "method": "FEL",
                         "alpha": rec.get("alpha"), "beta": rec.get("beta"),
                         "dnds": (rec.get("beta") / rec["alpha"]) if rec.get("alpha") else np.nan,
                         "p": rec.get("p-value")})
    if "MEME" in results:
        d = results["MEME"]["MLE"]["content"]["0"]
        _check_site_count(d, "MEME")
        hdr = [h[0] for h in results["MEME"]["MLE"]["headers"]]
        for i, r in enumerate(d):
            rec = dict(zip(hdr, r))
            rows.append({"match_col": i, "method": "MEME",
                         "alpha": rec.get("alpha"), "beta": rec.get("beta+"),
                         "dnds": np.nan, "p": rec.get("p-value")})

    df = pd.DataFrame(rows)
    if df.empty:
        pd.Series(stats).rename("value").rename_axis("metric").to_csv(a.out, sep="\t")
        return
    df["in_groove"] = df["match_col"].map(groove["in_groove"]).fillna(0).astype(int)
    df["is_triad"] = df["match_col"].isin(tri["match_columns"]).astype(int)
    from statsmodels.stats.multitest import multipletests
    df["q_bh"] = np.nan
    for m in df["method"].unique():
        sel = df["method"] == m
        df.loc[sel, "q_bh"] = multipletests(df.loc[sel, "p"].fillna(1),
                                            method="fdr_bh")[1]
    df.to_csv(a.out, sep="\t", index=False)

    # --- groove vs core, paired within one model fit -------------------------
    from scipy.stats import fisher_exact, mannwhitneyu
    summary = dict(stats)
    fel = df[(df["method"] == "FEL") & (df["is_triad"] == 0)]
    if len(fel) > 20 and fel["in_groove"].nunique() == 2:
        g = fel.loc[fel["in_groove"] == 1, "dnds"].dropna()
        c = fel.loc[fel["in_groove"] == 0, "dnds"].dropna()
        if len(g) > 2 and len(c) > 2:
            u, p = mannwhitneyu(g, c, alternative="greater")
            summary.update({"fel_median_dnds_groove": float(g.median()),
                            "fel_median_dnds_core": float(c.median()),
                            "fel_mannwhitney_p": float(p)})
    meme = df[(df["method"] == "MEME") & (df["is_triad"] == 0)]
    if len(meme) > 20:
        sig = meme["q_bh"] < 0.05
        tab = pd.crosstab(sig, meme["in_groove"])
        if tab.shape == (2, 2):
            orr, p = fisher_exact(tab.to_numpy(), alternative="greater")
            summary.update({"meme_sites": int(sig.sum()),
                            "meme_sites_in_groove": int((sig & (meme["in_groove"] == 1)).sum()),
                            "meme_groove_enrichment_or": float(orr),
                            "meme_groove_enrichment_p": float(p)})
    tri_fel = df[(df["method"] == "FEL") & (df["is_triad"] == 1)]
    if len(tri_fel):
        summary["fel_median_dnds_triad"] = float(tri_fel["dnds"].median())

    pd.Series(summary).rename("value").rename_axis("metric").to_csv(
        a.out.replace(".tsv", "_summary.tsv"), sep="\t")
    for k, v in summary.items():
        print(f"  {k}: {v}", file=sys.stderr)

    # --- figure --------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.4))
    ax = axes[0]
    if len(fel):
        data = [fel.loc[fel["in_groove"] == 0, "dnds"].dropna(),
                fel.loc[fel["in_groove"] == 1, "dnds"].dropna(),
                tri_fel["dnds"].dropna()]
        bp = ax.boxplot(data, labels=["core", "groove", "triad"], showfliers=False,
                        patch_artist=True, medianprops=dict(color="k", lw=0.8))
        for b, c in zip(bp["boxes"], [PALETTE[0], PALETTE[4], "0.5"]):
            b.set_facecolor(c); b.set_alpha(0.75); b.set_linewidth(0.5)
    ax.axhline(1.0, color="0.4", ls="--", lw=0.6)
    ax.set_ylabel("FEL $dN/dS$")
    ax.set_yscale("log")
    ax.set_title("A  Site-wise selection", loc="left")

    ax = axes[1]
    if len(meme):
        ax.scatter(meme["match_col"], -np.log10(meme["q_bh"].clip(lower=1e-300)),
                   s=6, c=["#D55E00" if x else "#0072B2" for x in meme["in_groove"]],
                   linewidths=0)
    ax.axhline(-np.log10(0.05), color="0.4", ls="--", lw=0.6)
    ax.set_xlabel("Match-state column")
    ax.set_ylabel("MEME $-\\log_{10}$ q")
    ax.set_title("B  Episodic diversifying selection", loc="left")
    fig.tight_layout()
    savefig(fig, a.figdir, "26_selection_groove_vs_core", fmts, dpi)


if __name__ == "__main__":
    main()
