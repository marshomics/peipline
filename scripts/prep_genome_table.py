#!/usr/bin/env python3
"""One row per searched genome: outcome, detection-opportunity covariates, tree tip.

This is where the sampling- and detection-bias correction is defined once, so
that every downstream analysis uses the same one.

Outcomes
  has_hit       any profile hit
  has_specific  a PF12386 hit at the curated gathering threshold
  has_c71       a triad-positive sequence  <- the outcome of record
  n_c71         count, for the offset-based prevalence model
  sg_<k>        among C71-positive genomes, does it carry subgroup k

Detection-opportunity covariates
  completeness    an incomplete MAG loses genes roughly in proportion to what is
                  missing, so a 60%-complete genome has ~60% of the chance of
                  showing a C71 it actually has
  contamination   inflates the apparent gene count and can import foreign genes
  log10_n50       fragmented assemblies split genes across contig boundaries;
                  the truncated fragments then fail the coverage filter
  log10_contigs   same phenomenon, other direction
  log10_n_proteins the raw number of chances to find a hit

Every genome that was searched appears here, hit or not. Genomes that were
skipped (missing/empty faa) are excluded, because "not searched" is not "no C71"
and a logistic model cannot tell the difference.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, resolve_taxonomy  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--assign", required=True)
    ap.add_argument("--idmap", required=True)
    ap.add_argument("--map-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    g = pd.read_csv(a.table, sep="\t", dtype={"sample": str}, low_memory=False)

    maps = pd.concat([pd.read_csv(m, sep="\t") for m in
                      sorted(glob.glob(os.path.join(a.map_dir, "batch_*.map.tsv.gz")))],
                     ignore_index=True)
    maps["sample"] = maps["sample"].astype(str)
    g = g.merge(maps[["sample", "status", "n_proteins", "n_residues"]], on="sample",
                how="left", validate="one_to_one")

    n0 = len(g)
    g = g[g["status"] == "ok"].copy()
    print(f"[genome_table] {len(g)}/{n0} genomes were actually searched", file=sys.stderr)

    hits = pd.read_csv(a.combined, sep="\t", dtype={"sample": str, "protein_id": str})
    g["has_hit"] = g["sample"].isin(hits["sample"]).astype(int)
    g["has_specific"] = g["sample"].isin(
        hits.loc[hits["evidence"] == "specific", "sample"]).astype(int)

    idmap = pd.read_csv(a.idmap, sep="\t", dtype=str)
    assign = pd.read_csv(a.assign, sep="\t", dtype={"seq_id": str})
    assign = assign.merge(idmap[["seq_id", "sample"]], on="seq_id", how="left",
                          suffixes=("", "_id"))
    if "sample" not in assign.columns:
        sys.exit("[genome_table] subgroup table has no sample column")

    c71 = assign.groupby("sample").size().rename("n_c71")
    g = g.merge(c71, on="sample", how="left")
    g["n_c71"] = g["n_c71"].fillna(0).astype(int)
    g["has_c71"] = (g["n_c71"] > 0).astype(int)

    for k in sorted(assign["subgroup"].unique()):
        s = set(assign.loc[assign["subgroup"] == k, "sample"])
        g[f"sg_{k}"] = g["sample"].isin(s).astype(int)

    # --- covariates ---------------------------------------------------------
    for c in ("completeness", "contamination", "n50", "contigs", "n_proteins"):
        g[c] = pd.to_numeric(g.get(c), errors="coerce")
    g["log10_n50"] = np.log10(g["n50"].clip(lower=1))
    g["log10_contigs"] = np.log10(g["contigs"].clip(lower=1))
    g["log10_n_proteins"] = np.log10(g["n_proteins"].clip(lower=1))

    tax = resolve_taxonomy(g, cfg["plots"]["taxonomy_col"], cfg["plots"]["taxonomy_rank"])
    g["taxon"] = tax if tax is not None else pd.NA

    cov = cfg["phyloglm"]["covariates"]
    miss = {c: int(g[c].isna().sum()) for c in cov if c in g.columns}
    print("[genome_table] missing covariate values:", miss, file=sys.stderr)
    absent = [c for c in cov if c not in g.columns]
    if absent:
        sys.exit(f"[genome_table] covariates absent from the table: {absent}. "
                 f"Either the QC metadata lacks them or `qc:` in config.yaml names "
                 f"the wrong columns.")

    keep = ["sample", "domain", "source", "prodigal_mode", "tree", "tree_tip",
            "species", "classification", "taxon", "status", "n_proteins", "n_residues",
            "completeness", "contamination", "n50", "contigs",
            "log10_n50", "log10_contigs", "log10_n_proteins",
            "has_hit", "has_specific", "has_c71", "n_c71"] + \
           [c for c in g.columns if c.startswith("sg_")]
    keep = [c for c in keep if c in g.columns]
    g[keep].to_csv(a.out, sep="\t", index=False)

    print(f"[genome_table] prevalence: has_hit {g['has_hit'].mean():.4g}, "
          f"has_specific {g['has_specific'].mean():.4g}, "
          f"has_c71 {g['has_c71'].mean():.4g}", file=sys.stderr)
    for d in g["domain"].dropna().unique():
        m = g["domain"] == d
        tt = int(g.loc[m, "tree_tip"].notna().sum())
        ntip = int(g.loc[m, "tree_tip"].nunique())
        print(f"  {d}: n={int(m.sum()):,}, on tree {tt:,} across {ntip:,} species tips, "
              f"C71+ {int(g.loc[m, 'has_c71'].sum()):,}", file=sys.stderr)


if __name__ == "__main__":
    main()
