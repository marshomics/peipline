#!/usr/bin/env python3
"""Combine every per-batch, per-profile domtblout into one filtered table.

Evidence tiers, kept distinct throughout:

  PF12386  (specific)    passed the curated Pfam gathering threshold. This is
                         evidence that the protein is a pseudomurein
                         endo-isopeptidase catalytic domain.
  SSF54001 (sensitivity) passed a bit score of 25 against a SCOP superfamily
                         model of *all* cysteine proteinases. This is evidence
                         that the protein has a papain- or transglutaminase-like
                         fold, and nothing more.

A protein that hits only SSF54001 is retained (you asked for that), but it is
labelled `ssf_only` in the `evidence` column and it never contributes to the
catalytic-triad column call.

The decoy table gives the empirical FDR: reversed sequences have the same length
and composition as the targets, so the rate at which they clear a bit score is
the rate at which chance clears it.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import DOMTBL_DTYPES, load_config, parse_domtblout  # noqa: E402

NUMERIC = list(DOMTBL_DTYPES)

OUT_COLS = [
    "sample", "protein_id", "faa", "profile", "role", "evidence",
    "hmm_name", "hmm_acc", "protein_len", "profile_len",
    "full_evalue", "full_score", "full_bias",
    "n_domains", "dom_score", "dom_bias", "c_evalue", "i_evalue",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "acc", "profile_coverage", "n_profiles_hit", "profiles_hit", "batch",
]


def load_domtbl(path: str, label: str, min_score, batch: str) -> pd.DataFrame:
    rows = []
    for r in parse_domtblout(path):
        try:
            fs, ds = float(r["full_score"]), float(r["dom_score"])
        except ValueError:
            continue
        if min_score is not None and (fs < min_score or ds < min_score):
            continue
        desc = r["description"].split()
        rows.append({
            "protein_id": desc[0] if desc else r["target_name"],
            "sample": desc[1] if len(desc) > 1 else pd.NA,
            "profile": label,
            "hmm_name": r["query_name"],
            "hmm_acc": r["query_accession"],
            "batch": batch,
            **{k: r[k] for k in NUMERIC},
        })
    if not rows:
        return pd.DataFrame(columns=["protein_id", "sample", "profile", "hmm_name",
                                     "hmm_acc", "batch"] + NUMERIC)
    df = pd.DataFrame(rows)
    for c, t in DOMTBL_DTYPES.items():
        df[c] = pd.to_numeric(df[c], errors="coerce").astype(t)
    return df


def scores_only(path: str, min_score) -> np.ndarray:
    out = []
    for r in parse_domtblout(path):
        try:
            fs = float(r["full_score"])
        except ValueError:
            continue
        if min_score is None or fs >= min_score:
            out.append(fs)
    return np.asarray(out, dtype=float)


def decoy_fdr(hmm_dir, profiles, thresholds, out_path):
    rows = []
    for label in profiles:
        tgt = np.concatenate([scores_only(p, None) for p in
                              sorted(glob.glob(os.path.join(hmm_dir, f"batch_*.{label}.domtblout")))]
                             or [np.array([])])
        dec = np.concatenate([scores_only(p, None) for p in
                              sorted(glob.glob(os.path.join(hmm_dir, "decoy",
                                                            f"batch_*.{label}.domtblout")))]
                             or [np.array([])])
        if tgt.size == 0:
            continue
        lo = float(np.floor(min(tgt.min(), dec.min() if dec.size else tgt.min())))
        hi = float(np.ceil(np.percentile(tgt, 99.5)))
        for s in np.linspace(lo, hi, 60):
            nt = int((tgt >= s).sum())
            nd = int((dec >= s).sum())
            rows.append({"profile": label, "bit_score": s, "n_target": nt,
                         "n_decoy": nd,
                         "fdr": (nd / nt) if nt else np.nan,
                         "applied_threshold": str(thresholds[label])})
    df = pd.DataFrame(rows)
    df.to_csv(out_path, sep="\t", index=False)
    for label in df["profile"].unique() if len(df) else []:
        d = df[df["profile"] == label]
        at = d.iloc[0]
        print(f"[fdr] {label}: at the lowest reported score {at['bit_score']:.1f}, "
              f"{at['n_decoy']} decoy vs {at['n_target']} target hits "
              f"(FDR ~ {at['fdr']:.3g})", file=sys.stderr)
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmm-dir", required=True)
    ap.add_argument("--map-dir", required=True)
    ap.add_argument("--profile-map", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-combined", required=True)
    ap.add_argument("--out-all", required=True)
    ap.add_argument("--out-stats", required=True)
    ap.add_argument("--out-fdr", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    # A disabled profile was never searched, so it must not appear in the funnel,
    # the FDR table or the evidence tiers. See config.yaml: PF03412 is declared
    # and disabled until the per-family alignment path exists.
    PROF = {k: v for k, v in cfg["profiles"].items() if v.get("enabled", True)}
    specific = cfg["specific_profile"]
    min_cov = float(cfg.get("min_profile_coverage", 0.0))

    for p in (a.out_combined, a.out_all, a.out_stats, a.out_fdr):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    pmap = pd.read_csv(a.profile_map, sep="\t", dtype=str).set_index("label")
    prot = ["sample", "protein_id"]

    # hmmsearch already applied every threshold. Re-applying the numeric ones is
    # defensive: it makes this script correct if you re-run it over tables
    # produced with a looser cutoff. `cut_ga` cannot be re-applied here (the
    # threshold lives in the model), so those rows pass through as-is.
    thresholds, reapply = {}, {}
    for label, d in PROF.items():
        thresholds[label] = d["threshold"]
        reapply[label] = None if str(d["threshold"]).startswith("cut_") else float(d["threshold"])

    parts = []
    for label in sorted(PROF):
        files = sorted(glob.glob(os.path.join(a.hmm_dir, f"batch_*.{label}.domtblout")))
        if not files:
            sys.exit(f"[combine] no domtblout for profile {label} under {a.hmm_dir}")
        for p in files:
            batch = os.path.basename(p).split(".")[0].replace("batch_", "")
            d = load_domtbl(p, label, reapply[label], batch)
            if len(d):
                parts.append(d)

    all_dom = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=OUT_COLS)
    print(f"[combine] {len(all_dom)} domain rows pass their profile's threshold",
          file=sys.stderr)

    map_files = sorted(glob.glob(os.path.join(a.map_dir, "batch_*.map.tsv.gz")))
    smap = pd.concat([pd.read_csv(m, sep="\t") for m in map_files], ignore_index=True)
    smap["sample"] = smap["sample"].astype(str)

    if len(all_dom):
        # A duplicate sample in the batch maps would multiply hit rows. Assert the
        # 1:1 that build_sample_table.py already guarantees, rather than trust it.
        all_dom = all_dom.merge(smap[["sample", "faa"]], on="sample", how="left",
                                validate="many_to_one")
        if all_dom["faa"].isna().any():
            n = int(all_dom["faa"].isna().sum())
            sys.exit(f"[combine] {n} hits could not be traced to a faa path. The batch "
                     f"faa description field was not '<protein_id> <sample>'.")

        all_dom["role"] = all_dom["profile"].map(lambda p: PROF[p]["role"])
        all_dom["profile_coverage"] = (
            (all_dom["hmm_to"] - all_dom["hmm_from"] + 1) / all_dom["qlen"]).astype("float32")
        all_dom.to_csv(a.out_all, sep="\t", index=False, compression="gzip")

        if min_cov > 0:
            before = len(all_dom)
            all_dom = all_dom[all_dom["profile_coverage"] >= min_cov]
            print(f"[combine] coverage filter >= {min_cov}: {before} -> {len(all_dom)} rows",
                  file=sys.stderr)

        key = prot + ["profile"]
        counts = all_dom.groupby(key, sort=False).size().rename("n_domains").reset_index()
        all_dom = all_dom.sort_values("dom_score", ascending=False, kind="mergesort")
        best = all_dom.drop_duplicates(subset=key, keep="first").merge(counts, on=key, how="left")

        best["protein_len"] = best["tlen"]
        best["profile_len"] = best["qlen"]

        g = best.groupby(prot, sort=False)["profile"]
        best["n_profiles_hit"] = g.transform("nunique").astype("int8")
        best["profiles_hit"] = g.transform(lambda s: ",".join(sorted(set(s))))
        # str.contains is a REGEX SUBSTRING test. It happens to work because
        # "PF12386" is not a substring of "SSF54001" and has no metacharacters.
        # Add a profile named e.g. "PF123" and every PF12386 hit silently becomes
        # "specific". Test membership in the comma-separated list instead.
        best["evidence"] = np.where(
            best["profiles_hit"].str.split(",").map(lambda xs: specific in xs),
            "specific", "ssf_only")

        best = best[OUT_COLS].sort_values(["sample", "protein_id", "profile"])
        best.to_csv(a.out_combined, sep="\t", index=False)
    else:
        pd.DataFrame(columns=OUT_COLS).to_csv(a.out_all, sep="\t", index=False,
                                              compression="gzip")
        pd.DataFrame(columns=OUT_COLS).to_csv(a.out_combined, sep="\t", index=False)
        best = pd.DataFrame(columns=OUT_COLS)

    if cfg.get("decoy_fdr", False):
        decoy_fdr(a.hmm_dir, sorted(PROF), thresholds, a.out_fdr)
    else:
        pd.DataFrame(columns=["profile", "bit_score", "n_target", "n_decoy", "fdr",
                              "applied_threshold"]).to_csv(a.out_fdr, sep="\t", index=False)

    n_ok = int((smap["status"] == "ok").sum())
    stats = {
        "samples_in_table": len(smap),
        "samples_searched": n_ok,
        "samples_skipped": len(smap) - n_ok,
        "proteins_searched": int(smap["n_proteins"].sum()),
        "residues_searched": int(smap["n_residues"].sum()),
        "min_profile_coverage": min_cov,
        "domain_rows_passing": len(all_dom),
        "protein_profile_pairs": len(best),
        "unique_proteins_with_hit": int(best[prot].drop_duplicates().shape[0]) if len(best) else 0,
        "samples_with_hit": int(best["sample"].nunique()) if len(best) else 0,
    }
    for label in sorted(PROF):
        stats[f"threshold_{label}"] = str(thresholds[label])
    if len(best):
        uniq = best.drop_duplicates(prot)
        for lab, n in best["profile"].value_counts().items():
            stats[f"proteins_hit_{lab}"] = int(n)
        stats["proteins_specific_evidence"] = int((uniq["evidence"] == "specific").sum())
        stats["proteins_ssf_only"] = int((uniq["evidence"] == "ssf_only").sum())
        stats["proteins_hit_both_profiles"] = int((uniq["n_profiles_hit"] == 2).sum())

    pd.Series(stats).rename("value").rename_axis("metric").to_csv(a.out_stats, sep="\t")
    for k, v in stats.items():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
