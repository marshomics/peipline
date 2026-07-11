#!/usr/bin/env python3
"""Build one sample table across the provided proteomes and the Prodigal-called
archaeal genomes, attach genome-quality covariates and GTDB tree tip labels,
then split it into hmmsearch batches.

Three things here are load-bearing, and all three were wrong in the first draft.

Domain comes from `gtdb_domain`, not from which file a genome arrived in.
`faa_sample_table_90percent.tsv` contains 3,215 archaea among its 342,759
proteomes, and Methanobacteriota is precisely where pseudomurein endopeptidases
live. Hardcoding "Bacteria" would have sent the most interesting genomes to the
wrong tree.

Tree tips are species, not genomes. `gtdbtk.rooted.speciesnames.tree` has ~9.8k
tips labelled `s__Genus species`. Matching a genome ID against them matches
nothing. Every genome is routed to its species' tip, so many genomes share a
tip and the downstream regression must aggregate (phyloglm.R does).

Archaeal reference trees carry species representatives only. A non-representative
genome inherits its species representative's tip, resolved through
`gtdb_representative` in the archaeal metadata.

Everything else fails loudly: duplicate sample IDs, unmatched tips below a
configured fraction, missing QC columns. A quiet NaN here becomes a wrong
p-value four rules later.
"""
from __future__ import annotations

import argparse
import glob
import math
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, parse_gtdb_lineage  # noqa: E402

QC_KEYS = ["completeness", "contamination", "n50", "contigs"]
TAX_CANDIDATES = ["classification", "gtdb_taxonomy", "taxonomy", "lineage",
                  "ncbi_taxonomy", "gtdb_classification"]
DOMAIN_CANDIDATES = ["gtdb_domain"]
SPECIES_CANDIDATES = ["gtdb_species"]


# ---------------------------------------------------------------------------
def newick_tip_labels(path: str) -> set:
    """Linear scan for tip labels. A label immediately following ')' is an
    internal node label (support value or name) and is skipped. Avoids pulling a
    10k-tip tree through a recursive parser just to read its leaves."""
    with open(path) as fh:
        s = fh.read()
    s = re.sub(r"\[[^\]]*\]", "", s)          # newick comments
    tips, buf, after_close, in_quote = set(), [], False, None

    def flush():
        if buf:
            lab = "".join(buf).strip()
            if lab and not after_close:
                tips.add(lab.split(":")[0].strip("'\" "))
        buf.clear()

    for ch in s:
        if in_quote:
            if ch == in_quote:
                in_quote = None
            else:
                buf.append(ch)
            continue
        if ch in "'\"":
            in_quote = ch
            continue
        if ch in "(,);":
            flush()
            after_close = (ch == ")")
            continue
        buf.append(ch)
    flush()
    return {t for t in tips if t}


def gtdb_variants(acc: str):
    """GTDB accessions appear as GB_GCA_000008665.1, RS_GCF_..., or bare
    GCA_000008665.1 depending on which file you are reading."""
    acc = str(acc).strip()
    out = {acc}
    m = re.match(r"^(GB|RS)_(.*)$", acc)
    if m:
        out.add(m.group(2))
    else:
        out.add("GB_" + acc)
        out.add("RS_" + acc)
    out.add(re.sub(r"\.\d+$", "", acc))
    return out


def find_col(df, candidates, what, required=False, path=""):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    msg = f"[sample_table] no '{what}' column in {path}; tried {candidates}"
    if required:
        sys.exit(msg)
    print(msg + " -- leaving it null", file=sys.stderr)
    return None


def autodetect_key(meta, ids, path):
    ids = set(map(str, ids))
    best, best_n = None, 0
    for c in meta.columns:
        if meta[c].dtype.kind not in "OU":
            continue
        n = len(ids & set(meta[c].astype(str)))
        if n > best_n:
            best, best_n = c, n
    if not best or best_n == 0:
        sys.exit(f"[sample_table] no column in {path} matches any sample ID. "
                 f"Set inputs.bacteria_metadata_key explicitly.")
    print(f"[sample_table] {path}: joining on '{best}' "
          f"({best_n}/{len(ids)} samples matched)", file=sys.stderr)
    return best


def attach_meta(df, meta, key_left, key_right, qc_cfg, path, extra=()):
    """Merge QC, taxonomy, domain, species and every gtdb_* column."""
    ren = {}
    missing_qc = []
    for k in QC_KEYS:
        col = find_col(meta, qc_cfg[k], k, path=path)
        if col:
            ren[col] = k
        else:
            missing_qc.append(k)
    if missing_qc:
        # The module docstring promises QC problems fail loudly. A missing QC
        # column becomes an all-NaN covariate that degrades the phyloglm
        # detection-bias model four rules downstream, so name it here where it is
        # visible rather than let it surface as a silent NaN. Not a hard error:
        # a metadata source may legitimately lack one, and the regression drops
        # that covariate; but it must never be silent.
        print(f"[sample_table] WARNING: QC column(s) {missing_qc} not found in "
              f"{path}. The matching covariate(s) will be NaN and any phyloglm "
              f"term on them is uninformative. Add the column or drop the covariate "
              f"from config.phyloglm.", file=sys.stderr)
    for canon, cands in extra:
        col = find_col(meta, cands, canon, path=path)
        if col:
            ren[col] = canon
    tax = find_col(meta, TAX_CANDIDATES, "taxonomy", path=path)
    if tax:
        ren[tax] = "classification"

    gtdb_cols = [c for c in meta.columns if c.lower().startswith("gtdb_")
                 and c not in ren]
    keep = [key_right] + list(ren) + gtdb_cols
    keep = list(dict.fromkeys(keep))
    m = meta[keep].rename(columns=ren).drop_duplicates(subset=[key_right])
    out = df.merge(m, left_on=key_left, right_on=key_right, how="left")
    if key_right != key_left and key_right in out.columns:
        out = out.drop(columns=[key_right])
    return out


def norm_domain(s):
    if not isinstance(s, str) or not s.strip():
        return pd.NA
    s = s.strip()
    s = re.sub(r"^d__", "", s)
    return s.capitalize() if s.lower() in ("bacteria", "archaea") else s


# ---------------------------------------------------------------------------
def species_rep_map(cfg):
    """species (with s__ prefix) -> representative accession, from the archaeal
    metadata's `gtdb_representative` flag."""
    path = cfg["inputs"]["archaea_metadata"]
    meta = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    acc = cfg["inputs"]["archaea_metadata_key"]
    if acc not in meta.columns:
        sys.exit(f"[sample_table] '{acc}' not in {path}")
    lin = find_col(meta, TAX_CANDIDATES, "taxonomy", required=True, path=path)
    rep = find_col(meta, ["gtdb_representative"], "gtdb_representative", path=path)
    meta["_species"] = parse_gtdb_lineage(meta[lin], "species")
    if not rep:
        print("[sample_table] no gtdb_representative column; species-rep tip "
              "matching will be unavailable", file=sys.stderr)
        return {}, meta
    isrep = meta[rep].astype(str).str.lower().isin(("t", "true", "1"))
    m = (meta.loc[isrep & meta["_species"].notna()]
             .drop_duplicates("_species")
             .set_index("_species")[acc].to_dict())
    print(f"[sample_table] {len(m):,} species representatives in {os.path.basename(path)}",
          file=sys.stderr)
    return m, meta


def _valid_species(sp):
    """`s__` with nothing after it means "unclassified", not a species.

    GTDB-Tk trees really do carry a tip labelled `s__`. Matching unclassified
    genomes to it would collapse every one of them onto a single tip and model
    them as one clade with one shared ancestor. That is a fabricated result, not
    a missing-data problem.
    """
    return (isinstance(sp, str) and len(sp.strip()) > 3
            and not re.fullmatch(r"[a-z]__", sp.strip()))


def resolve_tips(ids, species, tips, mode, reps):
    """Return one tip label (or None) per genome."""
    tips = {t for t in tips if not re.fullmatch(r"[a-z]__", t.strip())}
    out = []
    for x, sp in zip(ids, species):
        if not _valid_species(sp):
            sp = None
        hit = None
        if mode == "identity":
            hit = x if isinstance(x, str) and x in tips else None
        elif mode == "gtdb_accession":
            if isinstance(x, str):
                hit = next((v for v in gtdb_variants(x) if v in tips), None)
        elif mode == "gtdb_species":
            hit = sp if isinstance(sp, str) and sp in tips else None
        elif mode == "gtdb_species_rep":
            if isinstance(x, str):
                hit = next((v for v in gtdb_variants(x) if v in tips), None)
            if hit is None and isinstance(sp, str):
                r = reps.get(sp)
                if r:
                    hit = next((v for v in gtdb_variants(r) if v in tips), None)
        else:
            sys.exit(f"[sample_table] unknown tip_matching mode '{mode}'")
        out.append(hit)
    return out


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--manifests", nargs="*", default=[])
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-tips", required=True)
    ap.add_argument("--batchdir", required=True)
    ap.add_argument("--batch-size", type=int, required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    IN, QC = cfg["inputs"], cfg["qc"]
    reps, ameta = species_rep_map(cfg)

    # ---- provided proteomes ------------------------------------------------
    bac = pd.read_csv(IN["sample_table"], sep="\t", dtype=str, low_memory=False)
    for c in (IN["sample_col"], IN["faa_col"]):
        if c not in bac.columns:
            sys.exit(f"[sample_table] '{c}' not in {IN['sample_table']}")
    bac = bac.rename(columns={IN["sample_col"]: "sample", IN["faa_col"]: "faa"})
    bac["source"] = "provided_faa"
    bac["prodigal_mode"] = pd.NA
    # Optional per-genome GFF path for the pseudomurein synteny check. Only read
    # for out-of-order candidates, so a missing value is harmless (synteny is
    # then "not_evaluable" for that genome, never a false negative).
    gff_col = IN.get("gff_col")
    if gff_col and gff_col in bac.columns:
        bac["gff"] = bac[gff_col]
        print(f"[sample_table] bacterial GFF paths from column '{gff_col}'",
              file=sys.stderr)
    else:
        bac["gff"] = pd.NA
        if gff_col:
            print(f"[sample_table] inputs.gff_col='{gff_col}' not found in the "
                  f"sample table; bacterial synteny will be not_evaluable",
                  file=sys.stderr)

    if IN.get("bacteria_metadata"):
        bmeta = pd.read_csv(IN["bacteria_metadata"], sep="\t", dtype=str, low_memory=False)
        key = IN.get("bacteria_metadata_key") or autodetect_key(
            bmeta, bac["sample"], IN["bacteria_metadata"])
        bac = attach_meta(bac, bmeta, "sample", key, QC, IN["bacteria_metadata"],
                          extra=[("domain_raw", DOMAIN_CANDIDATES),
                                 ("species", SPECIES_CANDIDATES)])
    bac["domain"] = bac.get("domain_raw", pd.Series(pd.NA, index=bac.index)).map(norm_domain)
    n_nodom = int(bac["domain"].isna().sum())
    if n_nodom:
        print(f"[sample_table] {n_nodom} provided proteomes have no gtdb_domain; "
              f"assuming Bacteria", file=sys.stderr)
        bac["domain"] = bac["domain"].fillna("Bacteria")
    print("[sample_table] provided proteomes by domain: "
          f"{bac['domain'].value_counts().to_dict()}", file=sys.stderr)

    # ---- Prodigal archaea ---------------------------------------------------
    man = [pd.read_csv(m, sep="\t", dtype=str) for m in a.manifests]
    arc = pd.concat(man, ignore_index=True) if man else pd.DataFrame(
        columns=["stem", "fna", "faa", "prodigal_mode", "n_proteins", "genome_size"])
    n_failed = int((arc["prodigal_mode"] == "failed").sum()) if len(arc) else 0
    arc = arc[arc["prodigal_mode"] != "failed"].copy()
    if n_failed:
        print(f"[sample_table] {n_failed} archaeal genomes failed gene calling",
              file=sys.stderr)

    if len(arc):
        arc = arc.rename(columns={"stem": "sample"})
        arc["domain"] = "Archaea"
        arc["source"] = "prodigal"
        akey = IN["archaea_metadata_key"]
        acc_index = {}
        for acc in ameta[akey].astype(str):
            for v in gtdb_variants(acc):
                acc_index.setdefault(v, acc)
        arc["_join"] = arc["sample"].map(
            lambda s: next((acc_index[v] for v in gtdb_variants(s) if v in acc_index), None))
        n_nomatch = int(arc["_join"].isna().sum())
        if n_nomatch:
            print(f"[sample_table] {n_nomatch}/{len(arc)} archaeal genomes have no row in "
                  f"{IN['archaea_metadata']}; e.g. "
                  f"{arc.loc[arc['_join'].isna(), 'sample'].head(5).tolist()}",
                  file=sys.stderr)
        arc = attach_meta(arc, ameta, "_join", akey, QC, IN["archaea_metadata"])
        arc["species"] = parse_gtdb_lineage(arc.get("classification", pd.Series(dtype=str)),
                                            "species")
        arc["sample"] = arc["_join"].fillna(arc["sample"])
        arc = arc.drop(columns=["_join", "fna"], errors="ignore")

    # Archaea: the GFF Prodigal now writes beside each .faa (same stem).
    if len(arc):
        arc["gff"] = arc["faa"].astype(str).str.replace(r"\.faa$", ".gff", regex=True)

    cols = ["sample", "faa", "gff", "domain", "source", "prodigal_mode",
            "classification", "species", "completeness", "contamination", "n50",
            "contigs", "genome_size"]
    gtdb_extra = sorted({c for d in (bac, arc) for c in d.columns
                         if c.lower().startswith("gtdb_") and c != "gtdb_representative"})
    for d in (bac, arc):
        for c in cols + gtdb_extra:
            if c not in d.columns:
                d[c] = pd.NA
    df = pd.concat([bac[cols + gtdb_extra], arc[cols + gtdb_extra]], ignore_index=True)

    df = df[df["faa"].notna() & (df["faa"].astype(str) != "")]
    dup = df["sample"].duplicated()
    if dup.any():
        sys.exit(f"[sample_table] {int(dup.sum())} duplicate sample IDs across the two "
                 f"sources, e.g. {df.loc[dup, 'sample'].head(10).tolist()}")

    for c in ("completeness", "contamination", "n50", "contigs", "genome_size"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # ---- tree tips ---------------------------------------------------------
    tipcfg = cfg["tip_matching"]
    tip_rows = []
    df["tree"] = pd.NA
    df["tree_tip"] = pd.NA
    for domain, tree_key in (("Bacteria", "bacteria"), ("Archaea", "archaea")):
        m = df["domain"] == domain
        if not m.any():
            continue
        path = cfg["trees"][tree_key]
        if not os.path.exists(path):
            print(f"[sample_table] tree missing: {path}; {domain} rows get no tip",
                  file=sys.stderr)
            continue
        tips = newick_tip_labels(path)
        mode = tipcfg[tree_key]
        matched = resolve_tips(df.loc[m, "sample"], df.loc[m, "species"], tips, mode, reps)
        df.loc[m, "tree"] = tree_key
        df.loc[m, "tree_tip"] = matched

        n_ok = sum(x is not None for x in matched)
        frac = n_ok / max(m.sum(), 1)
        unmatched = [s for s, t in zip(df.loc[m, "sample"], matched) if t is None][:20]
        n_tip_used = len({t for t in matched if t})
        tip_rows.append({"domain": domain, "tree": path, "mode": mode,
                         "n_tips_in_tree": len(tips), "n_samples": int(m.sum()),
                         "n_matched": n_ok, "frac_matched": round(frac, 4),
                         "n_distinct_tips_used": n_tip_used,
                         "genomes_per_tip": round(n_ok / max(n_tip_used, 1), 2),
                         "examples_unmatched": ";".join(unmatched)})
        print(f"[sample_table] {domain}: {n_ok:,}/{int(m.sum()):,} samples "
              f"({100 * frac:.1f}%) matched to {n_tip_used:,} of {len(tips):,} tips "
              f"via '{mode}' ({n_ok / max(n_tip_used, 1):.1f} genomes per tip)",
              file=sys.stderr)
        if frac < float(tipcfg["min_match_fraction"]):
            sys.exit(f"[sample_table] only {100 * frac:.1f}% of {domain} samples matched "
                     f"tips in {path} using mode '{mode}'. That is a key problem, not a "
                     f"biological result. Examples: {unmatched[:5]}")

    pd.DataFrame(tip_rows).to_csv(a.out_tips, sep="\t", index=False)
    os.makedirs(os.path.dirname(os.path.abspath(a.out_table)), exist_ok=True)
    df.to_csv(a.out_table, sep="\t", index=False)

    # ---- batches -----------------------------------------------------------
    os.makedirs(a.batchdir, exist_ok=True)
    for f in glob.glob(os.path.join(a.batchdir, "batch_*.tsv")):
        os.remove(f)
    n_b = max(1, math.ceil(len(df) / a.batch_size))
    for i in range(n_b):
        chunk = df.iloc[i * a.batch_size:(i + 1) * a.batch_size][["sample", "faa"]]
        chunk.to_csv(os.path.join(a.batchdir, f"batch_{i:05d}.tsv"), sep="\t", index=False)

    print(f"[sample_table] {len(df):,} genomes "
          f"({int((df['domain'] == 'Bacteria').sum()):,} bacteria, "
          f"{int((df['domain'] == 'Archaea').sum()):,} archaea) -> {n_b} batches",
          file=sys.stderr)
    for c in QC_KEYS:
        print(f"  {c}: {int(df[c].notna().sum()):,}/{len(df):,} non-null", file=sys.stderr)


if __name__ == "__main__":
    main()
