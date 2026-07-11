#!/usr/bin/env python3
"""Check every input the pipeline will touch, before submitting 700 array jobs.

Reads nothing it does not have to and writes nothing except a report. Run it on
the cluster, where the data actually lives:

    python scripts/preflight.py --config config.yaml
    python scripts/preflight.py --config config.yaml --sample-faa 2000 --strict

Exit status is 0 if there are no ERRORs, 1 otherwise. WARNINGs never fail unless
--strict. Everything is written to preflight_report.tsv as well as stdout.

What gets checked
  config          parses, required keys present
  sample table    columns, duplicate sample IDs, whitespace in sample IDs
  faa paths       a random sample of them exist, are non-empty, are readable,
                  and actually look like protein FASTA
  bacteria QC     join key resolves, completeness/contamination/n50/contigs
                  present and numeric, coverage of the sample set
  archaea .fna    glob matches something, files are DNA, stems resolve to
                  accessions in ar53_metadata
  archaea QC      same column checks
  trees           parse, tip count, tip-label overlap with the sample IDs, with
                  20 examples of what failed to match
  HMMs            exist, HMMER3 format, NAME/LENG present, and -- the one that
                  wastes the most time when wrong -- a GA line if the config
                  asks for --cut_ga
  outputs         parent directories exist or can be created, and are writable
  tools           on PATH, with versions
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import random
import re
import shutil
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_sample_table import (QC_KEYS, TAX_CANDIDATES, autodetect_key,  # noqa: E402
                                gtdb_variants, newick_tip_labels, resolve_tips)
from utils import load_config, parse_gtdb_lineage, read_fasta, seq_id  # noqa: E402
import synteny  # noqa: E402

ROWS: list[dict] = []
COLOR = {"OK": "\033[32m", "WARN": "\033[33m", "ERROR": "\033[31m"}
RESET = "\033[0m"


def say(level, check, detail=""):
    ROWS.append({"level": level, "check": check, "detail": detail})
    c = COLOR.get(level, "") if sys.stdout.isatty() else ""
    r = RESET if c else ""
    print(f"{c}[{level:5s}]{r} {check}" + (f"\n         {detail}" if detail else ""),
          flush=True)


def ok(check, detail=""):
    say("OK", check, detail)


def warn(check, detail=""):
    say("WARN", check, detail)


def err(check, detail=""):
    say("ERROR", check, detail)


def head_lines(path, n=3):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        return [next(fh, "") for _ in range(n)]


# ---------------------------------------------------------------------------
def check_readable(path, what):
    if not os.path.exists(path):
        err(f"{what} missing", path)
        return False
    if not os.access(path, os.R_OK):
        err(f"{what} unreadable", path)
        return False
    if os.path.getsize(path) == 0:
        err(f"{what} empty", path)
        return False
    return True


def check_faa(path):
    """Returns (status, message)."""
    if not os.path.exists(path):
        return "missing", ""
    if os.path.getsize(path) == 0:
        return "empty", ""
    try:
        lines = head_lines(path, 4)
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        return "unreadable", f"{type(e).__name__}: {e}"
    if not lines[0].startswith(">"):
        return "not_fasta", lines[0][:60].rstrip()
    seq = "".join(l.strip() for l in lines[1:] if l and not l.startswith(">"))
    if seq and set(seq.upper()) <= set("ACGTUN"):
        return "looks_like_dna", seq[:40]
    return "ok", ""


def check_sample_table(cfg, n_sample):
    IN = cfg["inputs"]
    p = IN["sample_table"]
    if not check_readable(p, "sample_table"):
        return None

    try:
        hdr = pd.read_csv(p, sep="\t", nrows=0)
    except Exception as e:  # noqa: BLE001
        err("sample_table unparseable as TSV", str(e))
        return None

    for c in (IN["sample_col"], IN["faa_col"]):
        if c not in hdr.columns:
            err(f"sample_table has no '{c}' column",
                f"present: {list(hdr.columns)[:12]}")
            return None
    ok("sample_table columns", f"{len(hdr.columns)} columns, "
                               f"'{IN['sample_col']}' and '{IN['faa_col']}' present")

    df = pd.read_csv(p, sep="\t", usecols=[IN["sample_col"], IN["faa_col"]], dtype=str)
    df.columns = ["sample", "faa"]
    ok("sample_table rows", f"{len(df):,}")

    n_null = int(df.isna().any(axis=1).sum())
    (warn if n_null else ok)("sample_table null sample/faa", f"{n_null} rows")

    dup = df["sample"].duplicated()
    if dup.any():
        err("duplicate sample IDs",
            f"{int(dup.sum())} duplicates, e.g. {df.loc[dup, 'sample'].head(5).tolist()}. "
            f"Every left merge downstream would fan out rows.")
    else:
        ok("sample IDs unique")

    ws = df["sample"].dropna().map(lambda s: any(c.isspace() for c in s))
    if ws.any():
        err("whitespace in sample IDs",
            f"{int(ws.sum())} affected, e.g. {df.loc[ws, 'sample'].head(3).tolist()}. "
            f"Provenance rides in the whitespace-delimited FASTA description field.")
    else:
        ok("sample IDs whitespace-free")

    # sample the faa paths rather than stat 350k of them
    sub = df["faa"].dropna()
    take = min(n_sample, len(sub))
    picked = random.Random(0).sample(list(sub), take)
    counts: dict[str, int] = {}
    examples: dict[str, str] = {}
    for f in picked:
        st, msg = check_faa(f)
        counts[st] = counts.get(st, 0) + 1
        examples.setdefault(st, f"{f} {msg}".strip())

    good = counts.get("ok", 0)
    detail = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    if good == take:
        ok(f"faa spot check ({take} sampled)", detail)
    elif good == 0:
        err(f"faa spot check ({take} sampled)", f"{detail}\n         e.g. "
            + "\n         ".join(f"{k}: {v}" for k, v in examples.items() if k != "ok"))
    else:
        lvl = warn if cfg["allow_missing_faa"] else err
        lvl(f"faa spot check ({take} sampled)",
            f"{detail} -- extrapolates to ~{(take - good) / take * len(sub):,.0f} bad "
            f"of {len(sub):,}\n         e.g. "
            + "\n         ".join(f"{k}: {v}" for k, v in examples.items() if k != "ok"))
    return df


def check_gff_join(cfg, n_sample):
    """Do the bacterial GFFs exist and do their feature IDs join to the .faa?

    Synteny is corroboration for the out-of-order PM candidates, so a broken GFF
    warns rather than aborts. But a GFF whose IDs do not match the proteome is
    silently useless, so verify the join on a sample.
    """
    IN = cfg["inputs"]
    gff_col = IN.get("gff_col")
    if not gff_col:
        ok("inputs.gff_col not set", "no bacterial synteny; out-of-order PM "
           "candidates will be synteny_unknown (not an error)")
        return
    p = IN["sample_table"]
    hdr = pd.read_csv(p, sep="\t", nrows=0)
    if gff_col not in hdr.columns:
        warn(f"inputs.gff_col='{gff_col}' absent from the sample table",
             "bacterial synteny disabled")
        return
    df = pd.read_csv(p, sep="\t", usecols=[IN["faa_col"], gff_col], dtype=str).dropna()
    if df.empty:
        warn("no rows have both a faa and a gff", "bacterial synteny disabled")
        return
    take = min(n_sample, len(df))
    picked = df.sample(take, random_state=0)
    ok_join, bad, unreadable = 0, [], 0
    for _, row in picked.iterrows():
        faa, gff = row[IN["faa_col"]], row[gff_col]
        if not (os.path.exists(faa) and os.path.exists(gff)):
            unreadable += 1
            continue
        try:
            pid = next(seq_id(h) for h, _ in read_fasta(faa))
            coords = synteny.parse_gff(gff)
        except (OSError, StopIteration, Exception) as e:  # noqa: BLE001
            bad.append(f"{gff}: {type(e).__name__}")
            continue
        if pid in coords:
            ok_join += 1
        else:
            bad.append(f"{gff}: .faa id {pid!r} not found among {len(coords)} GFF features")
    if ok_join == take:
        ok(f"GFF join spot check ({take} sampled)", "faa protein IDs found in GFFs")
    elif ok_join == 0:
        warn(f"GFF join spot check ({take} sampled)",
             "no faa ID joined to its GFF -- IDs may differ, so bacterial synteny "
             "will be not_evaluable. e.g. " + (bad[0] if bad else f"{unreadable} unreadable"))
    else:
        warn(f"GFF join spot check ({take} sampled)",
             f"{ok_join}/{take} joined; {unreadable} unreadable. e.g. "
             + (bad[0] if bad else ""))


def check_metadata(path, label, keys, sample_ids=None, key_col=None):
    if not check_readable(path, label):
        return None
    try:
        meta = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, nrows=5)
    except Exception as e:  # noqa: BLE001
        err(f"{label} unparseable as TSV", str(e))
        return None
    if meta.shape[1] < 2:
        err(f"{label} has {meta.shape[1]} column(s)",
            f"delimiter is probably not a tab. First line: {head_lines(path,1)[0][:90]!r}")
        return None

    meta = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    ok(f"{label} rows", f"{len(meta):,} x {meta.shape[1]} columns")

    lower = {c.lower(): c for c in meta.columns}
    found = {}
    for k in QC_KEYS:
        hit = next((lower[c.lower()] for c in keys[k] if c.lower() in lower), None)
        if hit:
            found[k] = hit
            n_num = pd.to_numeric(meta[hit], errors="coerce").notna().sum()
            frac = n_num / len(meta)
            (ok if frac > 0.95 else warn)(
                f"{label}.{k} -> '{hit}'", f"{100*frac:.1f}% parse as numeric")
        else:
            err(f"{label}.{k} not found",
                f"tried {keys[k]}. Add the real column name to `qc:` in config.yaml.\n"
                f"         available: {sorted(meta.columns)[:20]}")

    # Either a lineage string, or per-rank gtdb_* columns. A bare `phylum`
    # column does NOT count: in these tables it is the host animal's phylum and
    # is empty for every microbial row.
    tax = next((lower[c.lower()] for c in TAX_CANDIDATES if c.lower() in lower), None)
    rank = "phylum"
    gtdb_rank = lower.get(f"gtdb_{rank}")
    bare_rank = lower.get(rank)
    def _populated(col):
        # astype(str) would turn NaN into the string "nan" and report 100%
        s = meta[col]
        return s.notna().mean() if s.dtype.kind not in "OU" else \
            s.fillna("").astype(str).str.strip().ne("").mean()

    if tax:
        ok(f"{label} taxonomy", f"lineage column '{tax}'")
    elif gtdb_rank:
        ok(f"{label} taxonomy",
           f"rank column '{gtdb_rank}' ({100*_populated(gtdb_rank):.1f}% populated)")
    else:
        warn(f"{label} taxonomy column",
             f"no lineage string and no gtdb_{rank}; taxon plots will be skipped")
    if bare_rank and gtdb_rank:
        warn(f"{label} has both '{bare_rank}' and '{gtdb_rank}'",
             f"'{bare_rank}' is {100*_populated(bare_rank):.1f}% populated and is a "
             f"HOST-organism field. resolve_taxonomy() prefers gtdb_{rank}; do not set "
             f"plots.taxonomy_col to '{bare_rank}'.")

    if sample_ids is not None:
        if key_col:
            if key_col not in meta.columns:
                err(f"{label} key '{key_col}' absent", f"available: {sorted(meta.columns)[:15]}")
                return meta
            k = key_col
        else:
            try:
                k = autodetect_key(meta, sample_ids, path)
            except SystemExit as e:
                err(f"{label} join key", str(e))
                return meta
        overlap = len(set(map(str, sample_ids)) & set(meta[k].astype(str)))
        frac = overlap / max(len(sample_ids), 1)
        lvl = ok if frac > 0.95 else warn if frac > 0.5 else err
        lvl(f"{label} joins on '{k}'", f"{overlap:,}/{len(sample_ids):,} samples "
                                       f"({100*frac:.1f}%) matched")
    return meta


def check_archaea(cfg):
    IN = cfg["inputs"]
    fna = sorted(glob.glob(IN["archaea_fna_glob"]))
    if not fna:
        err("archaea_fna_glob matched nothing", IN["archaea_fna_glob"])
        return None
    ok("archaea assemblies", f"{len(fna):,} files")

    bad = []
    for f in random.Random(1).sample(fna, min(20, len(fna))):
        try:
            l = head_lines(f, 2)
        except Exception as e:  # noqa: BLE001
            bad.append(f"{os.path.basename(f)}: {e}")
            continue
        if not l[0].startswith(">"):
            bad.append(f"{os.path.basename(f)}: no '>' on line 1")
        elif l[1] and not set(l[1].strip().upper()) <= set("ACGTUNRYKMSWBDHVN-"):
            bad.append(f"{os.path.basename(f)}: line 2 is not DNA")
    (ok if not bad else err)("archaea .fna spot check (20 sampled)", "\n         ".join(bad))

    meta = check_metadata(IN["archaea_metadata"], "archaea_metadata", cfg["qc"],
                          key_col=IN["archaea_metadata_key"])
    if meta is None or IN["archaea_metadata_key"] not in meta.columns:
        return fna

    idx = {}
    for acc in meta[IN["archaea_metadata_key"]].astype(str):
        for v in gtdb_variants(acc):
            idx.setdefault(v, acc)

    stems, unmatched = [], []
    for f in fna:
        s = os.path.basename(f)
        for ext in (".fna.gz", ".fna", ".fa", ".fasta"):
            if s.endswith(ext):
                s = s[: -len(ext)]
                break
        hit = next((idx[v] for v in gtdb_variants(s) if v in idx), None)
        stems.append(hit or s)
        if hit is None:
            unmatched.append(s)
    frac = 1 - len(unmatched) / len(fna)
    lvl = ok if frac > 0.95 else warn if frac > 0.5 else err
    lvl("archaea stems -> ar53_metadata accessions",
        f"{len(fna)-len(unmatched):,}/{len(fna):,} ({100*frac:.1f}%)"
        + (f"\n         unmatched e.g. {unmatched[:5]}" if unmatched else ""))
    return stems


def check_tree(path, label, ids, species, mode, min_frac, reps, meta_accessions=None):
    if not check_readable(path, f"{label} tree"):
        return
    try:
        tips = newick_tip_labels(path)
    except Exception as e:  # noqa: BLE001
        err(f"{label} tree unparseable", str(e))
        return
    if not tips:
        err(f"{label} tree has no tip labels", path)
        return

    looks_species = sum(t.startswith("s__") for t in tips) > 0.5 * len(tips)
    ok(f"{label} tree", f"{len(tips):,} tips, "
                        f"{'species names' if looks_species else 'genome/accession IDs'}"
                        f"; e.g. {sorted(tips)[:2]}")
    if looks_species and mode not in ("gtdb_species",):
        err(f"{label} tree tips are species names but tip_matching.{label}='{mode}'",
            "Set it to gtdb_species, or nothing will match.")
    if not looks_species and mode == "gtdb_species":
        err(f"{label} tip_matching='gtdb_species' but the tips are not species names",
            f"e.g. {sorted(tips)[:3]}")

    # GTDB release skew: tips that no longer exist in the metadata
    if meta_accessions is not None:
        gone = [t for t in tips if not (set(gtdb_variants(t)) & meta_accessions)]
        if gone:
            warn(f"{label} tree / metadata release skew",
                 f"{len(gone):,}/{len(tips):,} tips ({100*len(gone)/len(tips):.1f}%) "
                 f"are absent from the metadata; they will be unusable. "
                 f"e.g. {sorted(gone)[:3]}")
        else:
            ok(f"{label} tree tips all present in metadata")

    if not ids:
        warn(f"{label} tip matching skipped", "no genomes assigned to this domain")
        return

    matched_tips = resolve_tips(ids, species, tips, mode, reps)
    matched = sum(t is not None for t in matched_tips)
    used = len({t for t in matched_tips if t})
    unmatched = [x for x, t in zip(ids, matched_tips) if t is None][:20]
    frac = matched / len(ids)
    lvl = ok if frac >= min_frac else err
    lvl(f"{label} tip matching (mode={mode})",
        f"{matched:,}/{len(ids):,} genomes ({100*frac:.1f}%) -> {used:,} distinct tips "
        f"({matched/max(used,1):.1f} genomes per tip); threshold {100*min_frac:.0f}%"
        + (f"\n         unmatched e.g. {unmatched[:5]}" if unmatched else ""))
    if frac < min_frac:
        err(f"{label} tip matching below tip_matching.min_match_fraction",
            "The pipeline will abort here. Fix the key, not the threshold.")
    if used and matched / used > 1.5:
        ok(f"{label} tree is species-level",
           f"{matched/used:.1f} genomes per tip -- phyloglm.R aggregates to the tip")


def check_hmms(cfg):
    for label, d in cfg["profiles"].items():
        if not d.get("enabled", True):
            warn(f"HMM {label} declared but disabled", "not searched; PeiR-type C39 proteins will not be screened")
            continue
        p = d["path"]
        if not check_readable(p, f"HMM {label}"):
            continue
        hdr = {"NAME": None, "LENG": None, "GA": None, "TC": None, "NC": None, "fmt": None}
        with open(p, errors="replace") as fh:
            for line in fh:
                if line.startswith("HMMER3"):
                    hdr["fmt"] = line.strip()
                elif line.startswith("HMMER2"):
                    hdr["fmt"] = line.strip()
                if line.startswith("HMM "):
                    break
                k = line.split(None, 1)[0] if line.strip() else ""
                if k in hdr and hdr[k] is None and k not in ("fmt",):
                    hdr[k] = line[len(k):].strip().rstrip(";")

        if hdr["fmt"] is None:
            err(f"HMM {label} has no HMMER version line", p)
        elif hdr["fmt"].startswith("HMMER2"):
            warn(f"HMM {label} is HMMER2 format",
                 "prepare_hmms.py will run hmmconvert; make sure hmmconvert is on PATH")
        else:
            ok(f"HMM {label} format", hdr["fmt"])

        if not hdr["NAME"]:
            err(f"HMM {label} has no NAME line", p)
        else:
            ok(f"HMM {label}", f"NAME={hdr['NAME']} LENG={hdr['LENG']} "
                               f"GA={hdr['GA'] or '-'} threshold={d['threshold']}")

        thr = str(d["threshold"])
        if thr.startswith("cut_"):
            line = {"cut_ga": "GA", "cut_nc": "NC", "cut_tc": "TC"}[thr]
            if not hdr[line]:
                err(f"HMM {label}: config asks for --{thr} but there is no {line} line",
                    f"{p}\n         hmmsearch would abort on every batch. Either use the "
                    f"Pfam-A.hmm copy of this model (it carries GA/TC/NC; the InterPro "
                    f"single-model download often does not), or set "
                    f"profiles.{label}.threshold to a bit score.")
            else:
                ok(f"HMM {label} --{thr} available", f"{line} = {hdr[line]}")
        else:
            try:
                float(thr)
                ok(f"HMM {label} bit-score threshold", thr)
            except ValueError:
                err(f"HMM {label} threshold '{thr}' is neither cut_ga/nc/tc nor a number")

    names = []
    for label, d in cfg["profiles"].items():
        if not d.get("enabled", True):
            continue
        if os.path.exists(d["path"]):
            for line in open(d["path"], errors="replace"):
                if line.startswith("NAME"):
                    names.append(line[5:].strip())
                    break
    if len(names) != len(set(names)):
        err("HMM NAME collision", f"{names}: domtblout query names would be ambiguous")

    if cfg["align_profile"] not in cfg["profiles"]:
        err("align_profile not in profiles", cfg["align_profile"])
    if cfg["specific_profile"] not in cfg["profiles"]:
        err("specific_profile not in profiles", cfg["specific_profile"])
    if cfg["profiles"].get(cfg["specific_profile"], {}).get("role") != "specific":
        warn("specific_profile role is not 'specific'",
             "the evidence tiers will still work, but the label is misleading")


def check_outputs(cfg, min_free_gib=100):
    seen = set()
    for k, p in cfg["outputs"].items():
        d = p if k.endswith("_dir") or k in ("outdir", "workdir") else os.path.dirname(p)
        try:
            os.makedirs(d, exist_ok=True)
        except OSError as e:
            err(f"outputs.{k} not creatable", f"{d}: {e}")
            continue
        if not os.access(d, os.W_OK):
            err(f"outputs.{k} not writable", d)
            continue
        ok(f"outputs.{k} writable", d)
        mount = os.stat(d).st_dev
        if mount in seen:
            continue
        seen.add(mount)
        free = shutil.disk_usage(d).free / 2**30
        # 350k proteomes -> the per-batch faa are temp(), but hmm_output and the
        # all-domains table are not. Budget generously.
        (ok if free > min_free_gib else warn)(
            f"free space on {d}", f"{free:.0f} GiB (want > {min_free_gib})")


def check_specificity(cfg):
    """The specificity block needs four external things this cluster cannot
    download. Each one is checked, not assumed."""
    s = cfg.get("specificity")
    if not s or not s.get("enabled"):
        warn("specificity block disabled", "no target-specificity analyses will run")
        return

    # --- Pfam ---------------------------------------------------------------
    pfam = s["pfam_hmm"]
    if check_readable(pfam, "Pfam-A.hmm"):
        missing = [e for e in (".h3f", ".h3i", ".h3m", ".h3p")
                   if not os.path.exists(pfam + e)]
        if missing:
            err("Pfam-A.hmm is not pressed", f"missing {missing}. Run: hmmpress {pfam}")
        else:
            ok("Pfam-A.hmm is pressed")
        wanted = {str(s["pmbr_accession"]), str(s["catalytic_accession"])} | \
                 {str(x) for x in s["accessory_binding_domains"]}
        found = set()
        with open(pfam, errors="replace") as fh:
            for line in fh:
                if line.startswith("ACC "):
                    a = line[4:].strip().split(".")[0]
                    if a in {w.split(".")[0] for w in wanted}:
                        found.add(a)
        miss = {w.split(".")[0] for w in wanted} - found
        (ok if not miss else err)(
            "Pfam accessions present",
            f"{len(found)}/{len(wanted)} found"
            + (f"; MISSING {sorted(miss)} -- wrong Pfam release?" if miss else ""))

    # --- structure ----------------------------------------------------------
    st = s["structure"]
    num = s.get("structure_numbering", "reference")
    if check_readable(st, f"{num} structure ({os.path.basename(st)})"):
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from groove_map import parse_structure
            # parse_structure returns (residues, ions); preflight only needs the
            # residues. Unpacking is not optional -- `x not in (dict, list)` is
            # always True, which would silently report every seed as missing.
            res, _ions = parse_structure(st, s["structure_chain"])
            seeds = [int(x) + int(s.get("structure_offset", 0))
                     for x in s["groove_seed_residues"]]
            miss = [x for x in seeds if x not in res]
            if miss:
                err(f"structure chain {s['structure_chain']}: seed residues absent",
                    f"{miss} not in {min(res)}..{max(res)}. Check "
                    f"specificity.structure_offset and that the numbering is auth.")
            else:
                from groove_map import THREE_TO_ONE
                names_ = "".join(THREE_TO_ONE[res[x]["name"]] for x in seeds)
                ok(f"structure chain {s['structure_chain']}",
                   f"{len(res)} residues; seeds resolve to {names_} at {seeds}")
        except SystemExit as e:
            err("structure unusable", str(e))
        except Exception as e:  # noqa: BLE001
            err("structure could not be parsed", f"{type(e).__name__}: {e}")

    # --- Pmur HMMs ----------------------------------------------------------
    d = s["pmur_hmm_dir"]
    hmms = sorted(glob.glob(os.path.join(d, "*.hmm"))) if os.path.isdir(d) else []
    if not hmms:
        err("no Pmur marker HMMs", f"{d} has no *.hmm; cell-wall genotyping "
                                   f"cannot run and the P1 hypothesis is untestable")
    else:
        n_ga = 0
        for h in hmms:
            with open(h, errors="replace") as fh:
                for line in fh:
                    if line.startswith("HMM "):
                        break
                    if line.startswith("GA "):
                        n_ga += 1
                        break
        ok("Pmur marker HMMs", f"{len(hmms)} models, {n_ga} carry a GA line "
                               f"(the rest use pmur_score_threshold="
                               f"{s['pmur_score_threshold']})")
        if len(hmms) < int(s["pmur_min_markers"]):
            err("too few Pmur markers",
                f"{len(hmms)} models but pmur_min_markers={s['pmur_min_markers']}; "
                f"no genome could ever be called pseudomurein-positive")

    # --- selection ----------------------------------------------------------
    sel = s.get("selection", {})
    if sel.get("enabled"):
        if not shutil.which("hyphy"):
            warn("hyphy not on PATH", "fine if the `selection` conda env will be built")
        ok("selection enabled",
           f"methods {sel['hyphy_methods']}, max {sel['max_seqs']} tips")

    # --- geNomad ------------------------------------------------------------
    if s.get("genomad_enabled"):
        db = s.get("genomad_db")
        if db and os.path.isdir(db):
            ok("geNomad database", db)
        else:
            warn("geNomad database not found", f"{db}; prophage context skipped")

    # --- structure search (divergent-lineage confirmation) ------------------
    ss = s.get("structure_search") or {}
    if ss.get("enabled"):
        checks = [("esmfold_weights", ss.get("esmfold_weights"), False),
                  ("foldseek_pm_db", ss.get("foldseek_pm_db"), ".dbtype"),
                  ("hhsuite_pm_db", ss.get("hhsuite_pm_db"), "_hhm.ffdata")]
        for name, path, marker in checks:
            probe = (str(path) + marker) if marker else path
            if path and os.path.exists(probe):
                ok(f"structure_search {name}", str(path))
            else:
                warn(f"structure_search {name} not staged", f"{path}; that stage "
                     f"will be skipped (recorded, not fabricated)")
        for tool in ("foldseek", "hhsearch"):
            if not shutil.which(tool):
                warn(f"{tool} not on PATH",
                     "structure_search needs the optional `structure` env")

    # --- wall-chemistry reference -------------------------------------------
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from cellwall_reference import (CITATION, DISPUTED, GENUS_CHEMISTRIES,
                                        genus_is_homogeneous)
        het = [g for g in GENUS_CHEMISTRIES if not genus_is_homogeneous(g)]
        ok("wall-chemistry reference", f"{CITATION}; "
           f"{len(GENUS_CHEMISTRIES)} characterised genera, "
           f"{len(het)} heterogeneous ({het})")
        if DISPUTED:
            warn("wall-chemistry claims this pipeline will not act on",
                 f"{sorted(DISPUTED)}. They are called 'unsupported', not guessed. "
                 f"Add a primary citation to cellwall_reference.REFERENCE to use them.")
    except Exception as e:  # noqa: BLE001
        err("cellwall_reference unusable", f"{type(e).__name__}: {e}")

    # --- the circularity guard ----------------------------------------------
    parts = s["sdp_partitions"]
    if "active_site_subgroup" in parts:
        err("SDP partition list contains the active-site subgroups",
            "those are k-means on the very columns an SDP test examines. "
            "The result would be circular by construction.")
    if int(s["sdp_min_partitions"]) < 2:
        warn("sdp_min_partitions < 2",
             "single-partition SDPs are partition artefacts; 2 is the minimum "
             "that means anything")
    ok("SDP partitions", f"{parts}, replication >= {s['sdp_min_partitions']}")


def check_tools():
    need = {"prodigal": ["prodigal", "-v"], "hmmsearch": ["hmmsearch", "-h"],
            "hmmalign": ["hmmalign", "-h"], "hmmconvert": ["hmmconvert", "-h"],
            "trimal": ["trimal", "--version"], "iqtree2": ["iqtree2", "--version"],
            "mmseqs": ["mmseqs", "version"], "diamond": ["diamond", "version"],
            "seqkit": ["seqkit", "version"], "Rscript": ["Rscript", "--version"]}
    missing = []
    for name, cmd in need.items():
        if not shutil.which(name):
            missing.append(name)
            continue
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
            lines = [l for l in (r.stdout + r.stderr).strip().splitlines() if l.strip()]
            v = lines[0][:60] if lines else "present"
            if "Traceback" in v or "cannot stat" in v:
                v = "present (version string unavailable)"
        except Exception:  # noqa: BLE001
            v = "present"
        ok(f"tool {name}", v)
    if missing:
        warn("tools not on PATH", f"{missing}\n         Fine if snakemake will create the "
                                  f"conda envs (--software-deployment-method conda); "
                                  f"fatal otherwise.")

    if shutil.which("Rscript"):
        r = subprocess.run(["Rscript", "-e",
                            'cat(paste(sapply(c("ape","phylolm","caper","data.table","yaml"),'
                            'function(p) paste0(p,"=",requireNamespace(p,quietly=TRUE))),'
                            'collapse=" "))'],
                           capture_output=True, text=True)
        out = r.stdout.strip()
        (ok if "FALSE" not in out else warn)("R packages for phyloglm.R", out or r.stderr[:80])


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--sample-faa", type=int, default=500,
                    help="how many faa paths to spot-check (0 = all)")
    ap.add_argument("--strict", action="store_true", help="treat WARN as failure")
    ap.add_argument("--report", default="preflight_report.tsv")
    a = ap.parse_args()

    if not os.path.exists(a.config):
        sys.exit(f"no such config: {a.config}")
    cfg = load_config(a.config)
    print(f"preflight: {os.path.abspath(a.config)}\n" + "=" * 72)

    print("\n--- config ---")
    for k in ("inputs", "outputs", "profiles", "trees", "qc", "tip_matching",
              "triad", "phyloglm"):
        if k not in cfg:
            err(f"config missing top-level key '{k}'")
    if not [r for r in ROWS if r["level"] == "ERROR"]:
        ok("config structure")

    print("\n--- sample table and proteomes ---")
    n = a.sample_faa or 10**9
    bac = check_sample_table(cfg, n)
    check_gff_join(cfg, min(n, 20))

    print("\n--- bacteria metadata ---")
    bmeta = check_metadata(cfg["inputs"]["bacteria_metadata"], "bacteria_metadata", cfg["qc"],
                           sample_ids=(bac["sample"].dropna().tolist() if bac is not None else None),
                           key_col=cfg["inputs"].get("bacteria_metadata_key"))

    print("\n--- archaea assemblies and metadata ---")
    arc_ids = check_archaea(cfg)

    print("\n--- domain, species and tree tips ---")
    reps, ameta, meta_acc = {}, None, None
    try:
        from build_sample_table import species_rep_map
        reps, ameta = species_rep_map(cfg)
        ok("species representatives", f"{len(reps):,} from archaea_metadata")
        meta_acc = set()
        for x in ameta[cfg["inputs"]["archaea_metadata_key"]].astype(str):
            meta_acc |= gtdb_variants(x)
    except SystemExit as e:
        err("species representative map", str(e))

    bac_ids = bac_sp = []
    arc_sp = []
    if bac is not None and bmeta is not None:
        key = cfg["inputs"].get("bacteria_metadata_key") or autodetect_key(
            bmeta, bac["sample"], cfg["inputs"]["bacteria_metadata"])
        lower = {c.lower(): c for c in bmeta.columns}
        dom_c, sp_c = lower.get("gtdb_domain"), lower.get("gtdb_species")
        if not dom_c:
            warn("bacteria_metadata has no gtdb_domain",
                 "every provided proteome will be treated as Bacteria")
        if not sp_c:
            err("bacteria_metadata has no gtdb_species",
                "species-level tip matching is impossible without it")
        keep = [c for c in (key, dom_c, sp_c) if c]
        m = bmeta[keep].drop_duplicates(subset=[key]).set_index(key)
        j = bac.set_index("sample").join(m)
        dom = (j[dom_c].astype(str).str.replace("^d__", "", regex=True).str.capitalize()
               if dom_c else pd.Series("Bacteria", index=j.index))
        counts = dom.value_counts().to_dict()
        ok("provided proteomes by gtdb_domain", str(counts))
        if counts.get("Archaea", 0):
            warn(f"{counts['Archaea']:,} archaea inside the 'bacterial' sample table",
                 "they are routed to the archaeal tree, not the bacterial one")
        sel_b = (dom == "Bacteria")
        bac_ids = j.index[sel_b].tolist()
        bac_sp = (j[sp_c][sel_b].tolist() if sp_c else [None] * len(bac_ids))
        stray = j.index[~sel_b].tolist()
        stray_sp = (j[sp_c][~sel_b].tolist() if sp_c else [None] * len(stray))
    else:
        stray, stray_sp = [], []

    if ameta is not None and arc_ids:
        akey = cfg["inputs"]["archaea_metadata_key"]
        lin = next((c for c in TAX_CANDIDATES if c in
                    {x.lower(): x for x in ameta.columns}), None)
        sp_map = dict(zip(ameta[akey].astype(str),
                          parse_gtdb_lineage(ameta[lin], "species"))) if lin else {}
        arc_sp = [sp_map.get(x) for x in arc_ids]

    print("\n--- trees ---")
    tm = cfg["tip_matching"]
    check_tree(cfg["trees"]["bacteria"], "bacteria", bac_ids, bac_sp,
               tm["bacteria"], float(tm["min_match_fraction"]), reps)
    ar_ids = list(arc_ids or []) + list(stray)
    ar_sp = list(arc_sp or []) + list(stray_sp)
    check_tree(cfg["trees"]["archaea"], "archaea", ar_ids, ar_sp,
               tm["archaea"], float(tm["min_match_fraction"]), reps, meta_acc)

    print("\n--- HMM profiles ---")
    check_hmms(cfg)

    print("\n--- output directories ---")
    check_outputs(cfg)

    print("\n--- specificity block ---")
    check_specificity(cfg)

    print("\n--- tools ---")
    check_tools()

    df = pd.DataFrame(ROWS)
    df.to_csv(a.report, sep="\t", index=False)
    n_err = int((df["level"] == "ERROR").sum())
    n_warn = int((df["level"] == "WARN").sum())
    print("\n" + "=" * 72)
    print(f"{int((df['level']=='OK').sum())} ok, {n_warn} warnings, {n_err} errors "
          f"-> {a.report}")
    if n_err:
        print("\nERRORS:")
        for _, r in df[df["level"] == "ERROR"].iterrows():
            print(f"  - {r['check']}")
    sys.exit(1 if n_err or (a.strict and n_warn) else 0)


if __name__ == "__main__":
    main()
