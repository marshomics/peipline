#!/usr/bin/env python3
"""Pull full-length sequences back out of the original faa files.

Two modes:

  --hits <hmm_output_combine.txt>       -> writes hits.faa + hits_idmap.tsv.gz
      Every protein with >=1 passing hit. Sequences are given stable, clean IDs
      (`c71_000001`) because Stockholm and Newick both choke on the punctuation
      that shows up in real protein accessions.

  --keep-ids <ids.txt> --idmap <map>    -> writes c71.faa
      The subset that survived the catalytic-triad filter, re-extracted from the
      original faa files (not sliced out of hits.faa), so the deliverable is
      demonstrably byte-identical to the source proteomes.

Both modes group work by faa file and read each file exactly once, in parallel.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import read_fasta, seq_id, write_fasta  # noqa: E402


def _extract_one(job):
    """job = (faa_path, {protein_id: out_id}) -> [(out_id, seq), ...]"""
    faa, wanted = job
    found = []
    try:
        for header, seq in read_fasta(faa):
            pid = seq_id(header)
            if pid in wanted:
                found.append((wanted[pid], seq))
                if len(found) == len(wanted):
                    break
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        return faa, [], f"{type(e).__name__}: {e}"
    missing = len(wanted) - len(found)
    return faa, found, (f"{missing} protein_id(s) not found" if missing else None)


def run(jobs, out_faa, threads):
    os.makedirs(os.path.dirname(os.path.abspath(out_faa)), exist_ok=True)
    n, problems = 0, []
    with Pool(threads) as pool, open(out_faa, "w") as fo:
        for faa, found, err in pool.imap_unordered(_extract_one, jobs, chunksize=8):
            if err:
                problems.append(f"{faa}\t{err}")
            for oid, seq in found:
                write_fasta(fo, oid, seq)
                n += 1
    if problems:
        p = out_faa + ".problems.txt"
        with open(p, "w") as fh:
            fh.write("\n".join(problems) + "\n")
        print(f"[extract] {len(problems)} faa files had problems -> {p}", file=sys.stderr)
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits")
    ap.add_argument("--keep-ids")
    ap.add_argument("--idmap")
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-idmap")
    ap.add_argument("--out-evidence",
                    help="tier of every sequence written to the output FASTA")
    ap.add_argument("--family", default=None,
                    help="restrict --hits extraction to one family arm. With --config, "
                         "keeps only proteins whose profiles_hit intersects that "
                         "family's specific + sensitivity profiles, relabels evidence "
                         "relative to THAT family's specific profile, and prefixes "
                         "seq_ids with the family name. Omit for the historical "
                         "single-arm behaviour (all hits, c71_ prefix).")
    ap.add_argument("--config",
                    help="config.yaml; required with --family to resolve the arm's "
                         "profile set.")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    if a.hits:
        df = pd.read_csv(a.hits, sep="\t",
                         usecols=["sample", "protein_id", "faa", "profiles_hit", "evidence"],
                         dtype=str).drop_duplicates(subset=["sample", "protein_id"])
        df = df.reset_index(drop=True)

        # --- family scoping ---------------------------------------------------
        # Why this is here and not in combine_filter: a protein that hit ONLY
        # PF03412 is a C39 candidate, not a C71 `ssf_only` sequence. If it reached
        # the C71 arm it would be aligned to the PF12386 scaffold, score as
        # low_coverage against the C71 triad columns, and be silently miscounted.
        # Each arm takes the hits that belong to it: its specific model, plus the
        # shared SSF54001 fold net. A papain-fold protein (SSF54001-only) is a
        # candidate for BOTH arms and is tested against each family's own columns.
        prefix = "c71"
        if a.family:
            prefix = a.family
            if a.config:
                from utils import load_config  # noqa: E402
                fam = (load_config(a.config).get("families") or {}).get(a.family) or {}
                spec = fam.get("specific_profile")
                famset = set(([spec] if spec else []) +
                             list(fam.get("sensitivity_profiles") or []))
                if famset:
                    ph = df["profiles_hit"].fillna("")
                    keep = ph.map(lambda s: bool(set(str(s).split(",")) & famset))
                    df = df[keep.to_numpy()].reset_index(drop=True)
                    # evidence is now relative to THIS family's specific profile
                    df["evidence"] = ph[keep.to_numpy()].reset_index(drop=True).map(
                        lambda s: "specific" if spec in str(s).split(",") else "ssf_only")
                    print(f"[extract] family={a.family}: kept {len(df)} proteins whose "
                          f"profiles_hit intersects {sorted(famset)}", file=sys.stderr)
        df["seq_id"] = [f"{prefix}_{i:07d}" for i in range(len(df))]

        by_faa = defaultdict(dict)
        for faa, pid, sid, smp in zip(df["faa"], df["protein_id"], df["seq_id"], df["sample"]):
            by_faa[faa][pid] = f"{sid} {smp} {pid}"

        n = run(list(by_faa.items()), a.out_faa, a.threads)
        if a.out_idmap:
            df.to_csv(a.out_idmap, sep="\t", index=False, compression="gzip")
        print(f"[extract] wrote {n}/{len(df)} sequences -> {a.out_faa}", file=sys.stderr)
        if n != len(df):
            sys.exit("[extract] some hit sequences could not be recovered; see .problems.txt")

    elif a.keep_ids and a.idmap:
        keep = {l.strip() for l in open(a.keep_ids) if l.strip()}
        idmap = pd.read_csv(a.idmap, sep="\t", dtype=str)
        idmap = idmap[idmap["seq_id"].isin(keep)]
        if len(idmap) != len(keep):
            sys.exit(f"[extract] {len(keep) - len(idmap)} kept IDs absent from idmap")

        # c71.faa is the UNION of two evidence tiers. Nothing downstream of the
        # triad filter distinguishes a sequence that cleared PF12386's curated
        # gathering threshold from one that cleared only the SSF54001 fold model
        # and then happened to carry C/H/D at the right columns. Write the tier
        # out beside the FASTA so the fact is visible, and say the numbers.
        if "evidence" in idmap.columns:
            comp = idmap["evidence"].value_counts()
            n_ssf = int(comp.get("ssf_only", 0))
            frac = n_ssf / max(len(idmap), 1)
            print("[extract] c71.faa composition: "
                  + ", ".join(f"{k}={v:,}" for k, v in comp.items()), file=sys.stderr)
            if n_ssf:
                print(f"[extract] {frac:.1%} of c71.faa cleared only the SSF54001 "
                      f"fold model. Their triad is real; their family is not "
                      f"established, because SCOP 54001 spans ~22 families. NO "
                      f"downstream analysis conditions on this. Read the "
                      f"`evidence:ssf_only` row of convergence.tsv: if those tips "
                      f"form a clade, c71.faa contains two families.",
                      file=sys.stderr)
            if a.out_evidence:
                os.makedirs(os.path.dirname(os.path.abspath(a.out_evidence)),
                            exist_ok=True)
                cols = [c for c in ("seq_id", "sample", "protein_id", "evidence",
                                    "profiles_hit") if c in idmap.columns]
                idmap[cols].to_csv(a.out_evidence, sep="\t", index=False)

        by_faa = defaultdict(dict)
        for faa, pid, sid, smp in zip(idmap["faa"], idmap["protein_id"],
                                      idmap["seq_id"], idmap["sample"]):
            by_faa[faa][pid] = f"{sid} {smp} {pid}"

        n = run(list(by_faa.items()), a.out_faa, a.threads)
        print(f"[extract] wrote {n}/{len(keep)} triad-positive sequences -> {a.out_faa}",
              file=sys.stderr)
        if n != len(keep):
            sys.exit("[extract] some sequences could not be recovered; see .problems.txt")
    else:
        sys.exit("need --hits, or --keep-ids with --idmap")


if __name__ == "__main__":
    main()
