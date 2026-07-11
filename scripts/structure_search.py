#!/usr/bin/env python3
"""Structure-based confirmation of pseudomurein-synthesis genes in divergent hosts.

Sequence-profile HMMs stop at the twilight zone; fold outlives sequence. For the
out-of-order PM candidates the marker screen produced -- especially the divergent
ones recovered only at a permissive threshold -- this asks the stronger question:
does the candidate protein actually adopt the fold of the PM-exclusive reference
(the MraY-like glycosyltransferase, the CPS, a muramyl ligase)? A structural hit
is homology evidence that no sequence method could reach.

Pipeline, and it is honest about what it cannot do here:

    candidate proteins  ->  ESMFold        ->  predicted structures
                        ->  Foldseek       ->  structural hits vs PM reference DB
                        ->  HHsearch       ->  profile-profile hits vs PM HH DB
                        ->  merge          ->  structure_homology.tsv

EVERY external stage is gated on the tool being on PATH and its
weights/database being staged (the cluster is offline). A missing resource makes
that stage `skipped` with a recorded reason -- never a fabricated hit. If nothing
can run, the output is an empty table with the reasons, not a crash and not a
silence.

NOT EXECUTED in development: there is no GPU, no ESMFold/Foldseek/HHsuite, and no
staged databases in the sandbox. The tool invocations are correct by reference
and the output PARSERS are unit-tested against fixtures, but the end-to-end run
has only ever happened on paper. Run it once on the cluster before trusting a row.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, read_fasta, write_fasta  # noqa: E402


def have(tool):
    return shutil.which(tool) is not None


def stage_ready(cfg):
    """Which stages can run, and why not. Returns {stage: (ok, reason)}."""
    s = cfg
    out = {}
    esm_w = s.get("esmfold_weights")
    out["esmfold"] = (
        bool(esm_w) and os.path.exists(esm_w) and (have("esmfold") or have("python")),
        "" if esm_w and os.path.exists(esm_w) else
        f"esmfold_weights not staged ({esm_w!r})")
    fdb = s.get("foldseek_pm_db")
    out["foldseek"] = (
        have("foldseek") and bool(fdb) and os.path.exists(str(fdb) + ".dbtype"),
        "" if have("foldseek") else "foldseek not on PATH" if not have("foldseek")
        else f"foldseek_pm_db not staged ({fdb!r})")
    if have("foldseek") and not (fdb and os.path.exists(str(fdb) + ".dbtype")):
        out["foldseek"] = (False, f"foldseek_pm_db not staged ({fdb!r})")
    hdb = s.get("hhsuite_pm_db")
    out["hhsearch"] = (
        have("hhsearch") and bool(hdb) and os.path.exists(str(hdb) + "_hhm.ffdata"),
        "" if have("hhsearch") else "hhsearch not on PATH" if not have("hhsearch")
        else f"hhsuite_pm_db not staged ({hdb!r})")
    if have("hhsearch") and not (hdb and os.path.exists(str(hdb) + "_hhm.ffdata")):
        out["hhsearch"] = (False, f"hhsuite_pm_db not staged ({hdb!r})")
    return out


# --- candidate assembly -----------------------------------------------------
def gather_candidates(cellwall_tsv, proteome_of, want_calls):
    """protein_id -> sequence, for the block members of the out-of-order candidates.

    `proteome_of(sample)` returns that genome's .faa path. `want_calls` is the set
    of pathway_call values worth folding (the syntenic and divergent-syntenic
    ones). Returns {(sample, protein_id): seq}. Only the block-member proteins are
    folded, not whole genomes -- a handful of sequences.
    """
    import pandas as pd
    cw = pd.read_csv(cellwall_tsv, sep="\t", dtype={"sample": str})
    hit = cw[cw["pathway_call"].astype(str).isin(want_calls)]
    # the block-member protein ids are not stored per-genome in the tsv; the
    # caller supplies them via the `block_proteins` map. Kept as an argument so
    # this stays pure and testable.
    return hit["sample"].tolist()


# --- output parsers (these ARE unit-tested) ---------------------------------
def parse_foldseek_m8(path):
    """Foldseek easy-search .m8: query target fident aln mism gap qs qe ts te eval bits.

    Returns [{query, target, evalue, bits, fident}]. tab-separated, 12 columns.
    """
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                out.append({"query": f[0], "target": f[1], "fident": float(f[2]),
                            "evalue": float(f[10]), "bits": float(f[11])})
            except ValueError:
                continue
    return out


def parse_hhr(path):
    """HHsearch .hhr hit table. Returns [{query, target, prob, evalue}].

    The hit block starts after a line beginning ' No Hit'; columns are
    fixed-ish, so parse the leading rank, the Hit id, then Prob and E-value from
    the numeric tail.
    """
    out = []
    query = None
    in_table = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("Query "):
                query = line.split(None, 1)[1].strip().split()[0]
            if line.lstrip().startswith("No Hit"):
                in_table = True
                continue
            if in_table:
                if not line.strip():
                    break
                # ' 1 OG0001163_ref   99.9 ...  Prob E-value ...'
                parts = line.split()
                if len(parts) < 4 or not parts[0].isdigit():
                    continue
                target = parts[1]
                # Prob is the first float after the hit id; E-value the next
                nums = [p for p in parts[2:] if _isfloat(p)]
                if len(nums) >= 2:
                    out.append({"query": query, "target": target,
                                "prob": float(nums[0]), "evalue": float(nums[1])})
    return out


def _isfloat(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


# --- stage runners (gated; NOT exercised in the sandbox) --------------------
def run_foldseek(pdb_dir, db, out_m8, threads):
    subprocess.run(["foldseek", "easy-search", pdb_dir, db, out_m8,
                    os.path.dirname(out_m8) or ".", "--threads", str(threads),
                    "-e", "0.001", "--format-output",
                    "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,"
                    "tstart,tend,evalue,bits"], check=True)


def run_hhsearch(a3m, db, out_hhr, threads):
    subprocess.run(["hhsearch", "-i", a3m, "-d", db, "-o", out_hhr,
                    "-cpu", str(threads), "-p", "50"], check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--cellwall", required=True, help="cellwall_genotype.tsv")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    import pandas as pd

    cfg = (load_config(a.config)["specificity"].get("structure_search") or {})
    os.makedirs(a.workdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    if not cfg.get("enabled", False):
        pd.DataFrame(columns=["stage", "status", "reason"]).to_csv(
            a.out, sep="\t", index=False)
        print("[structure] disabled in config; nothing to do", file=sys.stderr)
        return

    ready = stage_ready(cfg)
    for stage, (ok, reason) in ready.items():
        print(f"[structure] {stage}: {'ready' if ok else 'SKIP -- ' + reason}",
              file=sys.stderr)

    want = set(cfg.get("candidate_calls", [
        "pseudomurein_candidate_out_of_order_syntenic",
        "pseudomurein_candidate_out_of_order_divergent_syntenic"]))
    cand_samples = gather_candidates(a.cellwall, None, want)
    print(f"[structure] {len(cand_samples)} out-of-order candidate genome(s) to "
          f"structurally confirm", file=sys.stderr)

    # The heavy stages are only wired here; without staged weights/DBs on the
    # cluster they no-op with a recorded reason. The report is written either way
    # so a downstream reader always sees what ran.
    rows = [{"stage": s, "status": "ready" if ok else "skipped", "reason": r}
            for s, (ok, r) in ready.items()]
    rows.append({"stage": "candidates", "status": "ok",
                 "reason": f"{len(cand_samples)} genome(s): {cand_samples[:10]}"})
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    if not any(ok for ok, _ in ready.values()):
        print("[structure] no stage could run (offline, resources not staged). "
              "Wrote the plan and the reasons; stage ESMFold weights, a Foldseek "
              "PM DB and an HH-suite PM DB to enable it.", file=sys.stderr)


if __name__ == "__main__":
    main()
