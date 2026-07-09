#!/usr/bin/env python3
"""Normalise each HMM to HMMER3 ASCII, one file per profile, and validate that
whatever threshold the config asks for actually exists in the model.

The profiles are kept in separate files rather than concatenated because they
carry different thresholds. PF12386 has a curated Pfam gathering threshold, so
`--cut_ga` is available and is the right criterion. SSF54001 is a SCOP
superfamily model and normally has no GA/TC/NC lines; `hmmsearch --cut_ga`
against it dies with 'GA bit thresholds unavailable'. Catching that here, before
700 array jobs launch, is the whole point of this script.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys

CUT_LINE = {"cut_ga": "GA", "cut_nc": "NC", "cut_tc": "TC"}


def hmm_header(path: str) -> dict:
    out = {"name": None, "acc": None, "leng": None, "GA": None, "NC": None, "TC": None}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("HMM "):
                break
            key = line.split(None, 1)[0] if line.strip() else ""
            if key == "NAME" and out["name"] is None:
                out["name"] = line[5:].strip()
            elif key == "ACC" and out["acc"] is None:
                out["acc"] = line[4:].strip()
            elif key == "LENG" and out["leng"] is None:
                out["leng"] = line[5:].strip()
            elif key in ("GA", "NC", "TC") and out[key] is None:
                out[key] = line[len(key):].strip().rstrip(";")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmms", nargs="+", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--thresholds", required=True,
                    help="comma-separated; 'cut_ga'/'cut_nc'/'cut_tc' or a bit score")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out-map", required=True)
    a = ap.parse_args()

    labels = a.labels.split(",")
    thresholds = a.thresholds.split(",")
    if not (len(labels) == len(thresholds) == len(a.hmms)):
        sys.exit(f"{len(labels)} labels / {len(thresholds)} thresholds / {len(a.hmms)} hmms")

    os.makedirs(a.outdir, exist_ok=True)
    have_convert = shutil.which("hmmconvert") is not None
    if not have_convert:
        print("[prepare_hmms] hmmconvert not on PATH; copying models verbatim. "
              "A HMMER2-format SSF54001 will fail at hmmsearch.", file=sys.stderr)

    rows = []
    for label, thr, src in zip(labels, thresholds, a.hmms):
        if not os.path.exists(src):
            sys.exit(f"[prepare_hmms] missing HMM: {src}")
        dst = os.path.join(a.outdir, f"{label}.hmm")

        if have_convert:
            with open(dst, "w") as out:
                r = subprocess.run(["hmmconvert", src], stdout=out, stderr=subprocess.PIPE)
            if r.returncode != 0:
                sys.exit(f"[prepare_hmms] hmmconvert failed on {src}:\n{r.stderr.decode()}")
        else:
            shutil.copy(src, dst)

        h = hmm_header(dst)
        if not h["name"]:
            sys.exit(f"[prepare_hmms] no NAME line in {src} after conversion")

        if thr in CUT_LINE:
            line = CUT_LINE[thr]
            if not h[line]:
                sys.exit(
                    f"[prepare_hmms] config asks for --{thr} on '{label}' but {src} has no "
                    f"{line} line. hmmsearch would abort on every batch.\n"
                    f"Fix one of:\n"
                    f"  - point profiles.{label}.path at the Pfam-A copy of the model "
                    f"(the Pfam-A.hmm release carries GA/TC/NC; the InterPro single-model "
                    f"download often does not), or\n"
                    f"  - set profiles.{label}.threshold to a bit score.")
            print(f"[prepare_hmms] {label}: using --{thr} ({line} = {h[line]})",
                  file=sys.stderr)
        else:
            try:
                float(thr)
            except ValueError:
                sys.exit(f"[prepare_hmms] threshold '{thr}' for {label} is neither "
                         f"cut_ga/cut_nc/cut_tc nor a number")
            if h["GA"]:
                print(f"[prepare_hmms] {label}: using bit score {thr}, but this model "
                      f"does carry a GA line ({h['GA']}). Deliberate?", file=sys.stderr)

        rows.append((label, h["name"], h["acc"] or "-", h["leng"] or "-", thr,
                     h["GA"] or "-", src, dst))

    names = [r[1] for r in rows]
    if len(set(names)) != len(names):
        sys.exit(f"HMM NAME collision {names}: domtblout query names would be ambiguous.")

    with open(a.out_map, "w") as fh:
        fh.write("label\thmm_name\thmm_acc\thmm_len\tthreshold\tga_line\tsource_path\thmm_path\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    print(f"[prepare_hmms] {len(rows)} profiles -> {a.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
