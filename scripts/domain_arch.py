#!/usr/bin/env python3
"""Domain architecture of every C71-carrying protein.

Pei is a two-module enzyme: four N-terminal pseudomurein-binding repeats
(PMBR, PF09373) and a C-terminal C71 catalytic domain (PF12386).

The repeat count is a phenotype, and the threshold is a cliff
------------------------------------------------------------
Visweswaran et al. 2011 fused one, two or three PMB motifs to GFP:

    3 motifs -> binds the pseudomurein sacculus, and bacterial spheroplasts
    2 motifs -> binds bacterial spheroplasts ONLY. No pseudomurein binding.
    1 motif  -> binds nothing

A C71 domain carried on fewer than three PMB motifs is therefore predicted unable
to engage an intact sacculus, whatever its active site looks like. That is a
falsifiable statement about every hit, not a covariate to regress out.

An earlier version of this file recorded the threshold as needed for "optimal"
binding. That is the abstract's wording. The Results, the section heading, the
title and the conclusion all say three motifs are required for *any* binding to
methanogen cells. See `pmbr_reference.py`.

Because the threshold is a cliff, the count is load-bearing. A PMB motif is 30-35
residues; at a strict domain E-value a weak repeat is dropped, the count falls
from 3 to 2, and the functional call inverts. So this module scans ONCE at a
permissive E-value and filters at two, then flags every protein whose class
depends on which threshold you chose. `pmbr_count_fragile` is not a diagnostic,
it is a result: it says the architecture of that protein is not determined by the
data.

What PMBR is not
----------------
It is not evidence that the substrate is pseudomurein. Two- and three-motif
constructs both bind lysozyme-treated *L. lactis* and *E. coli* spheroplasts, and
the signal survives 150 mM NaCl. NAG is the only sugar shared by murein and
pseudomurein, so the module is read as an NAG-binding module. It is also not
Pei-specific: MTH719, an S-layer protein with no catalytic domain, carries three
motifs.

So a bacterial C71 protein with a PMBR array is not automatically a binning
artefact. It may bind exposed murein. A C71 catalytic domain fused to a
*murein*-binding module (LysM, SH3b, PG_binding_1, a choline-binding repeat) is
flagged for the same reason. The architecture tells you which hypothesis to test.
It does not settle it, and neither does this script.

Overlapping domain hits are resolved by keeping the best-scoring non-overlapping
set (a greedy interval selection on the domain i-Evalue), which is what
`hmmscan --domtblout` leaves you to do yourself. Repeats of the same accession are
counted after that resolution, so a tandem array survives and a double-reported
single hit does not.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import Counter, defaultdict

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pmbr_reference import (  # noqa: E402
    CITATION as PMBR_CITATION, MOTIFS_FOR_PSEUDOMUREIN, MOTIFS_FOR_SPHEROPLAST,
    assay_ph_advice, isoelectric_point, predict_binding, rule_has_jurisdiction)
from utils import load_config, parse_domtblout, read_fasta, seq_id  # noqa: E402


def run_hmmscan(pfam, faa, out, evalue, threads):
    if not os.path.exists(pfam):
        sys.exit(f"[domain_arch] Pfam-A.hmm not found: {pfam}")
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        if not os.path.exists(pfam + ext):
            sys.exit(f"[domain_arch] {pfam} is not pressed (missing {ext}). "
                     f"Run: hmmpress {pfam}")
    # --domE alongside -E: without a domain E-value cutoff, a protein whose
    # FULL-sequence E-value exceeds `evalue` is dropped entirely even if it carries
    # one strong PMBR domain. The PMBR count is read at the 30-35 aa repeat level
    # and the 3-motif rule is a cliff, so a per-domain floor matters right where it
    # is load-bearing.
    cmd = ["hmmscan", "--cpu", str(threads), "--noali", "-E", str(evalue),
           "--domE", str(evalue),
           "--domtblout", out, "-o", os.devnull, pfam, faa]
    print("+ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def base_acc(a):
    return str(a).split(".")[0]


def resolve_overlaps(hits, max_overlap=0.4):
    """Greedy: take domains best-first, drop any that overlap an accepted one by
    more than `max_overlap` of the shorter envelope. hmmscan reports every model
    that clears the threshold, including nested and competing families."""
    hits = sorted(hits, key=lambda h: h["i_evalue"])
    kept = []
    for h in hits:
        ok = True
        for k in kept:
            lo = max(h["env_from"], k["env_from"])
            hi = min(h["env_to"], k["env_to"])
            ov = max(0, hi - lo + 1)
            shorter = min(h["env_to"] - h["env_from"] + 1,
                          k["env_to"] - k["env_from"] + 1)
            if shorter > 0 and ov / shorter > max_overlap:
                ok = False
                break
        if ok:
            kept.append(h)
    return sorted(kept, key=lambda h: h["env_from"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--domtbl", required=True, help="write/reuse hmmscan output here")
    ap.add_argument("--out-arch", required=True)
    ap.add_argument("--out-domains", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skip-scan", action="store_true")
    a = ap.parse_args()

    cfg = load_config(a.config)["specificity"]
    pmbr = base_acc(cfg["pmbr_accession"])
    cat = base_acc(cfg["catalytic_accession"])
    accessory = {base_acc(x) for x in cfg["accessory_binding_domains"]}
    min_cov = float(cfg["min_domain_coverage"])

    pcfg = cfg.get("pmbr") or {}
    e_strict = float(pcfg.get("count_evalue_strict", cfg["hmmscan_evalue"]))
    e_loose = float(pcfg.get("count_evalue_permissive", 1e-2))
    if e_loose < e_strict:
        sys.exit(f"[domain_arch] pmbr.count_evalue_permissive ({e_loose}) must be "
                 f"looser (larger) than count_evalue_strict ({e_strict})")

    for p in (a.domtbl, a.out_arch, a.out_domains):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    # Scan once, at the PERMISSIVE threshold, and filter twice in Python. Two
    # scans would double the cost and, worse, would let the two counts diverge
    # for reasons other than the threshold.
    if not a.skip_scan:
        run_hmmscan(cfg["pfam_hmm"], a.faa, a.domtbl, e_loose, a.threads)

    seqs = dict(read_fasta(a.faa))
    seqs = {seq_id(h): s for h, s in seqs.items()}
    lengths = {k: len(v) for k, v in seqs.items()}

    # hmmscan swaps the roles: `target` is the model, `query` is the protein.
    by_prot = defaultdict(list)
    for r in parse_domtblout(a.domtbl):
        try:
            cov = (int(r["hmm_to"]) - int(r["hmm_from"]) + 1) / int(r["tlen"])
            iev = float(r["i_evalue"])
        except (ValueError, ZeroDivisionError):
            continue
        if cov < min_cov or iev > e_loose:
            continue
        by_prot[r["query_name"]].append({
            "acc": base_acc(r["target_accession"]) or r["target_name"],
            "name": r["target_name"],
            "i_evalue": iev,
            "score": float(r["dom_score"]),
            "env_from": int(r["env_from"]),
            "env_to": int(r["env_to"]),
            "model_cov": cov,
        })

    rows, dom_rows = [], []
    for prot, hits in by_prot.items():
        strict_hits = [h for h in hits if h["i_evalue"] <= e_strict]
        kept = resolve_overlaps(strict_hits)
        kept_loose = resolve_overlaps(hits)
        for h in kept_loose:
            dom_rows.append({"seq_id": prot, "passes_strict": int(h["i_evalue"] <= e_strict),
                             **h})

        accs = [h["acc"] for h in kept]
        n_pmbr = sum(a_ == pmbr for a_ in accs)
        n_pmbr_loose = sum(h["acc"] == pmbr for h in kept_loose)
        n_cat = sum(a_ == cat for a_ in accs)
        acc_set = set(accs)
        acc_hits = sorted(acc_set & accessory)

        # N->C string, collapsing tandem repeats of the same accession
        arch = []
        for h in kept:
            if arch and arch[-1][0] == h["acc"]:
                arch[-1][1] += 1
            else:
                arch.append([h["acc"], 1])
        arch_str = "-".join(f"{x}x{n}" if n > 1 else x for x, n in arch)

        # is the catalytic domain C-terminal, as in PeiW/PeiP?
        cat_env = [h for h in kept if h["acc"] == cat]
        cat_cterm = bool(cat_env) and cat_env[-1] is kept[-1]

        # Three motifs bind the sacculus; two bind only lysozyme-exposed murein.
        # The call comes from pmbr_reference so it cannot drift from the paper.
        binding, can_bind, why = predict_binding(n_pmbr)
        binding_loose, can_bind_loose, _ = predict_binding(n_pmbr_loose)
        # If the strict and permissive counts land on opposite sides of the
        # cliff, the architecture of this protein is undetermined. Saying
        # "reduced_pmbr" would be asserting an absence that the E-value chose.
        fragile = int(can_bind != can_bind_loose)

        if n_pmbr == 0 and acc_hits:
            cls = "accessory_binding"     # murein-binding module: different wall?
        elif n_pmbr == 0:
            cls = "catalytic_only"
        elif n_pmbr >= MOTIFS_FOR_PSEUDOMUREIN:
            cls = "canonical_pei"
        else:
            cls = "reduced_pmbr"          # cannot dock on an intact sacculus
        if fragile:
            cls = "pmbr_count_ambiguous"

        # The binding module works near its own pI, and PMB pIs span 3-10, so
        # this is computed per protein rather than borrowed from MTH719.
        pmbr_seq = "".join(seqs.get(prot, "")[h["env_from"] - 1:h["env_to"]]
                           for h in kept if h["acc"] == pmbr)
        pi = isoelectric_point(pmbr_seq) if pmbr_seq else None

        rows.append({
            "seq_id": prot,
            "protein_len": lengths.get(prot, pd.NA),
            "n_pmbr": n_pmbr,
            "n_pmbr_permissive": n_pmbr_loose,
            "pmbr_count_fragile": fragile,
            "pmbr_binding_competent": int(can_bind),
            # Without a PF09373 array the 3-motif rule has no jurisdiction. PeiR
            # has none and lyses M. ruminantium M1. Anything downstream that reads
            # pmbr_binding_competent must check this column first.
            "pmbr_rule_applies": int(rule_has_jurisdiction(n_pmbr)),
            "predicted_binding": binding,
            "predicted_binding_permissive": binding_loose,
            "pmbr_pi": pi if pi is not None else pd.NA,
            "n_catalytic": n_cat,
            "catalytic_is_cterminal": int(cat_cterm),
            "accessory_binding_domains": ",".join(acc_hits),
            "n_domains": len(kept),
            "architecture": arch_str,
            "architecture_class": cls,
            "pmbr_span": (";".join(f"{h['env_from']}-{h['env_to']}"
                                   for h in kept if h["acc"] == pmbr)),
            "catalytic_span": (";".join(f"{h['env_from']}-{h['env_to']}"
                                        for h in kept if h["acc"] == cat)),
            "binding_note": why,
            "assay_ph": assay_ph_advice(pi),
        })

    # proteins with no Pfam hit at all still belong in the table
    for prot in lengths:
        if prot not in by_prot:
            binding, can_bind, why = predict_binding(0)
            rows.append({"seq_id": prot, "protein_len": lengths[prot], "n_pmbr": 0,
                         "n_pmbr_permissive": 0, "pmbr_count_fragile": 0,
                         "pmbr_binding_competent": int(can_bind),
                         "pmbr_rule_applies": 0,
                         "predicted_binding": binding,
                         "predicted_binding_permissive": binding,
                         "pmbr_pi": pd.NA,
                         "n_catalytic": 0, "catalytic_is_cterminal": 0,
                         "accessory_binding_domains": "", "n_domains": 0,
                         "architecture": "", "architecture_class": "no_pfam_hit",
                         "pmbr_span": "", "catalytic_span": "",
                         "binding_note": why, "assay_ph": assay_ph_advice(None)})

    arch = pd.DataFrame(rows).sort_values("seq_id")
    arch.to_csv(a.out_arch, sep="\t", index=False)
    pd.DataFrame(dom_rows).sort_values(["seq_id", "env_from"]).to_csv(
        a.out_domains, sep="\t", index=False)

    print(f"[domain_arch] {len(arch):,} proteins", file=sys.stderr)
    print(arch["architecture_class"].value_counts().to_string(), file=sys.stderr)
    print(f"\n[domain_arch] PMBR repeat count distribution (strict E<={e_strict:g}):",
          file=sys.stderr)
    print(arch["n_pmbr"].value_counts().sort_index().to_string(), file=sys.stderr)

    # --- the binding call ----------------------------------------------------
    n_ok = int(arch["pmbr_binding_competent"].sum())
    n_red = int((arch["architecture_class"] == "reduced_pmbr").sum())
    print(f"\n[domain_arch] {n_ok:,} of {len(arch):,} proteins carry "
          f">= {MOTIFS_FOR_PSEUDOMUREIN} PMB motifs and are predicted able to bind "
          f"an intact pseudomurein sacculus.", file=sys.stderr)
    if n_red:
        print(f"[domain_arch] {n_red:,} carry 1-2 motifs. The two-motif construct "
              f"of Visweswaran et al. 2011 bound lysozyme-treated bacterial "
              f"spheroplasts and did NOT bind pseudomurein, so their PMB arrays "
              f"are predicted unable to dock on an intact sacculus. That is a "
              f"claim about {n_red:,} proteins, and it is testable.",
              file=sys.stderr)

    n_nojur = int((arch["pmbr_rule_applies"] == 0).sum())
    if n_nojur:
        print(f"[domain_arch] {n_nojur:,} proteins carry NO PMB motif, so the "
              f"3-motif rule has no jurisdiction over them. Do not read their "
              f"pmbr_binding_competent=0 as 'cannot lyse': PeiR (D3DZZ6) has no "
              f"PMB motif and lyses Methanobrevibacter ruminantium M1, the one "
              f"host PeiW and PeiP cannot touch. Whatever PeiR docks with, it is "
              f"not a PF09373 array.", file=sys.stderr)

    n_frag = int(arch["pmbr_count_fragile"].sum())
    if n_frag:
        print(f"\n[domain_arch] {n_frag:,} proteins change binding class between "
              f"E<={e_strict:g} and E<={e_loose:g}. A PMB motif is 30-35 residues "
              f"and the threshold at {MOTIFS_FOR_PSEUDOMUREIN} motifs is a cliff, "
              f"so for these the architecture is not determined by the data. They "
              f"are classed `pmbr_count_ambiguous` rather than assigned to "
              f"whichever side the E-value happened to pick.", file=sys.stderr)

    acc = arch.loc[arch["accessory_binding_domains"] != "", "accessory_binding_domains"]
    if len(acc):
        print(f"\n[domain_arch] {len(acc)} proteins carry a non-PMBR binding module:",
              file=sys.stderr)
        print(Counter(acc).most_common(10), file=sys.stderr)
        print("  A separate hypothesis, not outliers. Note that PMBR itself is not "
              "a pseudomurein marker: the domain binds NAG and sticks to "
              "lysozyme-treated bacterial spheroplasts (Visweswaran et al. 2011), "
              "so a bacterial C71+PMBR protein may be binding exposed murein "
              "rather than being a binning artefact.", file=sys.stderr)

    if int((arch["n_pmbr"] > 0).sum()) and arch["pmbr_pi"].notna().any():
        pis = arch["pmbr_pi"].dropna()
        print(f"\n[domain_arch] PMB pI: median {pis.median():.1f}, range "
              f"{pis.min():.1f}-{pis.max():.1f}. The MTH719 domain binds "
              f"pseudomurein completely at pH 9.0 (its pI is 9.2), partially at "
              f"6.5, not at all at 4.0, and aggregates at 7.0. Every published Pei "
              f"lysis assay runs at pH 7.0-7.85.", file=sys.stderr)

    print(f"\n[domain_arch] binding rules from {PMBR_CITATION}", file=sys.stderr)
    if (arch["n_catalytic"] == 0).all():
        sys.exit("[domain_arch] no protein hit the catalytic accession "
                 f"{cat}. Wrong Pfam release, or wrong accession in config.")


if __name__ == "__main__":
    main()
