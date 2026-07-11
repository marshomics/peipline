#!/usr/bin/env python3
"""Genotype each genome's pseudomurein biosynthesis pathway.

Pei cleaves the epsilon-isopeptide bond between alanine and lysine. The alanine
is P1, and P1 is not fixed across hosts: it is threonine in Methanobrevibacter
ruminantium, serine in M. stadtmanae, and M. smithii / M. kandleri carry an
ornithine modification of the cross-link. If Pei subgroups mean anything
functional, they should track the substrate the host actually builds.

Taxonomy first, then the block, then never "confirmed"
------------------------------------------------------
Lupo et al. 2025 make two points that reorganise this step (see pmur_reference.py):

  * PM is restricted to two orders, Methanobacteriales and Methanopyrales, so
    "does this genome make PM?" is near-equivalent to "is it in one of those two
    orders?" GTDB placement answers that more reliably than an OG screen built
    from five genomes. So taxonomy is the PRIMARY signal:
        in an order that contains PM  -> pseudomurein_expected_by_taxonomy
    and the marker screen is corroboration, not the call.

  * The OG screen earns its keep only OUTSIDE the two orders -- testing whether
    PM turns up somewhere unexpected. For that, a bag of markers is worthless: 15
    of the 49 PM OGs are widespread across Bacteria and Archaea, and the Mur
    domain in particular (MurT, CapB, CphA, FPGS) is everywhere. The signal is a
    small PM-EXCLUSIVE block -- a muramyl ligase (Muralpha..delta) co-present
    with the MraY-like GT (OG0001163) and the CPS (OG0001014). Only an
    out-of-order genome carrying that whole block is flagged, as a candidate:
        out of order + block co-present -> pseudomurein_candidate_out_of_order

Neither call confirms PM. Genomic presence is a hypothesis. Confirmation is the
wall chemistry itself: N-acetyl-talosaminuronic acid (TalNAc) and the beta-(1->3)
glycan linkage. And because the HMMs come from five genomes, they miss divergent
lineages, so an out-of-order NEGATIVE is weak evidence, not absence.

The marker HMMs are still whatever you put in `pmur_hmm_dir`; this module only
knows which of them are block members (by OG id or name in the filename, or the
`specificity.pmur_block` config map).

Three further things this deliberately does not do.

It does not call a genome "Ala-type" because a marker is missing. A marker can
be absent because the pathway is absent, because the MAG is 70% complete, or
because the HMM is bad. Genomes below `pmur_min_markers` are called
`no_pathway_detected`, and the completeness of every genome is carried forward
so the regression can condition on it.

It does not infer P1 chemistry from the marker set, because no marker in the
published set is known to determine P1. The P1 call comes from
`cellwall_reference.py`, which encodes Kandler & Koenig 1978 at the level the
paper actually supports: species, and in two cases a single type strain.

It does not assign P1 from a genus. Methanobrevibacter contains M. ruminantium
(Thr at P1), M. smithii (Ala at P1, ~1/4 Orn at P1') and M. arboriphilus (Ala,
Lys). A genus-level call would put every gut Methanobrevibacter in one bucket
and the bucket would be wrong.

It does not silently use a bit-score cutoff when a model carries a GA line.
"""
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cellwall_reference import CITATION, annotate, genus_of  # noqa: E402
import synteny  # noqa: E402
from utils import load_config, parse_domtblout, read_fasta, seq_id  # noqa: E402
from pmur_reference import (  # noqa: E402
    BLOCK_ANCHORS, CITATION as PMUR_CITATION, MURAMYL_LIGASES, MUR_DOMAIN_TRAP,
    PM_EXCLUSIVE, PM_ORDERS, parse_order, pm_expected_by_taxonomy)

# The same alphabet batch_faa.py enforces. HMMER rejects "*", which Prodigal
# appends to every protein it calls -- and Prodigal calls the archaea, the only
# genomes that can make pseudomurein.
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUO")


def classify_markers(markers, block_cfg=None):
    """Map each marker (HMM filename stem) to a PM-block role.

    Roles: 'muramyl_ligase', 'mray_like', 'cps', 'trap', or 'other'. An explicit
    config map (specificity.pmur_block) wins; otherwise match the OG id or the
    canonical name inside the filename, case-insensitively.
    """
    # canonical name -> role, from the Lupo catalogue
    name_role = {}
    for og, d in PM_EXCLUSIVE.items():
        role = ("mray_like" if og == "OG0001163"
                else "cps" if og == "OG0001014"
                else "muramyl_ligase" if og in MURAMYL_LIGASES else "other")
        name_role[og.lower()] = role
        name_role[d["name"].lower()] = role
    trap_names = {t.lower() for t in MUR_DOMAIN_TRAP}

    override = {}
    if block_cfg:
        for role, names in block_cfg.items():
            for nm in (names or []):
                override[str(nm)] = ("muramyl_ligase" if role.startswith("muramyl")
                                     else "mray_like" if "mray" in role
                                     else "cps" if role == "cps" else "other")

    out = {}
    for m in markers:
        if m in override:
            out[m] = override[m]
            continue
        ml = m.lower()
        if ml in trap_names or any(t in ml for t in trap_names):
            out[m] = "trap"
            continue
        role = "other"
        for key, r in name_role.items():
            if key in ml:
                role = r
                break
        out[m] = role
    return out


def annotate_pathway_calls(g, lig, mray, cps, block_ok, min_markers, cls_col):
    """Add the taxonomy-first PM call and the out-of-order block flag to `g`.

    `g` already carries one `pmur_<marker>` 0/1 column per marker. `lig`, `mray`
    and `cps` are the marker names of the muramyl ligases, the MraY-like anchor
    and the CPS anchor. Pure: no I/O, so the test can drive the real code.

    Every value it writes is a GENOMIC hypothesis, not confirmed pseudomurein.
    """
    import pandas as pd

    marker_cols = [c for c in g.columns if c.startswith("pmur_")]
    g["n_pmur_markers"] = (g[marker_cols].sum(axis=1) if marker_cols
                           else pd.Series(0, index=g.index))

    # taxonomy is primary
    if cls_col and cls_col in g.columns:
        g["gtdb_order"] = g[cls_col].map(parse_order)
        g["pm_expected_by_taxonomy"] = g[cls_col].map(pm_expected_by_taxonomy).astype(int)
    else:
        g["gtdb_order"] = pd.NA
        g["pm_expected_by_taxonomy"] = 0

    def any_of(names):
        cols = [f"pmur_{m}" for m in names if f"pmur_{m}" in g.columns]
        return (g[cols].max(axis=1) if cols else pd.Series(0, index=g.index)).astype(int)

    lig_cols = [f"pmur_{m}" for m in lig if f"pmur_{m}" in g.columns]
    g["n_muramyl_ligases"] = (g[lig_cols].sum(axis=1) if lig_cols
                              else pd.Series(0, index=g.index))
    g["has_mray_like"] = any_of(mray)
    g["has_cps"] = any_of(cps)
    if block_ok:
        g["pm_block_present"] = (((g["n_muramyl_ligases"] >= 1) &
                                  (g["has_mray_like"] == 1) &
                                  (g["has_cps"] == 1)).astype(int))
    else:
        g["pm_block_present"] = 0
    g["pmur_count_pathway"] = (g["n_pmur_markers"] >= int(min_markers)).astype(int)

    g["pathway_call"] = "no_pathway_detected"
    g.loc[g["pm_expected_by_taxonomy"] == 1,
          "pathway_call"] = "pseudomurein_expected_by_taxonomy"
    ooo = (g["pm_expected_by_taxonomy"] == 0) & (g["pm_block_present"] == 1)
    g.loc[ooo, "pathway_call"] = "pseudomurein_candidate_out_of_order"
    weak = ((g["pm_expected_by_taxonomy"] == 0) & (g["pm_block_present"] == 0) &
            (g["n_pmur_markers"] > 0))
    g.loc[weak, "pathway_call"] = "markers_without_block"
    low = (g["pathway_call"] == "no_pathway_detected") & (g["completeness"] < 90)
    g.loc[low, "pathway_call"] = "indeterminate_low_completeness"
    return g


def refine_synteny(g, block_hits, block_hits_perm, scfg, verbose=False):
    """Run positional synteny for out-of-order candidates and refine the call.

    Two tiers feed it, and they are treated differently:

      strict     genome already flagged `..._candidate_out_of_order` by the strict
                 block. Synteny only REFINES it: _syntenic / _dispersed /
                 _synteny_unknown. It stays a candidate either way.
      divergent  out of order, NO strict block, but a permissive (below-GA) block
                 exists. Synteny is the ONLY thing that can elevate it. Syntenic
                 -> `..._divergent_syntenic`; anything else -> left as its base
                 call, because a permissive hit without synteny is noise. This is
                 the false-positive control for the divergent-lineage mode.

    Pure except for the optional log. Both main() and the test call it.
    """
    import pandas as pd

    def _has_perm_block(sample):
        bh = block_hits_perm.get(sample)
        return bool(bh and bh["muramyl_ligase"] and bh["mray_like"] and bh["cps"])
    g["pm_block_present_permissive"] = g["sample"].map(_has_perm_block).astype(int)

    for c in ("synteny_status", "synteny_detail", "synteny_tier"):
        g[c] = pd.NA

    strict_ooo = g["pathway_call"] == "pseudomurein_candidate_out_of_order"
    divergent = ((g["pm_expected_by_taxonomy"] == 0) &
                 (g["pm_block_present"] == 0) &
                 (g["pm_block_present_permissive"] == 1))
    if not scfg.get("enabled", True):
        return g

    gi = g.set_index("sample")
    for sample in g.loc[strict_ooo | divergent, "sample"]:
        row = gi.loc[sample]
        is_strict = bool(gi.loc[sample, "pm_block_present"])
        bh = block_hits.get(sample) if is_strict else block_hits_perm.get(sample)
        tier = "strict" if is_strict else "divergent"
        ncont = row.get("contigs")
        ncont = int(ncont) if pd.notna(ncont) else None
        coords, origin = synteny.load_coords(row.get("source"), row.get("faa"),
                                             row.get("gff"))
        if not coords:
            res = {"status": "not_evaluable", "reason": origin}
        else:
            res = synteny.block_synteny(
                bh, coords, n_contigs=ncont,
                window_genes=int(scfg.get("window_genes", 12)),
                window_bp=int(scfg.get("window_bp", 15000)),
                max_contigs=int(scfg.get("max_contigs_for_synteny", 200)))
        idx = g.index[g["sample"] == sample]
        g.loc[idx, "synteny_status"] = res["status"]
        g.loc[idx, "synteny_detail"] = res.get("reason", "")
        g.loc[idx, "synteny_tier"] = tier
        if is_strict:
            suffix = {"syntenic": "_syntenic", "dispersed": "_dispersed",
                      "not_evaluable": "_synteny_unknown"}[res["status"]]
            g.loc[idx, "pathway_call"] = "pseudomurein_candidate_out_of_order" + suffix
        elif res["status"] == "syntenic":
            g.loc[idx, "pathway_call"] = \
                "pseudomurein_candidate_out_of_order_divergent_syntenic"
        if verbose:
            print(f"[cellwall]   {sample} ({tier}): {res['status']} -- "
                  f"{res.get('reason', '')}", file=sys.stderr)
    return g


def hmm_has_ga(path):
    return hmm_ga_value(path) is not None


def hmm_ga_value(path):
    """The sequence GA (gathering) score, or None. GA line: `GA  25.00 22.00;`."""
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("HMM "):
                return None
            if line.startswith("GA "):
                try:
                    return float(line.split()[1].rstrip(";"))
                except (IndexError, ValueError):
                    return None
    return None


def scan_scores(hmm, faa, out, threshold, threads):
    """Run ONE permissive hmmsearch and return {protein_id: best full-seq score}.

    Searching once at a low threshold and splitting the hits in Python (strict vs
    permissive-only, by each model's own GA) is how the divergent-lineage mode
    stays cheap: it costs one search, not two. The strict set reproduces the old
    behaviour exactly; the permissive-only set is offered to the out-of-order
    synteny gate and nowhere else.
    """
    subprocess.run(["hmmsearch", "--cpu", str(threads), "--noali",
                    "-T", str(threshold), "--domT", str(threshold),
                    "--incT", str(threshold), "--incdomT", str(threshold),
                    "--domtblout", out, "-o", os.devnull, hmm, faa], check=True)
    best = {}
    for r in parse_domtblout(out):
        try:
            s = float(r["full_score"])
        except (KeyError, ValueError):
            continue
        t = r["target_name"]
        if s > best.get(t, float("-inf")):
            best[t] = s
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="unified sample table")
    ap.add_argument("--genomes", required=True, help="genome_level_table.tsv")
    ap.add_argument("--config", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-markers", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--only-with-c71", action="store_true",
                    help="screen only genomes carrying a C71 (much cheaper)")
    a = ap.parse_args()

    cfg = load_config(a.config)["specificity"]
    hmms = sorted(glob.glob(os.path.join(cfg["pmur_hmm_dir"], "*.hmm")))
    if not hmms:
        sys.exit(f"[cellwall] no .hmm in {cfg['pmur_hmm_dir']}")
    markers = [os.path.splitext(os.path.basename(h))[0] for h in hmms]
    print(f"[cellwall] {len(markers)} markers: {markers}", file=sys.stderr)

    # Classify each user-supplied marker against the Lupo 2025 PM-exclusive OGs,
    # so the "block" test knows which markers are muramyl ligases, which is the
    # MraY-like anchor, and which is the CPS anchor. Matching is by OG id or
    # canonical name appearing in the filename (case-insensitive), or an explicit
    # config map. A marker that matches nothing is still counted, but it cannot
    # contribute to the block -- and if it is really one of the Mur-domain trap
    # families, it should not be here at all.
    role_of = classify_markers(markers, cfg.get("pmur_block"))
    lig = [m for m, r in role_of.items() if r == "muramyl_ligase"]
    mray = [m for m, r in role_of.items() if r == "mray_like"]
    cps = [m for m, r in role_of.items() if r == "cps"]
    trap = [m for m, r in role_of.items() if r == "trap"]
    if trap:
        print(f"[cellwall] WARNING: {trap} look like Mur-domain TRAP families "
              f"(MurT/CapB/CphA/FPGS). Those carry a Mur domain but are scattered "
              f"across Bacteria and Archaea; a hit is not PM evidence. Remove them.",
              file=sys.stderr)
    block_ok = bool(lig) and bool(mray) and bool(cps)
    if not block_ok:
        print(f"[cellwall] NOTE: the marker set does not contain the full "
              f"PM-exclusive block (>=1 muramyl ligase + MraY-like + CPS). "
              f"Found ligases={lig}, MraY-like={mray}, CPS={cps}. Out-of-order "
              f"discovery falls back to marker count, which is weaker. See "
              f"pmur_reference.py for the OGs to build HMMs for.", file=sys.stderr)

    os.makedirs(a.workdir, exist_ok=True)
    for p in (a.out, a.out_markers):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    tab = pd.read_csv(a.table, sep="\t", dtype={"sample": str}, low_memory=False)
    gen = pd.read_csv(a.genomes, sep="\t", dtype={"sample": str}, low_memory=False)
    gen = gen[["sample", "has_c71", "completeness", "contamination", "domain"]]
    tab = tab.merge(gen, on="sample", how="inner", suffixes=("", "_g"))

    if a.only_with_c71:
        tab = tab[tab["has_c71"] == 1]
    # Only archaea can make pseudomurein. Screening 340k bacteria for Pmur is a
    # negative control worth running once, not every time.
    print(f"[cellwall] screening {len(tab):,} genomes "
          f"({int((tab['domain'] == 'Archaea').sum()):,} archaea)", file=sys.stderr)
    if tab.empty:
        sys.exit("[cellwall] no genomes to screen")

    # One concatenated faa, sequences prefixed by sample, so we scan once.
    #
    # This used to be a plain `open(faa, errors="replace")` line loop. Two silent
    # failures, both on inputs this pipeline actually produces:
    #
    #   * a gzipped proteome was opened as text, decode errors were SWALLOWED by
    #     errors="replace", no '>' was ever seen, and the gzip bytes were written
    #     out as "sequence". The genotype for that genome was quietly wrong.
    #   * Prodigal's `-a` output appends '*' to every protein. batch_faa.py strips
    #     it because HMMER rejects it; this did not. Prodigal calls the ARCHAEA,
    #     which are the only genomes that can make pseudomurein.
    #
    # Use the same gzip-aware reader and the same alphabet filter as the screen.
    cat = os.path.join(a.workdir, "genomes.faa")
    n_seq, n_missing = 0, 0
    with open(cat, "w") as out:
        for sample, faa in zip(tab["sample"], tab["faa"]):
            if not isinstance(faa, str) or not os.path.exists(faa):
                n_missing += 1
                continue
            if "|" in str(sample):
                sys.exit(f"[cellwall] sample id {sample!r} contains '|', which is "
                         f"the provenance delimiter. Hits would be misattributed.")
            for _hdr, seq in read_fasta(faa):
                clean = "".join(c for c in seq.upper() if c in VALID_AA)
                if not clean:
                    continue
                # Preserve the protein id (was a running counter). The synteny
                # check needs it to look the protein up in the genome's GFF / the
                # Prodigal header. sample carries no "|" (guarded above), so
                # split("|", 1) recovers (sample, protein_id) even if the id has a
                # pipe of its own.
                out.write(f">{sample}|{seq_id(_hdr)}\n{clean}\n")
                n_seq += 1
    if n_missing:
        print(f"[cellwall] {n_missing} genomes have no readable faa and were skipped",
              file=sys.stderr)
    n_screened = len(tab) - n_missing
    print(f"[cellwall] {n_seq:,} proteins from {n_screened:,} genomes", file=sys.stderr)

    strict_thr = float(cfg["pmur_score_threshold"])
    perm_thr = float(cfg.get("pmur_score_threshold_permissive", strict_thr))
    permissive_ooo = bool((cfg.get("synteny") or {}).get("permissive_out_of_order",
                                                         perm_thr < strict_thr))
    base_thr = min(perm_thr, strict_thr) if permissive_ooo else strict_thr
    # A model's strict bar is its own GA line (below). If any GA sits BELOW
    # base_thr, a hit scoring between that GA and base_thr passes the model's
    # gathering cutoff yet would never be reported by a single search floored at
    # base_thr -- silently truncating strict presence for that marker. Floor the
    # one search at min(base_thr, every model's GA) so nothing above a model's own
    # bar is lost. GA values are read once here and reused in the loop.
    gas = {hmm: hmm_ga_value(hmm) for hmm in hmms}
    ga_vals = [g for g in gas.values() if g is not None]
    search_thr = min([base_thr] + ga_vals)
    if ga_vals and min(ga_vals) < base_thr:
        print(f"[cellwall] lowered marker search floor to {search_thr:g} to cover a "
              f"model GA below the permissive threshold {base_thr:g}", file=sys.stderr)

    hits = defaultdict(set)              # genome presence, STRICT (unchanged)
    # per-genome, per-role protein ids for the block members, for synteny.
    # `_strict` drives the standard block; `_perm` (strict + permissive-only) is
    # offered ONLY to the out-of-order synteny gate, so a divergent homolog that
    # scores below GA can still be a candidate -- but only if the arrangement
    # holds. Permissive without synteny is never elevated.
    def _blk():
        return {"muramyl_ligase": set(), "mray_like": set(), "cps": set()}
    block_hits = defaultdict(_blk)
    block_hits_perm = defaultdict(_blk)
    thresholds = {}
    for hmm, m in zip(hmms, markers):
        with tempfile.NamedTemporaryFile(suffix=".domtbl", dir=a.workdir,
                                         delete=False) as tf:
            dom = tf.name
        # each model's strict bar is its GA if it has one, else the strict number
        ga = gas[hmm]
        strict_cut = ga if ga is not None else strict_thr
        thresholds[m] = f"cut_ga({ga})" if ga is not None else str(strict_thr)
        scores = scan_scores(hmm, cat, dom, search_thr, a.threads)
        os.remove(dom)
        role = role_of.get(m)
        n_strict = n_perm = 0
        for target, score in scores.items():
            sample_id, _, protein_id = target.partition("|")
            is_strict = score >= strict_cut
            if is_strict:
                hits[m].add(sample_id)
                n_strict += 1
            else:
                n_perm += 1
            if role in ("muramyl_ligase", "mray_like", "cps") and protein_id:
                block_hits_perm[sample_id][role].add(protein_id)
                if is_strict:
                    block_hits[sample_id][role].add(protein_id)
        print(f"[cellwall] {m}: {n_strict:,} strict hits in {len(hits[m]):,} "
              f"genomes (>= {strict_cut}); {n_perm:,} permissive-only "
              f"(>= {search_thr}, for out-of-order synteny only)", file=sys.stderr)

    pd.DataFrame({"marker": markers,
                  "threshold": [thresholds[m] for m in markers],
                  "n_genomes_hit": [len(hits[m]) for m in markers],
                  "frac_genomes_hit": [len(hits[m]) / len(tab) for m in markers],
                  "has_ga_line": [int(hmm_has_ga(h)) for h in hmms]}
                 ).to_csv(a.out_markers, sep="\t", index=False)

    cls_col = "classification" if "classification" in tab.columns else None
    keep = ["sample", "domain", "species", "completeness", "contamination", "has_c71"]
    for extra in (cls_col, "source", "faa", "gff", "contigs"):
        if extra and extra in tab.columns and extra not in keep:
            keep.append(extra)
    g = tab[keep].copy()
    for m in markers:
        g[f"pmur_{m}"] = g["sample"].isin(hits[m]).astype(int)

    g = annotate_pathway_calls(g, lig, mray, cps, block_ok,
                               int(cfg["pmur_min_markers"]), cls_col)
    if cls_col is None:
        print("[cellwall] no classification column: cannot use taxonomy, which is "
              "the strongest PM signal. Falling back to markers only.", file=sys.stderr)

    scfg = cfg.get("synteny") or {}
    g = refine_synteny(g, block_hits, block_hits_perm, scfg, verbose=True)
    strict_ooo = g["pathway_call"] == "pseudomurein_candidate_out_of_order_syntenic"
    low = g["pathway_call"] == "indeterminate_low_completeness"
    n_ooo = int(g["pathway_call"].astype(str).str.startswith(
        "pseudomurein_candidate_out_of_order").sum())
    n_div = int(g.get("synteny_tier", pd.Series(dtype=object)).eq("divergent").sum())

    n_div_syn = int((g["pathway_call"] ==
                     "pseudomurein_candidate_out_of_order_divergent_syntenic").sum())
    if n_ooo or n_div:
        syn = int((g["synteny_status"] == "syntenic").sum())
        print(f"[cellwall] *** out-of-order: {n_ooo} with the strict block, {n_div} "
              f"divergent (permissive block only). {syn} are SYNTENIC "
              f"({n_div_syn} of them divergent -- weak homologs rescued by the "
              f"cluster arrangement). Every one is a HYPOTHESIS: confirm at the "
              f"bench (TalNAc, beta-1,3), and consider the structure search.",
              file=sys.stderr)
    print(f"[cellwall] {int(low.sum()):,} genomes indeterminate (<90% complete, "
          f"no positive signal); not called negative.", file=sys.stderr)
    print(f"[cellwall] PM markers from {PMUR_CITATION}; taxonomy is the primary "
          f"signal inside {sorted(PM_ORDERS)}.", file=sys.stderr)

    # Literature wall chemistry, species-level, never inferred from markers and
    # never inherited from a heterogeneous genus. See cellwall_reference.py.
    g = pd.concat([g, annotate(g["species"])], axis=1)
    g["genus"] = g["species"].map(genus_of)
    g["p1_citation"] = CITATION

    g.to_csv(a.out, sep="\t", index=False)
    print("\n[cellwall] pathway calls:", file=sys.stderr)
    print(g["pathway_call"].value_counts().to_string(), file=sys.stderr)
    print(f"\n[cellwall] wall chemistry, from {CITATION}", file=sys.stderr)
    print("\n  P1 (acyl donor; the residue Pei's S1 pocket reads):", file=sys.stderr)
    print(g["p1_residue"].value_counts().to_string(), file=sys.stderr)
    print("\n  P1' (acyl acceptor; Lys epsilon-amine, or Orn delta-amine):",
          file=sys.stderr)
    print(g["p1_prime_residue"].value_counts().to_string(), file=sys.stderr)
    print("\n  provenance:", file=sys.stderr)
    print(g["p1_source"].value_counts().to_string(), file=sys.stderr)

    n_het = int((g["p1_source"] == "genus_heterogeneous").sum())
    if n_het:
        print(f"\n[cellwall] {n_het:,} genomes are in a genus whose characterised "
              f"species disagree about the wall (Methanobrevibacter has Ala/Lys, "
              f"Ala/Orn and Thr/Lys members). They are 'unknown', not guessed.",
              file=sys.stderr)
    disp = g.loc[g["p1_source"] == "disputed", "species"].value_counts()
    if len(disp):
        print(f"\n[cellwall] {int(disp.sum()):,} genomes belong to species whose "
              f"wall chemistry is asserted in secondary sources but not in the "
              f"primary reference this pipeline uses:", file=sys.stderr)
        print(disp.to_string(), file=sys.stderr)
        print("  Supply the primary citation and add them to "
              "cellwall_reference.REFERENCE before relying on them.", file=sys.stderr)
    n_nops = int((g["p1_residue"] == "no_pseudomurein").sum())
    if n_nops:
        print(f"\n[cellwall] {n_nops:,} genomes have no pseudomurein sacculus at "
              f"all (protein sheath). They are a negative control, not Ala-type.",
              file=sys.stderr)

    bac = g[(g["domain"] == "Bacteria") & (g["pmur_count_pathway"] == 1)]
    if len(bac):
        print(f"\n[cellwall] WARNING: {len(bac)} BACTERIA appear to carry the "
              f"pseudomurein pathway. Bacteria do not make pseudomurein. Either "
              f"the markers cross-react with MurC/MurE (they are homologous), or "
              f"those genomes are contaminated. Check before using this column.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
