#!/usr/bin/env python3
"""Is the pseudomurein block actually a syntenic cluster, or scattered hits?

Co-presence of the PM-exclusive block (a muramyl ligase + the MraY-like GT + the
CPS) is a much stronger signal than a bag of markers, but it is still weaker than
what Lupo et al. 2025 actually use to define the pathway: the genes sit together,
arranged as clusters A and B. A protein family transferred once, or two unrelated
Mur-domain hits on different replicons, gives co-presence without synteny.

This module tests the stronger claim: do the block members co-localize on ONE
contig, within a window? It is only ever run on the handful of OUT-OF-ORDER
candidates (genomes carrying the block but not in Methanobacteriales or
Methanopyrales), because inside those two orders taxonomy already answers the
question and synteny would be busywork.

The honest failure mode is a fragmented assembly. If the block members land on
different contigs of a 400-contig MAG, that is very likely assembly breakage, not
a biological absence of synteny. So the result is conditioned on contiguity, and
a broken assembly returns `not_evaluable`, never `dispersed`.

Coordinates come from a GFF3 per genome. Bacteria: the path the user added to the
sample table. Archaea: the GFF Prodigal now emits. The join is by protein id: the
.faa protein id must appear in the GFF (as ID=, locus_tag=, protein_id= or
Name=). Lupo-style co-localization needs order on a contig, not exact distances,
so a gene-index window is the primary test and a bp window is the fallback.
"""
from __future__ import annotations

import gzip
import re

# window, in genes, within which the whole block must fall on one contig.
# Cluster A is 5 genes, cluster B is 3, and the floating ligases sit between or
# just outside, so the block can legitimately span ~10-12 genes. Default is
# generous; tighten in config if you want a stricter call.
DEFAULT_WINDOW_GENES = 12
DEFAULT_WINDOW_BP = 15000
# above this many contigs, "different contigs" is uninformative: the assembly is
# too fragmented to conclude anything about synteny.
DEFAULT_MAX_CONTIGS_FOR_SYNTENY = 200

_ATTR_KEYS = ("ID", "locus_tag", "protein_id", "Name", "gene", "old_locus_tag")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_gff(path, feature_types=("CDS", "gene", "mRNA")):
    """protein_id -> (contig, start, end, strand, order_index_on_contig).

    Every attribute that could hold the .faa protein id points at the same coord
    record, because different annotators put the joinable id in different fields
    and the caller only knows it as "the .faa id". Order index is assigned by
    sorted start position per contig, so it survives an unsorted GFF.
    """
    rows = []  # (contig, start, end, strand, {ids})
    with _open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in feature_types:
                continue
            try:
                contig, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
            except ValueError:
                continue
            ids = set()
            for m in re.finditer(r"(\w+)=([^;]+)", attr):
                key, val = m.group(1), m.group(2)
                if key in _ATTR_KEYS:
                    ids.add(val)
                    # NCBI prefixes: cds-XXX, gene-XXX -> also index the bare form
                    if "-" in val:
                        ids.add(val.split("-", 1)[1])
            rows.append((contig, start, end, strand, ids))

    # order index per contig, by start
    by_contig = {}
    for r in rows:
        by_contig.setdefault(r[0], []).append(r)
    out = {}
    for contig, recs in by_contig.items():
        recs.sort(key=lambda r: r[1])
        for i, (c, s, e, strand, ids) in enumerate(recs):
            for pid in ids:
                out[pid] = (c, s, e, strand, i)
    return out


def coords_from_prodigal_faa(path):
    """protein_id -> (contig, start, end, strand, order) from Prodigal .faa headers.

    Used for the Prodigal-called archaea instead of the GFF. The real Prodigal C
    binary writes GFF feature IDs as `{seqnum}_{genenum}` (e.g. 1_1) while its
    .faa header is `{contig}_{genenum}` (e.g. contigX_1), so a GFF join would
    silently fail. The .faa header carries the coordinates directly and its id is
    exactly the protein id the marker hit uses, so it cannot mismatch.

    Header: `>contigX_7 # 1049 # 1174 # 1 # ID=1_7;...`  (strand 1 or -1)
    """
    rows = []
    with _open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:].rstrip("\n")
            parts = [p.strip() for p in head.split("#")]
            if len(parts) < 4:
                continue
            pid = parts[0].split()[0]
            try:
                start, end = int(parts[1]), int(parts[2])
                strand = "+" if parts[3].strip() in ("1", "+1", "+") else "-"
            except ValueError:
                continue
            contig = pid.rsplit("_", 1)[0]
            rows.append((contig, start, end, strand, pid))
    by_contig = {}
    for r in rows:
        by_contig.setdefault(r[0], []).append(r)
    out = {}
    for contig, recs in by_contig.items():
        recs.sort(key=lambda r: r[1])
        for i, (c, s, e, strand, pid) in enumerate(recs):
            out[pid] = (c, s, e, strand, i)
    return out


def load_coords(source, faa_path, gff_path):
    """Pick the right coordinate source: Prodigal .faa header for archaea, GFF for
    provided bacterial proteomes. Returns (coords, origin) or ({}, reason)."""
    import os
    if source == "prodigal" and faa_path and os.path.exists(faa_path):
        return coords_from_prodigal_faa(faa_path), "prodigal_faa_header"
    if gff_path and isinstance(gff_path, str) and os.path.exists(gff_path):
        return parse_gff(gff_path), "gff"
    return {}, ("no coordinates: "
                + ("prodigal .faa missing" if source == "prodigal"
                   else f"gff path {gff_path!r} missing or not set"))


def block_synteny(hits_by_role, coords, n_contigs=None,
                  window_genes=DEFAULT_WINDOW_GENES, window_bp=DEFAULT_WINDOW_BP,
                  max_contigs=DEFAULT_MAX_CONTIGS_FOR_SYNTENY):
    """Are the block members a syntenic cluster?

    hits_by_role: {"muramyl_ligase": {pid, ...}, "mray_like": {...}, "cps": {...}}
    coords:       protein_id -> (contig, start, end, strand, order) from parse_gff
    n_contigs:    contig count for the genome (for the fragmentation guard)

    Returns a dict: status in {syntenic, dispersed, not_evaluable}, plus the
    contig, the gene-index span, and which block members were placed.
    """
    placed, missing = {}, {}
    for role, pids in hits_by_role.items():
        placed[role] = {p: coords[p] for p in pids if p in coords}
        missing[role] = sorted(p for p in pids if p not in coords)

    # need at least one of each block role to have coordinates
    have_all_roles = all(placed.get(r) for r in ("muramyl_ligase", "mray_like", "cps"))
    n_missing = sum(len(v) for v in missing.values())
    if not have_all_roles:
        return {
            "status": "not_evaluable",
            "reason": ("a block role has no placeable coordinate "
                       f"({n_missing} hit(s) not found in the GFF; the id join "
                       f"may be wrong, or the hit is not a CDS)"),
            "contig": None, "span_genes": None,
            "roles_placed": {r: len(placed.get(r, {})) for r in placed},
            "missing_ids": {r: v for r, v in missing.items() if v},
        }

    # For every contig that carries at least one ligase, one MraY-like and one
    # CPS, take the tightest window over one representative of each role.
    contigs = {}
    for role, d in placed.items():
        for pid, (c, s, e, strand, idx) in d.items():
            contigs.setdefault(c, {}).setdefault(role, []).append((idx, s, e, pid))

    best = None
    for c, byrole in contigs.items():
        if not all(byrole.get(r) for r in ("muramyl_ligase", "mray_like", "cps")):
            continue
        idxs = [min(v)[0] for v in byrole.values()] + [max(v)[0] for v in byrole.values()]
        span_g = max(idxs) - min(idxs)
        starts = [x[1] for v in byrole.values() for x in v]
        ends = [x[2] for v in byrole.values() for x in v]
        span_bp = max(ends) - min(starts)
        if best is None or span_g < best["span_genes"]:
            best = {"contig": c, "span_genes": span_g, "span_bp": span_bp}

    # Gene-order span is the criterion when we have it (always, from a GFF): it
    # is robust to intergenic distance. bp span is only a fallback for the
    # degenerate case of no usable order. `A or B` on both would let a huge bp gap
    # pass on a tiny gene span; that is wrong.
    within = (best is not None and (
        best["span_genes"] <= window_genes if best.get("span_genes") is not None
        else best["span_bp"] <= window_bp))
    if within:
        return {
            "status": "syntenic",
            "reason": (f"block on {best['contig']} within {best['span_genes']} "
                       f"genes / {best['span_bp']} bp"),
            "contig": best["contig"], "span_genes": best["span_genes"],
            "span_bp": best["span_bp"],
            "roles_placed": {r: len(placed[r]) for r in placed},
            "missing_ids": {},
        }

    # block members exist and are placed, but not together on one contig within
    # the window. On a fragmented assembly that is uninformative.
    if n_contigs is not None and n_contigs > max_contigs:
        return {
            "status": "not_evaluable",
            "reason": (f"block members are on different contigs, but the assembly "
                       f"has {n_contigs} contigs (> {max_contigs}); this is "
                       f"probably assembly breakage, not biology"),
            "contig": None, "span_genes": (best or {}).get("span_genes"),
            "roles_placed": {r: len(placed[r]) for r in placed},
            "missing_ids": {},
        }
    # Contiguity unknown (no contig count for this genome). We cannot tell a
    # dispersed block from an assembly artefact, so we must NOT call it 'dispersed'
    # -- that reads as 'present but not co-localized', a biological statement we
    # have not earned. Unknown contiguity is not_evaluable, same as an over-
    # fragmented assembly.
    if n_contigs is None:
        return {
            "status": "not_evaluable",
            "reason": ("block members are on different contigs, but the contig "
                       "count for this genome is unknown; cannot distinguish "
                       "dispersal from assembly breakage"),
            "contig": None, "span_genes": (best or {}).get("span_genes"),
            "roles_placed": {r: len(placed[r]) for r in placed},
            "missing_ids": {},
        }
    return {
        "status": "dispersed",
        "reason": ("block members are present but not co-localized on one contig "
                   "within the window"
                   + (f" (best same-contig span {best['span_genes']} genes)"
                      if best else " (never all on the same contig)")),
        "contig": (best or {}).get("contig"), "span_genes": (best or {}).get("span_genes"),
        "roles_placed": {r: len(placed[r]) for r in placed},
        "missing_ids": {},
    }
