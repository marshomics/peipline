#!/usr/bin/env python3
"""The synteny test must separate a real cluster from scattered hits, and must
refuse to call a fragmented assembly.

Scenarios, all on a hand-built GFF:

  1. block members adjacent on one contig            -> syntenic
  2. block members present but spread across the      -> dispersed
     same contig beyond the window
  3. block members on DIFFERENT contigs, few contigs  -> dispersed
     (a near-complete genome: not co-localized is real)
  4. block members on different contigs, MANY contigs -> not_evaluable
     (a shredded MAG: absence of synteny is assembly breakage)
  5. a hit id absent from the GFF (bad join)          -> not_evaluable
  6. GFF join by locus_tag / NCBI cds- prefix works

Run:  python test/test_synteny.py
"""
from __future__ import annotations

import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

FAILURES = []


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def gff_line(contig, start, end, pid, key="ID", strand="+"):
    return f"{contig}\tprodigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{key}={pid}\n"


def main() -> None:
    import synteny as S

    # A GFF with:
    #   contig1: a tight cluster (lig, mray, cps within 3 genes) + filler
    #   contig1 also: a lone extra CPS far downstream (tests "tightest window")
    #   contig2: an isolated ligase (for the multi-contig scenarios)
    lines = []
    # contig1 tight block at gene order 5,6,7 (starts 5000,6000,7000)
    for i in range(1, 11):
        pid = {5: "lig1", 6: "mray1", 7: "cps1"}.get(i, f"c1_g{i}")
        lines.append(gff_line("contig1", i * 1000, i * 1000 + 800, pid))
    # contig2: a ligase far away
    lines.append(gff_line("contig2", 2000, 2800, "lig2"))
    # contig3: an mray + cps together (for a same-genome different-contig case)
    lines.append(gff_line("contig3", 3000, 3800, "mray2"))
    lines.append(gff_line("contig3", 3900, 4700, "cps2"))
    # a locus_tag-keyed feature and an NCBI cds- prefixed one
    lines.append(gff_line("contig1", 20000, 20800, "LOCUS_0042", key="locus_tag"))
    lines.append(gff_line("contig1", 21000, 21800, "cds-WP_999", key="ID"))

    with tempfile.TemporaryDirectory() as td:
        gff = os.path.join(td, "g.gff")
        open(gff, "w").write("##gff-version 3\n" + "".join(lines))
        coords = S.parse_gff(gff)

    # --- 0. parse + join ----------------------------------------------------
    check("lig1" in coords and coords["lig1"][0] == "contig1",
          "lig1 must parse onto contig1")
    check("LOCUS_0042" in coords, "locus_tag join failed")
    check("WP_999" in coords, "NCBI cds- prefix should also index the bare id")
    print(f"1. parse/join:       {len(coords)} features; ID, locus_tag and cds- "
          f"prefixes all joinable")

    # --- 1. syntenic --------------------------------------------------------
    hits = {"muramyl_ligase": {"lig1"}, "mray_like": {"mray1"}, "cps": {"cps1"}}
    r = S.block_synteny(hits, coords, n_contigs=3)
    check(r["status"] == "syntenic", f"tight block must be syntenic, got {r['status']}")
    check(r["contig"] == "contig1" and r["span_genes"] <= 3,
          f"span wrong: {r}")
    print(f"2. syntenic:         adjacent block on one contig -> syntenic "
          f"(span {r['span_genes']} genes)")

    # --- 2. dispersed on the same contig (beyond the window) ----------------
    hits2 = {"muramyl_ligase": {"lig1"}, "mray_like": {"mray1"}, "cps": {"cps2"}}
    # cps2 is on contig3; lig/mray on contig1 -> no single contig has all three
    r2 = S.block_synteny(hits2, coords, n_contigs=3)
    check(r2["status"] == "dispersed",
          f"block split across contigs on a 3-contig genome -> dispersed, got {r2['status']}")
    print(f"3. dispersed:        block split across contigs (few contigs) -> "
          f"dispersed")

    # --- 3. many contigs -> not_evaluable -----------------------------------
    r3 = S.block_synteny(hits2, coords, n_contigs=500)
    check(r3["status"] == "not_evaluable",
          f"split block on a 500-contig MAG must be not_evaluable, got {r3['status']}")
    check("breakage" in r3["reason"], "the reason must name assembly breakage")
    print(f"4. fragmented:       same split, 500 contigs -> not_evaluable "
          f"(assembly breakage, not biology)")

    # --- 4. missing id (bad join) -> not_evaluable --------------------------
    hits4 = {"muramyl_ligase": {"ghost"}, "mray_like": {"mray1"}, "cps": {"cps1"}}
    r4 = S.block_synteny(hits4, coords, n_contigs=3)
    check(r4["status"] == "not_evaluable",
          f"a hit absent from the GFF must be not_evaluable, got {r4['status']}")
    check(r4["missing_ids"].get("muramyl_ligase") == ["ghost"],
          "the unjoinable id must be reported")
    print(f"5. bad join:         a hit id absent from the GFF -> not_evaluable, "
          f"names the id")

    # --- 5. the window is honoured ------------------------------------------
    # put the three roles 20 genes apart on one contig; default window is 12
    lines_far = [gff_line("k", 1000, 1800, "L"), gff_line("k", 25000, 25800, "M"),
                 gff_line("k", 50000, 50800, "C")]
    with tempfile.TemporaryDirectory() as td:
        gff2 = os.path.join(td, "far.gff")
        open(gff2, "w").write("".join(lines_far))
        far = S.parse_gff(gff2)
    r5 = S.block_synteny({"muramyl_ligase": {"L"}, "mray_like": {"M"}, "cps": {"C"}},
                         far, n_contigs=1, window_genes=1, window_bp=1000)
    check(r5["status"] == "dispersed",
          f"three roles 2 genes apart but window=1 -> dispersed, got {r5['status']}")
    print(f"6. window:           roles beyond the window on one contig -> dispersed")

    # --- 6b. same contig, many filler genes between block members -----------
    # The `or window_bp` bug let a huge gene span pass on a small bp span. Build a
    # contig with 80 genes and the block at 5, 40, 78 (gene span 73).
    with tempfile.TemporaryDirectory() as td:
        gp = os.path.join(td, "spread.gff")
        mark = {5: "L", 40: "M", 78: "C"}
        open(gp, "w").write("".join(
            gff_line("one", i * 1000, i * 1000 + 800, mark.get(i, f"f{i}"))
            for i in range(1, 81)))
        spread = S.parse_gff(gp)
    rr = S.block_synteny({"muramyl_ligase": {"L"}, "mray_like": {"M"}, "cps": {"C"}},
                         spread, n_contigs=3, window_genes=12, window_bp=15000)
    check(rr["status"] == "dispersed" and rr["span_genes"] == 73,
          f"block spread 73 genes on one contig must be dispersed, got {rr}")
    print("6b. same-contig span: block 73 genes apart on one contig -> dispersed "
          "(gene order, not bp)")

    # --- 6. Prodigal .faa header coord source (archaea) ---------------------
    # Real Prodigal writes GFF ids that DON'T match its .faa ids, so archaea use
    # the header. Header: >contigX_7 # 1049 # 1174 # 1 # ID=1_7;...
    with tempfile.TemporaryDirectory() as td:
        faa = os.path.join(td, "p.faa")
        with open(faa, "w") as fh:
            fh.write(">ctg1_1 # 100 # 400 # 1 # ID=1_1;partial=00\nMAAA\n")
            fh.write(">ctg1_2 # 500 # 900 # -1 # ID=1_2;partial=00\nMBBB\n")
            fh.write(">ctg2_1 # 50 # 300 # 1 # ID=2_1;partial=00\nMCCC\n")
        pc = S.coords_from_prodigal_faa(faa)
    check(pc["ctg1_1"] == ("ctg1", 100, 400, "+", 0),
          f"prodigal header parse wrong: {pc.get('ctg1_1')}")
    check(pc["ctg1_2"][3] == "-", "strand -1 must parse to '-'")
    check(pc["ctg2_1"][0] == "ctg2", "contig must be the id minus the last _N")
    print("7. prodigal header: coords + strand parse; id == faa protein id")

    # --- 7. load_coords routes by source ------------------------------------
    with tempfile.TemporaryDirectory() as td:
        faa = os.path.join(td, "a.faa")
        open(faa, "w").write(">c_1 # 1 # 300 # 1 #\nM\n")
        gpath = os.path.join(td, "b.gff")
        open(gpath, "w").write(gff_line("c", 1, 300, "prot1"))
        c_arch, o_arch = S.load_coords("prodigal", faa, None)
        check(o_arch == "prodigal_faa_header" and "c_1" in c_arch,
              "an archaeal (prodigal) genome must use the .faa header")
        c_bact, o_bact = S.load_coords("provided_faa", None, gpath)
        check(o_bact == "gff" and "prot1" in c_bact,
              "a bacterial genome must use the GFF")
        c_none, o_none = S.load_coords("provided_faa", None, None)
        check(not c_none and "missing" in o_none,
              "no coordinate source -> empty + a reason")
    print("8. load_coords:     archaea->faa header, bacteria->gff, neither->reason")

    # --- 9. the divergent gate, through the REAL refine_synteny -------------
    # Permissive (below-GA) block hits may only elevate an out-of-order genome if
    # synteny holds. This is the false-positive control for divergent-lineage
    # mode: without it, permissive markers are noise.
    import pandas as pd
    from cellwall_genotype import annotate_pathway_calls, refine_synteny
    with tempfile.TemporaryDirectory() as td:
        def mkg(name, pos):
            p = os.path.join(td, name)
            open(p, "w").write("".join(
                gff_line("c1", i * 1000, i * 1000 + 800, pos.get(i, f"f{i}"))
                for i in range(1, 81)))
            return p
        tight = mkg("t.gff", {5: "lig", 6: "mray", 7: "cps"})
        spread = mkg("s.gff", {5: "lig", 40: "mray", 78: "cps"})
        BACT = "d__Bacteria;o__Lactobacillales"
        g = pd.DataFrame([
            # strict block, out of order, tight -> _syntenic
            {"sample": "S", "classification": BACT, "completeness": 100,
             "source": "provided_faa", "faa": None, "gff": tight, "contigs": 3,
             "pmur_lig": 1, "pmur_mray": 1, "pmur_cps": 1},
            # divergent (no strict block), tight permissive -> _divergent_syntenic
            {"sample": "D1", "classification": BACT, "completeness": 100,
             "source": "provided_faa", "faa": None, "gff": tight, "contigs": 3,
             "pmur_lig": 0, "pmur_mray": 0, "pmur_cps": 0},
            # divergent, dispersed permissive -> NOT elevated
            {"sample": "D2", "classification": BACT, "completeness": 100,
             "source": "provided_faa", "faa": None, "gff": spread, "contigs": 3,
             "pmur_lig": 0, "pmur_mray": 0, "pmur_cps": 0},
        ])
        g = annotate_pathway_calls(g, ["lig"], ["mray"], ["cps"], block_ok=True,
                                   min_markers=3, cls_col="classification")
        role = {"muramyl_ligase": {"lig"}, "mray_like": {"mray"}, "cps": {"cps"}}
        bh = {"S": role}
        bhp = {k: role for k in ("S", "D1", "D2")}
        g = refine_synteny(g, bh, bhp,
                           {"enabled": True, "window_genes": 12, "window_bp": 15000,
                            "max_contigs_for_synteny": 200})
    call = dict(zip(g["sample"], g["pathway_call"]))
    check(call["S"] == "pseudomurein_candidate_out_of_order_syntenic",
          f"strict + tight must be _syntenic, got {call['S']}")
    check(call["D1"] == "pseudomurein_candidate_out_of_order_divergent_syntenic",
          f"divergent + tight must be _divergent_syntenic, got {call['D1']}")
    check("candidate" not in call["D2"],
          f"divergent + dispersed must NOT be elevated (FP control), got {call['D2']}")
    print("9. divergent gate:   strict->syntenic; divergent+syntenic->elevated; "
          "divergent+dispersed->NOT elevated")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
