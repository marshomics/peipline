#!/usr/bin/env python3
"""The codon alignment must index residues, never consume them sequentially.

`triad_pass_matchcols.afa` holds only PF12386 match states. Insert residues and
anything `hmmalign --trim` clipped from the termini are gone. Two consequences,
and the old selection.py had both:

  1. `aligned[sid].replace("-","")` is NOT the protein, so the exact translated
     CDS match could never succeed. `stats["ok"]` was ~0 on real data and the
     module bailed with a message blaming "different annotation runs". The whole
     dN/dS analysis was a no-op that looked like a diagnosis.

  2. Threading a full CDS through a match-column row by taking the next codon at
     each non-gap column misassigns codons for any protein with an insertion:
     in frame, wrong position. The exact-match guard only "protected" this by
     rejecting every input.

This test builds a protein with a genuine insertion and trimmed termini, and
asserts that (a) the match-column -> residue-index map is right, (b) every codon
in the codon alignment translates to the residue in that column, and (c) the old
sequential-consumption approach would have produced a DIFFERENT, wrong answer.
If (c) ever stops being true the test is not testing anything.

Run:  python test/test_codon_mapping.py
"""
from __future__ import annotations

import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

FAILURES = []

# Reverse the standard code so we can build a CDS for any protein.
import selection as SEL  # noqa: E402

AA2CODON = {}
for cod, aa in SEL.GENCODE.items():
    AA2CODON.setdefault(aa, cod)


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def cds_for(prot):
    return "".join(AA2CODON[a] for a in prot)


def main() -> None:
    # Full protein. `head` and `tail` get trimmed away by hmmalign --trim.
    # The insertion sits between core[4] and core[5]. Match column 8 is a GAP in
    # this sequence, which means the protein genuinely LACKS core[8] -- a deletion
    # relative to the profile. Giving the protein that residue while gapping the
    # column would be an incoherent fixture, and the substring search says so.
    head, core, tail = "MKQ", "ACDEFGHIKL", "WYV"
    insertion = "ppp"                     # 3 insert residues, lowercase in Stockholm
    prot = head + core[:5] + insertion.upper() + core[5:8] + core[9] + tail
    assert len(prot) == 3 + 5 + 3 + 3 + 1 + 3

    # Stockholm: RF marks match states, '.' marks inserts, and hmmalign writes
    # insert residues in lowercase.
    seq_line = core[:5] + insertion + core[5:8] + "-" + core[9]
    rf_line = "x" * 5 + "..." + "x" * 5
    assert len(seq_line) == len(rf_line) == 13

    with tempfile.TemporaryDirectory() as td:
        sto = os.path.join(td, "hits.sto")
        with open(sto, "w") as fh:
            fh.write("# STOCKHOLM 1.0\n")
            fh.write(f"{'s1':<12}{seq_line}\n")
            fh.write(f"{'#=GC RF':<12}{rf_line}\n//\n")

        colmap, L = SEL.matchcol_residue_index(sto, {"s1": prot})

    check(L == 10, f"expected 10 match columns, got {L}")
    check("s1" in colmap, "the sequence was not mapped (substring search failed?)")
    idx = colmap["s1"]

    # Expected: match cols 0..4 -> residues 3..7 (after the 3-residue head).
    # Insert residues 8,9,10 occupy NO match column but DO consume the protein.
    # Match cols 5,6,7 -> residues 11,12,13. Col 8 is a gap. Col 9 -> residue 14.
    want = [3, 4, 5, 6, 7, 11, 12, 13, None, 14]
    check(idx == want, f"residue index map wrong:\n  got  {idx}\n  want {want}")
    print(f"1. residue map:      {idx}")
    print(f"   (head trimmed by 3, insertion of 3 skips match cols, col 8 is a gap)")

    # every mapped column must hold the residue the protein actually has there
    expect_res = [prot[i] if i is not None else "-" for i in idx]
    check("".join(expect_res) == core[:8] + "-" + core[9],
          f"columns do not carry the core residues: {''.join(expect_res)}")
    print(f"2. residues in cols: {''.join(expect_res)}")

    # --- codon alignment -----------------------------------------------------
    nt = cds_for(prot)
    ca = SEL.codon_align({"s1": idx}, {"s1": nt}, L)
    check("s1" in ca, "codon_align dropped the sequence")
    row = ca["s1"]
    check(len(row) == 3 * L, f"codon row is {len(row)} nt, expected {3 * L}")

    got = "".join("-" if row[3 * c:3 * c + 3] == "---"
                  else SEL.translate(row[3 * c:3 * c + 3]) for c in range(L))
    check(got == "".join(expect_res),
          f"codons translate to {got!r}, expected {''.join(expect_res)!r}")
    print(f"3. codons translate: {got}  (matches the protein, column for column)")

    # --- the old algorithm must be WRONG here, or this test proves nothing ----
    def old_sequential(aligned_matchcols, nt):
        out, k = [], 0
        for ch in aligned_matchcols:
            if ch == "-":
                out.append("---")
            else:
                out.append(nt[k:k + 3])
                k += 3
        return "".join(out)

    # what the old code fed it: the match-column row, ungapped == "protein"
    old_row = old_sequential("".join(expect_res), nt)
    old_aa = "".join("-" if old_row[3 * c:3 * c + 3] == "---"
                     else SEL.translate(old_row[3 * c:3 * c + 3]) for c in range(L))
    check(old_aa != got,
          "the old sequential-consumption algorithm produced the SAME answer, so "
          "this test does not demonstrate the bug it exists to prevent")
    print(f"4. old algorithm:    {old_aa}  <- wrong; it consumed the head and the "
          f"insertion as if they were match columns")

    # and the old exact-match target could never equal the protein
    old_target = "".join(expect_res).replace("-", "")
    check(old_target != prot,
          "the ungapped match-column row equals the protein, so the old exact "
          "match would have worked and there was no bug")
    print(f"5. old CDS target:   {old_target!r} != protein {prot!r}")
    print(f"   -> recover_cds() could never accept a single CDS. stats['ok'] = 0.")

    print()
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)})")
        sys.exit(1)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
