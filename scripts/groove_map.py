#!/usr/bin/env python3
"""Define the substrate groove from the reference structure and map it onto the
PF12386 match-state coordinate system that the rest of the pipeline uses.

Everything downstream that claims a residue is "in the groove" traces back to
this file. So it does three things and refuses to guess at any of them.

1. Reads the structure (mmCIF or PDB), extracts the requested chain's residues
   and their heavy-atom coordinates, and reconstructs the chain's one-letter
   sequence from the residue names. If the sequence has gaps in the auth
   numbering, they are reported, because a gap between the seed residues and
   the groove shell silently truncates the shell.

2. Takes the seed residues (the catalytic triad, plus Y174/V252/C265, which are
   PeiW numbering and belong to 8JX4) and collects every residue with a heavy
   atom within `groove_radius_a` of any seed heavy atom. This is a shell around
   the site, not a pocket-detection heuristic: it makes no claim about cavity
   volume, which is exactly the kind of number a predicted structure cannot
   support.

3. Aligns the structure's own sequence to PF12386 with hmmalign, so every
   structural residue lands on a profile match state or on an insert. Inserts
   cannot be compared across sequences and are dropped, with a count.

It also looks for a divalent cation. PeiW and PeiP are dead after EDTA and are
rescued by a metal -- PeiP only by Ca, PeiW also by Mn, Mg, Ba and Ni (Schofield et
al. 2015). That difference is the sharpest thing distinguishing the two
characterised enzymes, and the groove does not explain it. If a cation is in the
coordinates its first and second shells are emitted as separate columns, so
sdp.py can test them independently of the groove. If no cation is modelled, the
columns are zeros and the log says why. The site is not guessed at.

Output is a per-match-column table with `in_groove`, `in_metal_shell`,
`distance_to_seed_a`, and the structure's residue identity, plus the pairwise Ca
distance matrix over match columns, which coupling.py uses to test whether
coupled pairs are adjacent.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, match_columns, read_stockholm, write_fasta  # noqa: E402

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def parse_structure(path, chain, ion_names=()):
    """Return ({auth_seq_id: {name, atoms, ca}}, [ion records]).

    A minimal mmCIF/PDB reader. Biopython would do this, but its mmCIF parser
    pulls in a lot and silently renumbers; here the auth numbering is the whole
    point, because the literature residue numbers (Y174, V252, C265) are auth
    numbers.

    Ions are collected from EVERY chain, not just `chain`. A cation on a
    crystallographic interface still says the protein has a site, and both Pei
    enzymes are monomeric by gel filtration, so there is no biological interface
    for one to sit on.
    """
    res = {}
    ions = []
    ion_set = {s.upper() for s in ion_names}
    if path.endswith((".cif", ".cif.gz", ".mmcif")):
        cols, rows, in_loop = [], [], False
        with _open(path) as fh:
            for line in fh:
                s = line.strip()
                if s.startswith("_atom_site."):
                    cols.append(s.split(".", 1)[1])
                    in_loop = True
                    continue
                if in_loop:
                    if s.startswith("#") or not s:
                        break
                    rows.append(s.split())
        if not cols:
            sys.exit(f"[groove] no _atom_site loop in {path}")
        ix = {c: i for i, c in enumerate(cols)}
        need = ["group_PDB", "label_atom_id", "label_comp_id", "auth_asym_id",
                "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol"]
        for c in need:
            if c not in ix:
                sys.exit(f"[groove] mmCIF is missing _atom_site.{c}")
        for f in rows:
            if len(f) < len(cols) or f[ix["group_PDB"]] not in ("ATOM", "HETATM"):
                continue
            if f[ix["type_symbol"]] == "H":
                continue
            comp = f[ix["label_comp_id"]].upper()
            try:
                rid = int(f[ix["auth_seq_id"]])
                xyz = np.array([float(f[ix["Cartn_x"]]), float(f[ix["Cartn_y"]]),
                                float(f[ix["Cartn_z"]])])
            except ValueError:
                continue
            # A CA *atom* in an amino acid is the alpha carbon. A CA *residue*
            # is calcium. label_comp_id disambiguates; label_atom_id does not.
            if comp in ion_set and f[ix["group_PDB"]] == "HETATM":
                ions.append({"name": comp, "chain": f[ix["auth_asym_id"]],
                             "resnum": rid, "xyz": xyz})
                continue
            if f[ix["auth_asym_id"]] != chain or comp not in THREE_TO_ONE:
                continue
            e = res.setdefault(rid, {"name": comp, "atoms": [], "ca": None})
            e["atoms"].append(xyz)
            if f[ix["label_atom_id"]] == "CA":
                e["ca"] = xyz
    else:
        with _open(path) as fh:
            for line in fh:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                if line[76:78].strip() == "H":
                    continue
                comp = line[17:20].strip().upper()
                try:
                    rid = int(line[22:26])
                    xyz = np.array([float(line[30:38]), float(line[38:46]),
                                    float(line[46:54])])
                except ValueError:
                    continue
                if comp in ion_set and line.startswith("HETATM"):
                    ions.append({"name": comp, "chain": line[21],
                                 "resnum": rid, "xyz": xyz})
                    continue
                if line[21] != chain or comp not in THREE_TO_ONE:
                    continue
                e = res.setdefault(rid, {"name": comp, "atoms": [], "ca": None})
                e["atoms"].append(xyz)
                if line[12:16].strip() == "CA":
                    e["ca"] = xyz

    if not res:
        sys.exit(f"[groove] chain '{chain}' has no residues in {path}")
    for e in res.values():
        e["atoms"] = np.vstack(e["atoms"])
    return res, ions


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--align-hmm", required=True, help="PF12386 hmm (same one hmmalign used)")
    ap.add_argument("--chosen", required=True, help="triad_columns.json")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-columns", required=True)
    ap.add_argument("--out-contacts", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)["specificity"]
    os.makedirs(a.workdir, exist_ok=True)
    for p in (a.out_columns, a.out_contacts, a.out_json):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    path = cfg["structure"]
    if not os.path.exists(path):
        sys.exit(f"[groove] structure not found: {path}. This cluster has no "
                 f"internet; stage 8JX4 (PeiW-CD) and set "
                 f"specificity.structure. NOT 8Z4F: Y174/V252/C265 are "
                 f"PeiW numbering.")

    mcfg = cfg.get("metal_site") or {}
    res, ions = parse_structure(path, cfg["structure_chain"],
                                mcfg.get("ions", []) if mcfg.get("enabled") else ())
    rids = sorted(res)
    offset = int(cfg.get("structure_offset", 0))
    print(f"[groove] chain {cfg['structure_chain']}: {len(rids)} residues, "
          f"auth numbering {rids[0]}..{rids[-1]}", file=sys.stderr)

    gaps = [(rids[i], rids[i + 1]) for i in range(len(rids) - 1)
            if rids[i + 1] != rids[i] + 1]
    if gaps:
        print(f"[groove] {len(gaps)} chain breaks in the auth numbering, e.g. "
              f"{gaps[:3]}. Residues missing from the model cannot enter the "
              f"groove shell; the shell is therefore a lower bound.",
              file=sys.stderr)

    # --- groove shell --------------------------------------------------------
    seeds = [int(x) + offset for x in cfg["groove_seed_residues"]]
    with open(a.chosen) as fh:
        tri = json.load(fh)
    missing = [s for s in seeds if s not in res]
    if missing:
        sys.exit(f"[groove] seed residues {missing} are not in chain "
                 f"{cfg['structure_chain']}. Check specificity.structure_offset "
                 f"and that the numbering is auth, not label.")
    for s in seeds:
        print(f"[groove] seed {THREE_TO_ONE[res[s]['name']]}{s}", file=sys.stderr)

    # --- the guard that catches a structure/numbering mismatch --------------
    # Y174, V252 and C265 are PeiW numbering (8JX4). PeiP (8Z4F) has a different
    # triad (C213/H248/D272) and different residues at 252/265. Pointing at the
    # wrong file selects the wrong residues and every downstream analysis
    # inherits it. Check the identities, not just that the numbers exist.
    expect = {int(k): str(v).upper() for k, v in (cfg.get("structure_expect") or {}).items()}
    if expect:
        wrong = []
        for rid, want in sorted(expect.items()):
            rid_off = rid + offset
            got = THREE_TO_ONE[res[rid_off]["name"]] if rid_off in res else "-"
            if got != want:
                wrong.append(f"{rid}: expected {want}, found {got}")
        if wrong:
            sys.exit(
                "[groove] structure_expect mismatch:\n         "
                + "\n         ".join(wrong)
                + f"\n\n         {os.path.basename(path)} chain "
                  f"{cfg['structure_chain']} does not carry the residues this "
                  f"config names.\n"
                  f"         Wang et al. deposited PeiW-CD as 8JX4 (triad "
                  f"C198/H233/D250) and full-length PeiP as 8Z4F (triad "
                  f"C213/H248/D272).\n"
                  f"         Y174 / V252 / C265 are PeiW numbering. Use 8JX4, or "
                  f"renumber for PeiP.")
        ok_ids = ", ".join(f"{v}{k}" for k, v in sorted(expect.items()))
        print(f"[groove] structure_expect satisfied: {ok_ids} "
              f"({cfg.get('structure_numbering', '?')} numbering)", file=sys.stderr)

    R = float(cfg["groove_radius_a"])
    seed_atoms = np.vstack([res[s]["atoms"] for s in seeds])
    dist_to_seed = {}
    for rid in rids:
        d = np.sqrt(((res[rid]["atoms"][:, None, :] - seed_atoms[None, :, :]) ** 2)
                    .sum(-1)).min()
        dist_to_seed[rid] = float(d)
    groove = {rid for rid, d in dist_to_seed.items() if d <= R}
    print(f"[groove] {len(groove)} residues within {R} A of a seed heavy atom",
          file=sys.stderr)

    # --- divalent metal site -------------------------------------------------
    # Both enzymes are dead after EDTA and both are rescued by Ca, but only PeiW
    # is rescued by Mn/Mg/Ba/Ni (Schofield et al. 2015). If the coordinates hold a
    # cation, the residues around it are a second, independent hypothesis space
    # for specificity-determining positions -- one the groove analysis cannot see.
    #
    # If the coordinates hold no cation, that is reported as an absence, not
    # patched over with a prediction. A site inferred from sequence alone would
    # then propagate into the SDP enrichment as if it were structural.
    metal_coord, metal_shell, dist_to_metal = set(), set(), {}
    metal_records = []
    if mcfg.get("enabled"):
        r_coord = float(mcfg.get("coordination_radius_a", 3.2))
        r_shell = float(mcfg.get("shell_radius_a", 8.0))
        coordinating = set(mcfg.get("coordinating_residues") or [])
        if not ions:
            print(f"[metal] no {'/'.join(mcfg.get('ions', []))} ion in "
                  f"{os.path.basename(path)}. PeiW and PeiP require a divalent "
                  f"cation (<1% activity after EDTA; Schofield et al. 2015), and "
                  f"PeiP is Ca-specific while PeiW is not. No deposited Pei "
                  f"structure resolves the site. No metal shell will be tested; "
                  f"the columns are emitted as 0 so downstream code is uniform.",
                  file=sys.stderr)
        else:
            P = np.vstack([i["xyz"] for i in ions])
            for rid in rids:
                d = np.sqrt(((res[rid]["atoms"][:, None, :] - P[None, :, :]) ** 2)
                            .sum(-1)).min()
                dist_to_metal[rid] = float(d)
                if d <= r_coord:
                    metal_coord.add(rid)
                if d <= r_shell:
                    metal_shell.add(rid)
            for i in ions:
                near = [rid for rid in metal_coord
                        if np.sqrt(((res[rid]["atoms"] - i["xyz"]) ** 2)
                                   .sum(-1)).min() <= r_coord]
                aas = [THREE_TO_ONE[res[r]["name"]] for r in sorted(near)]
                n_ok = sum(1 for aa in aas if aa in coordinating)
                metal_records.append({
                    "ion": i["name"], "chain": i["chain"], "resnum": i["resnum"],
                    "n_coordinating_residues": len(near),
                    "coordinating": ",".join(f"{a}{r}" for a, r in
                                             zip(aas, sorted(near))),
                    "n_plausible_ligands": n_ok,
                })
                # An ion with no O/N donor within 3.2 A is a modelling artefact
                # or a cryoprotectant, and calling its neighbourhood a "site"
                # would launder noise into the SDP enrichment.
                flag = "" if n_ok >= 2 else "  <- fewer than 2 plausible ligands; " \
                                            "probably not a real site"
                print(f"[metal] {i['name']} {i['chain']}{i['resnum']}: "
                      f"{len(near)} residues within {r_coord} A "
                      f"({n_ok} of them O/N donors){flag}", file=sys.stderr)
            print(f"[metal] {len(metal_coord)} coordinating, {len(metal_shell)} "
                  f"in the {r_shell} A shell", file=sys.stderr)
            overlap = len(metal_shell & groove)
            print(f"[metal] {overlap} residues are in both the metal shell and "
                  f"the substrate groove; the two hypotheses are "
                  f"{'not ' if overlap > 0.5 * len(metal_shell) else ''}"
                  f"independent", file=sys.stderr)

    # --- map structure residues to PF12386 match states ----------------------
    seq = "".join(THREE_TO_ONE[res[r]["name"]] for r in rids)
    faa = os.path.join(a.workdir, "peip.faa")
    with open(faa, "w") as fh:
        write_fasta(fh, f"{cfg.get('structure_numbering', 'reference')}_structure", seq)
    sto = os.path.join(a.workdir, "peip.sto")
    subprocess.run(["hmmalign", "--trim", "--amino", "--outformat", "Stockholm",
                    "-o", sto, a.align_hmm, faa], check=True)

    order, seqs, rf = read_stockholm(sto)
    cols = match_columns(rf)
    aligned = seqs[order[0]]

    # `--trim` removes residues outside the profile envelope from the alignment.
    # They are removed from the ENDS, so the first aligned residue is not
    # necessarily the first residue of the chain. Indexing structure residues
    # from zero would shift the entire groove by the number of trimmed residues
    # and every downstream analysis would inherit the shift, silently.
    #
    # Locate the aligned fragment inside the chain sequence instead of assuming.
    ungapped = "".join(ch for ch in aligned if ch not in "-.").upper()
    off = seq.find(ungapped)
    if off < 0:
        sys.exit("[groove] the hmmalign output does not occur as a substring of the "
                 "chain sequence. Either the structure has non-standard residues "
                 "that were skipped, or hmmalign rewrote the sequence. Cannot map "
                 "structure residues to match states safely.")
    n_trim = len(rids) - len(ungapped)
    if n_trim:
        print(f"[groove] hmmalign --trim removed {n_trim} residues outside the "
              f"profile envelope; the aligned fragment starts at chain offset "
              f"{off} (auth {rids[off]})", file=sys.stderr)

    col_of_res, res_of_col = {}, {}
    k = off
    matchset = set(cols)
    colidx = {c: i for i, c in enumerate(cols)}
    for pos, ch in enumerate(aligned):
        if ch in "-.":
            continue
        rid = rids[k]
        k += 1
        if pos in matchset:
            col_of_res[rid] = colidx[pos]
            res_of_col[colidx[pos]] = rid
    assert k == off + len(ungapped)

    n_ins = len(groove) - sum(1 for r in groove if r in col_of_res)
    print(f"[groove] {len(col_of_res)} structure residues land on match states; "
          f"{n_ins} groove residues fall in insert columns and are dropped",
          file=sys.stderr)

    # sanity: the catalytic triad columns must be C/H/D in the structure
    ok = []
    for r, c in zip(tri["residues"], tri["match_columns"]):
        rid = res_of_col.get(c)
        aa = THREE_TO_ONE[res[rid]["name"]] if rid else "?"
        ok.append(aa == r)
        print(f"[groove] triad {r} at match col {c} -> "
              f"{cfg.get('structure_numbering', 'ref')} {aa}{rid}",
              file=sys.stderr)
    L = tri["n_match_columns"]
    cover = len(col_of_res) / L

    # A groove built on a bad alignment is worse than no groove: every SDP,
    # every coupling contact and every selection contrast downstream inherits it
    # silently. So this is fatal, not a warning.
    if not all(ok):
        sys.exit("[groove] the triad columns do not map to C/H/D in the structure. "
                 "Either the structure aligned badly to the profile, the chain is "
                 "wrong, structure_offset is wrong, or triad.override_columns is "
                 "wrong. Refusing to emit a groove that would corrupt every "
                 "downstream analysis.")
    if cover < 0.5:
        sys.exit(f"[groove] only {100 * cover:.1f}% of the {L} match columns have a "
                 f"structural residue. hmmalign --trim dropped most of the chain, "
                 f"which means the structure sequence and the profile disagree. "
                 f"Check the chain and the structure file.")
    print(f"[groove] {100 * cover:.1f}% of match columns have a structural residue",
          file=sys.stderr)

    lit = {int(x) + offset for x in (cfg.get("groove_literature_residues") or [])}
    rows = []
    for c in range(L):
        rid = res_of_col.get(c)
        rows.append({
            "match_col": c,
            "structure_resnum": rid if rid else pd.NA,
            "structure_residue": THREE_TO_ONE[res[rid]["name"]] if rid else pd.NA,
            "distance_to_seed_a": round(dist_to_seed[rid], 2) if rid else pd.NA,
            "in_groove": int(rid in groove) if rid else 0,
            "in_groove_literature": int(rid in lit) if rid else 0,
            "is_seed": int(rid in seeds) if rid else 0,
            "is_triad": int(c in tri["match_columns"]),
            "modelled": int(rid is not None),
            "metal_coordinating": int(rid in metal_coord) if rid else 0,
            "in_metal_shell": int(rid in metal_shell) if rid else 0,
            "distance_to_metal_a": (round(dist_to_metal[rid], 2)
                                    if rid and rid in dist_to_metal else pd.NA),
        })
    cols_df = pd.DataFrame(rows)
    # geometry OR the residues the paper names; the union is what downstream uses
    cols_df["in_groove_any"] = ((cols_df["in_groove"] == 1) |
                                (cols_df["in_groove_literature"] == 1)).astype(int)
    cols_df.to_csv(a.out_columns, sep="\t", index=False)
    print(f"[groove] {int(cols_df['in_groove'].sum())} of {L} columns in the "
          f"{R} A shell; {int(cols_df['in_groove_literature'].sum())} named in the "
          f"paper's conserved motif; {int(cols_df['in_groove_any'].sum())} in the union",
          file=sys.stderr)

    # class positions, resolved once here so pei_class.py never re-derives them
    pc = cfg.get("pei_class") or {}
    class_cols = {}
    if pc.get("enabled"):
        for key in ("position_1", "position_2"):
            rid = int(pc[key]) + offset
            col = col_of_res.get(rid)
            if col is None and key == "position_2" and pc.get("position_2_alt"):
                alt = int(pc["position_2_alt"]) + offset
                if col_of_res.get(alt) is not None:
                    print(f"[groove] class {key} residue {pc[key]} is not on a match "
                          f"state; falling back to the alternative numbering "
                          f"{pc['position_2_alt']} the paper's motif line uses",
                          file=sys.stderr)
                    rid, col = alt, col_of_res[alt]
            if col is None:
                sys.exit(f"[groove] class {key} (residue {pc[key]}) does not land on a "
                         f"profile match state; the four-class partition cannot be applied")
            class_cols[key] = {"auth": int(rid), "match_col": int(col),
                               "residue": THREE_TO_ONE[res[rid]["name"]]}
            print(f"[groove] class {key}: {class_cols[key]['residue']}{rid} "
                  f"-> match column {col}", file=sys.stderr)

    # --- Ca distance matrix over match columns -------------------------------
    ca = {c: res[r]["ca"] for c, r in res_of_col.items() if res[r]["ca"] is not None}
    cs = sorted(ca)
    D = np.full((L, L), np.nan)
    P = np.vstack([ca[c] for c in cs])
    M = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
    for i, ci in enumerate(cs):
        for j, cj in enumerate(cs):
            D[ci, cj] = M[i, j]
    np.savetxt(a.out_contacts, D, fmt="%.2f", delimiter="\t")

    with open(a.out_json, "w") as fh:
        json.dump({
            "structure": path,
            "chain": cfg["structure_chain"],
            "n_residues": len(rids),
            "chain_breaks": len(gaps),
            "seeds_auth": seeds,
            "groove_radius_a": R,
            "n_groove_residues": len(groove),
            "n_match_columns": L,
            "groove_columns": sorted(int(c) for c, r in res_of_col.items() if r in groove),
            "triad_maps_to_CHD": bool(all(ok)),
            "n_columns_with_structure": len(res_of_col),
            "match_column_coverage": round(cover, 4),
            "structure_numbering": cfg.get("structure_numbering"),
            "literature_groove_columns": sorted(
                int(c) for c, r in res_of_col.items() if r in lit),
            "class_positions": class_cols,
            "metal": {
                "enabled": bool(mcfg.get("enabled")),
                "ions_found": metal_records,
                "n_coordinating": len(metal_coord),
                "n_in_shell": len(metal_shell),
                "coordinating_columns": sorted(
                    int(c) for c, r in res_of_col.items() if r in metal_coord),
                "shell_columns": sorted(
                    int(c) for c, r in res_of_col.items() if r in metal_shell),
                "shell_overlaps_groove": len(metal_shell & groove),
                "note": ("PeiW and PeiP retain <1% activity after EDTA and are "
                         "rescued by a divalent cation; PeiP by Ca alone, PeiW "
                         "also by Mn/Mg/Ba/Ni (Schofield et al. 2015, "
                         "doi:10.1155/2015/828693). If ions_found is empty the "
                         "site is unresolved in this structure, not absent from "
                         "the enzyme."),
            },
        }, fh, indent=2)


if __name__ == "__main__":
    main()
