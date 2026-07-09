#!/usr/bin/env python3
"""Locate the catalytic-triad columns, then keep only sequences that carry them.

Coordinate system
-----------------
hmmalign emits `#=GC RF`, marking the PF12386 match states. Projecting every
sequence onto those columns gives a fixed-width matrix whose column index *is*
the profile position. Insert columns are dropped: they are not homologous across
sequences and cannot host a conserved triad.

Which sequences define the columns
----------------------------------
Only the PF12386 hits (`evidence == specific`). SSF54001 is a SCOP superfamily
model covering every papain- and transglutaminase-like fold; letting those
sequences vote on where the triad sits would let unrelated cysteine proteases
pull the columns off the Pei active site. The columns are then applied
*unchanged* to the SSF54001-only sequences, which is exactly the test we want:
does a protein that merely has a cysteine-protease fold also have the Pei triad
in the Pei positions?

The per-tier pass rates go to `triad_filter_by_tier.tsv`. If ssf_only sequences
pass at anything like the specific-tier rate, either SSF54001 is more specific
than advertised or the columns are not diagnostic. Both are worth knowing.

Weighting
---------
Residue frequencies are redundancy-weighted (1 / cluster size at 90% identity).
An unweighted frequency would tell you what has been sequenced, not what is
conserved.

Residue identity
----------------
The PeiP structure (PDB 8Z4F) shows a transglutaminase-like Cys-His-Asp triad,
so CHD is the hypothesis. CHN is scored alongside and written out, because a
heuristic you cannot audit is not a result.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from itertools import product

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, match_columns, read_stockholm, to_match_string  # noqa: E402

GAP = ord("-")
AA20 = "ACDEFGHIKLMNPQRSTVWY"


def load_matrix(sto_path: str):
    order, seqs, rf = read_stockholm(sto_path)
    if not rf:
        sys.exit("[triad] no '#=GC RF' line in the Stockholm file.")
    cols = match_columns(rf)
    if not cols:
        sys.exit("[triad] RF line contains no match states.")
    rows = [to_match_string(seqs[n], cols) for n in order]
    L = len(cols)
    if any(len(r) != L for r in rows):
        sys.exit("[triad] ragged alignment: Stockholm blocks did not concatenate cleanly.")
    arr = np.frombuffer("".join(rows).encode("ascii"), dtype=np.uint8).reshape(len(rows), L)
    return order, arr


def weighted_stats(arr: np.ndarray, w: np.ndarray):
    """Occupancy, per-residue frequency (over non-gap mass), Shannon entropy."""
    nongap = (arr != GAP)
    tot = w.sum()
    occ = (nongap * w[:, None]).sum(axis=0) / tot
    denom = np.maximum((nongap * w[:, None]).sum(axis=0), 1e-12)
    freq = {aa: ((arr == ord(aa)) * w[:, None]).sum(axis=0) / denom for aa in AA20}
    P = np.stack([freq[aa] for aa in AA20])
    ent = -np.nansum(np.where(P > 0, P * np.log2(np.where(P > 0, P, 1)), 0.0), axis=0)
    return occ, freq, ent


def candidates(freq, occ, res, min_occ, min_freq, cap):
    f = freq[res]
    idx = np.where((f >= min_freq) & (occ >= min_occ))[0]
    if idx.size == 0:
        return np.array([], dtype=int)
    return idx[np.argsort(-f[idx])][:cap]


def score_triples(freq, c1, c2, c3, r1, r2, r3, g12, g23, tol, w):
    out = []
    for i, j, k in product(c1.tolist(), sorted(c2.tolist()), sorted(c3.tolist())):
        if not (i < j < k):
            continue
        pen = (abs((j - i) - g12) + abs((k - j) - g23)) / float(tol)
        s = freq[r1][i] + freq[r2][j] + freq[r3][k] - w * pen
        out.append((float(s), int(i), int(j), int(k), float(pen)))
    out.sort(key=lambda t: -t[0])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sto", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--idmap", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--out-chosen", required=True)
    ap.add_argument("--out-keep", required=True)
    ap.add_argument("--out-afa", required=True)
    ap.add_argument("--out-colstats", required=True)
    ap.add_argument("--out-tiers", required=True)
    a = ap.parse_args()

    full = load_config(a.config)
    cfg = full["triad"]
    for p in (a.out_candidates, a.out_chosen, a.out_keep, a.out_afa, a.out_colstats,
              a.out_tiers):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    names, arr = load_matrix(a.sto)
    n, L = arr.shape
    name_idx = {nm: i for i, nm in enumerate(names)}

    # --- evidence tier and weight per aligned sequence -----------------------
    idmap = pd.read_csv(a.idmap, sep="\t", dtype=str).set_index("seq_id")
    wt = pd.read_csv(a.weights, sep="\t").set_index("seq_id")["weight"]

    evidence = pd.Series(index=names, dtype=object)
    evidence[:] = idmap.reindex(names)["evidence"].to_numpy()
    if evidence.isna().any():
        sys.exit(f"[triad] {int(evidence.isna().sum())} aligned sequences missing from idmap")

    w = wt.reindex(names).to_numpy(dtype=float)
    if np.isnan(w).any():
        sys.exit(f"[triad] {int(np.isnan(w).sum())} aligned sequences missing a weight")

    specific_mask = (evidence.to_numpy() == "specific")
    n_spec = int(specific_mask.sum())
    print(f"[triad] {n} aligned sequences x {L} match columns; "
          f"{n_spec} specific ({full['specific_profile']}), {n - n_spec} ssf_only",
          file=sys.stderr)
    if n_spec == 0:
        sys.exit(f"[triad] no sequence cleared {full['specific_profile']}. "
                 f"There is nothing to learn the triad columns from.")

    learn = specific_mask if cfg.get("restrict_to_specific", True) else np.ones(n, bool)
    occ, freq, ent = weighted_stats(arr[learn], w[learn])

    with open(a.out_colstats, "w") as fh:
        fh.write("match_col\toccupancy\tentropy_bits\t" +
                 "\t".join(f"freq_{aa}" for aa in AA20) + "\n")
        for j in range(L):
            fh.write(f"{j}\t{occ[j]:.4f}\t{ent[j]:.4f}\t" +
                     "\t".join(f"{freq[aa][j]:.4f}" for aa in AA20) + "\n")

    r1, r2 = cfg["first_residue"], cfg["second_residue"]
    thirds = [cfg["third_residue"]] + ([cfg["alt_third_residue"]]
                                       if cfg.get("alt_third_residue") else [])

    c1 = candidates(freq, occ, r1, cfg["min_occupancy"], cfg["min_residue_freq"],
                    cfg["max_candidates_per_residue"])
    c2 = candidates(freq, occ, r2, cfg["min_occupancy"], cfg["min_residue_freq"],
                    cfg["max_candidates_per_residue"])
    if c1.size == 0 or c2.size == 0:
        sys.exit(f"[triad] no candidate columns for {r1} or {r2} at "
                 f"min_residue_freq={cfg['min_residue_freq']}.")

    ranked = {}
    with open(a.out_candidates, "w") as fh:
        fh.write("rank\thypothesis\tscore\tcol_1\tcol_2\tcol_3\tres_1\tres_2\tres_3\t"
                 "freq_1\tfreq_2\tfreq_3\tocc_1\tocc_2\tocc_3\tgap_1_2\tgap_2_3\t"
                 "spacing_penalty\tn_specific_with_triad\tfrac_specific\n")
        for r3 in thirds:
            c3 = candidates(freq, occ, r3, cfg["min_occupancy"], cfg["min_residue_freq"],
                            cfg["max_candidates_per_residue"])
            if c3.size == 0:
                print(f"[triad] no candidate {r3} columns; skipping that hypothesis",
                      file=sys.stderr)
                continue
            trips = score_triples(freq, c1, c2, c3, r1, r2, r3,
                                  cfg["expected_gap_1_2"], cfg["expected_gap_2_3"],
                                  cfg["spacing_tolerance"], cfg["spacing_weight"])
            if not trips:
                continue
            ranked[f"{r1}{r2}{r3}"] = trips[0]
            sub = arr[learn]
            for rank, (s, i, j, k, pen) in enumerate(trips[:20], 1):
                hits = int(((sub[:, i] == ord(r1)) & (sub[:, j] == ord(r2)) &
                            (sub[:, k] == ord(r3))).sum())
                fh.write(f"{rank}\t{r1}{r2}{r3}\t{s:.4f}\t{i}\t{j}\t{k}\t{r1}\t{r2}\t{r3}\t"
                         f"{freq[r1][i]:.4f}\t{freq[r2][j]:.4f}\t{freq[r3][k]:.4f}\t"
                         f"{occ[i]:.4f}\t{occ[j]:.4f}\t{occ[k]:.4f}\t{j - i}\t{k - j}\t"
                         f"{pen:.4f}\t{hits}\t{hits / len(sub):.4f}\n")

    if not ranked:
        sys.exit("[triad] no triple satisfied the candidate criteria.")

    override = cfg.get("override_columns")
    if override:
        i, j, k = (int(x) for x in override)
        r3, source = cfg["third_residue"], "config_override"
    else:
        hyp, best = max(ranked.items(), key=lambda kv: kv[1][0])
        _, i, j, k, _ = best
        r3, source = hyp[2], "auto"
        for h, b in sorted(ranked.items(), key=lambda kv: -kv[1][0]):
            print(f"[triad] hypothesis {h}: score={b[0]:.4f} cols={b[1:4]}", file=sys.stderr)

    # --- apply to EVERY sequence, both tiers --------------------------------
    mask = (arr[:, i] == ord(r1)) & (arr[:, j] == ord(r2)) & (arr[:, k] == ord(r3))
    n_keep = int(mask.sum())
    if n_keep == 0:
        sys.exit(f"[triad] chosen columns {(i, j, k)} retain zero sequences.")

    tiers = []
    for tier in ("specific", "ssf_only"):
        m = evidence.to_numpy() == tier
        if not m.any():
            continue
        tiers.append({
            "evidence": tier,
            "n_aligned": int(m.sum()),
            "n_triad_positive": int((m & mask).sum()),
            "frac_triad_positive": float((m & mask).sum() / m.sum()),
            "effective_n_aligned": float(w[m].sum()),
            "effective_n_triad_positive": float(w[m & mask].sum()),
        })
    tdf = pd.DataFrame(tiers)
    tdf.to_csv(a.out_tiers, sep="\t", index=False)
    print(tdf.to_string(index=False), file=sys.stderr)

    chosen = {
        "source": source,
        "residues": [r1, r2, r3],
        "match_columns": [i, j, k],
        "learned_from": full["specific_profile"] if cfg.get("restrict_to_specific", True) else "all",
        "n_learning_sequences": int(learn.sum()),
        "effective_n_learning": float(w[learn].sum()),
        "residue_frequencies": [float(freq[r1][i]), float(freq[r2][j]), float(freq[r3][k])],
        "occupancies": [float(occ[i]), float(occ[j]), float(occ[k])],
        "gaps": [int(j - i), int(k - j)],
        "n_input_sequences": int(n),
        "n_triad_positive": n_keep,
        "frac_triad_positive": n_keep / n,
        "n_match_columns": int(L),
        "hypotheses_scored": {h: {"score": b[0], "columns": list(b[1:4])}
                              for h, b in ranked.items()},
        "by_tier": tiers,
    }
    with open(a.out_chosen, "w") as fh:
        json.dump(chosen, fh, indent=2)

    keep_names = [nm for nm, m in zip(names, mask.tolist()) if m]
    with open(a.out_keep, "w") as fh:
        fh.write("\n".join(keep_names) + "\n")

    with open(a.out_afa, "w") as fh:
        for nm, row in zip(keep_names, arr[mask]):
            fh.write(f">{nm}\n{row.tobytes().decode('ascii')}\n")

    print(f"[triad] {r1}{i}-{r2}{j}-{r3}{k} (match-column numbering); "
          f"{n_keep}/{n} sequences ({100 * n_keep / n:.1f}%) carry the full triad",
          file=sys.stderr)


if __name__ == "__main__":
    main()
