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

What the ssf_only pass rate does NOT tell you
---------------------------------------------
SCOP 54001 contains roughly 22 families that differ from the papain catalytic core
by insertion and by CIRCULAR PERMUTATION -- the transglutaminase core among them.
A permuted core presents its Cys, His and Asp in a different sequential order. This
filter selects columns with i < j < k in PF12386 coordinates and then demands C, H
and D at exactly those columns. It cannot see a permuted triad, and hmmalign will
not place a remote homologue's catalytic residues on PF12386 match states anyway.

So a sequence can fail here for three different reasons, and the tier table now
separates them: `triad_negative` (the test ran, the answer was no),
`gapped_at_triad`, and `low_coverage` (the test never ran). Only the first is a
result. Reading the other two as "this protein lacks the Pei active site" is an
alignment failure wearing a biological conclusion.

Weighting
---------
Residue frequencies are redundancy-weighted (1 / cluster size at 90% identity).
An unweighted frequency would tell you what has been sequenced, not what is
conserved.

Residue identity
----------------
Both solved structures use a Cys-His-Asp triad: PeiW-CD (8JX4) C198/H233/D250
and PeiP (8Z4F) C213/H248/D272, with every alanine mutant inactive (Wang et al.
2025). CHD is therefore the hypothesis, not CHN. CHN is scored alongside and
written out anyway, because a heuristic you cannot audit is not a result.

Note the spacing: the Cys->His gap is 35 in both, but the His->Asp gap is 24 in
PeiP and 17 in PeiW. The spacing prior in config is PeiP's; the tolerance covers
PeiW. Do not tighten it.
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
    # stable sort: ties in weighted frequency must break deterministically (by
    # column index), or which columns survive `cap` is platform-dependent and the
    # chosen catalytic triad -- a published result -- is not reproducible.
    return idx[np.argsort(-f[idx], kind="stable")][:cap]


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


def _write_empty(a, specific_profile, family, learn_spacing, why):
    """Write valid, empty outputs when a family arm's specific tier is empty.

    The C39 net (PF03412) can legitimately return zero Pei-grade hits across a set
    of proteomes. That must not abort the whole run, and it must not leave half-
    written tables that the report then misreads. Every output gets its header and
    nothing else; `chosen.json` records WHY, so the report can say 'the C39 arm ran
    and found nothing' rather than staying silent.
    """
    for p in (a.out_candidates, a.out_colstats, a.out_tiers, a.out_outcomes,
              a.out_chosen, a.out_keep, a.out_afa):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    with open(a.out_candidates, "w") as fh:
        fh.write("rank\thypothesis\tscore\tcol_1\tcol_2\tcol_3\tres_1\tres_2\tres_3\t"
                 "freq_1\tfreq_2\tfreq_3\tocc_1\tocc_2\tocc_3\tgap_1_2\tgap_2_3\t"
                 "spacing_penalty\tn_specific_with_triad\tfrac_specific\n")
    with open(a.out_colstats, "w") as fh:
        fh.write("match_col\toccupancy\tentropy_bits\t" +
                 "\t".join(f"freq_{aa}" for aa in AA20) + "\n")
    pd.DataFrame(columns=[
        "evidence", "n_aligned", "n_testable", "n_low_coverage", "n_gapped_at_triad",
        "n_triad_negative", "n_negative_with_all_three_residues", "n_triad_positive",
        "frac_triad_positive_of_testable", "frac_testable", "effective_n_aligned",
        "effective_n_triad_positive"]).to_csv(a.out_tiers, sep="\t", index=False)
    if a.out_outcomes:
        pd.DataFrame(columns=["seq_id", "evidence", "match_coverage", "outcome"]
                     ).to_csv(a.out_outcomes, sep="\t", index=False)
    with open(a.out_chosen, "w") as fh:
        json.dump({"source": "empty", "family": family,
                   "specific_profile": specific_profile,
                   "spacing_mode": "learned" if learn_spacing else "prior",
                   "match_columns": None, "n_triad_positive": 0,
                   "n_input_sequences": 0, "note": why}, fh, indent=2)
    open(a.out_keep, "w").close()
    open(a.out_afa, "w").close()


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
    ap.add_argument("--out-outcomes", default=None,
                    help="per-sequence outcome: why each sequence was kept or dropped")
    ap.add_argument("--family", default="c71",
                    help="which family arm this alignment belongs to. Selects the "
                         "specific profile that defines the tier and decides whether "
                         "the Cys->His/His->Asp spacing is a fixed prior (c71) or is "
                         "learned from the hits (c39, where one sequence is not a "
                         "prior). Default c71 reproduces the single-arm behaviour.")
    ap.add_argument("--allow-empty-specific", action="store_true",
                    help="write empty outputs and exit 0 if no sequence cleared the "
                         "family's specific profile, instead of failing. For the C39 "
                         "arm, whose PF03412 net may return nothing but must not abort "
                         "the whole run.")
    a = ap.parse_args()

    full = load_config(a.config)
    cfg = full["triad"]

    # --- which family arm, and where its spacing prior comes from -------------
    # C71 owns a fixed Cys->His/His->Asp prior (PeiW/PeiP; see the `triad:` block).
    # C39 owns none: PeiR is the only C39 Pei with an assigned catalytic residue
    # (C90), and its His was never assigned, so one sequence cannot seed a prior.
    # The C39 arm therefore ranks columns by residue frequency alone (i<j<k) and
    # REPORTS the spacing it found rather than scoring against a target. Borrowing
    # the C71 gap of 35 would reject PeiR, whose gap is 72; `pei_check` refuses to
    # start such a run, and this is the code path that keeps that refusal honest.
    fam = (full.get("families") or {}).get(a.family, {})
    specific_profile = fam.get("specific_profile") or full["specific_profile"]
    learn_spacing = bool(fam) and fam.get("expected_gap_1_2") is None
    if learn_spacing:
        g12_target, g23_target, spacing_w = 0, 0, 0.0   # no prior; frequency only
        spacing_tol = float(cfg["spacing_tolerance"])
    else:
        g12_target = cfg["expected_gap_1_2"]
        g23_target = cfg["expected_gap_2_3"]
        spacing_tol = float(cfg["spacing_tolerance"])
        spacing_w = cfg["spacing_weight"]
    for p in (a.out_candidates, a.out_chosen, a.out_keep, a.out_afa, a.out_colstats,
              a.out_tiers, a.out_outcomes):
        if p is None:
            continue
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    # An empty family net (PF03412 returning nothing) produces an empty or
    # header-only Stockholm. hmmalign cannot emit an RF line for zero sequences,
    # so load_matrix would abort. Catch it here when the arm is allowed to be empty.
    if a.allow_empty_specific:
        try:
            _order, _seqs, _rf = read_stockholm(a.sto)
        except Exception:
            _order, _rf = [], ""
        if not _order or not _rf or not match_columns(_rf):
            why = f"family={a.family}: alignment has no match states/sequences"
            _write_empty(a, specific_profile, a.family, learn_spacing, why)
            print(f"[triad] {why}; --allow-empty-specific: wrote empty outputs, "
                  f"exiting 0.", file=sys.stderr)
            return

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
    print(f"[triad] family={a.family}: {n} aligned sequences x {L} match columns; "
          f"{n_spec} specific ({specific_profile}), {n - n_spec} ssf_only; "
          f"spacing={'learned' if learn_spacing else 'prior'}", file=sys.stderr)
    if n_spec == 0:
        msg = (f"[triad] no sequence cleared {specific_profile}. There is nothing "
               f"to learn the triad columns from.")
        if a.allow_empty_specific:
            _write_empty(a, specific_profile, a.family, learn_spacing, msg)
            print(msg + " --allow-empty-specific: wrote empty outputs, exiting 0.",
                  file=sys.stderr)
            return
        sys.exit(msg)

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

    def _maybe_empty(why):
        """A near-empty family net (a handful of PF03412 hits with no clean CHD)
        must not abort the whole run when --allow-empty-specific is set. Returns
        True if it wrote empty outputs and the caller should return."""
        if a.allow_empty_specific:
            _write_empty(a, specific_profile, a.family, learn_spacing, why)
            print(f"[triad] {why}; --allow-empty-specific: wrote empty outputs, "
                  f"exiting 0.", file=sys.stderr)
            return True
        return False

    c1 = candidates(freq, occ, r1, cfg["min_occupancy"], cfg["min_residue_freq"],
                    cfg["max_candidates_per_residue"])
    c2 = candidates(freq, occ, r2, cfg["min_occupancy"], cfg["min_residue_freq"],
                    cfg["max_candidates_per_residue"])
    if c1.size == 0 or c2.size == 0:
        why = (f"family={a.family}: no candidate columns for {r1} or {r2} at "
               f"min_residue_freq={cfg['min_residue_freq']}")
        if _maybe_empty(why):
            return
        sys.exit(f"[triad] {why}.")

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
                                  g12_target, g23_target, spacing_tol, spacing_w)
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
        why = f"family={a.family}: no C/H/D triple satisfied the candidate criteria"
        if _maybe_empty(why):
            return
        sys.exit(f"[triad] {why}.")

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
    #
    # A sequence that fails this test has failed it for one of three reasons, and
    # they are not the same reason.
    #
    #   triad_negative   the sequence occupies the triad columns and the residues
    #                    are not C/H/D. A real negative. The test ran and the
    #                    answer was no.
    #
    #   gapped_at_triad  at least one triad column is a gap. The test never ran.
    #
    #   low_coverage     the sequence occupies fewer than `min_match_coverage` of
    #                    the profile's match states, so hmmalign could not place
    #                    it. The test never ran.
    #
    # Collapsing the last two into "negative" is how an alignment failure becomes
    # a biological conclusion. It matters most for the ssf_only tier: SCOP 54001
    # contains ~22 families that differ from the papain core by INSERTION and
    # CIRCULAR PERMUTATION. A permuted catalytic core has its Cys, His and Asp in
    # a different sequential order, and this filter demands i < j < k in PF12386
    # coordinates. It cannot see a permuted triad even when the 3D active site is
    # identical. Reporting that as "no triad" would be a false statement about
    # transglutaminase-like proteins.
    coverage = (arr != GAP).sum(axis=1) / float(L)
    min_cov = float(cfg.get("min_match_coverage", 0.5))
    gapped = ((arr[:, i] == GAP) | (arr[:, j] == GAP) | (arr[:, k] == GAP))
    low_cov = coverage < min_cov

    # A confident positive needs the triad residues at the three columns AND
    # enough of the profile occupied to trust the placement. A fragment that lands
    # C/H/D on three columns but covers < min_match_coverage of the match states
    # was not really placed by hmmalign; it is `low_coverage`, not a Pei. Gating
    # the call on ~low_cov (rather than asserting the situation impossible, which
    # it is NOT -- coverage counts all L columns, the triad is only 3) is what
    # keeps a fragment out of c71.faa and stops the run-ending assertion below from
    # firing on real data.
    triad_residues = ((arr[:, i] == ord(r1)) & (arr[:, j] == ord(r2)) &
                      (arr[:, k] == ord(r3)))
    mask = triad_residues & ~low_cov & ~gapped
    n_keep = int(mask.sum())
    if n_keep == 0:
        why = f"family={a.family}: chosen columns {(i, j, k)} retain zero sequences"
        if _maybe_empty(why):
            return
        sys.exit(f"[triad] {why}.")

    # The residual hole, stated rather than papered over. A permuted core that
    # happens to occupy all three columns with the wrong residues is scored
    # `triad_negative` and is indistinguishable, by any column test, from a
    # protein that simply lacks the site. What CAN be reported cheaply is whether
    # the sequence carries a cysteine, a histidine and an aspartate anywhere at
    # all: if it does not, no permutation could rescue it; if it does, the
    # negative is not safe and needs a superposition or a profile-profile map.
    has_c = (arr == ord(r1)).any(axis=1)
    has_h = (arr == ord(r2)).any(axis=1)
    has_d = (arr == ord(r3)).any(axis=1)
    has_all_three = has_c & has_h & has_d

    outcome = np.full(n, "triad_negative", dtype=object)
    outcome[low_cov] = "low_coverage"
    outcome[gapped & ~low_cov] = "gapped_at_triad"
    outcome[mask] = "triad_positive"
    if (mask & (low_cov | gapped)).any():
        sys.exit("[triad] a sequence is triad-positive and simultaneously gapped "
                 "or low-coverage. The coordinate system is inconsistent.")
    testable = mask | (outcome == "triad_negative")

    tiers = []
    for tier in ("specific", "ssf_only"):
        m = evidence.to_numpy() == tier
        if not m.any():
            continue
        n_test = int((m & testable).sum())
        neg = m & (outcome == "triad_negative")
        tiers.append({
            "evidence": tier,
            "n_aligned": int(m.sum()),
            "n_testable": n_test,
            "n_low_coverage": int((m & (outcome == "low_coverage")).sum()),
            "n_gapped_at_triad": int((m & (outcome == "gapped_at_triad")).sum()),
            "n_triad_negative": int(neg.sum()),
            # negatives that still carry a C, an H and a D somewhere. A circular
            # permutation of the catalytic core would look exactly like this, so
            # these negatives are not safe.
            "n_negative_with_all_three_residues": int((neg & has_all_three).sum()),
            "n_triad_positive": int((m & mask).sum()),
            # rate over the sequences the test could actually run on. The old
            # denominator was n_aligned, which silently counted unalignable
            # sequences as evidence that the columns discriminate.
            "frac_triad_positive_of_testable": (float((m & mask).sum() / n_test)
                                                if n_test else float("nan")),
            "frac_testable": float(n_test / m.sum()),
            "effective_n_aligned": float(w[m].sum()),
            "effective_n_triad_positive": float(w[m & mask].sum()),
        })
    tdf = pd.DataFrame(tiers)
    tdf.to_csv(a.out_tiers, sep="\t", index=False)
    print(tdf.to_string(index=False), file=sys.stderr)

    if a.out_outcomes:
        pd.DataFrame({"seq_id": names, "evidence": evidence.to_numpy(),
                      "match_coverage": np.round(coverage, 4),
                      "outcome": outcome,
                      f"has_{r1}_{r2}_{r3}_anywhere": has_all_three.astype(int),
                      }).to_csv(a.out_outcomes, sep="\t", index=False)

    # The warning that stops an alignment failure being read as biology.
    s = tdf.set_index("evidence")
    if "ssf_only" in s.index:
        n_unsafe = int(s.loc["ssf_only", "n_negative_with_all_three_residues"])
        n_neg = int(s.loc["ssf_only", "n_triad_negative"])
        if n_unsafe:
            print(f"\n[triad] {n_unsafe} of {n_neg} SSF54001-only `triad_negative` "
                  f"sequences carry a {r1}, an {r2} and a {r3} somewhere in the "
                  f"aligned region, just not at the triad columns. A circular "
                  f"permutation of the catalytic core looks exactly like that. "
                  f"Those negatives are not safe and this filter cannot make them "
                  f"safe: superpose them on the structure, or map the profiles.",
                  file=sys.stderr)
        ft = float(s.loc["ssf_only", "frac_testable"])
        if ft < 0.5:
            print(
                f"\n[triad] WARNING: only {100 * ft:.1f}% of the {int(s.loc['ssf_only', 'n_aligned'])} "
                f"SSF54001-only sequences could be tested at all. The rest are "
                f"gapped at the triad columns or below {min_cov:.0%} match-state "
                f"coverage.\n"
                f"[triad]   SCOP 54001 contains ~22 families related to the papain "
                f"core by insertion and CIRCULAR PERMUTATION. A permuted triad has "
                f"its Cys/His/Asp in a different sequential order; this filter "
                f"requires i<j<k in PF12386 coordinates and cannot see one.\n"
                f"[triad]   Do NOT read the ssf_only pass rate as 'transglutaminase-"
                f"like proteins lack the Pei active site'. It says the PF12386 "
                f"scaffold cannot address the question. Answering it needs a "
                f"profile-profile map (HH-suite) or a structural superposition, "
                f"not a column triple.", file=sys.stderr)

    chosen = {
        "source": source,
        "family": a.family,
        "specific_profile": specific_profile,
        # 'prior' -> scored against the C71 Cys->His/His->Asp targets; 'learned' ->
        # ranked by residue frequency alone and the spacing below is what was FOUND,
        # not what was assumed. The C39 arm is always 'learned': PeiR (C90, His
        # unassigned) is one sequence, and one sequence is not a prior.
        "spacing_mode": "learned" if learn_spacing else "prior",
        "spacing_prior_gaps": (None if learn_spacing else [g12_target, g23_target]),
        "learned_gaps": ([int(j - i), int(k - j)] if learn_spacing else None),
        "residues": [r1, r2, r3],
        "match_columns": [i, j, k],
        "learned_from": specific_profile if cfg.get("restrict_to_specific", True) else "all",
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
        "min_match_coverage": min_cov,
        "outcome_census": {k: int(v) for k, v in
                           pd.Series(outcome).value_counts().items()},
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
