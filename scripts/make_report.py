#!/usr/bin/env python3
"""Assemble the provenance report next to the figures."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, read_fasta  # noqa: E402


def tool_version(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        return (r.stdout + r.stderr).strip().splitlines()[0]
    except Exception:  # noqa: BLE001
        return "not found"


def yaml_safe_deps(path):
    """Return the pinned dependency strings from a conda env YAML."""
    try:
        import yaml
        y = yaml.safe_load(open(path)) or {}
        return [d for d in (y.get("dependencies") or []) if isinstance(d, str)]
    except Exception:  # noqa: BLE001
        return []


def load(path, **kw):
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except Exception:  # noqa: BLE001
        return pd.DataFrame()


def md_table(df, cols=None, floatfmt="{:.4g}"):
    if df.empty:
        return ["_(empty)_", ""]
    d = df[cols] if cols else df
    out = ["| " + " | ".join(map(str, d.columns)) + " |",
           "|" + "---|" * len(d.columns)]
    for _, r in d.iterrows():
        cells = [floatfmt.format(v) if isinstance(v, float) else f"{v:,}"
                 if isinstance(v, int) else str(v) for v in r]
        out.append("| " + " | ".join(cells) + " |")
    return out + [""]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tabdir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--chosen", required=True)
    ap.add_argument("--c71", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cfg = load_config(a.config)
    T = a.tabdir
    stats = load(os.path.join(T, "hmm_search_stats.tsv"), index_col="metric")
    tri = json.load(open(a.chosen))
    assign = load(os.path.join(T, "subgroup_assignments.tsv"))
    tiers = load(os.path.join(T, "triad_filter_by_tier.tsv"))
    fdr = load(os.path.join(T, "decoy_fdr.tsv"))
    conv = load(os.path.join(T, "convergence.tsv"))
    coup = load(os.path.join(T, "barcode_coupling.tsv"))
    coefs = load(os.path.join(T, "phyloglm_coefficients.tsv"))
    dstat = load(os.path.join(T, "phylogenetic_signal_D.tsv"))
    tips = load(os.path.join(T, "tree_tip_matching.tsv"))
    ssn_agree = load(os.path.join(T, "ssn_vs_subgroup_agreement.tsv"), index_col="metric")
    lysis = load(os.path.join(T, "lysis_reference_check.tsv"))
    try:
        gdef = json.load(open(os.path.join(T, "groove_definition.json")))
    except Exception:  # noqa: BLE001
        gdef = {}
    n_c71 = sum(1 for _ in read_fasta(a.c71))

    r1, r2, r3 = tri["residues"]
    i, j, k = tri["match_columns"]
    L = [f"# C71 / pseudomurein endo-isopeptidase screen\n",
         f"Generated {datetime.now():%Y-%m-%d %H:%M}\n"]

    # ---- funnel -------------------------------------------------------------
    # This funnel is the C71 arm only. combine_filter scopes these counts to the
    # C71-family profiles (PF12386 + SSF54001); PF03412/C39 hits are reported in
    # the "C39 arm" section, never pooled here.
    L.append("## C71 screening funnel\n")
    L.append("Counts below are the **C71 arm** (PF12386 + the shared SSF54001 fold "
             "net). PF03412 / C39 hits are reported separately in the C39 arm "
             "section and are not included here.\n")
    L.append("| Step | Count |")
    L.append("|---|---|")
    for key in ["samples_in_table", "samples_searched", "samples_skipped",
                "proteins_searched", "residues_searched", "unique_proteins_with_hit",
                "proteins_specific_evidence", "proteins_ssf_only",
                "proteins_hit_PF12386_and_SSF54001", "samples_with_hit"]:
        if key in stats.index:
            L.append(f"| {key.replace('_', ' ')} | {int(float(stats.loc[key, 'value'])):,} |")
    L.append(f"| triad-positive sequences | {tri['n_triad_positive']:,} |")
    L.append(f"| sequences in c71.faa | {n_c71:,} |\n")

    # c71.faa is a union of two evidence tiers, and after this point nothing
    # distinguishes them. Say so where it cannot be missed.
    c71ev = load(os.path.join(T, "c71_evidence.tsv"))
    if not c71ev.empty and "evidence" in c71ev.columns:
        comp = c71ev["evidence"].value_counts()
        n_ssf = int(comp.get("ssf_only", 0))
        L.append("`c71.faa` is the **union of two evidence tiers**: "
                 + ", ".join(f"{v:,} `{k}`" for k, v in comp.items()) + ".\n")
        if n_ssf:
            L.append(f"> {n_ssf / max(len(c71ev), 1):.1%} of `c71.faa` cleared only "
                     f"the SSF54001 fold model and then passed the triad filter. "
                     f"Their triad is real; their **family is not established**, "
                     f"because SCOP 54001 spans ~22 families related by insertion "
                     f"and circular permutation. Every analysis after this point — "
                     f"the tree, the subgroups, the SSN, the coupling, dN/dS, the "
                     f"prevalence model — treats them identically to PF12386 hits. "
                     f"The `evidence:ssf_only` row of `convergence.tsv` is the "
                     f"check: if those tips form a clade, `c71.faa` contains two "
                     f"families and the pooling has to be defended.\n")

    L.append("Thresholds applied: " + ", ".join(
        f"**{p}** = `{d['threshold']}` ({d['role']})" for p, d in cfg["profiles"].items()) + ".\n")

    # ---- FDR ----------------------------------------------------------------
    if not fdr.empty:
        L.append("## Empirical false-discovery rate\n")
        L.append("Reversed decoy sequences preserve length and composition but not "
                 "homology, so the rate at which they clear a bit score is the rate at "
                 "which chance clears it.\n")
        for p in sorted(fdr["profile"].unique()):
            d = fdr[fdr["profile"] == p].sort_values("bit_score")
            at = d.iloc[0]
            L.append(f"- **{p}** (threshold `{at['applied_threshold']}`): at the lowest "
                     f"reported score {at['bit_score']:.1f}, {int(at['n_decoy'])} decoy vs "
                     f"{int(at['n_target'])} target hits, FDR {at['fdr']:.3g}.")
        L.append("")

    # ---- triad --------------------------------------------------------------
    L.append("## Catalytic triad\n")
    L.append(f"Called **{r1}{i} - {r2}{j} - {r3}{k}** in PF12386 match-state coordinates "
             f"(source: `{tri['source']}`), learned from {tri['n_learning_sequences']:,} "
             f"`{tri['learned_from']}` sequences "
             f"(effective {tri['effective_n_learning']:.0f} after redundancy weighting), "
             f"then applied unchanged to every aligned sequence.\n")
    L.append("| Residue | Match column | Weighted frequency | Occupancy |")
    L.append("|---|---|---|---|")
    for res, col, f, o in zip(tri["residues"], tri["match_columns"],
                              tri["residue_frequencies"], tri["occupancies"]):
        L.append(f"| {res} | {col} | {f:.3f} | {o:.3f} |")
    L.append("")
    L.append(f"Spacing: {tri['gaps'][0]} and {tri['gaps'][1]} match columns "
             f"(prior {cfg['triad']['expected_gap_1_2']} and "
             f"{cfg['triad']['expected_gap_2_3']}).\n")

    if len(tri["hypotheses_scored"]) > 1:
        L.append("| Triad hypothesis | Score | Columns |")
        L.append("|---|---|---|")
        for h, d in sorted(tri["hypotheses_scored"].items(), key=lambda kv: -kv[1]["score"]):
            L.append(f"| {h} | {d['score']:.3f} | {d['columns']} |")
        L.append("")
        L.append("> PeiW (8JX4, C198/H233/D250) and PeiP (8Z4F, C213/H248/D272) both "
                 "use a Cys-His-Asp triad, and every alanine mutant is inactive "
                 "(Wang et al. 2025). CHD is the expected answer; CHN is scored only "
                 "as a control. If CHN wins, something is wrong with the alignment, "
                 "not with the biology.\n")

    if not tiers.empty:
        L.append("### Triad filter by evidence tier\n")
        L += md_table(tiers)
        s = tiers.set_index("evidence")
        col = ("frac_triad_positive_of_testable"
               if "frac_triad_positive_of_testable" in tiers.columns
               else "frac_triad_positive")
        if "ssf_only" in s.index and "specific" in s.index:
            rs, rf = float(s.loc["specific", col]), float(s.loc["ssf_only", col])
            ft = float(s.loc["ssf_only", "frac_testable"]) if "frac_testable" in \
                tiers.columns else 1.0
            L.append(f"> Of the SSF54001-only sequences, **{ft:.1%} could be tested "
                     f"at all**: the rest are gapped at the triad columns or below "
                     f"the match-state coverage floor, meaning hmmalign never placed "
                     f"them on the PF12386 scaffold. The test did not run on those; "
                     f"they are not negatives.\n")
            if "n_negative_with_all_three_residues" in tiers.columns:
                nu = int(s.loc["ssf_only", "n_negative_with_all_three_residues"])
                nn = int(s.loc["ssf_only", "n_triad_negative"])
                if nu:
                    L.append(f"> Of the {nn:,} sequences that *were* testable and "
                             f"came out negative, **{nu:,} carry a Cys, a His and an "
                             f"Asp somewhere in the aligned region** -- just not at "
                             f"the triad columns. A circularly permuted catalytic "
                             f"core looks exactly like this. Those negatives are not "
                             f"safe, and no column test can make them safe.\n")
            if ft < 0.5:
                L.append("> **The pass rate below cannot answer the "
                         "transglutaminase question.** SCOP 54001 spans ~22 "
                         "families related to the papain core by insertion and "
                         "*circular permutation*, including the transglutaminase "
                         "core. A permuted core has a Cys-His-Asp triad in 3D but "
                         "in a different sequential order, and this filter requires "
                         "C, H and D at fixed columns with `i < j < k`. It is "
                         "structurally incapable of detecting one. Use a "
                         "profile-profile map (HH-suite) or a superposition against "
                         "8JX4 instead.\n")
            else:
                L.append(f"> Among testable sequences, PF12386 hits carry the triad "
                         f"{rs:.1%} of the time and SSF54001-only hits {rf:.1%}. "
                         + ("The columns are diagnostic.\n" if rf < 0.5 * rs else
                            "These rates are close, which means the columns are not "
                            "discriminating fold from function. Investigate before "
                            "trusting c71.faa.\n"))

    # ---- C39 arm (PeiR-class), reported separately, never pooled ------------
    # The C39 arm searches PF03412 on its own scaffold with a spacing prior it
    # LEARNS rather than borrows from C71. Its numbers live in their own tables and
    # must never be added to the C71 counts above: a PF03412 hit is not a C71 hit.
    c39_chosen_path = os.path.join(T, "triad_columns_c39.json")
    if os.path.exists(c39_chosen_path):
        try:
            c39 = json.load(open(c39_chosen_path))
        except Exception:  # noqa: BLE001
            c39 = {}
        c39tiers = load(os.path.join(T, "triad_filter_by_tier_c39.tsv"))
        c39faa = os.path.join(a.outdir, "c39.faa")
        n_c39 = sum(1 for _ in read_fasta(c39faa)) if os.path.exists(c39faa) else 0
        L.append("## C39 arm (PeiR-class, PF03412)\n")
        L.append("PeiR (D3DZZ6) and the CRISPRTarget viral Peis are MEROPS C39 "
                 "(PF03412), not C71. They carry no PF12386 and cannot be aligned to "
                 "the PF12386 scaffold; PeiR's catalytic Cys->His gap is 72, where "
                 "C71's is 35. This arm therefore runs in parallel: its own PF03412 "
                 "alignment, its own triad columns, its own decoy FDR. **None of "
                 "these numbers are pooled with the C71 arm above.**\n")
        L.append("> MEROPS C39 proper is the bacteriocin-processing peptidase domain "
                 "of ABC exporters (double-glycine leader cleavage). Across 350,000 "
                 "proteomes most PF03412 hits are transporters, not Pei. The learned "
                 "triad and the per-family FDR below are the only things separating "
                 "signal from that background; treat every hit as a candidate to "
                 "confirm at the bench, not a called Pei.\n")
        if c39.get("match_columns") and c39.get("spacing_mode") == "learned":
            r1c, r2c, r3c = c39["residues"]
            ci, cj, ck = c39["match_columns"]
            lg = c39.get("learned_gaps") or [None, None]
            L.append(f"Learned triad: **{r1c}{ci} - {r2c}{cj} - {r3c}{ck}** in PF03412 "
                     f"match-state coordinates, from {c39.get('n_learning_sequences', 0):,} "
                     f"PF03412 sequences. The Cys->His / His->Asp spacing was **learned, "
                     f"not assumed**: {lg[0]} and {lg[1]} match columns. The C71 prior "
                     f"(gap 35) was NOT applied — applying it would reject a gap-72 "
                     f"enzyme like PeiR.\n")
            L.append(f"Sequences in `c39.faa`: {n_c39:,}.\n")
            if not c39tiers.empty:
                L.append("### C39 triad filter by evidence tier\n")
                L += md_table(c39tiers)
        else:
            L.append(f"> The C39 net returned no Pei-grade hits (empty specific "
                     f"tier). PF03412 was searched; nothing cleared it and carried the "
                     f"learned triad. `c39.faa` holds {n_c39:,} sequences. This is "
                     f"recorded, not hidden: an empty arm is a result about these "
                     f"proteomes, not a failure.\n")
        c39fdr = fdr[fdr["profile"] == "PF03412"] if not fdr.empty and "profile" in fdr.columns else pd.DataFrame()
        if not c39fdr.empty:
            d = c39fdr.sort_values("bit_score").iloc[0]
            L.append(f"> C39 decoy FDR (PF03412, own reversed decoys): at the lowest "
                     f"reported score {d['bit_score']:.1f}, {int(d['n_decoy'])} decoy "
                     f"vs {int(d['n_target'])} target hits, FDR {d['fdr']:.3g}. This is "
                     f"PF03412's alone.\n")
        L.append("> `pei_check` now certifies PeiR as findable (family c39 -> "
                 "PF03412). With this arm off, it reports PeiR as unfindable and the "
                 "screen is honest that a C71-only run cannot see it.\n")

    # ---- subgroups ----------------------------------------------------------
    if not assign.empty:
        sg = assign["subgroup"].value_counts().sort_index()
        L.append("## Active-site subgroups\n")
        L.append(f"k = {len(sg)}, chosen by `{cfg['active_site']['k_criterion']}`.\n")
        eff = assign.groupby("subgroup")["weight"].sum() if "weight" in assign.columns else None
        L.append("| Subgroup | n | effective n |")
        L.append("|---|---|---|")
        for s_, n_ in sg.items():
            e_ = f"{eff.loc[s_]:.1f}" if eff is not None else "-"
            L.append(f"| {s_} | {n_:,} | {e_} |")
        L.append("")

    if not ssn_agree.empty:
        L.append("### Sequence similarity network\n")
        L.append(f"{int(ssn_agree.loc['n_ssn_clusters', 'value'])} SSN clusters at "
                 f"alignment score {float(ssn_agree.loc['threshold', 'value']):g} vs "
                 f"{int(ssn_agree.loc['n_active_site_subgroups', 'value'])} active-site "
                 f"subgroups. Adjusted Rand "
                 f"{float(ssn_agree.loc['adjusted_rand', 'value']):.3f}, adjusted mutual "
                 f"information {float(ssn_agree.loc['adjusted_mutual_info', 'value']):.3f}.\n")
        L.append("> The SSN assumes only that homologues align. The k-means assumes "
                 "convex clusters in a BLOSUM embedding. Agreement between them is "
                 "evidence; disagreement is the result.\n")

    if not conv.empty:
        L.append("### Convergence\n")
        L += md_table(conv[["subgroup", "n_tips", "parsimony_changes", "null_random_mean",
                            "null_brownian_mean", "clustering_index", "p_clustered",
                            "interpretation"]])
        L.append("> `parsimony_changes` is the Fitch gains+losses count (an **upper "
                 "bound** on independent origins, not the origin count itself). "
                 "`clustering_index` is a Brownian-anchored 0/1 scaled statistic "
                 "(0 = as clumped as Brownian, 1 = random); it is **not** the "
                 "Fritz-Purvis *D* reported in `phylogenetic_signal_D.tsv`, which "
                 "is caper's sister-clade-difference *D*. The permutation "
                 "`p_clustered` is conditional on the single ML tree.\n")

    if not coup.empty:
        nsig = int(coup["significant"].sum())
        ncross = int((coup["significant"] & ~coup["same_block"]).sum())
        L.append("### Barcode coupling\n")
        L.append(f"{nsig} of {len(coup)} column pairs are coupled at "
                 f"q < {cfg['coupling']['fdr']} (APC-corrected MI, "
                 f"{cfg['coupling']['n_permutations']} permutations on cluster "
                 f"representatives). {ncross} span different catalytic residues, which is "
                 f"the evidence that the flanks form one surface rather than three.\n")
        if nsig:
            L += md_table(coup[coup["significant"]].head(10)[
                ["col_i", "col_j", "residue_i", "residue_j", "mi_apc", "q_bh"]])

    # ---- phylogenetic regression -------------------------------------------
    L.append("## Detection bias and phylogenetic non-independence\n")
    if not tips.empty:
        L.append("Tree tip matching:\n")
        L += md_table(tips[["domain", "n_samples", "n_matched", "frac_matched"]])
    if not dstat.empty:
        L.append("Fritz & Purvis *D* for C71 presence (0 = as clumped as Brownian, "
                 "1 = random):\n")
        L += md_table(dstat)
    if not coefs.empty:
        L.append("Phylogenetic logistic regression, `has_c71`. `phyloglm` rows model the "
                 "residual correlation induced by the GTDB tree; `glm_no_phylogeny` rows "
                 "are the same model fitted as if genomes were independent, shown so the "
                 "cost of that assumption is visible.\n")
        sub = coefs[coefs["response"] == "has_c71"]
        cols = [c for c in ["domain", "model", "term", "estimate", "std_error", "p",
                            "q_bh", "odds_ratio", "alpha", "n_tips", "n_positive"]
                if c in sub.columns]
        L += md_table(sub[cols])
        if "q_bh" in sub.columns:
            L.append("> `q_bh` is the Benjamini-Hochberg q-value across the "
                     "phyloglm hypothesis coefficients (subgroup/response presence "
                     "terms, both domains); the detection-bias covariates and the "
                     "non-phylogenetic comparison glm are excluded from that family "
                     "and carry raw `p` only.\n")
        L.append("> Coefficients on `completeness`, `log10_n50` and `log10_n_proteins` "
                 "quantify detection opportunity, not biology. A positive completeness "
                 "coefficient means the screen misses C71 in incomplete genomes, which "
                 "is why the prevalence table is reported at reference genome quality "
                 "rather than raw.\n")

    # ---- binding module -----------------------------------------------------
    archdf = load(os.path.join(T, "domain_architecture.tsv"))
    if not archdf.empty and "pmbr_binding_competent" in archdf.columns:
        n = len(archdf)
        n_ok = int(archdf["pmbr_binding_competent"].sum())
        n_frag = int(archdf.get("pmbr_count_fragile", pd.Series(0, index=archdf.index)).sum())
        L.append("## Binding module\n")
        L.append("Visweswaran et al. 2011 (doi:10.1371/journal.pone.0021582) fused "
                 "one, two or three PMB motifs to GFP. Three bind the pseudomurein "
                 "sacculus. Two bind lysozyme-treated bacterial spheroplasts and "
                 "**not** pseudomurein. One binds nothing.\n")
        L.append(f"| Binding call | Proteins |")
        L.append("|---|---|")
        for k, v in archdf["predicted_binding"].value_counts().items():
            L.append(f"| {k} | {v:,} |")
        L.append("")
        L.append(f"> {n_ok:,} of {n:,} proteins carry at least three motifs. The "
                 f"remainder are predicted unable to dock on an intact sacculus "
                 f"whatever their active site looks like, so a whole-cell lysis "
                 f"negative for one of them tests the module, not the triad.\n")
        if n_frag:
            L.append(f"> **{n_frag:,} proteins change binding class between the "
                     f"strict and permissive domain E-values.** A PMB motif is "
                     f"30-35 residues and the threshold is a cliff, so for these "
                     f"the architecture is not determined by the data. They are "
                     f"classed `pmbr_count_ambiguous` and dropped from the "
                     f"`pmbr_architecture` partition rather than assigned to "
                     f"whichever side the E-value picked.\n")
        L.append("> PMBR is not a pseudomurein marker. The domain binds NAG, the "
                 "one sugar shared by murein and pseudomurein, and it sticks to "
                 "*L. lactis* and *E. coli* spheroplasts. It is not Pei-specific "
                 "either: the S-layer protein MTH719 carries three motifs and no "
                 "catalytic domain. A bacterial C71+PMBR protein may be binding "
                 "exposed murein rather than being a binning artefact.\n")
        if "pmbr_pi" in archdf.columns and archdf["pmbr_pi"].notna().any():
            pis = pd.to_numeric(archdf["pmbr_pi"], errors="coerce").dropna()
            L.append(f"> PMB pI: median {pis.median():.1f}, range {pis.min():.1f} to "
                     f"{pis.max():.1f}. The characterised domain binds completely at "
                     f"pH 9.0 (pI 9.2), partially at 6.5, not at all at 4.0, and "
                     f"aggregates at 7.0. Every published Pei lysis assay runs at "
                     f"pH 7.0-7.85, which is where the module is least competent.\n")

    # ---- pseudomurein host genotype, and the discovery signal ---------------
    cw = load(os.path.join(T, "cellwall_genotype.tsv"))
    if not cw.empty and "pathway_call" in cw.columns:
        L.append("## Pseudomurein host genotype\n")
        L.append("Pseudomurein is restricted to two archaeal orders "
                 "(Methanobacteriales, Methanopyrales), so the call is "
                 "**taxonomy-first**: a genome in one of those orders is "
                 "`pseudomurein_expected_by_taxonomy` regardless of markers. The "
                 "marker block earns its keep only *outside* those orders "
                 "(Lupo et al. 2025, doi:10.1128/msystems.01201-24). Every call "
                 "below is a genomic hypothesis; confirmation is the wall "
                 "chemistry itself (TalNAc, the beta-1,3 linkage).\n")
        vc = cw["pathway_call"].value_counts()
        L.append("| Call | Genomes |")
        L.append("|---|---|")
        for k, v in vc.items():
            L.append(f"| {k} | {v:,} |")
        L.append("")

        # the headline: candidate PM OUTSIDE the two orders
        ooo = cw[cw["pathway_call"].astype(str).str.startswith(
            "pseudomurein_candidate_out_of_order")]
        if len(ooo):
            syn = ooo[ooo.get("synteny_status", "").astype(str) == "syntenic"] \
                if "synteny_status" in ooo.columns else ooo.iloc[0:0]
            L.append(f"### Candidate pseudomurein outside the two known orders "
                     f"({len(ooo)})\n")
            L.append("> These carry the PM-exclusive block (a muramyl ligase + the "
                     "MraY-like GT + the CPS) but are **not** in Methanobacteriales "
                     "or Methanopyrales. This is the discovery signal the screen "
                     "exists to surface. It is a hypothesis, not a finding: the "
                     "HMMs come from five genomes and miss divergent lineages, and "
                     "genomic presence is not wall chemistry.\n")
            if len(syn):
                L.append(f"**{len(syn)} are SYNTENIC** — the block co-localizes on "
                         f"one contig, the strongest genomic claim available. "
                         f"These are the ones to take to the bench first.\n")
                show = [c for c in ("sample", "species", "gtdb_order",
                                    "n_muramyl_ligases", "contigs", "synteny_detail")
                        if c in syn.columns]
                L += md_table(syn[show].head(30))
            n_disp = int((ooo.get("synteny_status", "").astype(str) == "dispersed").sum()) \
                if "synteny_status" in ooo.columns else 0
            n_unk = int((ooo.get("synteny_status", "").astype(str)
                         == "not_evaluable").sum()) if "synteny_status" in ooo.columns else 0
            if n_disp or n_unk:
                L.append(f"> Of the rest, {n_disp} are `dispersed` (block present but "
                         f"not co-localized — likely an HGT fragment) and {n_unk} are "
                         f"`synteny_unknown` (no usable coordinates, or an assembly "
                         f"too fragmented to judge — add a GFF via `inputs.gff_col`, "
                         f"or read this as 'not evaluated', never as 'not PM').\n")
        else:
            L.append("No candidate pseudomurein was found outside the two known "
                     "orders. Given how strongly PM tracks taxonomy, that is the "
                     "expected result, not a null one.\n")

    # ---- the one falsifiable check ------------------------------------------
    if not lysis.empty:
        L.append("## Validation against a measured phenotype\n")
        L.append("Schofield et al. 2015 (doi:10.1155/2015/828693) plated purified "
                 "PeiW and PeiP on eleven methanogens. This pipeline predicts "
                 "susceptibility from host wall chemistry (Kandler & Koenig 1978) "
                 "using a rule derived from the chromogenic substrate series. "
                 "Neither table saw the other, so these rows could have come out "
                 "wrong.\n")
        nt = lysis[lysis["nontrivial"].astype(str).str.lower() == "true"] \
            if "nontrivial" in lysis.columns else lysis
        testable = nt[nt["agrees"].notna()]
        n_ok = int(testable["agrees"].astype(str).str.lower().eq("true").sum())
        L += md_table(lysis[["strain", "enzyme", "wall", "p1_residue",
                             "observed", "predicted", "agrees"]])
        L.append(f"> {n_ok} of {len(testable)} falsifiable predictions agree. "
                 f"Rows whose host has no pseudomurein are excluded from that "
                 f"count: predicting that an enzyme cannot cut a wall which does "
                 f"not exist is not a test.\n")
        na = lysis[lysis["agrees"].isna()]
        if len(na):
            L.append(f"> {len(na)} rows are left unpredicted because no wall "
                     f"chemistry has been published for the host. *Methanobrevibacter* "
                     f"sp. SM9 is the interesting one: PeiW lyses it and PeiP does "
                     f"not, which is the only differential in the panel and the "
                     f"only place a specificity model could be tested against a "
                     f"phenotype it has never seen.\n")

    # ---- metal --------------------------------------------------------------
    met = (gdef or {}).get("metal") or {}
    if met:
        L.append("## Divalent metal site\n")
        L.append("Both characterised enzymes retain under 1% activity after EDTA. "
                 "Ca restores both; Mn, Mg, Ba and Ni restore PeiW but leave PeiP "
                 "under 15%. That is the sharpest measured difference between them, "
                 "and neither the substrate groove nor the four-class partition "
                 "predicts it.\n")
        if met.get("ions_found"):
            L += md_table(pd.DataFrame(met["ions_found"]))
            L.append(f"> {met['n_coordinating']} residues coordinate an ion and "
                     f"{met['n_in_shell']} lie in the shell; "
                     f"{met['shell_overlaps_groove']} of those are also in the "
                     f"substrate groove. The metal shell is tested for "
                     f"specificity-determining positions separately, so a hit in "
                     f"one region is not evidence for the other.\n")
        else:
            L.append("> **No cation is present in the deposited coordinates.** "
                     "The requirement is biochemical, not structural: no Pei "
                     "structure resolves the site. No metal shell was tested, and "
                     "none was invented. Locating this site is the most obvious "
                     "experiment this analysis cannot do.\n")

    # ---- provenance ---------------------------------------------------------
    # Seeds and pinned versions, so the run is reproducible from the report alone.
    L.append("## Provenance\n")
    seeds = {
        "tree.seed": (cfg.get("tree") or {}).get("seed"),
        "active_site.random_state": (cfg.get("active_site") or {}).get("random_state"),
        "convergence.seed": (cfg.get("convergence") or {}).get("seed"),
        "phyloglm.seed": (cfg.get("phyloglm") or {}).get("seed"),
    }
    L.append("Seeds: " + ", ".join(f"`{k}` = {v}" for k, v in seeds.items()
                                   if v is not None) + ".\n")
    L.append(f"Config: `{os.path.abspath(a.config)}`. "
             f"C39 arm (PF03412): "
             f"{'enabled' if (cfg['profiles'].get('PF03412') or {}).get('enabled') else 'disabled'}. "
             f"Its counts are reported only in the C39 section, never pooled into "
             f"the C71 funnel above.\n")

    # Pinned tool versions, read from the env YAMLs. The report runs in the `py`
    # env, which does not have hmmsearch/iqtree2/Rscript on PATH, so shelling out
    # to `--version` returns "not found" for every one of them (it always did).
    # The authoritative record is the pinned conda spec, and dist/locks/*.txt (if
    # present) holds the exact solve.
    L.append("### Pinned tool versions (from envs/*.yaml)\n```")
    envdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "envs")
    want = ("hmmer", "prodigal", "trimal", "iqtree", "mmseqs2", "diamond",
            "mafft", "hyphy", "foldseek", "hhsuite", "r-base", "python")
    seen = {}
    try:
        for ef in sorted(os.listdir(envdir)):
            if not ef.endswith(".yaml"):
                continue
            for dep in (yaml_safe_deps(os.path.join(envdir, ef))):
                name = str(dep).split("=")[0].strip()
                if name in want and name not in seen:
                    seen[name] = str(dep)
        for name in want:
            if name in seen:
                L.append(seen[name])
    except Exception as e:  # noqa: BLE001
        L.append(f"(could not read envs/: {e})")
    if not seen:
        L.append("(no env YAMLs found)")
    L.append("```")
    L.append("> Versions above are the pinned conda specs, not live `--version` "
             "probes: the report runs in the `py` environment, which does not "
             "carry the bioinformatics binaries. `dist/locks/*.txt` (written by "
             "setup/build_envs.sh) records the exact solved builds.\n")

    L.append("## Files\n")
    for root, _, files in os.walk(a.outdir):
        for f in sorted(files):
            if f.endswith((".svg", ".done")) or f == "report.md":
                continue
            L.append(f"- `{os.path.relpath(os.path.join(root, f), a.outdir)}`")

    with open(a.out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"[report] {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
