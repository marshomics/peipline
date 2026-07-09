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
    n_c71 = sum(1 for _ in read_fasta(a.c71))

    r1, r2, r3 = tri["residues"]
    i, j, k = tri["match_columns"]
    L = [f"# C71 / pseudomurein endo-isopeptidase screen\n",
         f"Generated {datetime.now():%Y-%m-%d %H:%M}\n"]

    # ---- funnel -------------------------------------------------------------
    L.append("## Screening funnel\n")
    L.append("| Step | Count |")
    L.append("|---|---|")
    for key in ["samples_in_table", "samples_searched", "samples_skipped",
                "proteins_searched", "residues_searched", "unique_proteins_with_hit",
                "proteins_specific_evidence", "proteins_ssf_only",
                "proteins_hit_both_profiles", "samples_with_hit"]:
        if key in stats.index:
            L.append(f"| {key.replace('_', ' ')} | {int(float(stats.loc[key, 'value'])):,} |")
    L.append(f"| triad-positive sequences | {tri['n_triad_positive']:,} |")
    L.append(f"| sequences in c71.faa | {n_c71:,} |\n")

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
        L.append("> PeiP (PDB 8Z4F) is a transglutaminase-like isopeptidase with a "
                 "Cys-His-Asp triad, so CHD is the expected answer. CHN is scored only "
                 "as a control. If CHN wins, something is wrong with the alignment, not "
                 "with the biology.\n")

    if not tiers.empty:
        L.append("### Triad filter by evidence tier\n")
        L += md_table(tiers)
        if len(tiers) == 2:
            s = tiers.set_index("evidence")
            if "ssf_only" in s.index and "specific" in s.index:
                rs, rf = s.loc["specific", "frac_triad_positive"], s.loc["ssf_only", "frac_triad_positive"]
                L.append(f"> PF12386 hits carry the triad {rs:.1%} of the time; "
                         f"SSF54001-only hits {rf:.1%} of the time. "
                         + ("The columns are diagnostic.\n" if rf < 0.5 * rs else
                            "These rates are close, which means the columns are not "
                            "discriminating fold from function. Investigate before "
                            "trusting c71.faa.\n"))

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
        L += md_table(conv[["subgroup", "n_tips", "observed_origins", "null_random_mean",
                            "null_brownian_mean", "fritz_purvis_D", "p_clustered",
                            "interpretation"]])

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
                            "odds_ratio", "alpha", "n_tips", "n_positive"] if c in sub.columns]
        L += md_table(sub[cols])
        L.append("> Coefficients on `completeness`, `log10_n50` and `log10_n_proteins` "
                 "quantify detection opportunity, not biology. A positive completeness "
                 "coefficient means the screen misses C71 in incomplete genomes, which "
                 "is why the prevalence table is reported at reference genome quality "
                 "rather than raw.\n")

    # ---- versions & files ---------------------------------------------------
    L.append("## Versions\n```")
    for cmd in (["prodigal", "-v"], ["hmmsearch", "-h"], ["hmmalign", "-h"],
                ["trimal", "--version"], ["iqtree2", "--version"], ["mmseqs", "version"],
                ["diamond", "version"], ["Rscript", "--version"]):
        L.append(f"{cmd[0]}: {tool_version(cmd)}")
    L.append("```\n")

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
