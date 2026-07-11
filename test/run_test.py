#!/usr/bin/env python3
"""Run every pipeline script against the synthetic dataset.

pyhmmer stands in for the hmmsearch/hmmalign binaries (it is the HMMER3 code);
test/stubs.py stands in for prodigal, seqkit, mmseqs, diamond and trimal.

Assertions
  1. the triad columns are recovered exactly, learned from PF12386 hits only
  2. no fold-family (ssf_only) sequence and no broken-triad sequence survives
  3. no intact Pei sequence is dropped
  4. the three planted active-site subgroups are recovered (ARI ~ 1)
  5. the archaeal genomes go through Prodigal and appear in the genome table
     with tree tips and QC covariates attached

phyloglm.R is NOT executed here: no R in this sandbox. Its inputs are checked.

    python test/make_testdata.py /tmp/c71_test
    python test/run_test.py     /tmp/c71_test
"""
from __future__ import annotations

import glob
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, os.path.join(ROOT, "test"))

FAILURES = []


def check(cond, msg):
    if not cond:
        print(f"!! {msg}")
        FAILURES.append(msg)
    return cond


def run(cmd, env=None):
    print("\n$ " + " ".join(map(str, cmd)), flush=True)
    r = subprocess.run(list(map(str, cmd)), env=env)
    if r.returncode != 0:
        sys.exit(f"FAILED: {' '.join(map(str, cmd))}")


# --- pyhmmer stand-ins ------------------------------------------------------
def _abc():
    import pyhmmer
    return pyhmmer.easel.Alphabet.amino()


def hmmsearch(hmm_path, faa, dom_out, threshold):
    import pyhmmer
    from pyhmmer.easel import SequenceFile
    from pyhmmer.plan7 import HMMFile

    with HMMFile(hmm_path) as hf:
        hmms = list(hf)
    if str(threshold).startswith("cut_"):
        ga = hmms[0].cutoffs.gathering
        assert ga is not None, f"{hmm_path} has no GA line but config asked for {threshold}"
        T = float(ga[0])
    else:
        T = float(threshold)

    os.makedirs(os.path.dirname(dom_out), exist_ok=True)
    if os.path.getsize(faa) == 0:
        open(dom_out, "w").write("#\n")
        return
    with SequenceFile(faa, digital=True, alphabet=_abc()) as sf:
        seqs = sf.read_block()
    with open(dom_out, "wb") as fh:
        first = True
        for hits in pyhmmer.hmmsearch(hmms, seqs, T=T, domT=T, incT=T, incdomT=T, cpus=1):
            hits.write(fh, format="domains", header=first)
            first = False


def hmmalign(hmm_path, faa, sto_out):
    import pyhmmer
    from pyhmmer.easel import SequenceFile
    from pyhmmer.plan7 import HMMFile
    with HMMFile(hmm_path) as hf:
        hmm = next(iter(hf))
    with SequenceFile(faa, digital=True, alphabet=_abc()) as sf:
        seqs = sf.read_block()
    msa = pyhmmer.hmmalign(hmm, seqs, trim=True)
    msa.name = b"hits"
    with open(sto_out, "wb") as fh:
        msa.write(fh, "stockholm")


def stub_tree(aln, out):
    """UPGMA stand-in for IQ-TREE, with fake aLRT/UFBoot labels."""
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist
    from utils import read_fasta
    names, seqs = zip(*read_fasta(aln))
    B = np.array([list(s) for s in seqs]).view(np.uint32).reshape(len(seqs), -1)
    Z = linkage(pdist(B, metric="hamming"), method="average")
    nodes = {i: n for i, n in enumerate(names)}
    n = len(names)
    heights = {i: 0.0 for i in range(n)}
    rng = np.random.default_rng(1)
    for r, (x, y, d, _) in enumerate(Z):
        x, y = int(x), int(y)
        s = rng.integers(50, 100)
        nodes[n + r] = (f"({nodes[x]}:{d/2 - heights[x]:.6f},"
                        f"{nodes[y]}:{d/2 - heights[y]:.6f}){s}.0/{s}")
        heights[n + r] = d / 2
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w").write(nodes[n + len(Z) - 1] + ";\n")


# ---------------------------------------------------------------------------
def main() -> None:
    base = sys.argv[1] if len(sys.argv) > 1 else "/tmp/c71_test"
    work, hmmout = os.path.join(base, "work"), os.path.join(base, "hmm_output")
    outdir = os.path.join(base, "outputs")
    tabdir, figdir = os.path.join(outdir, "tables"), os.path.join(outdir, "figures")
    prod = os.path.join(base, "prodigal")
    for d in (work, hmmout, outdir, tabdir, figdir, prod):
        os.makedirs(d, exist_ok=True)

    import stubs
    bindir = stubs.install(os.path.join(base, "bin"))
    env = dict(os.environ, PATH=f"{bindir}:{os.environ['PATH']}")

    # --- config pointed at the sandbox --------------------------------------
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml")))
    cfg["inputs"].update(
        sample_table=os.path.join(base, "faa_sample_table.tsv"),
        bacteria_metadata=os.path.join(base, "bacteria_qc.tsv"),
        bacteria_metadata_key=None,
        archaea_fna_glob=os.path.join(base, "gtdb_archaea", "seqs", "*.fna"),
        archaea_metadata=os.path.join(base, "ar53_metadata.tsv"),
        archaea_metadata_key="accession")
    cfg["trees"] = {"bacteria": os.path.join(base, "bacteria.tree"),
                    "archaea": os.path.join(base, "ar53.tree")}
    cfg["tip_matching"] = {"bacteria": "gtdb_species",
                           "archaea": "gtdb_species_rep",
                           "min_match_fraction": 0.50}
    cfg["profiles"]["PF12386"]["path"] = os.path.join(base, "hmms", "PF12386.hmm")
    cfg["profiles"]["SSF54001"]["path"] = os.path.join(base, "hmms", "SSF54001.hmm")
    # This end-to-end fixture exercises the C71 single arm and builds only PF12386 +
    # SSF54001 models. The C39 arm has its own test (test_c39_arm.py); pin PF03412
    # off here so combine_filter does not look for domtblout the fixture never made.
    if "PF03412" in cfg["profiles"]:
        cfg["profiles"]["PF03412"]["enabled"] = False
    cfg["outputs"].update(hmm_output_dir=hmmout, outdir=outdir, workdir=work,
                          prodigal_dir=prod,
                          combined_table=os.path.join(base, "hmm_output_combine.txt"),
                          merged_table=os.path.join(base, "sample_table_taxonomy.tsv"))
    cfg["batch_size"] = 25
    cfg["prodigal"]["batch_size"] = 10
    cfg["active_site"].update(k_max=6, hdbscan_min_cluster_size=5)
    cfg["tree"]["max_seqs"] = 100000
    cfg["coupling"]["n_permutations"] = 120
    cfg["convergence"].update(n_permutations=120, n_brownian=60)
    cfg["ssn"].update(threshold_min=5, threshold_max=120, threshold_step=5)
    cfg["redundancy"]["cluster_min_seq_id"] = 0.95
    tcfg = os.path.join(base, "config.test.yaml")
    yaml.safe_dump(cfg, open(tcfg, "w"))

    C = cfg["outputs"]["combined_table"]
    MG = cfg["outputs"]["merged_table"]

    # --- 1. prodigal --------------------------------------------------------
    fna = sorted(glob.glob(cfg["inputs"]["archaea_fna_glob"]))
    print(f"\n[test] {len(fna)} archaeal assemblies")
    for c in range(0, len(fna), 10):
        run([sys.executable, f"{SCRIPTS}/prodigal_run.py", "--fna", *fna[c:c + 10],
             "--outdir", prod, "--manifest", f"{work}/prodigal/chunk_{c//10:04d}.tsv",
             "--config", tcfg], env=env)
    manifests = sorted(glob.glob(f"{work}/prodigal/chunk_*.tsv"))

    # --- 2. profiles --------------------------------------------------------
    run([sys.executable, f"{SCRIPTS}/prepare_hmms.py",
         "--hmms", cfg["profiles"]["PF12386"]["path"], cfg["profiles"]["SSF54001"]["path"],
         "--labels", "PF12386,SSF54001", "--thresholds", "cut_ga,25",
         "--outdir", f"{work}/profiles", "--out-map", f"{work}/profiles/profile_map.tsv"])

    # --- 3. unified table ---------------------------------------------------
    run([sys.executable, f"{SCRIPTS}/build_sample_table.py", "--config", tcfg,
         "--manifests", *manifests, "--out-table", f"{work}/sample_table_unified.tsv",
         "--out-tips", f"{tabdir}/tree_tip_matching.tsv",
         "--batchdir", f"{work}/batches", "--batch-size", 25])

    batches = sorted(os.path.basename(p)[6:-4]
                     for p in glob.glob(f"{work}/batches/batch_*.tsv"))

    # --- 4-5. the fused search job: concat + decoys + hmmsearch on scratch ----
    # This runs the real scripts/search_batch.sh, with a pyhmmer-backed hmmsearch
    # on PATH. It is the code path the cluster will take.
    senv = dict(env, ENVS_ROOT="", ALLOW_MISSING="--allow-missing",
                TMPDIR=os.path.join(base, "scratch"))
    os.makedirs(senv["TMPDIR"], exist_ok=True)
    for b in batches:
        run(["bash", f"{SCRIPTS}/search_batch.sh", f"{work}/batches/batch_{b}.tsv", b,
             f"{work}/profiles", hmmout, f"{work}/map/batch_{b}.map.tsv.gz",
             "2", tcfg, "1"], env=senv)

    # scratch must be empty: the batch faa and its decoy never touch shared storage
    leftover = glob.glob(os.path.join(senv["TMPDIR"], "*"))
    check(not leftover, f"search_batch.sh left files in TMPDIR: {leftover[:3]}")
    check(not os.path.exists(f"{work}/faa"),
          "batch faa was written to the shared work dir; it should live in $TMPDIR")
    check(not os.path.exists(f"{work}/decoy"),
          "decoy faa was written to the shared work dir; it should live in $TMPDIR")

    # --- 6-7 ----------------------------------------------------------------
    run([sys.executable, f"{SCRIPTS}/combine_filter.py", "--hmm-dir", hmmout,
         "--map-dir", f"{work}/map", "--profile-map", f"{work}/profiles/profile_map.tsv",
         "--config", tcfg, "--out-combined", C,
         "--out-all", f"{tabdir}/hmm_hits_all_domains.tsv.gz",
         "--out-stats", f"{tabdir}/hmm_search_stats.tsv", "--out-fdr", f"{tabdir}/decoy_fdr.tsv"])
    run([sys.executable, f"{SCRIPTS}/merge_metadata.py", "--combined", C,
         "--table", f"{work}/sample_table_unified.tsv", "--sample-col", "sample", "--out", MG])
    run([sys.executable, f"{SCRIPTS}/extract_seqs.py", "--hits", C,
         "--out-faa", f"{work}/hits.faa", "--out-idmap", f"{work}/hits_idmap.tsv.gz",
         "--threads", 4])
    run([sys.executable, f"{SCRIPTS}/seq_weights.py", "--faa", f"{work}/hits.faa",
         "--config", tcfg, "--tmpdir", f"{work}/weights_tmp", "--threads", 2,
         "--out", f"{tabdir}/sequence_weights.tsv"], env=env)

    # --- 8-10 ---------------------------------------------------------------
    hmmalign(cfg["profiles"]["PF12386"]["path"], f"{work}/hits.faa", f"{work}/hits.sto")
    run([sys.executable, f"{SCRIPTS}/triad_detect_filter.py", "--sto", f"{work}/hits.sto",
         "--config", tcfg, "--combined", C, "--idmap", f"{work}/hits_idmap.tsv.gz",
         "--weights", f"{tabdir}/sequence_weights.tsv",
         "--out-candidates", f"{tabdir}/triad_candidates.tsv",
         "--out-chosen", f"{work}/triad_columns.json", "--out-keep", f"{work}/triad_pass_ids.txt",
         "--out-afa", f"{work}/triad_pass_matchcols.afa",
         "--out-colstats", f"{tabdir}/alignment_column_stats.tsv",
         "--out-tiers", f"{tabdir}/triad_filter_by_tier.tsv"])
    run([sys.executable, f"{SCRIPTS}/extract_seqs.py", "--keep-ids", f"{work}/triad_pass_ids.txt",
         "--idmap", f"{work}/hits_idmap.tsv.gz", "--out-faa", f"{outdir}/c71.faa", "--threads", 4])

    # --- assertions on the triad -------------------------------------------
    tri = json.load(open(f"{work}/triad_columns.json"))
    idmap = pd.read_csv(f"{work}/hits_idmap.tsv.gz", sep="\t")
    truth = pd.read_csv(f"{base}/truth.tsv", sep="\t")
    keep = {l.strip() for l in open(f"{work}/triad_pass_ids.txt")}

    # protein_id is the prodigal-assigned name for archaea, so match bacteria only
    bact = idmap[idmap["protein_id"].str.contains("_pei_|_fold_", regex=True)]
    bact = bact.merge(truth, on=["sample", "protein_id"], how="left")
    bact["kept"] = bact["seq_id"].isin(keep)

    print("\n" + "=" * 74)
    print(f"triad called at {tri['match_columns']} ({''.join(tri['residues'])}); "
          f"planted at [112, 147, 171]")
    print(f"learned from {tri['n_learning_sequences']} '{tri['learned_from']}' sequences")
    check(tri["match_columns"] == [112, 147, 171], "triad columns wrong")
    check(tri["learned_from"] == "PF12386", "triad not learned from PF12386 only")

    tiers = pd.read_csv(f"{tabdir}/triad_filter_by_tier.tsv", sep="\t")
    print("\n" + tiers.to_string(index=False))

    intact = bact[(bact["family"] == "pei") & (bact["broken_residue"] == "none")]
    broken = bact[(bact["family"] == "pei") & (bact["broken_residue"] != "none")]
    fold = bact[bact["family"] == "fold"]
    print(f"\nintact Pei kept   : {int(intact['kept'].sum())}/{len(intact)}")
    print(f"broken Pei kept   : {int(broken['kept'].sum())}/{len(broken)}  (want 0)")
    print(f"fold-family kept  : {int(fold['kept'].sum())}/{len(fold)}  (want ~0)")
    ev = bact.groupby("family")["evidence"].value_counts().unstack(fill_value=0)
    print("\nevidence tier by planted family:\n" + ev.to_string())

    check(int(intact["kept"].sum()) == len(intact), "intact Pei sequences were dropped")
    check(int(broken["kept"].sum()) == 0, "broken-triad sequences survived")
    check(int(fold["kept"].sum()) == 0, "fold-family sequences survived the triad filter")
    if "fold" in ev.index and "specific" in ev.columns:
        check(ev.loc["fold", "specific"] == 0, "fold family cleared PF12386's GA threshold")
    if "pei" in ev.index and "specific" in ev.columns:
        check(ev.loc["pei", "specific"] == len(bact[bact.family == "pei"]),
              "some Pei proteins failed PF12386's GA threshold")
    print("=" * 74)

    # --- 11. tree -----------------------------------------------------------
    run([sys.executable, f"{SCRIPTS}/tree_input.py", "--afa", f"{work}/triad_pass_matchcols.afa",
         "--out-aln", f"{work}/tree_input.aln", "--out-reps", f"{tabdir}/tree_representatives.tsv",
         "--config", tcfg, "--tmpdir", f"{work}/tree_prep", "--threads", 2], env=env)
    stub_tree(f"{work}/tree_input.aln", f"{outdir}/tree/c71.treefile")

    # --- 12-13. active site, coupling ---------------------------------------
    run([sys.executable, f"{SCRIPTS}/active_site_analysis.py",
         "--afa", f"{work}/triad_pass_matchcols.afa", "--chosen", f"{work}/triad_columns.json",
         "--merged", MG, "--idmap", f"{work}/hits_idmap.tsv.gz",
         "--weights", f"{tabdir}/sequence_weights.tsv", "--figdir", figdir,
         "--tabdir", tabdir, "--config", tcfg, "--out-assign", f"{tabdir}/subgroup_assignments.tsv"])
    run([sys.executable, f"{SCRIPTS}/coupling.py", "--assign", f"{tabdir}/subgroup_assignments.tsv",
         "--chosen", f"{work}/triad_columns.json", "--weights", f"{tabdir}/sequence_weights.tsv",
         "--config", tcfg, "--figdir", figdir, "--out", f"{tabdir}/barcode_coupling.tsv"])

    # --- 14. SSN ------------------------------------------------------------
    run(["bash", f"{SCRIPTS}/ssn_align.sh", f"{outdir}/c71.faa", f"{work}/ssn_edges.m8",
         f"{work}/ssn_nodes.faa", tcfg, f"{work}/ssn_tmp", "2"], env=env)
    run([sys.executable, f"{SCRIPTS}/ssn.py", "--m8", f"{work}/ssn_edges.m8",
         "--nodes", f"{work}/ssn_nodes.faa", "--assign", f"{tabdir}/subgroup_assignments.tsv",
         "--merged", MG, "--config", tcfg, "--out-clusters", f"{tabdir}/ssn_clusters.tsv",
         "--out-sweep", f"{tabdir}/ssn_threshold_sweep.tsv",
         "--out-graphml", f"{outdir}/ssn.graphml", "--figdir", figdir, "--tabdir", tabdir])

    # --- 15. convergence ----------------------------------------------------
    run([sys.executable, f"{SCRIPTS}/convergence.py", "--tree", f"{outdir}/tree/c71.treefile",
         "--assign", f"{tabdir}/subgroup_assignments.tsv",
         "--reps", f"{tabdir}/tree_representatives.tsv", "--config", tcfg,
         "--figdir", figdir, "--out", f"{tabdir}/convergence.tsv"])

    # --- 16. genome table (phyloglm input) ----------------------------------
    run([sys.executable, f"{SCRIPTS}/prep_genome_table.py",
         "--table", f"{work}/sample_table_unified.tsv", "--combined", C,
         "--assign", f"{tabdir}/subgroup_assignments.tsv", "--idmap", f"{work}/hits_idmap.tsv.gz",
         "--map-dir", f"{work}/map", "--config", tcfg,
         "--out", f"{tabdir}/genome_level_table.tsv"])

    gt = pd.read_csv(f"{tabdir}/genome_level_table.tsv", sep="\t")
    print("\n" + "=" * 74)
    print("genome-level table (phyloglm input):")
    print(gt.groupby("domain")[["has_hit", "has_specific", "has_c71"]].sum().to_string())
    for c in cfg["phyloglm"]["covariates"]:
        check(c in gt.columns, f"covariate {c} missing from genome table")
        check(gt[c].notna().mean() > 0.9, f"covariate {c} mostly null")
    check((gt["domain"] == "Archaea").sum() > 0, "no archaeal genomes reached the table")
    check(gt.loc[gt["domain"] == "Archaea", "tree_tip"].notna().all(),
          "archaeal genomes did not match ar53 tree tips")
    check(gt.loc[gt["domain"] == "Bacteria", "tree_tip"].notna().all(),
          "bacterial genomes did not match tree tips")
    check((gt["domain"] == "Archaea").sum() > 20,
          "the archaea hidden in the bacterial sample table were not detected")
    check(not (gt["tree_tip"] == "s__").any(),
          "genomes were matched to the bare 's__' tip")
    ntip = gt.loc[gt["domain"] == "Bacteria", "tree_tip"].nunique()
    check(ntip < (gt["domain"] == "Bacteria").sum(),
          "bacterial tree is not species-level (one tip per genome)")
    print(f"  bacterial: {(gt['domain']=='Bacteria').sum()} genomes -> {ntip} species tips")
    check(gt["has_c71"].sum() > 0, "no genome carries a C71")
    check(gt.loc[gt["domain"] == "Archaea", "has_c71"].sum() > 0,
          "Prodigal-called archaeal genomes yielded no C71 hits")
    print("phyloglm.R NOT executed: no R in this sandbox. Its inputs are validated above.")
    print("=" * 74)

    # stub the R outputs so the plotting/report steps can be exercised
    pd.DataFrame(columns=["term", "estimate", "std_error", "z", "p", "odds_ratio",
                          "model", "response", "domain", "n_tips", "n_positive", "alpha"]
                 ).to_csv(f"{tabdir}/phyloglm_coefficients.tsv", sep="\t", index=False)
    pd.DataFrame(columns=["domain", "trait", "D", "p_random", "p_brownian"]
                 ).to_csv(f"{tabdir}/phylogenetic_signal_D.tsv", sep="\t", index=False)
    prev = gt.groupby("taxon").agg(n=("has_c71", "size"), raw_prevalence=("has_c71", "mean"))
    prev = prev.reset_index()
    prev["adjusted_prevalence"] = prev["raw_prevalence"].clip(0.01, 0.99)
    prev["lo"] = (prev["adjusted_prevalence"] - 0.05).clip(0, 1)
    prev["hi"] = (prev["adjusted_prevalence"] + 0.05).clip(0, 1)
    prev["domain"] = "mixed"
    prev.to_csv(f"{tabdir}/bias_adjusted_prevalence.tsv", sep="\t", index=False)

    # --- 17-19. figures + report --------------------------------------------
    run([sys.executable, f"{SCRIPTS}/plots_overview.py", "--combined", C, "--merged", MG,
         "--search-stats", f"{tabdir}/hmm_search_stats.tsv", "--fdr", f"{tabdir}/decoy_fdr.tsv",
         "--colstats", f"{tabdir}/alignment_column_stats.tsv",
         "--chosen", f"{work}/triad_columns.json", "--tiers", f"{tabdir}/triad_filter_by_tier.tsv",
         "--prevalence", f"{tabdir}/bias_adjusted_prevalence.tsv",
         "--figdir", figdir, "--tabdir", tabdir, "--config", tcfg])
    run([sys.executable, f"{SCRIPTS}/plot_tree.py", "--tree", f"{outdir}/tree/c71.treefile",
         "--assign", f"{tabdir}/subgroup_assignments.tsv", "--merged", MG,
         "--figdir", figdir, "--config", tcfg])
    run([sys.executable, f"{SCRIPTS}/make_report.py", "--outdir", outdir, "--tabdir", tabdir,
         "--config", tcfg, "--chosen", f"{work}/triad_columns.json", "--c71", f"{outdir}/c71.faa",
         "--out", f"{outdir}/report.md"])

    # --- subgroup recovery ---------------------------------------------------
    from sklearn.metrics import adjusted_rand_score
    assign = pd.read_csv(f"{tabdir}/subgroup_assignments.tsv", sep="\t")
    m = assign.merge(truth[truth["family"] == "pei"], on=["sample", "protein_id"],
                     suffixes=("_pred", "_true"))
    if len(m) > 10:
        ct = pd.crosstab(m["subgroup_true"], m["subgroup_pred"])
        ari = adjusted_rand_score(m["subgroup_true"], m["subgroup_pred"])
        print("\nplanted vs recovered subgroup:\n" + ct.to_string())
        print(f"adjusted Rand {ari:.3f}   k_pred={ct.shape[1]} k_true={ct.shape[0]}")
        check(ari > 0.9, f"subgroup recovery poor (ARI {ari:.3f})")
        check(ct.shape[1] == ct.shape[0], f"recovered {ct.shape[1]} subgroups, planted {ct.shape[0]}")

    # The five planted subgroup-defining positions are, by construction, perfectly
    # determined by subgroup and therefore perfectly coupled to one another.
    # Barcode layout: blocks of 11 around C, H, D -> offsets -3,-2,+2 in block 0,
    # -1 in block 1, +1 in block 2. Centres are at 5, 16, 27.
    planted = {2, 3, 7, 15, 28}
    coup = pd.read_csv(f"{tabdir}/barcode_coupling.tsv", sep="\t")
    conv = pd.read_csv(f"{tabdir}/convergence.tsv", sep="\t")
    ssn = pd.read_csv(f"{tabdir}/ssn_clusters.tsv", sep="\t")

    coup["planted_pair"] = coup["col_i"].isin(planted) & coup["col_j"].isin(planted)
    n_planted = int(coup["planted_pair"].sum())
    sig_planted = int((coup["planted_pair"] & coup["significant"]).sum())
    sig_other = int((~coup["planted_pair"] & coup["significant"]).sum())
    print(f"\ncoupling: {int(coup['significant'].sum())}/{len(coup)} pairs significant")
    print(f"  planted co-varying pairs: {sig_planted}/{n_planted} recovered")
    print(f"  other pairs called significant: {sig_other}/{len(coup) - n_planted}")
    print(coup.nlargest(5, "z")[["col_i", "col_j", "mi_apc", "z", "q_bh",
                                 "planted_pair"]].to_string(index=False))
    check(sig_planted >= 0.8 * n_planted,
          f"coupling missed the planted co-varying positions ({sig_planted}/{n_planted})")
    check(sig_other <= 0.05 * (len(coup) - n_planted),
          f"coupling called too many unplanted pairs ({sig_other})")

    check(ssn["ssn_cluster_size"].max() >= 2, "SSN degenerated to all singletons")
    print(f"convergence:\n{conv[['subgroup','parsimony_changes','null_random_mean','clustering_index']].to_string(index=False)}")
    print(f"ssn: {ssn['ssn_cluster'].nunique()} clusters over {len(ssn)} nodes")

    nsvg = len([f for f in os.listdir(figdir) if f.endswith(".svg")])
    print(f"\nfigures: {nsvg} svg / {len([f for f in os.listdir(figdir) if f.endswith('.png')])} png")
    print("\nRESULT:", "PASS" if not FAILURES else f"FAIL ({len(FAILURES)})")
    for f in FAILURES:
        print("  -", f)
    sys.exit(0 if not FAILURES else 1)


if __name__ == "__main__":
    main()
