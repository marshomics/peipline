# ============================================================================
# C71 / pseudomurein endo-isopeptidase Pei screening pipeline
#
#   prodigal (archaea) + existing faa (bacteria)
#     -> hmmsearch: PF12386 at --cut_ga, SSF54001 at score 25
#     -> combine, filter, left-merge onto metadata
#     -> redundancy weights -> hmmalign -> triad columns from PF12386 hits only,
#        then applied to the SSF54001-only sequences
#     -> c71.faa -> IQ-TREE
#     -> active-site subgroups, barcode coupling, SSN, convergence
#     -> phylogenetic logistic regression on the GTDB trees
#     -> figures
#
# Run: see run.sh
# ============================================================================
import glob as _glob
import math
import os
import sys

configfile: "config.yaml"

IN   = config["inputs"]
OUT  = config["outputs"]
PROF = config["profiles"]

HMM_DIR   = OUT["hmm_output_dir"]
COMBINED  = OUT["combined_table"]
MERGED    = OUT["merged_table"]
OUTDIR    = OUT["outdir"]
WORK      = OUT["workdir"]
PRODIGAL  = OUT["prodigal_dir"]

SCRIPTS = os.path.join(workflow.basedir, "scripts")
CFG     = os.path.join(workflow.basedir, "config.yaml")

FIGDIR = os.path.join(OUTDIR, "figures")
TABDIR = os.path.join(OUTDIR, "tables")
LOGDIR = os.path.join(WORK, "logs")


# --- environments -----------------------------------------------------------
# Two ways to get the tools, and the cluster forces the second.
#
#   envs_root: null   Snakemake solves and creates the conda envs itself from
#                     envs/*.yaml. Needs internet on whatever node runs the job.
#                     Fine on a laptop; impossible on an offline compute node.
#
#   envs_root: /path  Pre-built environments already exist at
#                     <envs_root>/{hmmer,py,phylo,network,prodigal,r}. Every
#                     shell command activates the one it needs. Snakemake's conda
#                     machinery is not used at all, so nothing is downloaded and
#                     no channel is contacted. See setup/ and README.
ENVS_ROOT = (config.get("envs_root") or "").rstrip("/")


def _act(name):
    if not ENVS_ROOT:
        return ""
    # `set +u` because conda's activate script reads unbound variables, and
    # snakemake runs shells with `set -euo pipefail`.
    return (f"set +u; source {ENVS_ROOT}/{name}/bin/activate; set -u; ")


ENV_HMMER, ENV_PY, ENV_PHYLO = _act("hmmer"), _act("py"), _act("phylo")
ENV_NET, ENV_PROD, ENV_R = _act("network"), _act("prodigal"), _act("r")

if ENVS_ROOT:
    for _n in ("hmmer", "py", "phylo", "network", "prodigal", "r"):
        _p = os.path.join(ENVS_ROOT, _n, "bin", "activate")
        if not os.path.exists(_p):
            sys.exit(f"[c71] envs_root is set but {_p} does not exist. "
                     f"Run setup/install_envs.sh first, or set envs_root: null.")
    sys.stderr.write(f"[c71] using pre-built environments under {ENVS_ROOT}\n")

PROFILE_NAMES = sorted(PROF)
SPECIFIC      = config["specific_profile"]
ALIGN_HMM     = os.path.join(WORK, "profiles", f"{config['align_profile']}.hmm")

# --- Prodigal chunks are knowable at DAG time: the .fna already exist --------
FNA = sorted(_glob.glob(IN["archaea_fna_glob"]))
if not FNA:
    sys.stderr.write(f"[c71] WARNING: no .fna matched {IN['archaea_fna_glob']}\n")
PCHUNK = max(1, math.ceil(len(FNA) / int(config["prodigal"]["batch_size"])))
PCHUNKS = [f"{i:04d}" for i in range(PCHUNK)]
sys.stderr.write(f"[c71] {len(FNA)} archaeal assemblies -> {PCHUNK} prodigal jobs\n")


def fna_chunk(wildcards):
    i = int(wildcards.c)
    n = int(config["prodigal"]["batch_size"])
    return FNA[i * n:(i + 1) * n]


# --- hmmsearch batches are NOT knowable until the unified table exists -------
def batch_ids():
    d = checkpoints.unified_table.get().output.batchdir
    return sorted(_glob.glob(os.path.join(d, "batch_*.tsv")))


def batches(wildcards):
    return [os.path.basename(p)[6:-4] for p in batch_ids()]


def all_domtbl(wildcards):
    return expand(os.path.join(HMM_DIR, "batch_{b}.{p}.domtblout"),
                  b=batches(wildcards), p=PROFILE_NAMES)


def all_maps(wildcards):
    return expand(os.path.join(WORK, "map", "batch_{b}.map.tsv.gz"), b=batches(wildcards))


def all_decoys(wildcards):
    if not config["decoy_fdr"]:
        return []
    return expand(os.path.join(HMM_DIR, "decoy", "batch_{b}.{p}.domtblout"),
                  b=batches(wildcards), p=PROFILE_NAMES)


# ============================================================================
rule all:
    input:
        COMBINED,
        MERGED,
        os.path.join(OUTDIR, "c71.faa"),
        os.path.join(OUTDIR, "tree", "c71.treefile"),
        os.path.join(TABDIR, "triad_candidates.tsv"),
        os.path.join(TABDIR, "subgroup_assignments.tsv"),
        os.path.join(TABDIR, "ssn_clusters.tsv"),
        os.path.join(TABDIR, "convergence.tsv"),
        os.path.join(TABDIR, "barcode_coupling.tsv"),
        os.path.join(TABDIR, "phyloglm_coefficients.tsv"),
        os.path.join(FIGDIR, ".overview.done"),
        os.path.join(FIGDIR, ".tree.done"),
        os.path.join(FIGDIR, ".activesite.done"),
        os.path.join(OUTDIR, "report.md"),


# ---------------------------------------------------------------------------
# 1. Prodigal on the archaeal assemblies.
# ---------------------------------------------------------------------------
rule prodigal:
    input:
        fna=fna_chunk,
    output:
        manifest=os.path.join(WORK, "prodigal", "chunk_{c}.tsv"),
    params:
        outdir=PRODIGAL, cfg=CFG,
    threads: 1
    conda: "envs/prodigal.yaml"
    log: os.path.join(LOGDIR, "prodigal", "chunk_{c}.log")
    shell:
        "{ENV_PROD}python {SCRIPTS}/prodigal_run.py --fna {input.fna} --outdir {params.outdir} "
        "--manifest {output.manifest} --config {params.cfg} &> {log}"


# ---------------------------------------------------------------------------
# 2. Normalise both HMMs, and refuse to proceed if PF12386 has no GA line.
# ---------------------------------------------------------------------------
rule prepare_hmms:
    input:
        hmms=[PROF[k]["path"] for k in PROFILE_NAMES],
    output:
        hmms=expand(os.path.join(WORK, "profiles", "{p}.hmm"), p=PROFILE_NAMES),
        map=os.path.join(WORK, "profiles", "profile_map.tsv"),
    params:
        labels=",".join(PROFILE_NAMES),
        thresholds=",".join(str(PROF[k]["threshold"]) for k in PROFILE_NAMES),
        outdir=os.path.join(WORK, "profiles"),
    conda: "envs/hmmer.yaml"
    log: os.path.join(LOGDIR, "prepare_hmms.log")
    shell:
        "{ENV_HMMER}python {SCRIPTS}/prepare_hmms.py --hmms {input.hmms} --labels {params.labels} "
        "--thresholds {params.thresholds} --outdir {params.outdir} "
        "--out-map {output.map} &> {log}"


# ---------------------------------------------------------------------------
# 3. Build one sample table across bacteria and archaea, attach QC metadata and
#    tree tip labels, and split it into hmmsearch batches. A checkpoint, because
#    the number of batches is not known until Prodigal has run.
# ---------------------------------------------------------------------------
checkpoint unified_table:
    input:
        manifests=expand(os.path.join(WORK, "prodigal", "chunk_{c}.tsv"), c=PCHUNKS),
        table=IN["sample_table"],
        bac_meta=IN["bacteria_metadata"],
        ar_meta=IN["archaea_metadata"],
        bac_tree=config["trees"]["bacteria"],
        ar_tree=config["trees"]["archaea"],
    output:
        table=os.path.join(WORK, "sample_table_unified.tsv"),
        tips=os.path.join(TABDIR, "tree_tip_matching.tsv"),
        batchdir=directory(os.path.join(WORK, "batches")),
    params:
        cfg=CFG, size=config["batch_size"],
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "unified_table.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/build_sample_table.py --config {params.cfg} "
        "--manifests {input.manifests} --out-table {output.table} "
        "--out-tips {output.tips} --batchdir {output.batchdir} "
        "--batch-size {params.size} &> {log}"


# ---------------------------------------------------------------------------
# 4+5. One job per batch: concatenate the proteomes, build the reversed decoys,
#      and run hmmsearch for every profile against both.
#
# The batch FASTA and its decoy are ~1 GB each and are read by nothing except
# the hmmsearch that immediately follows. Split across four rules they would be
# written to shared storage and read straight back: ~700 GB of pointless
# round-trip over ~365 batches, which on this cluster costs far more than the
# search itself (3x10^11 residues x 2 profiles x 2 is only tens of core-hours).
#
# Fused, the intermediates live in $TMPDIR on the execution host and only the
# domtblout files and the per-sample map reach shared storage.
#
# PF12386 uses its curated Pfam gathering threshold; SSF54001, a SCOP
# superfamily model with no GA line, uses the bit-score cutoff of 25.
# search_batch.sh reads both from config.
# ---------------------------------------------------------------------------
rule search_batch:
    input:
        batch=os.path.join(WORK, "batches", "batch_{b}.tsv"),
        hmms=expand(os.path.join(WORK, "profiles", "{p}.hmm"), p=PROFILE_NAMES),
    output:
        dom=expand(os.path.join(HMM_DIR, "batch_{{b}}.{p}.domtblout"), p=PROFILE_NAMES),
        tbl=expand(os.path.join(HMM_DIR, "batch_{{b}}.{p}.tblout"), p=PROFILE_NAMES),
        ddom=expand(os.path.join(HMM_DIR, "decoy", "batch_{{b}}.{p}.domtblout"),
                    p=PROFILE_NAMES) if config["decoy_fdr"] else [],
        dtbl=temp(expand(os.path.join(HMM_DIR, "decoy", "batch_{{b}}.{p}.tblout"),
                         p=PROFILE_NAMES)) if config["decoy_fdr"] else [],
        map=os.path.join(WORK, "map", "batch_{b}.map.tsv.gz"),
    params:
        profiles_dir=os.path.join(WORK, "profiles"),
        hmm_dir=HMM_DIR,
        cfg=CFG,
        decoy=1 if config["decoy_fdr"] else 0,
        envs_root=ENVS_ROOT,
        allow_missing="--allow-missing" if config["allow_missing_faa"] else "",
    threads: config["hmmsearch_cpu"]
    conda: "envs/hmmer.yaml"
    log: os.path.join(LOGDIR, "search_batch", "batch_{b}.log")
    shell:
        "ENVS_ROOT='{params.envs_root}' ALLOW_MISSING='{params.allow_missing}' "
        "bash {SCRIPTS}/search_batch.sh {input.batch} {wildcards.b} "
        "{params.profiles_dir} {params.hmm_dir} {output.map} {threads} "
        "{params.cfg} {params.decoy} &> {log}"


# ---------------------------------------------------------------------------
# 6. Combine, filter, tier the evidence, and estimate the decoy FDR.
# ---------------------------------------------------------------------------
rule combine_filter:
    input:
        dom=all_domtbl,
        map=all_maps,
        decoy=all_decoys,
        pmap=os.path.join(WORK, "profiles", "profile_map.tsv"),
    output:
        combined=COMBINED,
        all_domains=os.path.join(TABDIR, "hmm_hits_all_domains.tsv.gz"),
        stats=os.path.join(TABDIR, "hmm_search_stats.tsv"),
        fdr=os.path.join(TABDIR, "decoy_fdr.tsv"),
    params:
        hmm_dir=HMM_DIR, map_dir=os.path.join(WORK, "map"), cfg=CFG,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "combine_filter.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/combine_filter.py --hmm-dir {params.hmm_dir} "
        "--map-dir {params.map_dir} --profile-map {input.pmap} --config {params.cfg} "
        "--out-combined {output.combined} --out-all {output.all_domains} "
        "--out-stats {output.stats} --out-fdr {output.fdr} &> {log}"


rule merge_metadata:
    input:
        combined=COMBINED,
        table=os.path.join(WORK, "sample_table_unified.tsv"),
    output:
        merged=MERGED,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "merge_metadata.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/merge_metadata.py --combined {input.combined} "
        "--table {input.table} --sample-col sample --out {output.merged} &> {log}"


# ---------------------------------------------------------------------------
# 7. Recover full-length hit sequences, then compute redundancy weights before
#    anything counts residue frequencies.
# ---------------------------------------------------------------------------
rule extract_hits:
    input:
        combined=COMBINED,
    output:
        faa=os.path.join(WORK, "hits.faa"),
        idmap=os.path.join(WORK, "hits_idmap.tsv.gz"),
    threads: 8
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "extract_hits.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/extract_seqs.py --hits {input.combined} "
        "--out-faa {output.faa} --out-idmap {output.idmap} --threads {threads} &> {log}"


rule seq_weights:
    input:
        faa=os.path.join(WORK, "hits.faa"),
    output:
        weights=os.path.join(TABDIR, "sequence_weights.tsv"),
    params:
        cfg=CFG, tmp=os.path.join(WORK, "weights_tmp"),
    threads: 8
    conda: "envs/phylo.yaml"
    log: os.path.join(LOGDIR, "seq_weights.log")
    shell:
        "{ENV_PHYLO}python {SCRIPTS}/seq_weights.py --faa {input.faa} --config {params.cfg} "
        "--tmpdir {params.tmp} --threads {threads} --out {output.weights} &> {log}"


rule hmmalign:
    input:
        faa=os.path.join(WORK, "hits.faa"),
        hmm=ALIGN_HMM,
    output:
        sto=os.path.join(WORK, "hits.sto"),
    params:
        trim="--trim" if config["hmmalign_trim"] else "",
    conda: "envs/hmmer.yaml"
    log: os.path.join(LOGDIR, "hmmalign.log")
    shell:
        "{ENV_HMMER}hmmalign {params.trim} --amino --outformat Stockholm "
        "-o {output.sto} {input.hmm} {input.faa} &> {log}"


# ---------------------------------------------------------------------------
# 8. Triad columns are learned from the PF12386 hits alone, using redundancy-
#    weighted residue frequencies, then applied unchanged to every sequence in
#    the alignment including the SSF54001-only ones.
# ---------------------------------------------------------------------------
rule triad:
    input:
        sto=os.path.join(WORK, "hits.sto"),
        combined=COMBINED,
        idmap=os.path.join(WORK, "hits_idmap.tsv.gz"),
        weights=os.path.join(TABDIR, "sequence_weights.tsv"),
    output:
        cands=os.path.join(TABDIR, "triad_candidates.tsv"),
        chosen=os.path.join(WORK, "triad_columns.json"),
        keep=os.path.join(WORK, "triad_pass_ids.txt"),
        afa=os.path.join(WORK, "triad_pass_matchcols.afa"),
        colstats=os.path.join(TABDIR, "alignment_column_stats.tsv"),
        tiers=os.path.join(TABDIR, "triad_filter_by_tier.tsv"),
    params:
        cfg=CFG,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "triad.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/triad_detect_filter.py --sto {input.sto} --config {params.cfg} "
        "--combined {input.combined} --idmap {input.idmap} --weights {input.weights} "
        "--out-candidates {output.cands} --out-chosen {output.chosen} "
        "--out-keep {output.keep} --out-afa {output.afa} "
        "--out-colstats {output.colstats} --out-tiers {output.tiers} &> {log}"


rule c71_faa:
    input:
        keep=os.path.join(WORK, "triad_pass_ids.txt"),
        idmap=os.path.join(WORK, "hits_idmap.tsv.gz"),
    output:
        faa=os.path.join(OUTDIR, "c71.faa"),
    threads: 8
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "c71_faa.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/extract_seqs.py --keep-ids {input.keep} --idmap {input.idmap} "
        "--out-faa {output.faa} --threads {threads} &> {log}"


# ---------------------------------------------------------------------------
# 9. Phylogeny.
# ---------------------------------------------------------------------------
rule tree_input:
    input:
        afa=os.path.join(WORK, "triad_pass_matchcols.afa"),
    output:
        aln=os.path.join(WORK, "tree_input.aln"),
        reps=os.path.join(TABDIR, "tree_representatives.tsv"),
    params:
        cfg=CFG, tmp=os.path.join(WORK, "tree_prep"),
    threads: 8
    conda: "envs/phylo.yaml"
    log: os.path.join(LOGDIR, "tree_input.log")
    shell:
        "{ENV_PHYLO}python {SCRIPTS}/tree_input.py --afa {input.afa} --out-aln {output.aln} "
        "--out-reps {output.reps} --config {params.cfg} --tmpdir {params.tmp} "
        "--threads {threads} &> {log}"


rule iqtree:
    input:
        aln=os.path.join(WORK, "tree_input.aln"),
    output:
        tree=os.path.join(OUTDIR, "tree", "c71.treefile"),
    params:
        prefix=os.path.join(OUTDIR, "tree", "c71"),
        mode=config["tree"]["mode"], seed=config["tree"]["seed"],
    threads: config["tree"]["iqtree_threads"]
    conda: "envs/phylo.yaml"
    log: os.path.join(LOGDIR, "iqtree.log")
    shell:
        """
        {ENV_PHYLO}
        mkdir -p $(dirname {params.prefix})
        if [ "{params.mode}" = "fast" ]; then
            iqtree2 -s {input.aln} --prefix {params.prefix} --seqtype AA \
                -m LG+F+G4 --fast -T {threads} -seed {params.seed} -redo &> {log}
        else
            iqtree2 -s {input.aln} --prefix {params.prefix} --seqtype AA \
                -m MFP -B 1000 --alrt 1000 -T {threads} -seed {params.seed} -redo &> {log}
        fi
        """


# ---------------------------------------------------------------------------
# 10. Active-site subgroups (redundancy-weighted), then coupling.
# ---------------------------------------------------------------------------
rule active_site:
    input:
        afa=os.path.join(WORK, "triad_pass_matchcols.afa"),
        chosen=os.path.join(WORK, "triad_columns.json"),
        merged=MERGED,
        idmap=os.path.join(WORK, "hits_idmap.tsv.gz"),
        weights=os.path.join(TABDIR, "sequence_weights.tsv"),
    output:
        assign=os.path.join(TABDIR, "subgroup_assignments.tsv"),
        done=touch(os.path.join(FIGDIR, ".activesite.done")),
    params:
        figdir=FIGDIR, tabdir=TABDIR, cfg=CFG,
    threads: 8
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "active_site.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/active_site_analysis.py --afa {input.afa} --chosen {input.chosen} "
        "--merged {input.merged} --idmap {input.idmap} --weights {input.weights} "
        "--figdir {params.figdir} --tabdir {params.tabdir} --config {params.cfg} "
        "--out-assign {output.assign} &> {log}"


rule coupling:
    input:
        assign=os.path.join(TABDIR, "subgroup_assignments.tsv"),
        chosen=os.path.join(WORK, "triad_columns.json"),
        weights=os.path.join(TABDIR, "sequence_weights.tsv"),
    output:
        table=os.path.join(TABDIR, "barcode_coupling.tsv"),
    params:
        figdir=FIGDIR, cfg=CFG,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "coupling.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/coupling.py --assign {input.assign} --chosen {input.chosen} "
        "--weights {input.weights} --config {params.cfg} --figdir {params.figdir} "
        "--out {output.table} &> {log}"


# ---------------------------------------------------------------------------
# 11. Sequence similarity network.
# ---------------------------------------------------------------------------
rule ssn_align:
    input:
        faa=os.path.join(OUTDIR, "c71.faa"),
    output:
        m8=os.path.join(WORK, "ssn_edges.m8"),
        nodes=os.path.join(WORK, "ssn_nodes.faa"),
    params:
        cfg=CFG, tmp=os.path.join(WORK, "ssn_tmp"),
    threads: config["ssn"]["threads"]
    conda: "envs/network.yaml"
    log: os.path.join(LOGDIR, "ssn_align.log")
    shell:
        "{ENV_NET}bash {SCRIPTS}/ssn_align.sh {input.faa} {output.m8} {output.nodes} "
        "{params.cfg} {params.tmp} {threads} &> {log}"


rule ssn:
    input:
        m8=os.path.join(WORK, "ssn_edges.m8"),
        nodes=os.path.join(WORK, "ssn_nodes.faa"),
        assign=os.path.join(TABDIR, "subgroup_assignments.tsv"),
        merged=MERGED,
    output:
        clusters=os.path.join(TABDIR, "ssn_clusters.tsv"),
        sweep=os.path.join(TABDIR, "ssn_threshold_sweep.tsv"),
        graphml=os.path.join(OUTDIR, "ssn.graphml"),
    params:
        figdir=FIGDIR, tabdir=TABDIR, cfg=CFG,
    conda: "envs/network.yaml"
    log: os.path.join(LOGDIR, "ssn.log")
    shell:
        "{ENV_NET}python {SCRIPTS}/ssn.py --m8 {input.m8} --nodes {input.nodes} "
        "--assign {input.assign} --merged {input.merged} --config {params.cfg} "
        "--out-clusters {output.clusters} --out-sweep {output.sweep} "
        "--out-graphml {output.graphml} --figdir {params.figdir} "
        "--tabdir {params.tabdir} &> {log}"


# ---------------------------------------------------------------------------
# 12. Convergence: how many times did each subgroup arise?
# ---------------------------------------------------------------------------
rule convergence:
    input:
        tree=os.path.join(OUTDIR, "tree", "c71.treefile"),
        assign=os.path.join(TABDIR, "subgroup_assignments.tsv"),
        reps=os.path.join(TABDIR, "tree_representatives.tsv"),
    output:
        table=os.path.join(TABDIR, "convergence.tsv"),
    params:
        figdir=FIGDIR, cfg=CFG,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "convergence.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/convergence.py --tree {input.tree} --assign {input.assign} "
        "--reps {input.reps} --config {params.cfg} --figdir {params.figdir} "
        "--out {output.table} &> {log}"


# ---------------------------------------------------------------------------
# 13. Genome-level table, then phylogenetic logistic regression on the GTDB
#     trees, with the genome-quality covariates as the detection-bias model.
# ---------------------------------------------------------------------------
rule genome_table:
    input:
        table=os.path.join(WORK, "sample_table_unified.tsv"),
        combined=COMBINED,
        assign=os.path.join(TABDIR, "subgroup_assignments.tsv"),
        idmap=os.path.join(WORK, "hits_idmap.tsv.gz"),
        maps=all_maps,
    output:
        genomes=os.path.join(TABDIR, "genome_level_table.tsv"),
    params:
        map_dir=os.path.join(WORK, "map"), cfg=CFG,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "genome_table.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/prep_genome_table.py --table {input.table} "
        "--combined {input.combined} --assign {input.assign} --idmap {input.idmap} "
        "--map-dir {params.map_dir} --config {params.cfg} --out {output.genomes} &> {log}"


rule phyloglm:
    input:
        genomes=os.path.join(TABDIR, "genome_level_table.tsv"),
        bac_tree=config["trees"]["bacteria"],
        ar_tree=config["trees"]["archaea"],
    output:
        coefs=os.path.join(TABDIR, "phyloglm_coefficients.tsv"),
        prev=os.path.join(TABDIR, "bias_adjusted_prevalence.tsv"),
        dstat=os.path.join(TABDIR, "phylogenetic_signal_D.tsv"),
    params:
        cfg=CFG, tabdir=TABDIR,
    threads: 4
    conda: "envs/r.yaml"
    log: os.path.join(LOGDIR, "phyloglm.log")
    shell:
        "{ENV_R}Rscript {SCRIPTS}/phyloglm.R {input.genomes} {input.bac_tree} {input.ar_tree} "
        "{params.cfg} {output.coefs} {output.prev} {output.dstat} &> {log}"


# ---------------------------------------------------------------------------
# 14. Figures and report.
# ---------------------------------------------------------------------------
rule plots_overview:
    input:
        combined=COMBINED,
        merged=MERGED,
        stats=os.path.join(TABDIR, "hmm_search_stats.tsv"),
        fdr=os.path.join(TABDIR, "decoy_fdr.tsv"),
        colstats=os.path.join(TABDIR, "alignment_column_stats.tsv"),
        chosen=os.path.join(WORK, "triad_columns.json"),
        tiers=os.path.join(TABDIR, "triad_filter_by_tier.tsv"),
        prev=os.path.join(TABDIR, "bias_adjusted_prevalence.tsv"),
    output:
        done=touch(os.path.join(FIGDIR, ".overview.done")),
    params:
        figdir=FIGDIR, tabdir=TABDIR, cfg=CFG,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "plots_overview.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/plots_overview.py --combined {input.combined} --merged {input.merged} "
        "--search-stats {input.stats} --fdr {input.fdr} --colstats {input.colstats} "
        "--chosen {input.chosen} --tiers {input.tiers} --prevalence {input.prev} "
        "--figdir {params.figdir} --tabdir {params.tabdir} --config {params.cfg} &> {log}"


rule plot_tree:
    input:
        tree=os.path.join(OUTDIR, "tree", "c71.treefile"),
        assign=os.path.join(TABDIR, "subgroup_assignments.tsv"),
        merged=MERGED,
    output:
        done=touch(os.path.join(FIGDIR, ".tree.done")),
    params:
        figdir=FIGDIR, cfg=CFG,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "plot_tree.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/plot_tree.py --tree {input.tree} --assign {input.assign} "
        "--merged {input.merged} --figdir {params.figdir} --config {params.cfg} &> {log}"


rule report:
    input:
        os.path.join(FIGDIR, ".overview.done"),
        os.path.join(FIGDIR, ".tree.done"),
        os.path.join(FIGDIR, ".activesite.done"),
        stats=os.path.join(TABDIR, "hmm_search_stats.tsv"),
        chosen=os.path.join(WORK, "triad_columns.json"),
        assign=os.path.join(TABDIR, "subgroup_assignments.tsv"),
        c71=os.path.join(OUTDIR, "c71.faa"),
        tiers=os.path.join(TABDIR, "triad_filter_by_tier.tsv"),
        fdr=os.path.join(TABDIR, "decoy_fdr.tsv"),
        conv=os.path.join(TABDIR, "convergence.tsv"),
        ssn=os.path.join(TABDIR, "ssn_clusters.tsv"),
        coup=os.path.join(TABDIR, "barcode_coupling.tsv"),
        coefs=os.path.join(TABDIR, "phyloglm_coefficients.tsv"),
        dstat=os.path.join(TABDIR, "phylogenetic_signal_D.tsv"),
        tips=os.path.join(TABDIR, "tree_tip_matching.tsv"),
    output:
        md=os.path.join(OUTDIR, "report.md"),
    params:
        outdir=OUTDIR, tabdir=TABDIR, cfg=CFG,
    conda: "envs/py.yaml"
    log: os.path.join(LOGDIR, "report.log")
    shell:
        "{ENV_PY}python {SCRIPTS}/make_report.py --outdir {params.outdir} --tabdir {params.tabdir} "
        "--config {params.cfg} --chosen {input.chosen} --c71 {input.c71} "
        "--out {output.md} &> {log}"
