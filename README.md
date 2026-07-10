# C71 / pseudomurein endo-isopeptidase Pei screening pipeline

Calls genes on the GTDB archaeal assemblies with Prodigal, screens those plus
~350,000 existing bacterial proteomes for PF12386 and SSF54001, filters to
sequences with an intact Cys-His-Asp catalytic triad, builds a phylogeny and a
sequence similarity network, partitions the active sites into subgroups, tests
whether those subgroups are convergent, tests whether the barcode positions are
coupled, and regresses C71 presence on genome quality with a phylogenetic
logistic model on the GTDB trees.

```
./run.sh preflight  # validate every input before anything is submitted
./run.sh dry        # DAG only
./run.sh cluster    # submit to SGE
./run.sh local      # single node, CORES=32 ./run.sh local
```

Install first (below), then edit `config.yaml`, then run preflight.

---

## Installation

The compute nodes have no internet, and the `defaults` / `anaconda` channels are
off limits. So the environments are **built somewhere else and shipped in**.
Snakemake's own conda machinery is never used on the cluster: if `envs_root` is
set in `config.yaml`, every rule activates a pre-built environment and nothing
is ever downloaded or solved.

Seven environments: `hmmer`, `prodigal`, `phylo`, `network`, `py`, `r`,
`selection`. All are built from **conda-forge and bioconda only**. Every
`envs/*.yaml` carries `- nodefaults`, and `build_envs.sh` passes
`--override-channels` with `CONDA_CHANNEL_PRIORITY=strict`. Either alone can be
undone by a stray `~/.condarc`; both together cannot.

### Pick one route, not both

There are two ways to get the environments onto the cluster, and they are
mutually exclusive. Do **not** run Steps 1–2 and the login-node route both.

| | Route A — conda-pack (default) | Route B — login node has internet |
|---|---|---|
| When | The cluster, including the login node, cannot reach the internet | Your login node can reach conda-forge, and the compute nodes share its filesystem |
| Build | Step 1 `build_envs.sh` on a separate internet machine | `snakemake --conda-create-envs-only` on the login node |
| Ship | Step 2 `install_envs.sh` on the cluster | nothing to ship; compute nodes reuse the login node's prefix |
| `config.yaml` | `envs_root: /…/c71_envs` | `envs_root: null` |
| Profile | leave as shipped | add `conda-prefix: /…` |

Both routes then run the same **Step 3** (preflight → dry → cluster). Route A is
the default because it works even when nothing on the cluster can reach the
internet. If you take Route B, skip straight to
["If your login node has internet"](#if-your-login-node-does-have-internet) and
then do Step 3 with `envs_root` left `null`.

### Route A, Step 1 — on a machine with internet

```bash
git clone <repo> c71_pipeline && cd c71_pipeline
./setup/build_envs.sh ./dist
```

This solves the seven environments, packs each into a relocatable tarball with
`conda-pack`, downloads the Snakemake wheels for offline `pip`, and writes an
explicit package lock per environment into `dist/locks/`.

`conda-pack` tarballs are architecture- and glibc-specific. Build on Linux
x86-64, on an OS no newer than the cluster's. If your laptop is macOS or arm64,
build in a container:

```bash
docker run --rm -v "$PWD":/w -w /w --platform linux/amd64 \
    condaforge/miniforge3:latest ./setup/build_envs.sh /w/dist
```

Miniforge is the conda-forge-only installer. Using a stock Anaconda install as
the build base is what pulls in `defaults` behind your back.

### Route A, Step 2 — on the cluster

```bash
rsync -a dist/ user@cluster:/ebio/abt3_scratch/jmarsh/c71_dist/
# then, on the cluster:
cd c71_pipeline
./setup/install_envs.sh /ebio/abt3_scratch/jmarsh/c71_dist \
                        /ebio/abt3_projects/software/c71_envs
```

It untars each environment, runs `conda-unpack` (which rewrites the build
machine's absolute paths baked into shebangs, scripts and rpaths), installs
Snakemake from the vendored wheels, and then **verifies**: every binary the
pipeline calls is on the activated PATH, every Python module imports, and
`ape`, `phylolm`, `caper`, `data.table`, `yaml` all load in R. It exits non-zero
if anything is missing, before you have burned any queue time.

Put `envs_root` on a filesystem the compute nodes can read. `/tmp` on the submit
host is not one.

### Step 3 — both routes

Set `envs_root` to match the route you took:

```yaml
# config.yaml
# Route A (conda-pack):  the prefix install_envs.sh wrote to
envs_root: /ebio/abt3_projects/software/c71_envs
# Route B (login node):  leave it null
# envs_root: null
```

```bash
python scripts/check_snakefile.py   # static wiring check, no snakemake needed
./run.sh preflight                  # validates trees, tables, HMMs, faa paths
./run.sh dry
./run.sh cluster
```

`run.sh` reads `envs_root` and decides the deployment for you. Route A, with
`envs_root` set: it never passes `--software-deployment-method conda`, so no job
tries to download anything. Route B, with `envs_root` null: it passes that flag
automatically, and the compute nodes reuse the prefix the login node built. You
do **not** add the flag by hand. `profiles/sge/config.yaml` deliberately omits it,
because a profile setting would override `run.sh` and make all 700 array jobs
attempt a download.

### If your login node does have internet

This is **Route B**, and it replaces Steps 1 and 2 — do not run `build_envs.sh`
or `install_envs.sh` as well. Leave `envs_root: null`.

**Stage the external inputs first, then run preflight, then create the envs.**
`snakemake --conda-create-envs-only` builds the whole job graph before it makes a
single environment, and the graph will not build until every *external* input —
the files no rule produces — is present. A missing structure or HMM aborts env
creation with a `MissingInputException`, which looks like a conda problem and is
not. Preflight checks all of them in one pass and needs no envs, so run it first:

```bash
./run.sh preflight        # lists every missing input with the path it expected
```

The external inputs preflight validates (stage each before proceeding):

| Input | config key | how to get it |
|---|---|---|
| PeiW-CD structure `8jx4.cif` | `specificity.structure` | `wget https://files.rcsb.org/download/8JX4.cif -O <path>/8jx4.cif` |
| PF12386, SSF54001 HMMs | `profiles.*.path` | `hmmfetch` from Pfam / a SUPERFAMILY build |
| Pfam-A.hmm | `specificity.pfam_hmm` | Pfam release, `hmmpress`ed |
| Pmur marker HMMs | `specificity.pmur_hmm_dir` | your marker set |
| geNomad DB | `specificity.genomad_db` | `genomad download-database` |
| bacteria + archaea GTDB trees | `trees.*` | GTDB release |
| sample table, metadata, archaeal `.fna` | `inputs.*` | your data |

Then create the environments once, into a shared prefix:

```bash
snakemake --software-deployment-method conda \
          --conda-prefix /ebio/abt3_projects/software/c71_snakemake_envs \
          --conda-create-envs-only --cores 1
```

(Drop `--conda-frontend mamba`: recent Snakemake ignores it and warns, because
the classic conda solver already uses libmamba.)

Then add `conda-prefix:` (pointing at that same prefix) to
`profiles/sge/config.yaml`. That is the only profile edit; `run.sh` supplies
`--software-deployment-method conda` on its own because `envs_root` is null. The
compute nodes reuse what the login node built. Then continue with Step 3 above.

The Route A conda-pack path works even when nothing on the cluster can reach the
internet, which is why it is the default.

### Versions

Pinned in `envs/*.yaml`: HMMER 3.4, Prodigal 2.6.3, IQ-TREE 2.3.6, trimAl 1.4.1,
MMseqs2 15, DIAMOND 2.1.9, seqkit 2.8.2, Python 3.11, R 4.3 with `phylolm` 2.6.5
and `caper` 1.0.3. `dist/locks/` records the exact build that was solved.

---

Edit `config.yaml` first, then run preflight. It takes about a minute and it
catches the things that otherwise surface six hours into a 700-job array:
a PF12386 copy with no GA line, a tree whose tip labels don't match your sample
IDs, a QC column named something other than what `qc:` expects, duplicate sample
IDs, dead faa paths, a metadata file that is actually comma-separated. Each of
those is a tested failure case, not a hypothetical.

```
$ ./run.sh preflight
[OK   ] sample table rows            352,914
[ERROR] HMM PF12386: config asks for --cut_ga but there is no GA line
[ERROR] bacteria tip matching (mode=identity)
         41,203/352,914 (11.7%), threshold 50%
         unmatched e.g. ['SRR12345678', ...]
```

---

## What the biology is

PeiW and PeiP are prophage-encoded enzymes that cleave the ε(Ala)–Lys isopeptide
bond of archaeal pseudomurein, in Methanobacteriales and Methanopyrales.
PF12386 is the C-terminal catalytic cysteine-protease domain; the full enzymes
carry four N-terminal pseudomurein-binding repeats.

Both structures use a **Cys–His–Asp** triad, so CHD is the hypothesis this
pipeline tests, and the C→A, H→A and D→A mutants are inactive
([Wang et al. 2025](https://doi.org/10.1016/j.ijbiomac.2025.141813)). PeiW-CD is
8JX4 (C198/H233/D250); PeiP is 8Z4F (C213/H248/D272). Wang et al. place the
catalytic domain in the **papain superfamily** fold, not the transglutaminase
fold an earlier draft of this README asserted from a secondary summary. The
triad is CHD either way, and that is what the filter uses.

(A still earlier draft argued for Cys–His–**Asn** on clan-CA grounds. That was
also wrong.)

Two consequences run through the design:

- The two profiles are not equivalent evidence. PF12386 says "Pei catalytic
  domain". SSF54001 is the SCOP *cysteine proteinases* superfamily and says only
  "papain- or transglutaminase-like fold". They are kept in separate evidence
  tiers everywhere.
- These are phage genes in a narrow host clade, so genome-level prevalence is
  confounded by both phylogeny and assembly quality. Both are modelled, not
  assumed away.

---

## Probing target specificity

Pei cleaves the ε-isopeptide bond between **alanine** and lysine, so the alanine
is P1. [Kandler & König 1978](https://doi.org/10.1007/BF00415722) report the
molar ratio Lys : Ala : Glu as 1 : 1.2 : 2 — **one alanine per subunit** —
therefore the alanine they find completely replaced by threonine in
*Methanobrevibacter ruminantium* M1 is *the* P1 alanine. That closes the question
an earlier draft of this README left open.

The same paper reports a second substitution on the **other side** of the
scissile bond: about a quarter of the lysine in *Methanobrevibacter smithii* PS
is ornithine, whose δ-amine is one methylene shorter than a lysine ε-amine. So
the wall varies at P1 *and* at P1′, and these are different questions.

`scripts/cellwall_reference.py` encodes that paper, and only that paper. What
the enzymes do with those walls comes from Wang et al. 2025 (below).

**Genus-level assignment is not uncertain, it is wrong.** Both strains the 1978
paper calls "*M. ruminantium*" now sit in *Methanobrevibacter*: M1 is the type
strain of *M. ruminantium* (Thr/Lys), PS is the type strain of *M. smithii*
(Ala/Lys–Orn). A third species of the same genus, *M. arboriphilus*, is Ala/Lys.
One genus, three wall chemistries. `chemistry_for()` returns `unknown` rather
than guessing, and `genus_is_homogeneous()` is what decides.

**Two claims are refused.** Serine-for-alanine in *Methanosphaera stadtmanae* and
an ornithine modification in *Methanopyrus kandleri* both circulate in secondary
summaries; neither is in Kandler & König 1978 (those species were described in
1985 and 1991). They are labelled `unsupported`, not silently used. Supply the
primary reference and add them to `REFERENCE` to activate them.

**The published assay only probes one axis.** H-Glu-γ-Ala-pNA has a chromophore
where the acceptor residue should be, so it tests P1 and cannot test P1′.
`assay_panel.tsv` therefore emits `prediction_p1` and `prediction_p1_prime`
separately, and says which picks need real Ala–Lys / Ala–Orn isopeptide
substrates instead.

| Species | P1 | P1′ | Evidence |
|---|---|---|---|
| *M. ruminantium* | **Thr** | Lys | type strain M1, complete replacement |
| *M. smithii* | Ala | **Lys/Orn** | type strain PS, ~1/4 Orn |
| *M. arboriphilus* | Ala | Lys | species |
| *M. thermautotrophicus* | Ala | Lys | species (the PeiW/PeiP host genus) |
| *M. formicicum*, *M. bryantii* | Ala | Lys | species / strain M.o.H. |
| *Methanospirillum hungatei* | — | — | protein sheath, no sacculus |

The enzymes are two-module: four N-terminal pseudomurein-binding repeats (PMBR,
PF09373) and a C-terminal C71 catalytic domain. The 2010 review states that the
repeat array "probably serves as a determinant of substrate specificity." So
specificity has two axes — which bond the groove cuts, and which sacculus the
repeats bind — and the pipeline measures both.

`test/test_cellwall_reference.py` checks every row of that table, the genus
refusal, and the two disputed claims.

### What Wang et al. 2025 settles, and what it broke

[Wang et al., *Int J Biol Macromol*](https://doi.org/10.1016/j.ijbiomac.2025.141813)
solved both enzymes and assayed them on synthetic isopeptides. Four consequences.

**Two structures, two numberings.** PeiW-CD is **8JX4** (triad C198 / H233 /
D250); full-length PeiP is **8Z4F** (triad **C213 / H248 / D272** — exactly the
numbers you gave me at the outset, which were right). Y174, V252 and C265 are
quoted in **PeiW** numbering. An earlier version of this config pointed at 8Z4F
while using those numbers, silently selecting the wrong residues.
`specificity.structure_expect` now asserts the identity of every named residue
and `groove_map.py` exits non-zero on a mismatch, naming both PDB entries.

**The P1 axis is real, already measured, and the two papers disagree about it.**
Subedi et al. 2015 bought Glu-γ-**Thr**-pNA from JPT, assayed both enzymes, and
saw nothing; they also found *M. ruminantium* M1, the Thr-walled organism, is not
lysed by either. Wang et al. 2025 report PeiW cleaving the Glu-γ-Thr-ε-**Lys**
isopeptide. Both can be true: the pNA leaving group sits where the acyl acceptor
belongs, so the chromogenic series leaves S1′ empty. They reconcile if S1 and S1′
are coupled — PeiW takes the extra methyl of threonine only when a real ε-amino
acceptor is bound.

The consequence is operational, not philosophical: **a P1 claim cannot be tested
with a pNA substrate** unless the residue is Ala or Ser, the two the pNA series
does cleave. `assay_panel.py` emits a `substrate_format` column and routes every
other P1 to an isopeptide. An earlier version of this pipeline asserted that the
Thr chromogen had never been synthesised; `test_lysis_reference.py` greps the
whole repo and fails the build if that claim reappears.

Cleaving a soluble isopeptide is still not lysing a sacculus. PeiW cuts the Thr
isopeptide and cannot lyse the Thr-walled organism. That gap is unexplained and
the code says so rather than smoothing it.

**Serine is accepted at P1.** Both enzymes cleave Glu-γ-Ser-pNA. Aspartate is
rejected at P2 (Asp-β-Ala-pNA is a poor substrate), and the Glu–Ala bond need not
be a γ-isopeptide at all — only the Ala–Lys bond does.

**No published substrate tests P1′.** The pNA series has a chromophore where the
acceptor sits; the six isopeptides all carry Lys. The Orn acceptor of
*M. smithii* is untested, and the panel says so.

**Both enzymes need a divalent metal, and they disagree about which.** After
EDTA, PeiW and PeiP retain under 1% activity. Ca²⁺ restores both. Mn, Mg, Ba and
Ni restore PeiW; for PeiP they give under 15%. This is the sharpest measured
difference between the two characterised enzymes, and neither the groove nor the
four-class partition predicts it. No deposited Pei structure resolves the cation.
`groove_map.py` looks for one, builds a coordination shell if it finds one, and
reports the absence plainly if it does not — `sdp.py` then tests that shell for
specificity-determining positions separately from the groove.

**The prespecified hypothesis: Pei sequence tracks host wall chemistry.**
*M. ruminantium* M1 has Thr at P1, resists PeiW and PeiP, and its own prophage
Φmru encodes an endoisopeptidase, PeiR, that lyses it while showing little
homology to either. A divergent enzyme for a divergent wall, in the one host where
it can be checked. PeiR is registered as a sentinel in `config.yaml`: if PF12386's
curated threshold misses it, the screen is blind in exactly the direction it
cares about.

**Three PMB motifs are required to bind, and the threshold is a cliff.**
Visweswaran et al. 2011 fused one, two or three motifs to GFP. Three bind the
pseudomurein sacculus. Two bind lysozyme-treated *L. lactis* and *E. coli*
spheroplasts and **not** pseudomurein. One binds nothing. So a C71 domain carried
on fewer than three motifs cannot dock on an intact sacculus whatever its active
site looks like — a falsifiable claim about every hit, not a covariate to regress
out. `domain_arch.py` emits `pmbr_binding_competent` and `predicted_binding`.

That makes the repeat count load-bearing. A PMB motif is 30–35 residues; one
dropped repeat at a strict domain E-value flips the call. hmmscan now runs once
at a permissive E-value and `domain_arch.py` filters at two, classing anything
that changes side as `pmbr_count_ambiguous` and dropping it from the
`pmbr_architecture` partition. Assigning it to a threshold is not the same as
assigning it to a phenotype.

**PMBR is not a pseudomurein marker.** The domain binds NAG, the one sugar shared
by murein and pseudomurein, and it is not Pei-specific: the S-layer protein
MTH719 carries three motifs and no catalytic domain. A bacterial C71+PMBR protein
may be binding exposed murein rather than being a binning artefact. Usefully, this
is also why `pmbr_architecture` is a legitimate SDP partition against a catalytic
barcode: the module reads the glycan, the groove reads the peptide. (Wang et al.'s
claim that the PB repeats improve recognition of Glu-γ-Thr/Ser is in tension with
that, and `sdp.py` says so.)

`sdp.py` will not accept `pmbr_partition_mode: both`. The repeat-count partition
and the binarized one are nested, so replication across them would count one
observation twice and satisfy `sdp_min_partitions` on a single piece of evidence.

**The published assays run where the binding module works worst.** The domain
binds pseudomurein completely at pH 9.0 (its pI is 9.2), partially at 6.5, not at
all at 4.0, and aggregates at pH 7.0. Every published Pei lysis assay runs at
pH 7.0–7.85. PMB pIs span 3–10, so `domain_arch.py` computes each protein's own
pI from its own PMB span and `assay_panel.py` carries the advice per pick.

This reopens the M1 paradox. PeiW cleaves the Thr isopeptide but cannot lyse the
Thr-walled organism. If the failure were catalytic, the isopeptide would resist
too. It may be a docking failure. The discriminating experiment is cheap: PeiW
catalytic domain vs full-length PeiW, on M1 sacculi. The tension is that the PMB
domain is thought to read the glycan while Thr-for-Ala is a peptide substitution,
so a binding explanation needs the two coupled somehow. Unresolved, and recorded
as unresolved.

**The only experiment that can prove this pipeline wrong.** `rule lysis_check`
scores the host-chemistry rule against Subedi et al.'s eleven-strain plate-lysate
panel. It runs in milliseconds, before 350,000 proteomes are searched, and
`--strict` aborts the run on a disagreement. Four falsifiable rows, four
agreements. The rule declines to predict *Methanobrevibacter* sp. SM9 — lysed by
PeiW, not by PeiP — because no chemistry has been published for it. That is the
honest output, and the differential is a free hypothesis for the SDP analysis.

**The family has a published four-class partition — from two residues.**
V252 and C265 are the only non-conserved positions near the His/Asp centre, and
they set activity:

| Class | 252 | 265 | Measured activity |
|---|---|---|---|
| I | V | C | active (PeiW, PeiP) |
| II | V | V | ~50% of PeiW-CD |
| III | T/S | I | strongly reduced |
| IV | A | M/W | strongly reduced |

`scripts/pei_class.py` assigns every sequence, and reports the adjusted Rand
index against our own k-means subgroups and SSN clusters. This is an **external
prior for k = 4** that owes nothing to our clustering. If the k-means doesn't
recover it, that is the result, and the script says so rather than smoothing it.

It also creates a fresh circularity hazard, which is handled: **V252 sits two
residues from the catalytic D250**, so it falls inside the ±5 triad-flank
barcode. `sdp.py` excludes both class-defining columns from the SDP test
whenever the `pei_class` partition is used, and only for that partition. The
synthetic test asserts exactly this.

A large class III/IV population among triad-positive sequences would mean the
triad filter is retaining proteins that cannot cleave. That's a finding about
the filter, and `pei_class.py` prints it as one.

### The circularity problem, and how it is avoided

The active-site subgroups are k-means on the triad-flank barcode. Testing
specificity-determining positions on those same columns, against those same
groups, rediscovers them by construction. The result would be beautiful and
empty.

`sdp.py` therefore calls SDPs under four partitions that never saw the barcode —
tree clades, SSN clusters, PMBR architecture, host genus — and reports a column
only if it clears FDR under at least two. `sdp_concordance.tsv` gives the
pairwise Jaccard of the four SDP sets. `preflight.py` fails if
`active_site_subgroup` appears in `sdp_partitions`.

The permutation null shuffles group labels **within tree clades**, so it
preserves phylogeny instead of destroying it. A global shuffle makes every
phylogenetically clustered position look like an SDP.

### What each script contributes

| Script | Question |
|---|---|
| `domain_arch.py` | How many PMBR repeats? Any non-PMBR binding module (LysM, SH3b, PG_binding_1, choline-binding)? Those cannot bind pseudomurein. |
| `module_trees.py` | Do the binding and catalytic modules share a history? Mantel on patristic distances, normalised Robinson-Foulds, tanglegram. Incongruence means retargeting by module swap. |
| `cellwall_genotype.py` | What cross-link does the host actually build? Pmur marker screen plus the literature P1 call, by taxon. |
| `groove_map.py` | Which PF12386 match columns line the substrate groove in PeiW-CD (8JX4)? Refuses a structure whose residues contradict the config. |
| `pei_class.py` | Which of Wang et al.'s four classes is each sequence, and does our k-means agree? |
| `sdp.py` | Which columns determine specificity, non-circularly, and are they in the groove? |
| `coupling.py` | Are the coupled barcode pairs spatially adjacent, or is that phylogenetic signal APC missed? |
| `selection.py` | Purifying selection on the triad, episodic diversifying selection on the groove? FEL, MEME. |
| `assay_panel.py` | Which N proteins to synthesise, and what does each predict on EγA/EγT/EγS-pNA? |

### Three things it deliberately refuses to do

It does not call a genome "Ala-type" because a Pmur marker is missing. A marker
can be absent because the pathway is absent, because the MAG is 70% complete, or
because the HMM is bad. Genomes below `pmur_min_markers` are
`no_pathway_detected`, and completeness is carried forward so the regression can
condition on it.

It does not infer P1 chemistry from the marker set, because no published Pmur
marker is known to determine P1. The call comes from `cellwall_reference.py` at
the level the primary paper supports (species, and in two cases a single type
strain), and `p1_source` records which.

`groove_map.py` **exits non-zero** if the triad columns do not map to C/H/D in
the structure, or if under half the match columns have a structural residue. A
groove built on a bad alignment silently corrupts every SDP, every coupling
contact and every selection contrast downstream. Better to stop.

### Verification

`test/test_specificity.py` builds a world with a known answer: a synthetic PeiP
whose triad, two specificity seeds and five genus-driven SDP columns cluster
within 8 Å while everything else is 30 Å away; five decoy columns driven by the
barcode subgroup and nothing else; two trees, one star and one balanced.

```
groove columns recovered: [112, 147, 150, 168, 171, 145, 149, 152, 170, 173]  (exact)
--- star tree ---      null: global_fallback   (degenerate within-clade null, detected)
--- balanced tree ---  null: within_clade
replicated SDPs      : [145, 149, 152, 170, 173]   = the planted genus-driven columns
circular decoys      : [40, 60, 80, 200, 220]      = none called
groove enrichment    : OR = inf, p = 1.13e-09
RESULT: PASS
```

Two real bugs surfaced here and both are fixed. `hmmalign --trim` removes
terminal residues, and the first draft kept indexing the structure from residue
zero, shifting the entire groove by the number trimmed. And cutting tree clades
at a fixed relative depth produced 240 singleton clades on an unbalanced tree,
which makes a within-clade shuffle the identity permutation: the null becomes the
observation, every z is zero, and the test can never reject anything. It fails
silently and looks exactly like "no signal". Clades are now cut by subtree size,
and `null_is_degenerate()` catches the residual case and falls back loudly.

### What it needs that this cluster cannot download

`Pfam-A.hmm` (pressed), the PeiP structure (8Z4F), a directory of Pmur marker
HMMs, and the geNomad database. Paths are in `config.yaml` under `specificity:`;
`./run.sh preflight` checks every one, verifies the Pfam accessions are present
in your release, and resolves the groove seed residues against the structure
before anything is submitted.

---

## Tuning for this SGE cluster

From `sge_probe_report` (SGE 8.1.9, 30 hosts, 2026-07-06). Four things were
wrong for this cluster; the fourth was costing more than the science.

**`h_vmem` is consumable, so SGE reads it per slot.** The old profile passed
`-pe parallel {threads} -l h_vmem={mem_mb}M`. For IQ-TREE that requested 128 GB
**per slot** across 32 slots: 4 TB, more than any node has. The job would have
pended forever with no error. `profiles/sge/submit.sh` now takes the job's total
memory and divides by the slot count. It also drops `-pe` entirely for
single-threaded rules, and caps slots at 32 — the smallest node (node529) — so
a job asking for 64 doesn't restrict itself to nine hosts.

**`standard.q` caps at 24 h.** IQ-TREE asked for 168 h and hmmalign for 72 h.
Neither is schedulable there. Long rules (`iqtree`, `hmmalign`, `ssn_align`,
`phyloglm`) are routed to `long.q` (672 h, 1,265 slots free); everything else
stays on `standard.q` (865 slots free). `submit.sh` promotes anything
mis-routed and says so on stderr. `scripts/check_snakefile.py` fails if a rule's
`h_rt` exceeds its queue's cap.

**`openmpi` is a `$round_robin` PE.** It scatters slots across hosts. Every
threaded tool here — hmmsearch, IQ-TREE, DIAMOND, MMseqs2 — forks on one node.
The profile uses `parallel` (`$pe_slots`), and the checker rejects `openmpi`.

**The batch FASTA never should have touched shared storage.** `batch_faa`,
`decoy_faa` and two `hmmsearch` rules used to write ~1 GB of concatenated
proteome plus ~1 GB of reversed decoy to `workdir`, then read both straight
back. Across ~365 batches that is roughly **700 GB of pointless round-trip**,
for data no other job ever reads. The actual search is cheap by comparison:
3×10¹¹ residues × 2 profiles × 2 (target + decoy) is tens of core-hours.

Those four rules are now one, `search_batch`, which does the concatenation, the
decoys and all four hmmsearch passes on `$TMPDIR` on the execution host. Only
the domtblout files and the per-sample map reach shared storage. Job count drops
from ~2,920 to ~365. The synthetic test asserts that `work/faa` and `work/decoy`
are never created and that scratch is empty afterwards.

The trade is granularity: re-running the search now re-reads the proteomes. That
was always the expensive half, so there was nothing to save by keeping it split.

**Sizing.** 1,000 genomes per `search_batch` job, 8 slots, 16 GB total
(2 GB/slot), ~365 jobs, `jobs: 100` in the profile so the search occupies about
800 of `standard.q`'s 865 free slots. Prodigal runs 100 genomes per job at one
slot. Raise `batch_size` to cut job count; lower it if `$TMPDIR` on the exec
hosts is tight, since scratch use is about twice the proteome bytes in a batch.

**Before you run:** `SGE_ROOT` is unset in a fresh shell on this cluster. Source
the SGE settings file first, or `status.sh` will not find `qstat`/`qacct` and
every job will look like it is still running.

```bash
SGE_DRY_RUN=1 profiles/sge/submit.sh iqtree 32 128000 168:00:00 standard.q job.sh
# submit.sh: iqtree wants 168:00:00 but standard.q caps at 23:55:00; promoting to long.q
# qsub -cwd -V -terse -N c71.iqtree -q long.q -l h_rt=168:00:00 -l h_vmem=4000M ... -pe parallel 32 job.sh
```

`cryo-em.q` is idle (0/576 slots used) and also allows 672 h. If it is fair game
for you, set `queue: "cryo-em.q"` on the long rules, or
`SGE_LONG_QUEUE=cryo-em.q` to redirect promotions.

---

## What the real inputs turned out to be

Checked against the actual files, not assumed. Every one of these was wrong in
the first version of the config, and four of the five would have failed loudly.
The fourth would not have.

| File | Reality |
|---|---|
| `faa_sample_table_90percent.tsv` | 342,759 rows, `sample` + `faa`. **3,215 of them are `d__Archaea`.** |
| `total_metadata_qc_bbmap_animals_extra.txt` | 581,395 rows, key `genome`, covers 100% of samples. `Completeness`, `Contamination`, `n_contigs`. |
| `gtdbtk.rooted.speciesnames.tree` | 9,784 tips, labelled **`s__Genus species`**, plus one bare `s__`. Not genome IDs. |
| `ar53.tree` | 4,416 tips, GTDB accessions, **species representatives only**, from an older release than the metadata. |
| `ar53_metadata_r232.tsv` | 22,343 genomes, 10,122 species reps, `checkm2_*`, `n50_contigs`, `contig_count`, `gtdb_taxonomy`. |

**Both trees are species-level.** Genomes route to their species' tip:
304,568/339,544 bacterial genomes (89.7%) onto 9,783 tips, at 31.1 genomes per
tip. `phyloglm.R` therefore aggregates to the tip and models a *species*, with
`log10_n_genomes` as a sequencing-effort covariate. An earlier draft used
`!duplicated(tree_tip)`, which kept one arbitrary genome per species and threw
away 97% of the data.

**The sample table is not bacteria-only.** Domain comes from `gtdb_domain`,
never from the filename. The 3,215 archaea go to `ar53.tree` via their species
representative. Given that Pei is a Methanobacteriota enzyme, sending them to a
bacterial tree, or dropping them, would have been the expensive mistake.

**`ctg_L50` is an N50 length, not an L50 count.** Median 39,413 against a 2.1 Mb
median genome and 89 contigs. The column name lies; the data is what we want. It
is mapped to `n50` deliberately, and the config says so.

**The taxonomy trap, which would have been silent.** The bacterial metadata
carries a block of host-organism fields — `kingdom`, `phylum`, `class`, `order`,
`family`, `common_name_x` — that is **0% populated** for microbial rows, right
next to the populated `gtdb_phylum`. `resolve_taxonomy()` used to take the first
column literally named `phylum`. Every taxonomy figure and the entire
`subgroup × taxon` enrichment would have come out empty, without an error.
It now prefers `gtdb_<rank>`, and refuses a rank column that is under half
populated. `preflight.py` warns when both columns exist.

**GTDB release skew.** 359 of `ar53.tree`'s 4,416 tips (8.1%) are absent from
`ar53_metadata_r232.tsv`. Routing through `gtdb_representative` reaches only
1,446/3,215 (45%) of the stray archaea, below `min_match_fraction`, so the
pipeline stops. Supply the `ar53.tree` from the r232 release and it clears.

**`s__` is a real tip in a real GTDB-Tk tree.** Matching unclassified genomes to
it would collapse them onto one tip and model them as a single clade. Bare rank
prefixes are excluded from tip matching. No genome in your table is affected;
the guard is there for the next one.

---

## Deliverables

| Path | What |
|---|---|
| `prodigal/*.faa` | genes called on the archaeal assemblies |
| `hmm_output/batch_*.{PF12386,SSF54001}.domtblout` | raw HMMER output, per profile |
| `hmm_output_combine.txt` | filtered hits, one row per (protein, profile), with an `evidence` tier |
| `faa_sample_table_90percent_taxonomy.tsv` | the above, left-merged onto the unified sample table |
| `outputs/c71.faa` | full-length sequences with a complete catalytic triad |
| `outputs/tree/c71.treefile` | IQ-TREE ML phylogeny, UFBoot + SH-aLRT |
| `outputs/ssn.graphml` | sequence similarity network, for Cytoscape |
| `outputs/figures/*.{png,svg}` | 25 figures, 600 dpi PNG and SVG with live text |
| `outputs/tables/*.tsv` | every number behind every figure |
| `outputs/report.md` | funnel, FDR, triad call, subgroups, coupling, convergence, regression |

Key tables: `decoy_fdr.tsv`, `triad_candidates.tsv`, `triad_filter_by_tier.tsv`,
`sequence_weights.tsv`, `subgroup_assignments.tsv`, `barcode_coupling.tsv`,
`ssn_threshold_sweep.tsv`, `ssn_vs_subgroup_agreement.tsv`, `convergence.tsv`,
`genome_level_table.tsv`, `phyloglm_coefficients.tsv`,
`bias_adjusted_prevalence.tsv`, `phylogenetic_signal_D.tsv`.

---

## The decisions that matter

### Two profiles, two thresholds

PF12386 runs at `--cut_ga`, its curated Pfam gathering threshold. SSF54001 runs
at a bit score of 25. They cannot share an `hmmsearch` invocation, so each
profile gets its own job per batch.

`prepare_hmms.py` refuses to start if the config asks for `--cut_ga` on a model
with no GA line. The single-model downloads from InterPro often lack it while
the `Pfam-A.hmm` release carries it, and finding that out 700 array jobs later
is not fun.

A protein hitting only SSF54001 is retained and labelled `ssf_only`. It never
votes on where the triad columns are. Whether it clears the triad filter anyway
is reported in `triad_filter_by_tier.tsv` and is one of the more informative
numbers the pipeline produces: if `ssf_only` sequences pass at anything like the
`specific` rate, the columns are not discriminating fold from function.

### An empirical FDR, not a trusted E-value

Every batch is also searched against its own reversed sequences. Reversal keeps
length and composition, including the low-complexity and biased-composition
regions that inflate HMMER scores, and destroys homology. The rate at which
decoys clear a bit score is the rate at which chance clears it.
`decoy_fdr.tsv` and figure 20 are the defence of the cutoffs.

### 350k proteomes without 700,000 process start-ups

Samples are batched 500 at a time and each sequence renamed to a synthetic ID,
with its provenance written into the FASTA description:

```
>b00042_s00000123 NODE_7_length_4412_cov_9.2_11 SRR12345678
```

HMMER copies the description into the `description of target` column of the
domtblout, so a *hit* carries its own sample and protein ID and no billion-row
mapping table is ever written. Renaming rather than prefixing avoids the fact
that across 350k assemblies from mixed sources there is no delimiter guaranteed
absent from both sample and protein IDs.

The number of batches is not known until Prodigal has run, so `unified_table` is
a Snakemake **checkpoint**.

### Triad columns learned from PF12386 only, then applied to everything

`hmmalign` against PF12386 emits `#=GC RF`, marking profile match states. Those
columns are homologous across every sequence by construction, which turns "find
the catalytic triad" from an inference into a lookup.

Column conservation is measured **only on the PF12386 hits**, redundancy-weighted
(below). The winning C/H/D columns are then applied unchanged to every aligned
sequence, `ssf_only` included. Spacing (your C213/H248/D272, gaps of 35 and 24)
enters as a soft prior; conservation does the work. Alternatives are ranked in
`triad_candidates.tsv`, and `triad.override_columns` pins the answer once you
have looked.

### Redundancy weighting, applied consistently

Every sequence-level frequency in this pipeline — triad column conservation,
logos, entropy, consensus motifs, mutual information, subgroup sizes — is
weighted by `1 / (size of the sequence's 90%-identity cluster)`. A genus with
8,000 sequenced isolates otherwise dictates what the family's active site
"looks like". Weights sum to the number of clusters, so weighted counts read as
"effective number of independent sequences" and are reported as such.

Where a *test* needs exchangeable observations rather than a corrected point
estimate (the coupling permutation test, the Fisher taxon enrichment), the
pipeline dereplicates instead of weighting. Down-weighting fixes the estimate
but not the null distribution.

### Coupling: z-scores, not permutation ranks

With `B` permutations the smallest achievable rank-based p is `1/(B+1)`. After
Benjamini-Hochberg across ~435 column pairs, nothing can reach q = 0.05 until
`B` exceeds about 9,000. The pipeline therefore scores each pair as a z-score
against its own permutation null (mean and sd), which extracts resolution beyond
`1/B` from the same `B` permutations, and reports the raw empirical p alongside.
This is the standard treatment of APC-corrected MI in the coevolution
literature. The synthetic test has five positions perfectly determined by
subgroup; the rank-based version recovered 0 of their 10 pairs, the z-score
version recovers 10 of 10 with one false positive in the remaining 425.

### Convergence, measured rather than eyeballed

Minimum independent origins of each subgroup under Fitch parsimony on the C71
gene tree, against a tip-shuffle null (same tree, same prevalence, no phylogenetic
arrangement) and a Brownian-threshold null. Fritz & Purvis *D* rescales the
observed count between the two nulls so subgroups of different prevalence are
comparable.

*D* is a ratio of differences between nulls. On short, star-like trees the two
nulls coincide, the denominator collapses and *D* explodes. The pipeline reports
`D = NA, "Brownian and random nulls coincide"` rather than a spurious −16.

### SSN

EFI-EST convention: edge weight is `-log10(E)` from an all-vs-all DIAMOND search,
swept across thresholds, knee picked by the Kneedle rule. The knee search is
restricted to thresholds that leave at least one component with ≥ 2 nodes,
because a knee on a flat curve otherwise lands where every node is isolated. If
most nodes end up singletons anyway, the pipeline says so loudly.

`ssn_vs_subgroup_agreement.tsv` compares SSN clusters to the k-means subgroups
(ARI, AMI). The SSN assumes only that homologues align; the k-means assumes
convex clusters in a BLOSUM embedding; the tree assumes a substitution model.
Where all three agree the subgroup is real. Where they disagree, that is the
result.

### Detection and sampling bias, one model used everywhere

`prep_genome_table.py` builds one row per **searched** genome. Genomes that were
skipped are excluded, because "not searched" is not "no C71" and a logistic model
cannot tell the difference.

Covariates, the same set everywhere:

`completeness`, `contamination`, `log10_n50`, `log10_contigs`, `log10_n_proteins`.

An incomplete MAG loses genes roughly in proportion to what is missing. A
fragmented assembly splits genes across contig ends, so the fragments fail the
coverage filter. `n_proteins` is the raw number of chances to find a hit. These
quantify detection opportunity, not biology, which is why prevalence is reported
adjusted to a reference genome (completeness 100, contamination 0, median
assembly quality) rather than raw.

`phyloglm.R` fits `has_c71 ~ covariates` with `phylolm::phyloglm`
(`logistic_MPLE`, Ives & Garland) separately on the bacterial and archaeal GTDB
trees, plus per-subgroup models among C71-positive genomes. The matched
non-phylogenetic `glm()` is fitted and reported alongside so the cost of ignoring
the tree is visible. Fritz & Purvis *D* for C71 presence on the GTDB tree comes
from `caper::phylo.d`.

`phyloglm` builds the n × n phylogenetic covariance. Above `phyloglm.max_tips`
(default 8,000) the script draws outcome-stratified subsamples, fits each, and
pools with Rubin's rules — the between-replicate variance is added to the mean
within-replicate variance, so the pooled standard errors are wider than any
single fit's, as they should be. The fraction of missing information is
reported per coefficient.

Tip-label matching is checked before anything runs: `tree_tip_matching.tsv` gives
the matched fraction per domain and the pipeline aborts below
`tip_matching.min_match_fraction`. A bad join key is not a biological result.

---

## Verification

`test/` builds a miniature dataset with a known answer:

- 60 bacterial proteomes and 20 archaeal assemblies (real DNA, realistic codon
  usage and Shine-Dalgarno motifs, so Prodigal has something to train on)
- a **Pei family** with the triad planted at match columns 112 / 147 / 171 and
  three active-site subgroups with distinct flanking motifs
- a **fold family**: same structural template, 55% diverged, triad at 120 / 160 /
  200. It clears SSF54001 at score 25 and must fail PF12386's GA threshold
- Pei sequences with a deliberately broken catalytic residue (C→S, H→Y, D→E)
- a missing faa and an empty faa
- a PF12386 model **with** a GA line, an SSF54001 model **without** one
- GTDB-style trees and QC metadata for both domains

Profiles are built with `hmmbuild`; search and alignment run through `pyhmmer`,
which is the HMMER3 code. `test/stubs.py` stands in for prodigal, seqkit,
mmseqs, diamond and trimal.

```bash
pip install pyhmmer pyrodigal python-igraph
python test/make_testdata.py /tmp/c71_test
python test/run_test.py     /tmp/c71_test
```

Current result:

```
triad called at [112, 147, 171] (CHD); planted at [112, 147, 171]
learned from 83 'PF12386' sequences

evidence  n_aligned  n_triad_positive  frac_triad_positive
specific         83                76             0.916
ssf_only         56                 0             0.000

intact Pei kept  : 56/56
broken Pei kept  : 0/3
fold-family kept : 0/45

evidence tier by planted family:      fold: 0 specific, 45 ssf_only
                                      pei : 59 specific,  0 ssf_only
archaea: 20 genomes on the ar53 tree, 16 C71-positive after Prodigal

subgroups: k=3 (silhouette_cosine), gap statistic 3, adjusted Rand 1.000
coupling : 10/10 planted co-varying pairs recovered, 1 false positive in 425
RESULT: PASS
```

**`phyloglm.R` is not executed by the test** — there is no R in the test sandbox.
Its inputs are validated: every configured covariate is present and non-null,
every genome has a tree tip, and both domains reach the table with C71-positive
members. Run it once on a small subset before committing a cluster run.

---

## Layout

```
Snakefile                 rules, checkpoint, DAG, env activation
config.yaml               every knob, including envs_root
run.sh                    driver; picks conda vs pre-built envs from config
setup/
  build_envs.sh           STEP 1, on a machine with internet: solve + conda-pack
  install_envs.sh         STEP 2, on the cluster: unpack, conda-unpack, verify
profiles/sge/
  config.yaml             per-rule slots, memory, walltime, queue
  submit.sh               per-slot h_vmem, queue routing, no -pe when threads=1
  status.sh               qstat/qacct probe (needs SGE_ROOT sourced)
envs/                     hmmer, prodigal, phylo, network, py, r, selection (all nodefaults)
scripts/
  check_snakefile.py      static wiring + SGE profile check; no snakemake needed
  preflight.py            validate every input path, column, tip label, HMM
  search_batch.sh         fused concat + decoy + hmmsearch on node-local $TMPDIR
  prodigal_run.py         gene calling, with a meta-mode fallback
  build_sample_table.py   unify bacteria + archaea, QC metadata, tree tips, batches
  prepare_hmms.py         hmmconvert, per-profile files, GA-line validation
  batch_faa.py            concat + rename + provenance in the description
  make_decoys.py          reversed sequences for the empirical FDR
  combine_filter.py       parse domtblout, tier the evidence, decoy FDR, stats
  merge_metadata.py       left merge onto the sample table
  extract_seqs.py         pull sequences back from the original faa
  seq_weights.py          1 / cluster size at 90% identity
  triad_detect_filter.py  RF match states -> columns from PF12386 -> filter all
  tree_input.py           derep + MMseqs2 + trimAl
  active_site_analysis.py subgroups, logos, entropy, PCA (figures 10-15)
  coupling.py             APC-corrected MI, z-scores (figure 16)
  ssn_align.sh, ssn.py    all-vs-all + threshold sweep (figures 17-18)
  convergence.py          Fitch origins, tip-shuffle and Brownian nulls (figure 19)
  prep_genome_table.py    one row per searched genome, detection covariates
  phyloglm.R              phylogenetic logistic regression, adjusted prevalence, D
  plots_overview.py       figures 01-07, 20-22
  plot_tree.py            figures 08-09, rectangular and circular
  make_report.py          report.md
  utils.py                fasta/stockholm/domtblout I/O, plot style
  --- target specificity ---
  domain_arch.py          PMBR repeats + accessory binding modules (figure 23)
  module_trees.py         catalytic vs PMBR tree congruence, tanglegram
  cellwall_reference.py   Kandler & Koenig 1978, encoded; refuses genus guesses
  cellwall_genotype.py    Pmur markers, host P1 / P1' call
  groove_map.py           8JX4 (PeiW-CD) groove -> PF12386 match columns; identity guard
  pei_class.py            the published four-class partition (figure 28)
  sdp.py                  non-circular SDPs across 4 partitions (figures 24-25)
  selection.py            FEL/MEME, groove vs core (figure 26)
  assay_panel.py          which proteins to synthesise (figure 27)
test/
  make_testdata.py        synthetic screen dataset
  run_test.py             end-to-end screen check
  test_specificity.py     planted groove, planted SDPs, circularity trap
  test_cellwall_reference.py  every literature row, the genus refusal, disputed claims
  stubs.py                prodigal/seqkit/mmseqs/diamond/hmmsearch/trimal stubs
```

## Before you run it

1. `source` your SGE settings file so `SGE_ROOT` is set and `qstat`/`qacct` are
   on PATH. `profiles/sge/config.yaml` is already tuned from your probe report;
   check the queue names still match.
2. `chmod +x profiles/sge/{status,submit}.sh run.sh setup/*.sh scripts/*.sh`
3. `pip install snakemake-executor-plugin-cluster-generic`
4. `./run.sh preflight`, and fix every ERROR. It checks the GA line, the tree tip
   labels, the QC column names and a sample of the faa paths.
5. `./run.sh dry`, and read the job count.
6. After `combine_filter`, look at `decoy_fdr.tsv` and `hmm_search_stats.tsv`
   before letting the rest run.
7. After `triad`, read `triad_candidates.tsv` and `triad_filter_by_tier.tsv`.
   Pin the columns with `triad.override_columns` and re-run from that rule.
8. After `phyloglm`, check `tree_tip_matching.tsv` and the `fmi` column of
   `phyloglm_coefficients.tsv`. High fmi means the subsampling, not the data, is
   driving the standard error; raise `max_tips` or `n_replicates`.
