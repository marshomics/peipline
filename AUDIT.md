# Pipeline audit — pseudomurein endo-isopeptidase (Pei) screen

Date: 2026-07-11. Scope: the whole pipeline (40 scripts, ~12,200 LOC, plus the
Snakefile, config, SGE profile, conda envs, setup scripts, and README), audited
against six axes: robust, defensible, reproducible, informative, useful,
publication-grade.

Method: five independent domain audits (search→triad core; phylogeny/statistics;
specificity/structure/synteny; reference/literature; reproducibility/reporting),
then triage, fix, and re-verification. Every literature accession, PDB ID, residue
number, and citation was checked against UniProt and Crossref. Two citation
corrections were confirmed against the Crossref record for the DOI.

Bottom line: two guaranteed run-ending crashes and one core citation error were
found and fixed, along with a systematic C39-into-C71 pooling problem introduced
when the C39 arm was enabled by default. Thirty-odd defects were fixed in total;
a set of design-level statistical limitations is documented below as known
limitations to state in the methods or address before submission. After the fixes,
all 13 test suites (12 static + a new audit-regression suite) and both end-to-end
suites pass, and the Snakefile wiring check is clean.

---

## Severity summary

| Severity | Found | Fixed | Documented / deferred |
|---|---|---|---|
| Critical | 1 | 1 | 0 |
| High | 8 | 8 | 0 |
| Medium | 19 | 15 | 4 |
| Low | 15 | 8 | 7 |

"Documented / deferred" are design-level items (mostly statistical) where a code
change would be a research decision, not a bug fix. Several were already disclosed
in code comments; they are collected here so they can be stated as limitations or
scheduled deliberately.

---

## Critical

**C1 — Core experimental citation misattributed.** The substrate / lysis / kinetics
ground truth was cited as "Subedi et al. 2015" across eight modules. The Crossref
record for DOI 10.1155/2015/828693 gives the authors as **Schofield LR, Beattie AK,
Tootill CM, Dey D, Ronimus RS (2015)**, "Biochemical Characterisation of Phage
Pseudomurein Endoisopeptidases PeiW and PeiP Using Synthetic Peptides", *Archaea*
2015:828693. The title, journal, volume, and DOI were correct; the author list was
not. A referee checking the reference would find it does not exist under "Subedi".
**Fixed**: replaced throughout (lysis_reference, pei_reference, assay_panel,
lysis_check, groove_map, cellwall_reference, make_report, pmbr_reference, sdp), with
a note recording that earlier drafts miscited it. Please confirm the corrected
attribution against your own records — the change was driven by the Crossref data
for the DOI you supplied.

---

## High

**H1 — module_trees.py crashed on every real run.** The tanglegram title formatted
an undefined name `p`; it exists only inside two helper functions, never in `main()`.
Any run with ≥20 shared tips reached it and raised `NameError`, failing the rule
(the results table was written first, so numbers survived, the figure did not).
The intended value was the tip-shuffle RF p, not a Mantel p (the Mantel p is
deliberately not reported). **Fixed**: uses `p_rf`, relabelled so the Mantel r is
marked descriptive and the significance statistic is the RF null. Regression-tested.

**H2 — cellwall_genotype.py crashed at the end of every run.** The bacteria
cross-reaction check referenced `pmur_pathway`, a column that is never created (the
real name is `pmur_count_pathway`), raising `KeyError` after the output table was
written. Invisible to the test battery because `main()` was never exercised.
**Fixed**: correct column name. Regression-tested; a `main()`-level smoke path is
now covered.

**H3 — triad filter could abort the run on a low-coverage fragment, and could
also let one into c71.faa.** A sequence carrying C/H/D at the three triad columns
but covering < `min_match_coverage` of the profile satisfied the positive mask
(computed independently of coverage) and tripped an assertion that declared the
situation "impossible" — it is not, since coverage counts all match states and the
triad is only three. On real data (the SSF54001-only tier is exactly the population
at risk) this is a data-dependent hard stop; deleting the assertion would instead
have written a fragment into c71.faa. **Fixed**: a positive now requires the triad
residues AND adequate coverage AND no gap at the triad columns; the fragment is
labelled `low_coverage` and kept out of c71.faa. Regression-tested.

**H4 — C39 hits were pooled into the C71 screening funnel and overview figures.**
With PF03412 enabled by default, `combine_filter` counted every PF03412-only hit
(mostly bacteriocin exporters) as C71 `ssf_only`, inflating `unique_proteins_with_hit`,
`proteins_ssf_only`, and `samples_with_hit`, and the overview figures drew a PF03412
series inside the C71 panels (also silently dropping the two-profile concordance
panel). This undercut the careful C39/C71 separation elsewhere. **Fixed**: the
funnel is scoped to the C71-family profiles (PF12386 + shared SSF54001); PF03412 is
counted separately and reported only in the C39 section; the report funnel is
relabelled "C71 screening funnel"; `plots_overview` filters to the C71 family. All
no-ops when PF03412 is disabled.

**H5 — C39 arm rules would be killed on the cluster, blocking the whole run.** The
SGE profile predates the C39 arm, so `hmmalign_c39`, `seq_weights_c39`, etc. fell to
default resources (8 GB / 4 h / standard.q / 1 slot). `hmmalign_c39` streams every
PF03412 hit and would exceed the 24 h standard.q cap; because `report` depends on
`triad_c39`, a C39 failure fails `rule all` after potentially days of compute.
**Fixed**: added `set-threads` / `set-resources` for all five C39 rules mirroring
their C71 twins (`hmmalign_c39` routed to long.q).

**H6 — the report could silently drop whole sections and go stale.** `make_report`
reads several specificity tables by path that the `report` rule did not declare as
inputs. With `keep-going: true` and a load function that swallows exceptions, a
failed `domain_arch` / `cellwall_genotype` / `lysis_check` would still produce a
report — minus those sections, with exit 0 — and editing a marker would not rebuild
the report. **Fixed**: `report` now depends on the whole specificity block, so a
missing table is a hard error and a changed one triggers a rebuild.

**H7 — the report's version block was uniformly "not found".** `make_report` runs
in the `py` conda env, which has none of the bioinformatics binaries, so every
`--version` probe failed and the provenance block implied the tools were absent.
**Fixed**: versions are read from the pinned `envs/*.yaml`, and a provenance block
now records the seeds, the config path, and the C39 arm status. `dist/locks/*.txt`
holds the exact solved builds.

**H8 — phyloglm reported many coefficient p-values with no multiple-testing
correction.** The `has_c71`, `has_specific`, and every `sg_*` subgroup model, across
two domains, produced dozens of raw p-values with no correction, while the two
Python statistics scripts correct theirs. **Fixed**: a Benjamini-Hochberg `q_bh` is
computed within the phylogenetic-model hypothesis family (subgroup / response terms,
excluding intercepts and the detection-bias covariates), surfaced in the report,
and the family is stated. Raw p is kept visible.

---

## Medium — fixed

- **M1 PeiR "nearest histidine" claim was false.** The narrative said the nearest
  downstream histidine to the catalytic C90 is H162 (gap 72); the module's own
  `peir_his_candidates()` returns H104 (gap 14) as nearest. The valid point is that
  no histidine sits at the C71 distance of 35, and PeiR's catalytic His is
  unassigned. Reworded in pei_reference and pei_check to match the recomputed gaps.
- **M2 PeiW/PeiP citation wrong.** "Steenbakkers et al. 2001, FEMS Microbiol Lett
  208:47-53" is, per Crossref, **Luo et al. 2002, 208:47-51**. Corrected.
- **M3 active-site subgroups double-corrected for redundancy.** k-means and
  k-selection were weighted by 1/cluster_size on rows that were already dereplicated
  to one representative per cluster, so a large cluster's representative nearly
  dropped out and the clustering geometry disagreed with the (unweighted) PCA basis.
  Now uniform weight on dereplicated representatives; the 1/size weight remains only
  as the fallback when no cluster column is available.
- **M4 non-deterministic triad call.** Candidate-column selection used an unstable
  argsort, so ties near the `max_candidates` cap could change the chosen catalytic
  triad across platforms. Now a stable sort.
- **M5 `--allow-empty-specific` did not cover all abort paths.** Three downstream
  exits (no candidate columns, no triple, zero retained) could still crash the C39
  arm on a near-empty net. All now route through the empty-output writer.
- **M6 join hazard in combine_filter.** The batch maps were read without
  `dtype={"sample": str}`, so a numeric or leading-zero sample ID would mismatch the
  hit table and abort the merge. Pinned.
- **M7 synteny mislabelled unknown contiguity.** A split block on a genome with an
  unknown contig count returned `dispersed` (a biological statement) instead of
  `not_evaluable`. Unknown contiguity is now `not_evaluable`.
- **M8 marker search floor could truncate strict hits.** The cellwall marker search
  floored at `min(perm, strict)`; a model whose GA line sits below that floor would
  lose hits between its GA and the floor. Now floors at the min of the base
  threshold and every model's GA.
- **M9 prodigal cache overwrote a covariate.** A resumed run recorded
  `prodigal_mode = "cached"` instead of the real gene-call mode, making a
  detection-bias covariate depend on execution history. The real mode is now
  persisted in the `.done` marker and read back.
- **M10 convergence statistics were mislabelled.** `observed_origins` is the Fitch
  gains+losses count (an upper bound on origins, not the origin count), and
  `fritz_purvis_D` is a home-rolled scaled index, not caper's Fritz-Purvis D (which
  the pipeline computes separately). Renamed to `parsimony_changes` and
  `clustering_index`, with the distinction stated in the report.
- **M11 sdp had no NaN-weight guard.** A sequence missing a redundancy weight would
  inject NaN into every weighted MI. Now fails with a clear message.
- **M12 SSN figure layout unseeded.** The force-directed layout used igraph's
  unseeded RNG, so node coordinates changed run to run. Now seeded (numbers were
  already deterministic; only the figure was not).
- **M13 phantom subgroup column.** `prep_genome_table` iterated `subgroup.unique()`
  including NaN, which could create an `sg_nan` column that phyloglm then treated as
  a real subgroup. Now drops NaN.
- **M14 missing QC columns were silent.** Despite a "fails loudly" docstring, a
  missing QC column was nulled without warning, degrading a regression covariate.
  Now prints a prominent warning naming the columns and the affected covariates.
- **M15 assay_panel dropped a column and used dead logic.** `pmbr_rule_applies` was
  referenced but never mapped from the architecture table (which does emit it), so
  the declared preference was dead and the output column was dropped. Now mapped.

## Medium — documented (address or state as a limitation)

- **M16 Rubin pooling is anti-conservative for the bacterial headline effect.** With
  ~9,800 bacterial tips > `max_tips` = 8,000, the model is fit on 20 outcome-
  stratified subsamples that reuse every positive, so the between-replicate variance
  underestimates true sampling variance and the pooled SE for bacterial `has_c71` is
  too small. This is avoidable here: the ~9,800×9,800 covariance is ≈0.77 GB, so
  raising `max_tips` above the bacterial tip count fits the full tree once with no
  subsampling. Recommended before submission. (phyloglm.R)
- **M17 coupling permutation null assumes sequence exchangeability.** The column
  shuffle destroys all row structure; 90% dereplication and APC reduce but do not
  remove residual phylogenetic covariation between two columns, so MI p-values may be
  anti-conservative and the "co-adapted surface" claim can be inflated. The
  structural-contact enrichment is the intended guard and should be reported
  prominently; a phylogeny-aware sensitivity null (simulate columns on the tree, or
  DCA-style reweighting) would settle it. (coupling.py)
- **M18 coupling significance rests on Gamma-tail extrapolation.** At B = 1000
  permutations the empirical floor is ~1e-3, but BH over the barcode pairs needs
  ~1e-4, so every significant call comes from the moment-matched Gamma tail
  extrapolated below what the permutations can corroborate. Recommend a tail
  goodness-of-fit check, or raising B for the handful of survivors. (coupling.py)
- **M19 convergence and module-tree congruence condition on a single ML tree.** The
  tip-shuffle and RF nulls hold one point-estimate topology fixed, so p-values are
  valid conditional on the tree but ignore gene-tree and rooting uncertainty.
  Recommend recomputing origin counts / RF across the UFBoot replicates and
  reporting the distribution. (convergence.py, module_trees.py)

---

## Low — fixed

- **L1** HyPhy FEL/MEME site count is now asserted equal to the match-column count
  before the positional site→column mapping, closing a silent coordinate-drift path.
  (selection.py)
- **L2** `hmmscan` now passes `--domE` alongside `-E`, so a protein with one strong
  PMBR domain but a weak full-sequence E-value is not dropped at the load-bearing
  3-motif cliff. (domain_arch.py)
- **L3** FASTA and domtbl readers now decode with `errors="replace"`, so a stray
  non-ASCII byte in a header cannot crash a 350k-proteome read on a LANG=C node.
  (utils.py)
- **L4** Setup-script comments corrected from "six" to "seven" environments.
- **L5–L8** Report provenance additions: seeds, config path, C39 status, pinned
  versions (folded into H7).

## Low — documented

- **L9** The decoy FDR search is thresholded, so the FDR curve starts at the
  operating cutoff and cannot inform the choice of that cutoff (only report FDR at
  it). An unthresholded or low-`-T` decoy pass would let the curve span below the
  operating point. The decoy is also a single fixed permutation (n=1 per target),
  so low-count FDR estimates are noisy. (search_batch.sh / combine_filter.py)
- **L10** Reported triad columns are 0-based match-column indices; a reader
  comparing against 1-based Pfam match states or structure residues is off by one.
  Internally consistent (all downstream code uses the same indices); worth a label.
- **L11** `sdp` counts the global-null `tree_clade` partition toward the ≥2-partition
  replication rule alongside the correlated `ssn_cluster`, so a purely
  phylogenetically-clustered column can reach "replicated" with one uncontrolled and
  one correlated partition. Mitigated by the within-clade nulls and the concordance
  matrix; consider excluding the global-null partition from the count.
- **L12** `pei_class` vs subgroup ARI is partially circular: the subgroup barcode
  includes class position 252, so agreement on that axis is partly guaranteed. State
  it as an upper bound, or mask the class columns before comparing.
- **L13** The non-phylogenetic comparison `glm` has no separation guard and its odds
  ratios lack CIs; the `z` column holds a Wald z for glm rows and a Rubin t for
  pooled rows. Comparison-only, but worth splitting the column.
- **L14** The optional `structure` env is not built by the setup scripts (documented
  as optional); enabling `structure_search` under the pre-built-env route would fail
  late. Add it to the setup list if you intend to run the structure stack.
- **L15** Latent, dormant in the current config: the archaeal stem↔accession join,
  the metadata-key autodetect min-overlap, and the Newick underscore convention are
  not validated. They do not fire with the configured inputs but would be silent
  mis-joins if inputs change.

---

## Not defects (design choices confirmed sound)

The family-aware C71/C39 split does not regress the C71 default path; the evidence
tiering uses list membership (not a fragile substring); the dN/dS CDS→codon mapping
fix is real and regression-tested; the `structure_expect` identity guard prevents
the 8JX4/8Z4F swap; taxonomy-first PM calling, permissive-out-of-order gating (a
permissive marker never elevates without positive synteny), and PMBR partition
nesting are correct; conda channel hygiene (conda-forge/bioconda/nodefaults, strict
priority, --override-channels) and the per-slot h_vmem arithmetic are sound; every
stochastic step is seeded from config. phyloglm.R, selection.py, and
structure_search.py are not executed in development (no R / HyPhy / GPU / Foldseek in
the sandbox) — they are guarded statically and honestly labelled.

---

## Verification

- 12 static suites + 1 new audit-regression suite: **all pass**.
- Two end-to-end suites (C71 screen; specificity block) on synthetic data with
  stubbed HMMER/tree tools: **both pass**.
- Snakefile wiring check: **OK** (45 rules).
- The new `test/test_audit_regressions.py` locks the two crashes, the triad
  low-coverage gating, and the convergence renames.

Standing items that can only be run on the real cluster: phyloglm.R, selection.py
(HyPhy), and the structure stack (GPU / Foldseek / HHsuite). These were audited
statically.
