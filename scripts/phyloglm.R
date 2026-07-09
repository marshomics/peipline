#!/usr/bin/env Rscript
# ============================================================================
# Phylogenetic logistic regression of C71 presence on genome-quality covariates.
#
# Why phylogenetic. Genomes are not independent draws. Two Methanobacterium
# strains share C71 because they share an ancestor, not because they
# independently acquired it, and a plain glm() counts that as two pieces of
# evidence. Ives & Garland's phylogenetic logistic regression (phylolm::phyloglm,
# method = "logistic_MPLE") models the residual correlation with a single
# parameter alpha: large alpha means the phylogenetic signal has decayed and the
# model collapses towards ordinary logistic regression.
#
# The unit is a SPECIES, because both GTDB trees supplied are species-level.
# 342,759 genomes collapse onto ~9,800 bacterial tips and a few thousand archaeal
# ones. Taking one arbitrary genome per tip (which an earlier draft of this
# script did, via !duplicated) would have thrown away 97% of the data and let
# whichever genome sorted first decide the species' status. Instead:
#
#   has_c71      1 if ANY genome of the species carries a triad-positive protein
#   covariates   per-species means of completeness / contamination / log10 N50 /
#                log10 contigs / log10 n_proteins
#   log10_n_genomes  added automatically: a species represented by 400 genomes
#                has 400 chances to reveal a C71 that a singleton species does
#                not. Without this term, "prevalence" would largely measure how
#                hard each species has been sequenced.
#
# prevalence_within_species is written out so a proportion model can be fitted
# as a sensitivity analysis, but the primary model is Bernoulli, which is what
# phylolm::phyloglm actually fits.
#
# Why subsample. phyloglm builds the n x n phylogenetic covariance matrix.
# Above `phyloglm.max_tips` we draw outcome-stratified subsamples of the tree,
# fit each, and pool the coefficients with Rubin's rules (the between-replicate
# variance is added to the mean within-replicate variance). This is honest and
# it is reported: the pooled SEs are wider than any single fit's, as they should
# be. At ~9,800 bacterial species this barely bites.
#
# What comes out
#   phyloglm_coefficients.tsv    coefficient, SE, z, p, per model, plus the
#                                matched non-phylogenetic glm for comparison
#   bias_adjusted_prevalence.tsv per-taxon prevalence predicted at
#                                completeness = 100, contamination = 0 and median
#                                assembly quality: what you would have seen with
#                                perfect genomes
#   phylogenetic_signal_D.tsv    Fritz & Purvis D for C71 presence on the GTDB
#                                tree. D ~ 0: as clumped as Brownian. D ~ 1:
#                                random. Negative: more clumped than Brownian,
#                                which is what vertical inheritance of a
#                                prophage-borne gene in one order would look like.
# ============================================================================
suppressPackageStartupMessages({
  library(ape); library(phylolm); library(data.table); library(yaml)
})

args     <- commandArgs(trailingOnly = TRUE)
genomes  <- args[1]; bac_tree <- args[2]; ar_tree <- args[3]; cfgpath <- args[4]
out_coef <- args[5]; out_prev <- args[6]; out_d <- args[7]

cfg  <- yaml::read_yaml(cfgpath)
P    <- cfg$phyloglm
COV  <- unlist(P$covariates)
set.seed(P$seed)

dat <- fread(genomes, sep = "\t", na.strings = c("", "NA"))
message(sprintf("[phyloglm] %d searched genomes", nrow(dat)))

# ---------------------------------------------------------------------------
# Collapse genomes onto their tree tip (a species). Everything downstream works
# at this level.
aggregate_to_tips <- function(d) {
  sgcols <- grep("^sg_", names(d), value = TRUE)
  numcov <- intersect(COV, names(d))
  agg <- d[, c(
    list(n_genomes = .N,
         has_c71 = as.integer(any(has_c71 == 1)),
         has_specific = as.integer(any(has_specific == 1)),
         has_hit = as.integer(any(has_hit == 1)),
         prevalence_within_species = mean(has_c71),
         n_c71 = sum(n_c71),
         domain = domain[1],
         taxon = taxon[1]),
    lapply(setNames(numcov, numcov), function(x) mean(x, na.rm = TRUE)),
    lapply(setNames(sgcols, sgcols), function(x) as.integer(any(x == 1)))
  ), by = tree_tip]
  agg[, log10_n_genomes := log10(n_genomes)]
  agg[]
}

prep <- function(d, treefile) {
  if (!file.exists(treefile)) return(NULL)
  tr <- read.tree(treefile)
  d  <- d[!is.na(tree_tip)]
  if (!nrow(d)) return(NULL)

  n_gen <- nrow(d)
  d <- aggregate_to_tips(d)
  message(sprintf("  %d genomes -> %d species tips (%.1f genomes per tip)",
                  n_gen, nrow(d), n_gen / nrow(d)))

  d <- d[tree_tip %in% tr$tip.label]
  if (nrow(d) < 50) {
    message(sprintf("  only %d species on this tree; skipping", nrow(d)))
    return(NULL)
  }
  tr <- keep.tip(tr, d$tree_tip)
  # phylolm needs strictly positive, non-zero branch lengths
  if (is.null(tr$edge.length)) tr <- compute.brlen(tr)
  tr$edge.length[tr$edge.length <= 0] <- 1e-8

  # NaN, not NA, is what mean(na.rm=TRUE) returns for an all-missing species
  for (c in intersect(COV, names(d))) set(d, which(is.nan(d[[c]])), c, NA_real_)

  ok <- complete.cases(d[, ..COV])
  if (sum(!ok)) message(sprintf("  dropping %d species with missing covariates",
                                sum(!ok)))
  d <- d[ok]; tr <- keep.tip(tr, d$tree_tip)
  d <- d[match(tr$tip.label, tree_tip)]
  list(d = d, tr = tr)
}

# Rubin's rules across replicate fits.
pool <- function(fits) {
  terms <- Reduce(intersect, lapply(fits, function(f) rownames(f)))
  m <- length(fits)
  rbindlist(lapply(terms, function(tm) {
    est <- sapply(fits, function(f) f[tm, "Estimate"])
    se  <- sapply(fits, function(f) f[tm, "StdErr"])
    qbar <- mean(est); ubar <- mean(se^2); b <- if (m > 1) var(est) else 0
    tot  <- ubar + (1 + 1 / m) * b
    z    <- qbar / sqrt(tot)
    data.table(term = tm, estimate = qbar, std_error = sqrt(tot),
               z = z, p = 2 * pnorm(-abs(z)),
               odds_ratio = exp(qbar), n_replicates = m,
               fmi = ((1 + 1 / m) * b) / tot)
  }))
}

fit_one <- function(d, tr, response, label, domain) {
  y <- d[[response]]
  if (sum(y) < P$min_positives || sum(1 - y) < P$min_positives) {
    message(sprintf("  %s/%s: only %d positives; skipping", domain, label, sum(y)))
    return(NULL)
  }
  # Sampling effort per species is a detection-opportunity term, not a nuisance:
  # a species with 400 sequenced genomes has 400 chances to reveal a C71.
  allcov <- unique(c(COV, "log10_n_genomes"))
  allcov <- intersect(allcov, names(d))
  keep <- allcov[sapply(allcov, function(c) length(unique(d[[c]])) > 1)]
  if (!length(keep)) { message("  no varying covariates; skipping"); return(NULL) }
  form <- as.formula(paste(response, "~", paste(keep, collapse = " + ")))

  # plain glm, for the reader to see what ignoring the tree would have said
  g  <- summary(glm(form, data = d, family = binomial()))$coefficients
  gl <- data.table(term = rownames(g), estimate = g[, 1], std_error = g[, 2],
                   z = g[, 3], p = g[, 4], odds_ratio = exp(g[, 1]),
                   n_replicates = 1L, fmi = NA_real_,
                   model = "glm_no_phylogeny", response = label, domain = domain,
                   n_tips = nrow(d), n_positive = sum(y), alpha = NA_real_,
                   n_genomes = sum(d$n_genomes))

  ntip <- length(tr$tip.label)
  reps <- if (ntip > P$max_tips) P$n_replicates else 1L
  fits <- list(); alphas <- numeric(0)
  for (r in seq_len(reps)) {
    if (ntip > P$max_tips) {
      pos <- which(y == 1); neg <- which(y == 0)
      npos <- min(length(pos), floor(P$max_tips / 2))
      nneg <- min(length(neg), P$max_tips - npos)
      idx  <- c(sample(pos, npos), sample(neg, nneg))
      dd <- d[idx]; tt <- keep.tip(tr, dd$tree_tip); dd <- dd[match(tt$tip.label, tree_tip)]
    } else { dd <- d; tt <- tr }
    fit <- try(phyloglm(form, data = as.data.frame(dd), phy = tt,
                        method = "logistic_MPLE", btol = P$btol,
                        boot = 0), silent = TRUE)
    if (inherits(fit, "try-error")) { message("    phyloglm failed on a replicate"); next }
    s <- summary(fit)$coefficients
    colnames(s)[1:2] <- c("Estimate", "StdErr")
    fits[[length(fits) + 1]] <- s
    alphas <- c(alphas, fit$alpha)
  }
  if (!length(fits)) return(gl)

  pl <- pool(fits)
  pl[, `:=`(model = "phyloglm", response = label, domain = domain,
            n_tips = ntip, n_positive = sum(y), alpha = mean(alphas),
            n_genomes = sum(d$n_genomes))]
  message(sprintf("  %s/%s: %d species (%d positive, %d genomes), %d replicate(s), alpha=%.4g",
                  domain, label, ntip, sum(y), sum(d$n_genomes), length(fits), mean(alphas)))
  rbindlist(list(pl, gl), use.names = TRUE, fill = TRUE)
}

# ---------------------------------------------------------------------------
coefs <- list(); dstats <- list(); prevs <- list()

for (dom in c("Bacteria", "Archaea")) {
  tf <- if (dom == "Bacteria") bac_tree else ar_tree
  sub <- dat[domain == dom]
  if (!nrow(sub)) next
  message(sprintf("[phyloglm] %s: %d searched genomes", dom, nrow(sub)))
  pp <- prep(sub, tf)
  if (is.null(pp)) { message(sprintf("[phyloglm] %s: no usable tree overlap", dom)); next }
  d <- pp$d; tr <- pp$tr
  message(sprintf("[phyloglm] %s: %d species on the tree (%d genomes), %d C71-positive",
                  dom, nrow(d), sum(d$n_genomes), sum(d$has_c71)))
  fwrite(d, file.path(dirname(out_coef),
                      sprintf("species_level_table_%s.tsv", tolower(dom))), sep = "\t")

  r <- fit_one(d, tr, "has_c71", "has_c71", dom); if (!is.null(r)) coefs[[length(coefs)+1]] <- r
  r <- fit_one(d, tr, "has_specific", "has_specific", dom); if (!is.null(r)) coefs[[length(coefs)+1]] <- r

  if (isTRUE(P$run_subgroup_models)) {
    sgcols <- grep("^sg_", names(d), value = TRUE)
    dp <- d[has_c71 == 1]     # subgroup identity is only defined where C71 exists
    if (nrow(dp) >= 2 * P$min_positives) {
      trp <- keep.tip(tr, dp$tree_tip); dp <- dp[match(trp$tip.label, tree_tip)]
      for (sc in sgcols) {
        r <- fit_one(dp, trp, sc, sc, dom)
        if (!is.null(r)) coefs[[length(coefs) + 1]] <- r
      }
    } else message(sprintf("  %s: too few C71-positive genomes for subgroup models", dom))
  }

  # --- Fritz & Purvis D --------------------------------------------------
  if (requireNamespace("caper", quietly = TRUE) && sum(d$has_c71) >= 5) {
    dd <- as.data.frame(d[, .(tree_tip, has_c71)])
    cd <- try(caper::comparative.data(tr, dd, names.col = "tree_tip"), silent = TRUE)
    if (!inherits(cd, "try-error")) {
      pd <- try(caper::phylo.d(cd, binvar = has_c71, permut = 500), silent = TRUE)
      if (!inherits(pd, "try-error")) {
        dstats[[length(dstats) + 1]] <- data.table(
          domain = dom, trait = "has_c71", D = pd$DEstimate,
          p_random = pd$Pval1, p_brownian = pd$Pval0,
          n_tips = length(tr$tip.label), n_positive = sum(d$has_c71))
        message(sprintf("  D = %.3f (p vs random %.3g, p vs Brownian %.3g)",
                        pd$DEstimate, pd$Pval1, pd$Pval0))
      }
    }
  }

  # --- bias-adjusted prevalence ------------------------------------------
  # Refit with taxon as a fixed effect, then predict at reference genome
  # quality. This is the prevalence you would have observed had every genome
  # been complete, uncontaminated and contiguous.
  if ("taxon" %in% names(d) && sum(!is.na(d$taxon)) > 50) {
    dt <- d[!is.na(taxon)]
    tab <- dt[, .N, by = taxon][N >= 20]
    dt <- dt[taxon %in% tab$taxon]
    if (nrow(dt) > 50 && length(unique(dt$taxon)) > 1 && sum(dt$has_c71) >= 5) {
      allcov <- intersect(unique(c(COV, "log10_n_genomes")), names(dt))
      keep <- allcov[sapply(allcov, function(c) length(unique(dt[[c]])) > 1)]
      f <- as.formula(paste("has_c71 ~ taxon +", paste(keep, collapse = " + ")))
      m <- glm(f, data = dt, family = binomial())
      nd <- unique(dt[, .(taxon)])
      # Predict at a reference genome: complete, uncontaminated, median assembly
      # quality, and the median amount of sequencing effort.
      for (c in keep) nd[[c]] <- if (c %in% names(P$reference)) P$reference[[c]]
                                 else median(dt[[c]], na.rm = TRUE)
      pr <- predict(m, newdata = as.data.frame(nd), type = "link", se.fit = TRUE)
      raw <- dt[, .(n_species = .N, n_genomes = sum(n_genomes),
                    raw_prevalence = mean(has_c71)), by = taxon]
      nd <- as.data.table(nd)
      nd[, `:=`(adjusted_prevalence = plogis(pr$fit),
                lo = plogis(pr$fit - 1.96 * pr$se.fit),
                hi = plogis(pr$fit + 1.96 * pr$se.fit),
                domain = dom)]
      prevs[[length(prevs) + 1]] <- merge(nd, raw, by = "taxon")
    }
  }
}

emptydt <- function(cols) setNames(data.table(matrix(nrow = 0, ncol = length(cols))), cols)

fwrite(if (length(coefs)) rbindlist(coefs, use.names = TRUE, fill = TRUE)
       else emptydt(c("term","estimate","std_error","z","p","odds_ratio","model",
                      "response","domain","n_tips","n_positive","alpha","n_genomes")),
       out_coef, sep = "\t")
fwrite(if (length(prevs)) rbindlist(prevs, use.names = TRUE, fill = TRUE)
       else emptydt(c("taxon","adjusted_prevalence","lo","hi","domain","n_species",
                      "n_genomes","raw_prevalence")),
       out_prev, sep = "\t")
fwrite(if (length(dstats)) rbindlist(dstats, use.names = TRUE, fill = TRUE)
       else emptydt(c("domain","trait","D","p_random","p_brownian","n_tips","n_positive")),
       out_d, sep = "\t")

message("[phyloglm] done")
