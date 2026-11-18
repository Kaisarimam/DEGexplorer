# =============================================================================
# dtu.R -- differential transcript usage (DTU)
#
# Does the share of a gene's expression carried by each of its transcripts
# change between conditions? Five methods, each run the way its authors
# recommend:
#
#   DEXSeq            Love, Soneson & Patro 2018 (F1000Research 7:952): scaledTPM
#                     counts, DRIMSeq's filter, ~sample + exon + condition:exon,
#                     per-gene q-values (perGeneQValue), stageR.
#   DRIMSeq           same workflow: Dirichlet-multinomial, dmPrecision / dmFit /
#                     dmTest, the post-hoc proportion-SD filter, stageR.
#   satuRn            Gilis et al. 2022 (vignette): filterByExpr, singletons removed,
#                     quasi-binomial fitDTU / testDTU, empirical p-values,
#                     per-gene q-values, stageR.
#   edgeR diffSplice  Baldoni et al. 2025 (NAR, gkaf1305): counts divided by their
#   limma diffSplice  read-to-transcript ambiguity, filterByExpr, normLibSizes,
#                     glmQLFit (edgeR) or voomLmFit (limma), diffSplice, gene-level
#                     Simes (or F) p-values.
#
# The dtu_* functions below the "core" line depend only on base R and the
# method packages: the generated reproduce.R copies them verbatim, so the
# script and the app run the same code.
# =============================================================================

DTU_ENGINES <- c("edgeR diffSplice" = "edger_dtu", "limma diffSplice" = "limma_dtu",
                 "DEXSeq" = "dexseq", "DRIMSeq" = "drimseq", "satuRn" = "saturn")

DTU_ENGINE_PACKAGES <- list(edger_dtu = "edgeR", limma_dtu = c("limma", "edgeR"),
                            dexseq = "DEXSeq", drimseq = "DRIMSeq",
                            saturn = c("satuRn", "SummarizedExperiment", "BiocParallel"))

#' Display name of any engine key.
engine_label <- function(engine) {
  all <- c(DE_ENGINES, DTU_ENGINES)
  hit <- names(all)[match(engine, all)]
  if (is.na(hit)) engine else hit
}

is_dtu_engine <- function(engine) engine %in% DTU_ENGINES

#' Count type each method's workflow uses ("auto" resolved), given what is loaded.
#' @return list(mode, note)
dtu_resolve_counts <- function(engine, requested = "auto", has_lengths = FALSE, has_od = FALSE) {
  want <- if (is.null(requested) || identical(requested, "auto")) {
    if (engine %in% c("edger_dtu", "limma_dtu")) "divided" else "scaledTPM"
  } else requested
  note <- NULL
  if (want %in% c("scaledTPM", "dtuScaledTPM") && !has_lengths) {
    note <- sprintf(paste("%s counts need each sample's effective transcript lengths, which are not loaded,",
                          "so the estimated counts were used as they are. Read the Salmon / kallisto / RSEM folders",
                          "(or load an effective-lengths file) to scale them."), want)
    want <- "estimated"
  }
  if (identical(want, "divided") && !has_od) {
    note <- paste("No quantification resamples are loaded, so the counts could not be divided by their",
                  "read-to-transcript ambiguity and were used as they are. Read the Salmon / kallisto folders",
                  "quantified with Gibbs or bootstrap samples to add it (Baldoni et al. 2025).")
    want <- "estimated"
  }
  list(mode = want, note = note)
}

DTU_COUNT_LABELS <- c(
  divided = "estimated counts divided by their quantification overdispersion (edgeR catchSalmon / catchKallisto; Baldoni et al. 2025)",
  scaledTPM = "scaledTPM: TPM scaled up to each sample's library size (tximport countsFromAbundance = \"scaledTPM\"; Love et al. 2018)",
  dtuScaledTPM = "dtuScaledTPM: TPM times the median transcript length of the gene, scaled to library size (tximport)",
  estimated = "estimated counts as quantified")

#' Filter method each workflow uses ("auto" resolved).
dtu_resolve_filter <- function(engine, requested = "auto") {
  if (is.null(requested) || identical(requested, "auto"))
    return(if (engine %in% c("dexseq", "drimseq")) "drimseq" else "edger")
  requested
}

# ---- core (copied into reproduce.R) ---------------------------------------------------

#' A function exported by the first of `pkgs` that has it (S4 generics such as
#' estimateSizeFactors are exported by different packages in different releases).
dtu_fn <- function(name, pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE) && name %in% getNamespaceExports(p))
      return(getExportedValue(p, name))
  }
  stop("None of ", paste(pkgs, collapse = ", "), " provides ", name, "().")
}

#' Counts for a DTU analysis.
#'
#' scaledTPM / dtuScaledTPM follow tximport's makeCountsFromAbundance: abundance
#' (counts / effective length; for dtuScaledTPM times the median length of the
#' gene's transcripts) rescaled to each sample's original count total.
#' @param lengths transcripts x samples effective lengths (aligned, no NA)
#' @param od      named overdispersion per transcript (>= 1)
dtu_counts <- function(counts, mode, lengths = NULL, od = NULL, gene = NULL) {
  if (identical(mode, "divided")) {
    o <- unname(od[rownames(counts)]); o[is.na(o)] <- 1
    return(counts / pmax(o, 1))
  }
  if (mode %in% c("scaledTPM", "dtuScaledTPM")) {
    len <- lengths[rownames(counts), colnames(counts), drop = FALSE]
    ab <- ifelse(len > 0, counts / len, 0)
    if (identical(mode, "dtuScaledTPM")) {
      g <- ifelse(is.na(gene) | gene == "", rownames(counts), gene)
      med <- tapply(rowMeans(len), g, stats::median)
      ab <- ab * as.numeric(med[g])
    }
    tot <- colSums(counts); s <- colSums(ab)
    out <- t(t(ab) * ifelse(s > 0, tot / s, 0))
    dimnames(out) <- dimnames(counts)
    return(out)
  }
  counts
}

#' DRIMSeq's dmFilter (dmDS_filter), as used by Love et al. 2018: a gene needs
#' `min_gene_expr` counts in `min_samps_gene_expr` samples; a transcript needs
#' `min_feature_expr` counts in `min_samps_feature_expr` samples and a share of at
#' least `min_feature_prop` of its gene in `min_samps_feature_prop` samples; genes
#' are kept only with two or more transcripts left.
#' @return logical keep per row of counts
dtu_filter_drimseq <- function(counts, gene, min_samps_gene_expr, min_gene_expr = 10,
                               min_samps_feature_expr, min_feature_expr = 10,
                               min_samps_feature_prop, min_feature_prop = 0.1) {
  keep <- stats::setNames(logical(nrow(counts)), rownames(counts))
  idx <- split(seq_len(nrow(counts)), gene)
  for (rows in idx[lengths(idx) > 1]) {
    ef <- counts[rows, , drop = FALSE]
    if (sum(ef, na.rm = TRUE) == 0) next
    if (!sum(colSums(ef) >= min_gene_expr, na.rm = TRUE) >= min_samps_gene_expr) next
    f2k <- rowSums(ef > 0, na.rm = TRUE) > 0
    if (sum(f2k) <= 1) next
    ef <- ef[f2k, , drop = FALSE]; rows <- rows[f2k]
    f2k <- rowSums(ef >= min_feature_expr, na.rm = TRUE) >= min_samps_feature_expr
    if (sum(f2k) <= 1) next
    ef <- ef[f2k, , drop = FALSE]; rows <- rows[f2k]
    s2k <- colSums(ef) > 0 & !is.na(ef[1, ])
    if (sum(s2k) < max(1, min_samps_feature_prop)) next
    prop <- prop.table(ef[, s2k, drop = FALSE], 2)
    f2k <- rowSums(prop >= min_feature_prop) >= min_samps_feature_prop
    if (sum(f2k) <= 1) next
    keep[rows[f2k]] <- TRUE
  }
  keep
}

#' Drop transcripts whose gene has fewer than two transcripts left.
dtu_drop_singletons <- function(keep, gene) {
  n <- table(gene[keep])
  keep & !is.na(gene) & gene != "" & gene %in% names(n)[n >= 2]
}

#' The filter each method's workflow uses.
#' @param method "drimseq" or "edger"; engine decides filterByExpr's arguments
#' @return list(keep, description)
dtu_filter <- function(counts, gene, group, design, method, engine, opts) {
  gene[is.na(gene) | gene == ""] <- rownames(counts)[is.na(gene) | gene == ""]
  n <- ncol(counts); n_small <- min(table(droplevels(as.factor(group))))
  if (identical(method, "drimseq")) {
    keep <- dtu_filter_drimseq(counts, gene, min_samps_gene_expr = n, min_gene_expr = opts$min_gene_expr,
                               min_samps_feature_expr = n_small, min_feature_expr = opts$min_feature_expr,
                               min_samps_feature_prop = n_small, min_feature_prop = opts$min_feature_prop)
    return(list(keep = keep, description = sprintf(paste(
      "DRIMSeq dmFilter: gene >= %g counts in all %d samples; transcript >= %g counts and >= %g%% of its",
      "gene in at least %d samples; genes keep two or more transcripts"),
      opts$min_gene_expr, n, opts$min_feature_expr, 100 * opts$min_feature_prop, n_small)))
  }
  if (identical(engine, "saturn")) {
    # satuRn vignette: single-isoform genes out, filterByExpr, single-isoform genes out again
    pre <- dtu_drop_singletons(rep(TRUE, nrow(counts)), gene)
    keep <- pre
    keep[pre] <- edgeR::filterByExpr(counts[pre, , drop = FALSE], design = NULL, group = group,
                                     lib.size = NULL, min.count = 10, min.total.count = 30,
                                     large.n = 20, min.prop = 0.7)
    keep <- dtu_drop_singletons(keep, gene)
    return(list(keep = keep, description = paste(
      "edgeR filterByExpr with satuRn's settings (min.count 10, min.total.count 30, large.n 20,",
      "min.prop 0.7); genes with a single transcript removed before and after")))
  }
  if (engine %in% c("edger_dtu", "limma_dtu")) {
    keep <- edgeR::filterByExpr(edgeR::DGEList(counts), design = design)
    return(list(keep = keep, description = paste(
      "edgeR filterByExpr (defaults, on the design); genes with a single transcript are skipped by diffSplice")))
  }
  keep <- edgeR::filterByExpr(counts, group = group)
  list(keep = dtu_drop_singletons(keep, gene), description = paste(
    "edgeR filterByExpr (defaults, by group); genes with a single transcript removed"))
}

#' Cell-means design ~ 0 + main + covariates and a target-vs-reference contrast.
dtu_design <- function(cd, main, covariates) {
  stats::model.matrix(stats::as.formula(paste("~ 0 +", paste(c(main, covariates), collapse = " + "))), cd)
}
dtu_contrast <- function(design, main, target, reference) {
  v <- stats::setNames(numeric(ncol(design)), colnames(design))
  v[paste0(main, target)] <- 1; v[paste0(main, reference)] <- -1
  v
}

#' DEXSeq's per-gene q-value (perGeneQValueExact): the FDR for calling a gene
#' with at least one differentially used transcript, from the minimum p-value of
#' its transcripts, accounting for the number of transcripts per gene.
#' @return named numeric, one q per gene
dtu_per_gene_qvalue <- function(p, gene) {
  ok <- !is.na(p)
  geneID <- factor(gene[ok]); pv <- p[ok]
  split_i <- split(seq_along(geneID), geneID)
  pGene <- vapply(split_i, function(i) min(pv[i]), numeric(1))
  theta <- unique(sort(pGene))
  tab <- tabulate(lengths(split_i)); nz <- which(tab > 0)
  numerator <- vapply(theta, function(t) sum(tab[nz] * (1 - (1 - t)^nz)), numeric(1))
  bins <- cut(pGene, breaks = c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  denom <- cumsum(tabulate(bins, nbins = nlevels(bins)))
  q <- pmin(1, (numerator / denom)[match(pGene, theta)])
  stats::setNames(q, names(split_i))
}

#' Mean share of each transcript within its gene (among the transcripts tested)
#' in the reference and target samples.
dtu_props <- function(counts, gene, ref_cols, tgt_cols) {
  tot <- rowsum(counts, gene, reorder = FALSE)[gene, , drop = FALSE]
  pr <- counts / tot
  pr[!is.finite(pr)] <- NA
  data.frame(prop_reference = rowMeans(pr[, ref_cols, drop = FALSE], na.rm = TRUE),
             prop_target = rowMeans(pr[, tgt_cols, drop = FALSE], na.rm = TRUE))
}

#' Mean count after library-size scaling.
dtu_base_mean <- function(counts, cols) {
  lib <- colSums(counts); sf <- lib / exp(mean(log(lib)))
  rowMeans(t(t(counts[, cols, drop = FALSE]) / sf[cols]))
}

#' Stage-wise DTU testing (stageR, method "dtu"): genes are screened on the
#' gene-level p-value (BH, or q-values already adjusted), then transcripts are
#' confirmed within significant genes. Adds stageR_gene_padj / stageR_tx_padj.
dtu_stagewise <- function(df, alpha, screen_adjusted) {
  g <- unique(df$parent_gene_id)
  sp <- if (screen_adjusted) df$gene_padj else df$gene_pvalue
  screen <- stats::setNames(sp[match(g, df$parent_gene_id)], g)
  screen[is.na(screen)] <- 1
  pc <- df$pvalue; pc[is.na(pc)] <- 1
  conf <- matrix(pc, ncol = 1, dimnames = list(df$gene_id, "transcript"))
  obj <- stageR::stageRTx(pScreen = screen, pConfirmation = conf, pScreenAdjusted = screen_adjusted,
                          tx2gene = data.frame(transcript = df$gene_id, gene = df$parent_gene_id,
                                               stringsAsFactors = FALSE))
  obj <- suppressMessages(suppressWarnings(stageR::stageWiseAdjustment(obj, method = "dtu", alpha = alpha,
                                                                       allowNA = TRUE)))
  adj <- as.data.frame(suppressMessages(suppressWarnings(
    stageR::getAdjustedPValues(obj, order = FALSE, onlySignificantGenes = FALSE))))
  i <- match(df$gene_id, adj$txID)
  df$stageR_gene_padj <- suppressWarnings(as.numeric(adj$gene[i]))
  df$stageR_tx_padj <- suppressWarnings(as.numeric(adj$transcript[i]))
  df
}

# -- engines: each returns list(tables = list(id = data.frame), qc, notes) ------------

#' DEXSeq on one comparison's samples (Love et al. 2018).
dtu_engine_dexseq <- function(counts, cd, main, covariates, gene, target, reference) {
  sd <- cd[, c(covariates, main), drop = FALSE]
  sd[[main]] <- factor(as.character(sd[[main]]), levels = c(reference, target))
  rownames(sd) <- colnames(counts)
  terms <- paste0(c(covariates, main), ":exon")
  full <- stats::as.formula(paste("~ sample + exon +", paste(terms, collapse = " + ")))
  reduced <- stats::as.formula(paste("~ sample + exon", if (length(covariates))
    paste("+", paste(paste0(covariates, ":exon"), collapse = " + ")) else ""))
  dxd <- DEXSeq::DEXSeqDataSet(countData = round(counts), sampleData = sd, design = full,
                               featureID = rownames(counts), groupID = gene)
  dxd <- dtu_fn("estimateSizeFactors", c("DEXSeq", "BiocGenerics", "DESeq2"))(dxd)
  dxd <- dtu_fn("estimateDispersions", c("DEXSeq", "BiocGenerics", "DESeq2"))(dxd, formula = full, quiet = TRUE)
  dxd <- DEXSeq::testForDEU(dxd, fullModel = full, reducedModel = reduced)
  dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = main, denominator = reference)
  dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
  q <- DEXSeq::perGeneQValue(dxr)
  lfc_col <- grep("^log2fold_", colnames(dxr), value = TRUE)[1]
  tab <- data.frame(gene_id = as.character(dxr$featureID), parent_gene_id = as.character(dxr$groupID),
                    baseMean = as.numeric(dxr$exonBaseMean),
                    log2FoldChange = as.numeric(dxr[[lfc_col]]), stat = as.numeric(dxr$stat),
                    pvalue = as.numeric(dxr$pvalue), padj = as.numeric(dxr$padj),
                    gene_pvalue = NA_real_, gene_padj = as.numeric(q[as.character(dxr$groupID)]),
                    dispersion = as.numeric(dxr$dispersion), stringsAsFactors = FALSE)
  qc <- tryCatch({
    mc <- S4Vectors::mcols(dxd)
    pick <- function(...) { for (n in c(...)) if (n %in% colnames(mc)) return(as.numeric(mc[[n]])); rep(NA_real_, nrow(mc)) }
    list(size_factors = stats::setNames(as.numeric(DESeq2::sizeFactors(dxd))[seq_len(ncol(counts))], colnames(counts)),
         dispersion = data.frame(baseMean = pick("exonBaseMean", "baseMean"), gene_est = pick("dispGeneEst", "dispBeforeSharing"),
                                 fit = pick("dispFit", "dispFitted"), final = pick("dispersion"),
                                 outlier = if ("dispOutlier" %in% colnames(mc)) mc$dispOutlier %in% TRUE else FALSE))
  }, error = function(e) list())
  list(table = tab, qc = qc, screen_adjusted = TRUE)
}

#' DRIMSeq on one comparison's samples (Love et al. 2018).
dtu_engine_drimseq <- function(counts, cd, main, covariates, gene, target, reference, sd_filter = 0.1, seed = 1) {
  smp <- data.frame(sample_id = colnames(counts), cd[, c(covariates, main), drop = FALSE],
                    stringsAsFactors = FALSE, check.names = FALSE)
  smp[[main]] <- factor(as.character(smp[[main]]), levels = c(reference, target))
  rownames(smp) <- NULL
  cts <- data.frame(gene_id = gene, feature_id = rownames(counts), counts, check.names = FALSE,
                    stringsAsFactors = FALSE)
  d <- DRIMSeq::dmDSdata(counts = cts, samples = smp)
  design <- stats::model.matrix(stats::as.formula(paste("~", paste(c(covariates, main), collapse = " + "))),
                                data = DRIMSeq::samples(d))
  set.seed(seed)
  d <- DRIMSeq::dmPrecision(d, design = design)
  d <- DRIMSeq::dmFit(d, design = design)
  d <- DRIMSeq::dmTest(d, coef = colnames(design)[ncol(design)])
  res_fun <- dtu_fn("results", "DRIMSeq")
  rg <- res_fun(d); rt <- res_fun(d, level = "feature")
  i <- match(rownames(counts), rt$feature_id)
  # Post-hoc filter of the workflow: transcripts whose share varies little across
  # samples (SD < sd_filter) get p = 1.
  p <- as.numeric(rt$pvalue[i]); pa <- as.numeric(rt$adj_pvalue[i])
  if (sd_filter > 0) {
    tot <- rowsum(counts, gene, reorder = FALSE)[gene, , drop = FALSE]
    pr <- counts / tot; pr[!is.finite(pr)] <- NA
    small <- apply(pr, 1, stats::sd, na.rm = TRUE) < sd_filter
    small[is.na(small)] <- FALSE
    p[small] <- 1; pa[small] <- 1
  }
  fp <- dtu_fn("proportions", "DRIMSeq")(d)
  fpm <- as.matrix(fp[match(rownames(counts), fp$feature_id), colnames(counts), drop = FALSE])
  grp <- as.character(smp[[main]])
  lfc <- log2(rowMeans(fpm[, grp == target, drop = FALSE]) / rowMeans(fpm[, grp == reference, drop = FALSE]))
  gi <- match(gene, rg$gene_id)
  tab <- data.frame(gene_id = rownames(counts), parent_gene_id = gene,
                    baseMean = dtu_base_mean(counts, seq_len(ncol(counts))),
                    log2FoldChange = as.numeric(lfc), stat = as.numeric(rt$lr[i]),
                    pvalue = p, padj = pa,
                    gene_pvalue = as.numeric(rg$pvalue[gi]), gene_padj = as.numeric(rg$adj_pvalue[gi]),
                    stringsAsFactors = FALSE)
  qc <- tryCatch({
    me <- dtu_fn("mean_expression", "DRIMSeq")(d); gp <- dtu_fn("genewise_precision", "DRIMSeq")(d)
    list(precision = data.frame(mean_expression = me$mean_expression[match(gp$gene_id, me$gene_id)],
                                precision = gp$genewise_precision,
                                common = dtu_fn("common_precision", "DRIMSeq")(d)))
  }, error = function(e) list())
  list(table = tab, qc = qc, screen_adjusted = FALSE)
}

#' satuRn on all samples, every comparison at once (Gilis et al. 2022).
dtu_engine_saturn <- function(counts, cd, main, covariates, gene, contrasts, pvalues = "empirical") {
  form <- stats::as.formula(paste("~ 0 +", paste(c(main, covariates), collapse = " + ")))
  cdf <- cd[, c(main, covariates), drop = FALSE]
  rownames(cdf) <- colnames(counts)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = S4Vectors::DataFrame(cdf),
    rowData = S4Vectors::DataFrame(isoform_id = rownames(counts), gene_id = gene))
  se <- satuRn::fitDTU(object = se, formula = form, parallel = FALSE,
                       BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  design <- stats::model.matrix(form, cdf)
  L <- matrix(0, ncol(design), nrow(contrasts), dimnames = list(colnames(design), contrasts$id))
  for (k in seq_len(nrow(contrasts)))
    L[, k] <- dtu_contrast(design, main, contrasts$target[k], contrasts$reference[k])
  se <- satuRn::testDTU(object = se, contrasts = L, diagplot1 = FALSE, diagplot2 = FALSE, sort = FALSE)
  rd <- SummarizedExperiment::rowData(se)
  tabs <- list(); notes <- character(0)
  for (k in seq_len(nrow(contrasts))) {
    r <- as.data.frame(rd[[paste0("fitDTUResult_", contrasts$id[k])]])
    use_emp <- identical(pvalues, "empirical") && any(is.finite(r$empirical_pval))
    if (identical(pvalues, "empirical") && !use_emp)
      notes <- c(notes, sprintf("satuRn could not compute empirical p-values for %s; the regular ones are used.", contrasts$id[k]))
    p <- if (use_emp) r$empirical_pval else r$pval
    q <- dtu_per_gene_qvalue(p, gene)
    cols <- which(cd[[main]] %in% c(contrasts$reference[k], contrasts$target[k]))
    tabs[[contrasts$id[k]]] <- data.frame(
      gene_id = rownames(counts), parent_gene_id = gene, baseMean = dtu_base_mean(counts, cols),
      log2FoldChange = r$estimates / log(2), stat = r$t, pvalue = p,
      padj = if (use_emp) r$empirical_FDR else r$regular_FDR,
      gene_pvalue = NA_real_, gene_padj = as.numeric(q[gene]),
      log_odds_ratio = r$estimates, se = r$se, stringsAsFactors = FALSE)
  }
  list(tables = tabs, qc = list(), notes = notes, screen_adjusted = TRUE)
}

#' edgeR quasi-likelihood DTU (Baldoni et al. 2025): diffSplice on a glmQLFit,
#' or diffSpliceDGE on edgeR releases before 4.6.
dtu_engine_edger <- function(counts, cd, main, covariates, gene, contrasts, gene_test = "simes",
                             robust = TRUE) {
  design <- dtu_design(cd, main, covariates)
  y <- edgeR::DGEList(counts, genes = data.frame(TranscriptID = rownames(counts), GeneID = gene,
                                                 stringsAsFactors = FALSE))
  norm_lib <- if ("normLibSizes" %in% getNamespaceExports("edgeR")) edgeR::normLibSizes else edgeR::calcNormFactors
  y <- norm_lib(y)
  fit <- edgeR::glmQLFit(y, design, robust = robust, legacy = FALSE)
  generic <- exists("diffSplice.DGEGLM", envir = asNamespace("edgeR"), inherits = FALSE)
  eff <- y$samples$lib.size * y$samples$norm.factors
  tabs <- list()
  for (k in seq_len(nrow(contrasts))) {
    cv <- dtu_contrast(design, main, contrasts$target[k], contrasts$reference[k])
    if (generic) {
      # edgeR >= 4.6: natural-log coefficients relative to the mean of the gene's other transcripts
      ds <- limma::diffSplice(fit, contrast = cv, geneid = "GeneID", exonid = "TranscriptID", verbose = FALSE)
      tx <- ds$genes$TranscriptID; lfc <- as.numeric(ds$coefficients) / log(2)
      stat <- as.numeric(ds$t); p <- as.numeric(ds$p.value)
      gid <- ds$gene.genes$GeneID
      gp <- as.numeric(if (identical(gene_test, "F")) ds$gene.F.p.value else ds$gene.simes.p.value)
    } else {
      # diffSpliceDGE: natural-log coefficients relative to the whole gene
      ds <- edgeR::diffSpliceDGE(fit, contrast = cv, geneid = "GeneID", exonid = "TranscriptID", verbose = FALSE)
      tx <- ds$genes$TranscriptID; lfc <- as.numeric(ds$coefficients) / log(2)
      stat <- as.numeric(ds$exon.F); p <- as.numeric(ds$exon.p.value)
      gid <- ds$gene.genes$GeneID
      gp <- as.numeric(if (identical(gene_test, "F")) ds$gene.p.value else ds$gene.Simes.p.value)
    }
    cols <- which(cd[[main]] %in% c(contrasts$reference[k], contrasts$target[k]))
    sf <- eff / exp(mean(log(eff)))
    base <- rowMeans(t(t(counts[tx, cols, drop = FALSE]) / sf[cols]))
    g <- gene[match(tx, rownames(counts))]
    gpa <- stats::setNames(stats::p.adjust(gp, "BH"), gid)
    tabs[[contrasts$id[k]]] <- data.frame(
      gene_id = tx, parent_gene_id = g, baseMean = base, log2FoldChange = lfc, stat = stat,
      pvalue = p, padj = stats::p.adjust(p, "BH"),
      gene_pvalue = as.numeric(stats::setNames(gp, gid)[g]), gene_padj = as.numeric(gpa[g]),
      stringsAsFactors = FALSE)
  }
  ql <- tryCatch(data.frame(AveLogCPM = fit$AveLogCPM,
                            raw = if (!is.null(fit$deviance.adj)) fit$deviance.adj / fit$df.residual.adj else fit$deviance / fit$df.residual.zeros,
                            squeezed = if (!is.null(fit$var.post)) fit$var.post else fit$s2.post,
                            trend = if (!is.null(fit$var.prior)) fit$var.prior else fit$s2.prior),
                 error = function(e) NULL)
  list(tables = tabs, qc = list(size_factors = stats::setNames(eff / exp(mean(log(eff))), colnames(counts)), ql = ql),
       notes = sprintf("edgeR %s on a quasi-likelihood fit (glmQLFit, robust = %s).",
                       if (generic) "diffSplice" else "diffSpliceDGE", robust),
       screen_adjusted = FALSE)
}

#' limma DTU (Baldoni et al. 2025): voomLmFit, then diffSplice per comparison.
dtu_engine_limma <- function(counts, cd, main, covariates, gene, contrasts, gene_test = "simes",
                             robust = TRUE) {
  design <- dtu_design(cd, main, covariates)
  y <- edgeR::DGEList(counts, genes = data.frame(TranscriptID = rownames(counts), GeneID = gene,
                                                 stringsAsFactors = FALSE))
  norm_lib <- if ("normLibSizes" %in% getNamespaceExports("edgeR")) edgeR::normLibSizes else edgeR::calcNormFactors
  y <- norm_lib(y)
  fit <- edgeR::voomLmFit(y, design, save.plot = TRUE)
  eff <- y$samples$lib.size * y$samples$norm.factors
  tabs <- list()
  for (k in seq_len(nrow(contrasts))) {
    cv <- dtu_contrast(design, main, contrasts$target[k], contrasts$reference[k])
    f2 <- limma::contrasts.fit(fit, cv)
    ds <- limma::diffSplice(f2, geneid = "GeneID", exonid = "TranscriptID", robust = robust, verbose = FALSE)
    tx <- ds$genes$TranscriptID; p <- as.numeric(ds$p.value)
    gid <- ds$gene.genes$GeneID
    gp <- as.numeric(if (identical(gene_test, "F")) ds$gene.F.p.value else ds$gene.simes.p.value)
    cols <- which(cd[[main]] %in% c(contrasts$reference[k], contrasts$target[k]))
    sf <- eff / exp(mean(log(eff)))
    g <- gene[match(tx, rownames(counts))]
    gpa <- stats::setNames(stats::p.adjust(gp, "BH"), gid)
    tabs[[contrasts$id[k]]] <- data.frame(
      gene_id = tx, parent_gene_id = g,
      baseMean = rowMeans(t(t(counts[tx, cols, drop = FALSE]) / sf[cols])),
      log2FoldChange = as.numeric(ds$coefficients), stat = as.numeric(ds$t),
      pvalue = p, padj = stats::p.adjust(p, "BH"),
      gene_pvalue = as.numeric(stats::setNames(gp, gid)[g]), gene_padj = as.numeric(gpa[g]),
      stringsAsFactors = FALSE)
  }
  vq <- tryCatch(list(x = fit$voom.xy$x, y = fit$voom.xy$y, line_x = fit$voom.line$x, line_y = fit$voom.line$y),
                 error = function(e) NULL)
  list(tables = tabs, qc = list(size_factors = stats::setNames(eff / exp(mean(log(eff))), colnames(counts)), voom = vq),
       notes = sprintf("limma diffSplice on voomLmFit (robust = %s).", robust), screen_adjusted = FALSE)
}

#' Run DTU for every comparison.
#'
#' @param counts  transcripts x samples, already in the method's count type
#' @param cd      design table (rows = columns of counts)
#' @param gene    parent gene of each row of counts
#' @param opts    list: filter, min_gene_expr, min_feature_expr, min_feature_prop,
#'                gene_test, robust, sd_filter, saturn_pvalues, seed, stager, stager_alpha
#' @return list(results, kept, filter, qc, notes)
dtu_run <- function(counts, cd, main, covariates, contrasts, engine, gene, opts, progress = NULL) {
  tick <- function(f, d) if (is.function(progress)) progress(f, d)
  gene <- ifelse(is.na(gene) | gene == "", rownames(counts), gene)
  results <- list(); qc <- list(); notes <- character(0); kept <- character(0); fdesc <- NULL
  finish <- function(tab, ref_cols, tgt_cols, cnt, screen_adjusted) {
    pr <- dtu_props(cnt[tab$gene_id, , drop = FALSE], tab$parent_gene_id, ref_cols, tgt_cols)
    tab$prop_reference <- pr$prop_reference; tab$prop_target <- pr$prop_target
    tab$delta_prop <- pr$prop_target - pr$prop_reference
    if (isTRUE(opts$stager)) tab <- dtu_stagewise(tab, opts$stager_alpha, screen_adjusted)
    tab[order(tab$gene_padj, tab$pvalue, na.last = TRUE), , drop = FALSE]
  }
  if (engine %in% c("dexseq", "drimseq")) {
    for (k in seq_len(nrow(contrasts))) {
      cr <- contrasts[k, ]
      tick(0.1 + 0.8 * (k - 1) / nrow(contrasts), paste(if (engine == "dexseq") "DEXSeq:" else "DRIMSeq:", cr$id))
      cols <- which(as.character(cd[[main]]) %in% c(cr$reference, cr$target))
      cnt <- counts[, cols, drop = FALSE]; cdk <- cd[cols, , drop = FALSE]
      for (v in c(main, covariates)) if (is.factor(cdk[[v]])) cdk[[v]] <- droplevels(cdk[[v]])
      fl <- dtu_filter(cnt, gene, cdk[[main]], NULL, opts$filter, engine, opts)
      if (sum(fl$keep) < 4) stop("Only ", sum(fl$keep), " transcripts pass the filter for ", cr$id, "; loosen it.")
      fdesc <- fl$description
      x <- cnt[fl$keep, , drop = FALSE]; gx <- gene[fl$keep]
      out <- if (engine == "dexseq") dtu_engine_dexseq(x, cdk, main, covariates, gx, cr$target, cr$reference)
             else dtu_engine_drimseq(x, cdk, main, covariates, gx, cr$target, cr$reference,
                                     opts$sd_filter, opts$seed)
      grp <- as.character(cdk[[main]])
      results[[cr$id]] <- finish(out$table, which(grp == cr$reference), which(grp == cr$target), x, out$screen_adjusted)
      kept <- union(kept, rownames(x))
      if (!length(qc)) qc <- out$qc
    }
  } else {
    tick(0.1, "Filtering transcripts")
    design <- dtu_design(cd, main, covariates)
    fl <- dtu_filter(counts, gene, cd[[main]], design, opts$filter, engine, opts)
    if (sum(fl$keep) < 4) stop("Only ", sum(fl$keep), " transcripts pass the filter; loosen it.")
    fdesc <- fl$description
    x <- counts[fl$keep, , drop = FALSE]; gx <- gene[fl$keep]
    tick(0.3, switch(engine, saturn = "satuRn: fitting quasi-binomial models",
                     edger_dtu = "edgeR: quasi-likelihood fit", "limma: voomLmFit"))
    out <- switch(engine,
      saturn = dtu_engine_saturn(x, cd, main, covariates, gx, contrasts, opts$saturn_pvalues),
      edger_dtu = dtu_engine_edger(x, cd, main, covariates, gx, contrasts, opts$gene_test, opts$robust),
      limma_dtu = dtu_engine_limma(x, cd, main, covariates, gx, contrasts, opts$gene_test, opts$robust))
    grp <- as.character(cd[[main]])
    for (k in seq_len(nrow(contrasts))) {
      cr <- contrasts[k, ]
      tick(0.7 + 0.25 * (k - 1) / nrow(contrasts), paste("Collecting", cr$id))
      results[[cr$id]] <- finish(out$tables[[cr$id]], which(grp == cr$reference), which(grp == cr$target),
                                 x, out$screen_adjusted)
    }
    kept <- unique(unlist(lapply(results, function(d) d$gene_id))); qc <- out$qc; notes <- c(notes, out$notes)
  }
  for (id in names(results)) rownames(results[[id]]) <- NULL
  list(results = results, kept = kept, filter = fdesc, qc = qc, notes = notes)
}

DTU_CORE_FUNS <- c("dtu_fn", "dtu_counts", "dtu_filter_drimseq", "dtu_drop_singletons", "dtu_filter",
                   "dtu_design", "dtu_contrast", "dtu_per_gene_qvalue", "dtu_props", "dtu_base_mean",
                   "dtu_stagewise", "dtu_engine_dexseq", "dtu_engine_drimseq", "dtu_engine_saturn",
                   "dtu_engine_edger", "dtu_engine_limma", "dtu_run")

# ---- app wrapper -------------------------------------------------------------------

#' Effective lengths aligned to the count matrix (missing: the transcript's mean
#' over samples, else the overall median).
dtu_align_lengths <- function(lengths, rows, samples) {
  if (is.null(lengths) || is.null(dim(lengths)) || !all(samples %in% colnames(lengths))) return(NULL)
  m <- matrix(NA_real_, length(rows), length(samples), dimnames = list(rows, samples))
  hit <- intersect(rows, rownames(lengths))
  if (!length(hit)) return(NULL)
  m[hit, ] <- suppressWarnings(as.numeric(as.matrix(lengths[hit, samples, drop = FALSE])))
  m[!is.finite(m) | m < 0] <- NA
  rm <- rowMeans(m, na.rm = TRUE); med <- stats::median(m, na.rm = TRUE)
  idx <- which(is.na(m), arr.ind = TRUE)
  if (nrow(idx)) m[idx] <- ifelse(is.finite(rm[idx[, 1]]), rm[idx[, 1]], med)
  m
}

#' DTU run from the Analysis tab (called by run_de_analysis()).
run_dtu_analysis <- function(counts, sheet, included, main, covariates = character(0), types = list(),
                             contrasts, engine = "edger_dtu", options = de_default_options(),
                             progress = NULL, lengths = NULL, transcripts = NULL) {
  tick <- function(f, d) if (is.function(progress)) progress(f, d)
  options <- utils::modifyList(de_default_options(), options)
  miss <- de_engine_missing(engine)
  if (!is.null(miss)) stop("The ", engine_label(engine), " package is not installed. Run: ", miss)
  covariates <- setdiff(if (is.null(covariates)) character(0) else covariates, main)
  chk <- validate_design(sheet, included, main, covariates, types, contrasts)
  chk$errors <- c(chk$errors, engine_design_errors(engine, c(covariates, main), types, options))
  if (length(chk$errors)) stop(paste(chk$errors, collapse = "\n"))
  vars <- c(covariates, main)
  cd <- design_coldata(sheet, included, vars, types)
  for (v in vars) if (is.factor(cd[[v]])) cd[[v]] <- droplevels(cd[[v]])
  notes <- chk$warnings

  mat <- counts[, cd$sample, drop = FALSE]
  storage.mode(mat) <- "double"
  if (any(is.na(mat))) { mat[is.na(mat)] <- 0; notes <- c(notes, "Missing counts were set to 0.") }
  t2g <- merge_tx2gene(rownames(counts), transcripts$tx2gene)
  if (all(is.na(t2g$gene_id)))
    stop("Transcript usage needs a transcript-to-gene map: load a GTF or a tx2gene table.")
  gene <- t2g$gene_id[match(rownames(mat), t2g$transcript_id)]
  if (any(is.na(gene))) notes <- c(notes, sprintf(
    "%s transcript(s) have no gene in the map and were left out (a gene needs two or more transcripts).",
    format(sum(is.na(gene)), big.mark = ",")))

  len <- dtu_align_lengths(lengths, rownames(mat), colnames(mat))
  od <- transcripts$overdispersion
  rc <- dtu_resolve_counts(engine, options$dtu_counts, !is.null(len), !is.null(od))
  notes <- c(notes, rc$note, sprintf("Counts: %s.", DTU_COUNT_LABELS[[rc$mode]]))
  tick(0.05, "Preparing counts")
  x <- dtu_counts(mat, rc$mode, len, od, gene)

  stager <- isTRUE(options$stager)
  if (stager && !requireNamespace("stageR", quietly = TRUE)) {
    stager <- FALSE
    notes <- c(notes, 'Stage-wise testing was skipped: the stageR package is not installed (BiocManager::install("stageR")).')
  }
  opts <- list(filter = dtu_resolve_filter(engine, options$dtu_filter),
               min_gene_expr = options$dtu_min_gene_expr, min_feature_expr = options$dtu_min_feature_expr,
               min_feature_prop = options$dtu_min_feature_prop,
               gene_test = options$dtu_gene_test, robust = isTRUE(options$dtu_robust),
               sd_filter = options$drimseq_sd_filter, saturn_pvalues = options$saturn_pvalues,
               seed = 1L, stager = stager, stager_alpha = options$stager_alpha)
  out <- dtu_run(x, cd, main, covariates, contrasts, engine, gene, opts, progress)

  # transcript annotation next to the ids
  tk <- t2g[match(rownames(mat), t2g$transcript_id), , drop = FALSE]
  for (id in names(out$results)) {
    df <- out$results[[id]]
    a <- tk[match(df$gene_id, tk$transcript_id), , drop = FALSE]
    out$results[[id]] <- data.frame(df[, 1:2], gene_name = a$gene_name, transcript_name = a$transcript_name,
                                    transcript_biotype = a$transcript_biotype, df[, -(1:2), drop = FALSE],
                                    stringsAsFactors = FALSE, check.names = FALSE)
  }
  n_g <- vapply(out$results, function(d) length(unique(d$parent_gene_id[!is.na(d$gene_padj) &
                                                                          d$gene_padj <= 0.05])), integer(1))
  notes <- c(notes, out$notes,
    sprintf(paste("Transcript usage: log2FoldChange is the change in a transcript's share of its gene",
                  "(%s); prop_reference / prop_target are its mean share in each group and delta_prop their",
                  "difference. gene_padj is the gene-level FDR (%s); padj is per transcript."),
            switch(engine, dexseq = "DEXSeq exon-usage fold change", drimseq = "ratio of DRIMSeq's fitted proportions",
                   saturn = "satuRn log-odds ratio / log(2)",
                   edger_dtu = if (exists("diffSplice.DGEGLM", envir = asNamespace("edgeR"), inherits = FALSE))
                     "diffSplice coefficient: the transcript's log2 fold change minus the mean of its gene's other transcripts"
                     else "diffSpliceDGE coefficient: the transcript's log2 fold change minus its gene's",
                   "diffSplice coefficient: the transcript's log2 fold change minus the mean of its gene's transcripts"),
            switch(engine, dexseq = , saturn = "per-gene q-value from the transcripts' p-values (DEXSeq perGeneQValue)",
                   drimseq = "DRIMSeq gene-level likelihood-ratio test, BH",
                   if (identical(opts$gene_test, "F")) "diffSplice gene-level F-test, BH" else "Simes combination of the transcripts' p-values, BH")),
    sprintf("Genes with DTU at gene-level FDR 5%%: %s.", paste(sprintf("%s %s", names(n_g), n_g), collapse = "; ")),
    if (stager) sprintf(paste("Stage-wise testing (stageR %s, method \"dtu\", alpha %g): stageR_gene_padj and",
                              "stageR_tx_padj are valid at this alpha only; a transcript is significant when both are <= alpha."),
                        as.character(utils::packageVersion("stageR")), opts$stager_alpha))
  if (engine == "drimseq" && opts$sd_filter > 0) notes <- c(notes, sprintf(
    "DRIMSeq post-hoc filter: transcripts whose share of their gene varies with SD < %g across samples were given p = 1 (Love et al. 2018).",
    opts$sd_filter))

  pkgs <- DTU_ENGINE_PACKAGES[[engine]]
  pkgs <- setdiff(pkgs, c("SummarizedExperiment", "BiocParallel"))
  vers <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), "")
  if (stager) vers <- c(vers, stageR = as.character(utils::packageVersion("stageR")))

  keep_rows <- out$kept
  xk <- x[keep_rows, , drop = FALSE]
  lib <- colSums(xk); sf <- lib / exp(mean(log(lib)))
  qc <- c(list(lib_size = colSums(x), detected = colSums(x > 0)), out$qc)
  if (is.null(qc$size_factors)) qc$size_factors <- sf
  tick(1, "Done")
  list(engine = engine_label(engine), engine_label = engine_label(engine), engine_key = engine,
       results = out$results,
       transformed = log2(t(t(xk) / lib * 1e6) + 1), transform_label = "log2(CPM + 1)",
       normalised = t(t(xk) / sf), size_factors = sf,
       notes = notes, versions = vers,
       kept = keep_rows, n_input = nrow(counts), counts = x,
       formula = paste(deparse(design_formula(main, covariates)), collapse = ""),
       filter = out$filter, coldata = cd, main = main, qc = qc,
       covariates = covariates, contrasts = contrasts, settings = options, dtu_opts = opts,
       length_offset = FALSE, lengths_used = NULL,
       level = "dtu", dtu_counts = rc$mode, dtu_lengths = if (rc$mode %in% c("scaledTPM", "dtuScaledTPM")) len else NULL,
       divided = identical(rc$mode, "divided"), tx2gene = t2g, ngvector = NULL,
       overdispersion = if (identical(rc$mode, "divided")) od else NULL,
       input_counts = counts[, cd$sample, drop = FALSE], filter_n = NA,
       sheet = sheet[match(cd$sample, sheet$sample), , drop = FALSE], types = types,
       run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

#' Genes with DTU (gene-level FDR) in one result table.
dtu_gene_count <- function(df, padj_cutoff = 0.05) {
  if (!"gene_padj" %in% colnames(df)) return(NA_integer_)
  length(unique(df$parent_gene_id[!is.na(df$gene_padj) & df$gene_padj <= padj_cutoff]))
}

# ---- transcript usage plot -----------------------------------------------------------

#' Genes of a DTU result for the plot picker, most significant first.
dtu_gene_choices <- function(res, id) {
  df <- res$results[[id]]
  if (is.null(df)) return(character(0))
  o <- df[order(df$gene_padj, df$pvalue, na.last = TRUE), , drop = FALSE]
  o <- o[!duplicated(o$parent_gene_id), , drop = FALSE]
  lab <- ifelse(!is.na(o$gene_name) & nzchar(o$gene_name), o$gene_name, o$parent_gene_id)
  stats::setNames(o$parent_gene_id, sprintf("%s (gene FDR %s)", lab, signif(o$gene_padj, 2)))
}

#' Share of each tested transcript of one gene, per sample, by group.
dtu_usage_plot <- function(res, id, gene, style = default_style(), padj_cutoff = 0.05) {
  st <- resolve_style(style)
  df <- res$results[[id]]
  tx <- df[df$parent_gene_id %in% gene, , drop = FALSE]
  if (!nrow(tx)) stop("That gene was not tested in this comparison.")
  cr <- res$contrasts[res$contrasts$id == id, ]
  smp <- res$coldata$sample[as.character(res$coldata[[res$main]]) %in% c(cr$reference, cr$target)]
  cnt <- res$counts[tx$gene_id, smp, drop = FALSE]
  tot <- colSums(cnt)
  pr <- t(t(cnt) / ifelse(tot > 0, tot, NA))
  lab <- ifelse(!is.na(tx$transcript_name) & nzchar(tx$transcript_name), tx$transcript_name, tx$gene_id)
  sig <- if ("stageR_tx_padj" %in% colnames(tx)) !is.na(tx$stageR_tx_padj) & tx$stageR_tx_padj <= padj_cutoff
         else !is.na(tx$padj) & tx$padj <= padj_cutoff
  lab <- paste0(lab, ifelse(sig, " *", ""))
  long <- data.frame(transcript = factor(rep(lab, ncol(pr)), levels = lab[order(-tx$prop_reference)]),
                     sample = rep(colnames(pr), each = nrow(pr)), share = as.vector(pr),
                     stringsAsFactors = FALSE)
  long$group <- factor(as.character(res$coldata[[res$main]])[match(long$sample, res$coldata$sample)],
                       levels = c(cr$reference, cr$target))
  means <- stats::aggregate(share ~ transcript + group, long, mean)
  gname <- tx$gene_name[1]; if (is.na(gname) || !nzchar(gname)) gname <- gene
  pal <- stats::setNames(group_palette(2), levels(long$group))
  ggplot2::ggplot(long, ggplot2::aes(x = transcript, y = share, colour = group)) +
    ggplot2::geom_col(data = means, ggplot2::aes(fill = group), colour = NA, alpha = 0.25,
                      position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 1),
                        size = st$point_size * 0.9) +
    ggplot2::scale_colour_manual(values = pal, name = res$main) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"), limits = c(0, 1)) +
    deg_theme(st) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    ggplot2::labs(title = sprintf("%s: transcript usage, %s vs %s", gname, cr$target, cr$reference),
                  subtitle = sprintf("%s; gene-level FDR %s; * transcript significant%s. Shares among the %d tested transcripts.",
                                     res$engine_label, signif(tx$gene_padj[1], 3),
                                     if ("stageR_tx_padj" %in% colnames(tx)) " (stageR)" else "", nrow(tx)),
                  x = NULL, y = "Share of the gene's counts")
}

# ---- reproduce script --------------------------------------------------------------

#' Standalone script for a DTU run: the same dtu_* functions, applied to the bundle's files.
dtu_reproduce_script <- function(res) {
  o <- res$dtu_opts; key <- res$engine_key; main <- res$main; cov <- res$covariates; cd <- res$coldata
  pkgs <- c(DTU_ENGINE_PACKAGES[[key]], if (o$filter == "edger" || key %in% c("edger_dtu", "limma_dtu", "saturn")) "edgeR",
            if (isTRUE(o$stager)) "stageR", if (key == "dexseq") c("DESeq2", "S4Vectors"))
  pkgs <- unique(pkgs)
  vers <- paste(names(res$versions), res$versions, collapse = ", ")
  L <- character(0); add <- function(...) L <<- c(L, ...)
  add("# =============================================================================",
      sprintf("# Reproduces the DEGexplorer transcript-usage (DTU) run of %s", res$run_at),
      sprintf("# Method: %s (%s, R %s.%s)", res$engine_label, vers, R.version$major, R.version$minor),
      sprintf("# Model:  %s", res$formula),
      "#",
      "# Reads the files exported with it, from the same folder:",
      "#   counts.tsv       estimated transcript counts of the included samples, before filtering",
      "#   samplesheet.csv  sample sheet",
      "#   contrasts.csv    comparisons (id, variable, reference, target, blocking)",
      "#   tx2gene.tsv      transcript -> gene map",
      if (!is.null(res$dtu_lengths)) "#   lengths.tsv      effective transcript lengths per sample (for scaledTPM)",
      if (isTRUE(res$divided)) "#   overdispersion.tsv  quantification overdispersion per transcript (catchSalmon / catchKallisto)",
      "# and writes reproduced/<comparison>.dtu.results.tsv.",
      "#",
      "# The dtu_* functions are the ones the app ran.",
      "# Run it from that folder:  Rscript reproduce.R",
      "# =============================================================================",
      "",
      "suppressPackageStartupMessages({",
      paste0("  library(", pkgs, ")"),
      "})",
      "",
      "counts <- as.matrix(read.delim(\"counts.tsv\", row.names = 1, check.names = FALSE))",
      "sheet <- read.csv(\"samplesheet.csv\", colClasses = \"character\", check.names = FALSE)",
      "contrasts <- read.csv(\"contrasts.csv\", colClasses = \"character\")",
      sprintf("main <- %s", .dq(main)),
      sprintf("covariates <- %s", .dq(if (length(cov)) cov else character(0))),
      "",
      "# Design table; factor levels are written out so they match the app exactly",
      "coldata <- data.frame(row.names = sheet$sample)")
  for (v in c(cov, main)) {
    if (is.numeric(cd[[v]])) add(sprintf("coldata[[%s]] <- as.numeric(sheet[[%s]])", .dq(v), .dq(v)))
    else add(sprintf("coldata[[%s]] <- factor(trimws(sheet[[%s]]), levels = %s)", .dq(v), .dq(v), .dq(levels(cd[[v]]))))
  }
  add("counts <- counts[, rownames(coldata), drop = FALSE]",
      "counts[is.na(counts)] <- 0",
      "tx2gene <- read.delim(\"tx2gene.tsv\", colClasses = \"character\", check.names = FALSE)",
      "gene <- tx2gene$gene_id[match(rownames(counts), tx2gene$transcript_id)]",
      if (!is.null(res$dtu_lengths)) c(
        "lengths <- as.matrix(read.delim(\"lengths.tsv\", row.names = 1, check.names = FALSE))",
        "lengths <- lengths[rownames(counts), colnames(counts)]") else "lengths <- NULL",
      if (isTRUE(res$divided)) c(
        "odt <- read.delim(\"overdispersion.tsv\", colClasses = c(\"character\", \"numeric\"))",
        "od <- setNames(odt$overdispersion, odt$transcript_id)") else "od <- NULL",
      "", "# ---- the app's DTU functions --------------------------------------------------")
  for (f in DTU_CORE_FUNS) {
    add("", paste(f, "<-", paste(deparse(get(f, mode = "function")), collapse = "\n")))
  }
  add("", "# ---- run -----------------------------------------------------------------------",
      sprintf("# Counts: %s", DTU_COUNT_LABELS[[res$dtu_counts]]),
      sprintf("x <- dtu_counts(counts, %s, lengths, od, gene)", .dq(res$dtu_counts)),
      sprintf("opts <- %s", .dq(o)),
      sprintf("out <- dtu_run(x, coldata, main, covariates, contrasts, %s, gene, opts)", .dq(key)),
      "cat(length(out$kept), \"transcripts tested\\n\")",
      "",
      "dir.create(\"reproduced\", showWarnings = FALSE)",
      "for (id in names(out$results)) {",
      "  df <- out$results[[id]]",
      "  a <- tx2gene[match(df$gene_id, tx2gene$transcript_id), , drop = FALSE]",
      "  df <- data.frame(df[, 1:2], gene_name = a$gene_name, transcript_name = a$transcript_name,",
      "                   transcript_biotype = a$transcript_biotype, df[, -(1:2), drop = FALSE], check.names = FALSE)",
      "  colnames(df)[1] <- \"transcript_id\"",
      "  write.table(df, file.path(\"reproduced\", paste0(id, \".dtu.results.tsv\")),",
      "              sep = \"\\t\", quote = FALSE, row.names = FALSE, na = \"NA\")",
      "}",
      "writeLines(capture.output(sessionInfo()), file.path(\"reproduced\", \"session_info.txt\"))",
      "cat(\"Results written to\", normalizePath(\"reproduced\"), \"\\n\")")
  L
}
