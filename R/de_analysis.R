# =============================================================================
# de_analysis.R -- differential expression run inside the app
#
# Phase 1: DESeq2. limma-voom, edgeR and EBSeq plug into the same interface
# (run_de_analysis() dispatches on `engine`).
#
# Everything here is plain R with no Shiny, so the smoke test can exercise it.
# The experimental design lives in a "sample sheet" data frame built in the
# app: one row per sample (column `sample`), then one character column per
# user-defined variable. `types` records whether each variable is a factor or
# a numeric covariate; `included` lists the samples taking part.
# =============================================================================

HAS_DESEQ2 <- requireNamespace("DESeq2", quietly = TRUE)

DE_ENGINES <- c("DESeq2" = "deseq2", "limma-voom" = "limma", "edgeR" = "edger", "EBSeq" = "ebseq")

# limma-voom normalises with edgeR's TMM, so it needs both packages.
DE_ENGINE_PACKAGES <- list(deseq2 = "DESeq2", limma = c("limma", "edgeR"), edger = "edgeR",
                           ebseq = "EBSeq")

#' Install line for a missing engine, or NULL when it is available.
de_engine_missing <- function(engine) {
  pkgs <- if (exists("DTU_ENGINE_PACKAGES")) c(DE_ENGINE_PACKAGES, DTU_ENGINE_PACKAGES) else DE_ENGINE_PACKAGES
  pkg <- pkgs[[engine]]
  miss <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(miss)) return(NULL)
  sprintf("BiocManager::install(c(%s))", paste0('"', miss, '"', collapse = ", "))
}

#' Design problems specific to one engine (on top of validate_design()).
engine_design_errors <- function(engine, covariates = character(0), types = list(),
                                 options = de_default_options()) {
  err <- character(0)
  covariates <- if (is.null(covariates)) character(0) else covariates
  if (identical(engine, "ebseq") && length(covariates)) {
    err <- c(err, paste("EBSeq compares two groups at a time and cannot adjust for covariates or",
                        "blocking factors. Untick them, or choose DESeq2, limma-voom or edgeR."))
  }
  b <- options$limma_block
  if (identical(engine, "limma") && !is.null(b) && nzchar(b) && b != "(none)") {
    if (!b %in% covariates) err <- c(err, sprintf("The random block '%s' must also be ticked as a covariate.", b))
    else if (identical(types[[b]], "numeric")) err <- c(err, sprintf("The random block '%s' must be a factor.", b))
  }
  if (engine %in% c("dexseq", "drimseq")) {
    bad <- intersect(covariates, c("exon", "sample_id", "gene_id", "feature_id"))
    if (length(bad)) err <- c(err, sprintf("%s uses the name(s) %s internally; rename that variable.",
                                           if (engine == "dexseq") "DEXSeq" else "DRIMSeq", paste(bad, collapse = ", ")))
  }
  sh <- options$shrink
  if (identical(engine, "deseq2") && !is.null(sh) && sh %in% c("apeglm", "ashr") &&
      !requireNamespace(sh, quietly = TRUE)) {
    err <- c(err, sprintf("Fold-change shrinkage with %s needs the '%s' package: %s, or choose another shrinkage type.",
                          sh, sh, if (sh == "apeglm") 'BiocManager::install("apeglm")' else 'install.packages("ashr")'))
  }
  err
}

# ---- sample sheet ------------------------------------------------------------

#' A fresh sample sheet from the count-matrix column names.
design_scaffold <- function(samples) {
  data.frame(sample = as.character(samples), stringsAsFactors = FALSE)
}

#' Turn a user-typed name into a valid R variable name that is not taken.
clean_variable_name <- function(name, existing = character(0)) {
  if (is.null(name) || !nzchar(trimws(name))) stop("Give the variable a name.")
  nm <- make.names(trimws(name))
  if (nm %in% c("sample", existing)) stop("A variable called '", nm, "' already exists.")
  nm
}

#' Guess a value for each sample from its name.
#'
#' @param sep      delimiter, e.g. "_"
#' @param position which piece to take: 1 = first, -1 = last
guess_from_names <- function(samples, sep = "_", position = 1) {
  if (is.null(sep) || !nzchar(sep)) stop("Give a delimiter.")
  position <- as.integer(position)
  parts <- strsplit(as.character(samples), sep, fixed = TRUE)
  vapply(parts, function(p) {
    i <- if (position < 0) length(p) + position + 1L else position
    if (!is.na(i) && i >= 1 && i <= length(p)) p[i] else NA_character_
  }, character(1))
}

#' Sample sheet restricted to included samples, with typed columns.
design_coldata <- function(sheet, included, variables, types) {
  cd <- sheet[sheet$sample %in% included, c("sample", variables), drop = FALSE]
  for (v in variables) {
    x <- trimws(as.character(cd[[v]]))
    x[!is.na(x) & x == ""] <- NA
    cd[[v]] <- if (identical(types[[v]], "numeric")) suppressWarnings(as.numeric(x)) else factor(x)
  }
  rownames(cd) <- cd$sample
  cd
}

#' Model formula: covariates first, the main variable last (DESeq2 convention).
design_formula <- function(main, covariates = character(0)) {
  stats::as.formula(paste("~", paste(c(covariates, main), collapse = " + ")))
}

#' Default identifier for a comparison.
contrast_id <- function(variable, target, reference) {
  make.names(sprintf("%s_%s_vs_%s", variable, target, reference))
}

empty_contrasts <- function() {
  data.frame(id = character(0), variable = character(0), reference = character(0),
             target = character(0), stringsAsFactors = FALSE)
}

# ---- design checks -------------------------------------------------------------

#' Check a design before running.
#'
#' @return list(errors = chr, warnings = chr)
validate_design <- function(sheet, included, main, covariates, types, contrasts) {
  err <- character(0); wrn <- character(0)
  covariates <- setdiff(if (is.null(covariates)) character(0) else covariates, main)
  if (length(included) < 2) err <- c(err, "Include at least 2 samples.")
  if (is.null(main) || !nzchar(main)) {
    return(list(errors = c(err, "Choose the main variable (the one your comparisons are between)."),
                warnings = wrn))
  }
  vars <- unique(c(covariates, main))
  missing_vars <- setdiff(vars, colnames(sheet))
  if (length(missing_vars)) {
    return(list(errors = c(err, paste("Unknown variable(s):", paste(missing_vars, collapse = ", "))),
                warnings = wrn))
  }
  if (identical(types[[main]], "numeric")) {
    err <- c(err, sprintf("The main variable '%s' is numeric; comparisons need a factor.", main))
  }
  if (length(err)) return(list(errors = err, warnings = wrn))
  cd <- design_coldata(sheet, included, vars, types)

  for (v in vars) {
    n_na <- sum(is.na(cd[[v]]))
    if (n_na) err <- c(err, sprintf("'%s' has %d included sample(s) with no value%s.", v, n_na,
                                    if (identical(types[[v]], "numeric")) " (or a value that is not a number)" else ""))
  }
  if (length(err)) return(list(errors = err, warnings = wrn))

  grp <- table(droplevels(cd[[main]]))
  if (length(grp) < 2) err <- c(err, sprintf("'%s' has only one level among the included samples.", main))
  one <- names(grp)[grp == 1]
  if (length(one)) wrn <- c(wrn, sprintf(
    "Level(s) of '%s' with a single sample: %s. Comparisons involving them are unreliable.",
    main, paste(one, collapse = ", ")))

  for (v in covariates) {
    if (identical(types[[v]], "numeric")) {
      if (length(unique(cd[[v]])) < 2) err <- c(err, sprintf("Covariate '%s' has the same value for every sample.", v))
      next
    }
    if (nlevels(droplevels(cd[[v]])) < 2) {
      err <- c(err, sprintf("Covariate '%s' has only one level among the included samples; remove it.", v))
      next
    }
    tab <- table(droplevels(cd[[main]]), droplevels(cd[[v]]))
    if (all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)) {
      err <- c(err, sprintf(paste(
        "'%s' is completely confounded with '%s': each level of one matches exactly one level of the other,",
        "so their effects cannot be separated. Remove '%s' from the covariates."), v, main, v))
    }
  }
  if (length(err)) return(list(errors = err, warnings = wrn))

  mm <- stats::model.matrix(design_formula(main, covariates), cd)
  if (qr(mm)$rank < ncol(mm)) {
    err <- c(err, paste("The design matrix is not full rank: some combination of variables is",
                        "confounded or has an empty group. Remove a covariate or include more samples."))
  } else if (nrow(mm) - ncol(mm) < 1) {
    err <- c(err, sprintf(paste("No residual degrees of freedom: %d samples for %d model terms.",
                                "Variability cannot be estimated without replicates."),
                          nrow(mm), ncol(mm)))
  }

  if (is.null(contrasts) || !nrow(contrasts)) {
    err <- c(err, "Add at least one comparison.")
  } else {
    lv <- levels(droplevels(cd[[main]]))
    for (i in seq_len(nrow(contrasts))) {
      cr <- contrasts[i, ]
      if (!identical(cr$variable, main)) {
        err <- c(err, sprintf("Comparison '%s' uses '%s', not the main variable '%s'.", cr$id, cr$variable, main))
      } else if (!all(c(cr$reference, cr$target) %in% lv)) {
        err <- c(err, sprintf("Comparison '%s' uses a level with no included samples.", cr$id))
      } else if (identical(cr$reference, cr$target)) {
        err <- c(err, sprintf("Comparison '%s' compares a level with itself.", cr$id))
      }
    }
    if (anyDuplicated(contrasts$id)) err <- c(err, "Two comparisons share the same name.")
  }
  list(errors = err, warnings = wrn)
}

# ---- filtering -------------------------------------------------------------------

#' Abundance filter: a minimum count (or CPM) in enough samples.
#'
#' A gene is kept when at least `n_required` samples reach `min_abundance`
#' (counts, or CPM when `use_cpm`). `n_required` comes from `mode`:
#'   "samples"    -- `min_samples`
#'   "proportion" -- ceiling(`min_proportion` x number of samples)
#'   "group"      -- size of the smallest group of `group`
#' @return list(keep = logical, n_required = int, description = chr)
filter_counts <- function(counts, min_abundance = 10, mode = "group", min_samples = 1,
                          min_proportion = 0.5, group = NULL, use_cpm = FALSE) {
  n <- ncol(counts)
  n_req <- switch(mode,
    samples    = min_samples,
    proportion = ceiling(min_proportion * n),
    group      = if (is.null(group)) min_samples else min(table(droplevels(as.factor(group)))),
    stop("Unknown filtering mode: ", mode))
  n_req <- max(1L, min(as.integer(n_req), n))
  vals <- if (use_cpm) t(t(counts) / pmax(colSums(counts, na.rm = TRUE), 1) * 1e6) else counts
  keep <- rowSums(vals >= min_abundance, na.rm = TRUE) >= n_req
  list(keep = keep, n_required = n_req,
       description = sprintf("%s >= %g in at least %d of %d samples",
                             if (use_cpm) "CPM" else "count", min_abundance, n_req, n))
}

# ---- engines ---------------------------------------------------------------------

#' Default options for each engine.
de_default_options <- function() {
  list(
    # filtering
    min_abundance = 10, filter_mode = "group", min_samples = 1, min_proportion = 0.5,
    use_cpm = FALSE,
    # DESeq2
    fit_type = "parametric", alpha = 0.1, lfc_threshold = 0,
    independent_filtering = TRUE, cooks_cutoff = TRUE,
    transform = "vst", blind = TRUE,
    deseq2_test = "wald",            # "wald" or "lrt" (reduced model drops the main variable)
    shrink = "none",                 # "none", "normal", "apeglm" or "ashr" (lfcShrink)
    # gene-length offsets (DESeq2 and edgeR; needs sample-specific lengths, e.g. RSEM / Salmon)
    length_offset = FALSE,
    # transcript (isoform)-level analysis
    feature_level = "gene",          # "gene" or "transcript"
    divide_overdispersion = TRUE,    # edgeR / limma: counts / RTA overdispersion (Baldoni et al. 2024)
    ebseq_ng = "isoforms",           # EBSeq Ng vector: "isoforms" (GetNg), "file" or "none"
    ebseq_trunthre = 3,              # GetNg TrunThre: Ng groups 1, 2, 3+ isoforms
    stager = TRUE,                   # stage-wise testing (stageR, method "dte" / "dtu")
    stager_alpha = 0.05,
    # differential transcript usage (feature_level = "dtu")
    dtu_counts = "auto",             # "auto", "divided", "scaledTPM", "dtuScaledTPM", "estimated"
    dtu_filter = "auto",             # "auto", "drimseq" (dmFilter) or "edger" (filterByExpr)
    dtu_min_gene_expr = 10, dtu_min_feature_expr = 10, dtu_min_feature_prop = 0.1,
    dtu_gene_test = "simes",         # edgeR / limma diffSplice gene-level test: "simes" or "F"
    dtu_robust = TRUE,               # edgeR glmQLFit / limma diffSplice robust
    drimseq_sd_filter = 0.1,         # DRIMSeq post-hoc filter (0 = off)
    saturn_pvalues = "empirical",    # satuRn: "empirical" or "regular"
    # limma-voom and edgeR
    norm_method = "TMM",
    limma_voom = "voom", limma_trend = FALSE, limma_robust = FALSE, limma_block = "(none)",
    limma_lfc = 0,
    edger_test = "ql", edger_robust = TRUE, edger_lfc = 0,
    # EBSeq
    ebseq_maxround = 5
  )
}

#' Stage-wise testing of differential transcript expression (stageR, Van den
#' Berge et al. 2017): genes are screened first (Simes combination of their
#' transcripts' p-values, BH across genes), then transcripts are confirmed
#' within significant genes (Holm, method "dte"). This controls the FDR at the
#' gene level and reports which transcripts drive each gene.
#' @return list(df with stageR_gene_padj / stageR_tx_padj, note)
stage_wise_dte <- function(df, alpha = 0.05) {
  if (!requireNamespace("stageR", quietly = TRUE)) {
    return(list(df = df, note = paste(
      "Stage-wise testing was skipped: the stageR package is not installed",
      '(BiocManager::install("stageR")).')))
  }
  ok <- !is.na(df$pvalue) & !is.na(df$parent_gene_id)
  p <- stats::setNames(df$pvalue[ok], df$gene_id[ok])
  g <- df$parent_gene_id[ok]
  sc <- tapply(p, g, function(x) { x <- sort(x); min(1, min(length(x) * x / seq_along(x))) })
  screen <- stats::setNames(as.numeric(sc), names(sc))   # stageR needs a plain numeric vector
  conf <- matrix(p, ncol = 1, dimnames = list(names(p), "transcript"))
  obj <- stageR::stageRTx(pScreen = screen, pConfirmation = conf, pScreenAdjusted = FALSE,
                          tx2gene = data.frame(transcript = names(p), gene = g, stringsAsFactors = FALSE))
  obj <- suppressMessages(suppressWarnings(stageR::stageWiseAdjustment(obj, method = "dte", alpha = alpha)))
  adj <- as.data.frame(suppressMessages(stageR::getAdjustedPValues(obj, order = FALSE, onlySignificantGenes = FALSE)))
  i <- match(df$gene_id, adj$txID)
  df$stageR_gene_padj <- suppressWarnings(as.numeric(adj$gene[i]))
  df$stageR_tx_padj <- suppressWarnings(as.numeric(adj$transcript[i]))
  n_g <- length(unique(df$parent_gene_id[!is.na(df$stageR_gene_padj) & df$stageR_gene_padj <= alpha]))
  n_t <- sum(!is.na(df$stageR_tx_padj) & df$stageR_tx_padj <= alpha)
  list(df = df, note = sprintf(paste(
    "Stage-wise testing (stageR %s, alpha = %g): %s gene(s) pass the gene-level screen (Simes, BH) and",
    "%s transcript(s) are confirmed within them (Holm). stageR_gene_padj and stageR_tx_padj are valid",
    "at this alpha only; a transcript is significant when both are <= alpha."),
    as.character(utils::packageVersion("stageR")), alpha, format(n_g, big.mark = ","), format(n_t, big.mark = ",")))
}

#' Recommended DESeq2 settings: a gene is kept with at least 10 counts in at least
#' 1 sample, gene-length offsets when sample-specific lengths are loaded, Wald test,
#' and ashr-shrunken fold changes.
recommended_options <- function() {
  utils::modifyList(de_default_options(), list(
    min_abundance = 10, filter_mode = "samples", min_samples = 1, use_cpm = FALSE,
    fit_type = "parametric", alpha = 0.1, lfc_threshold = 0, independent_filtering = TRUE,
    cooks_cutoff = TRUE, deseq2_test = "wald", shrink = "ashr", length_offset = TRUE,
    transform = "vst", blind = TRUE))
}

#' Run differential expression for every comparison.
#'
#' @param counts     gene x sample matrix (restricted here to `included`)
#' @param contrasts  data frame: id, variable, reference, target
#' @param progress   optional function(fraction, detail) for a progress bar
#' @return list(engine, results = named list of data frames, transformed,
#'              transform_label, normalised, kept, formula, filter, notes,
#'              coldata, main, covariates, contrasts, settings, versions, run_at)
run_de_analysis <- function(counts, sheet, included, main, covariates = character(0),
                            types = list(), contrasts, engine = "deseq2",
                            options = de_default_options(), progress = NULL, lengths = NULL,
                            transcripts = NULL) {
  # literal list, so gene-level runs work even when dtu.R is not loaded
  dtu_keys <- c("edger_dtu", "limma_dtu", "dexseq", "drimseq", "saturn")
  if (engine %in% dtu_keys || identical(options$feature_level, "dtu")) {
    if (!engine %in% dtu_keys) stop("Transcript usage needs a DTU method: ", paste(names(DTU_ENGINES), collapse = ", "), ".")
    return(run_dtu_analysis(counts, sheet, included, main, covariates, types, contrasts, engine,
                            options, progress, lengths, transcripts))
  }
  tick <- function(f, d) if (is.function(progress)) progress(f, d)
  options <- utils::modifyList(de_default_options(), options)
  miss <- de_engine_missing(engine)
  if (!is.null(miss)) {
    stop("The ", names(DE_ENGINES)[DE_ENGINES == engine], " package is not installed. Run: ", miss)
  }
  covariates <- setdiff(if (is.null(covariates)) character(0) else covariates, main)
  chk <- validate_design(sheet, included, main, covariates, types, contrasts)
  chk$errors <- c(chk$errors, engine_design_errors(engine, covariates, types, options))
  if (length(chk$errors)) stop(paste(chk$errors, collapse = "\n"))

  vars <- c(covariates, main)
  cd <- design_coldata(sheet, included, vars, types)
  for (v in vars) if (is.factor(cd[[v]])) cd[[v]] <- droplevels(cd[[v]])

  notes <- chk$warnings
  mat <- counts[, cd$sample, drop = FALSE]
  storage.mode(mat) <- "double"
  if (any(is.na(mat))) { mat[is.na(mat)] <- 0; notes <- c(notes, "Missing counts were set to 0.") }
  level <- if (identical(options$feature_level, "transcript")) "transcript" else "gene"
  feat <- if (level == "transcript") "transcripts" else "genes"

  # Transcript level, edgeR / limma: divide out the read-to-transcript ambiguity
  # (Baldoni et al. 2024). Divided counts are fractional and analysed as such.
  divided <- FALSE
  if (level == "transcript" && engine %in% c("edger", "limma") && isTRUE(options$divide_overdispersion)) {
    od <- transcripts$overdispersion
    if (!is.null(od)) {
      o <- unname(od[rownames(mat)]); o[is.na(o)] <- 1
      mat <- mat / pmax(o, 1)
      divided <- TRUE
      notes <- c(notes, sprintf(paste(
        "Counts were divided by each transcript's quantification overdispersion (read-to-transcript",
        "ambiguity, from the resamples; median %.2f, %s%% of transcripts above 2), as edgeR's",
        "divided-count method recommends (Baldoni et al. 2024). Divided counts are fractional and are not rounded."),
        stats::median(o), format(round(100 * mean(o > 2), 1))))
    } else {
      notes <- c(notes, paste(
        "No quantification-uncertainty estimates are loaded, so transcript counts were used as they are.",
        "Load the Salmon / kallisto folders with resamples (Gibbs or bootstrap) to divide them out."))
    }
  }
  if (!divided) {
    if (any(abs(mat - round(mat)) > 1e-8)) {
      notes <- c(notes, "Non-integer counts (e.g. RSEM or Salmon expected counts) were rounded for the count model.")
    }
    mat <- round(mat)
  }
  qc <- list(lib_size = colSums(mat), detected = colSums(mat > 0))

  tick(0.05, sprintf("Filtering low-abundance %s", feat))
  flt <- filter_counts(mat, options$min_abundance, options$filter_mode, options$min_samples,
                       options$min_proportion, group = cd[[main]], use_cpm = isTRUE(options$use_cpm))
  mat <- mat[flt$keep, , drop = FALSE]
  if (nrow(mat) < 10) stop("Only ", nrow(mat), " ", feat, " pass the filter; loosen it.")

  # Transcript -> gene map for every input transcript (EBSeq's Ng, stageR, annotation)
  t2g <- NULL; ng <- NULL
  if (level == "transcript") {
    t2g <- merge_tx2gene(rownames(counts), transcripts$tx2gene)
    if (all(is.na(t2g$gene_id))) notes <- c(notes, paste(
      "No transcript-to-gene map is loaded: results have no gene symbols, EBSeq cannot group isoforms",
      "and stage-wise testing is skipped. Load a GTF or a tx2gene table."))
    if (engine == "ebseq") {
      mode <- options$ebseq_ng %||% "isoforms"
      if (identical(mode, "file")) {
        if (is.null(transcripts$ngvec)) stop("Choose an Ng vector file for EBSeq, or use 'isoforms per gene'.")
        ng_all <- transcripts$ngvec
      } else if (identical(mode, "isoforms") && !all(is.na(t2g$gene_id))) {
        g <- ifelse(is.na(t2g$gene_id), t2g$transcript_id, t2g$gene_id)
        tt <- min(as.integer(options$ebseq_trunthre %||% 3), max(table(g)))
        ng_all <- EBSeq::GetNg(t2g$transcript_id, g, TrunThre = tt)$IsoformNgTrun
        notes <- c(notes, sprintf(paste(
          "EBSeq isoform mode: transcripts are grouped by the number of isoforms of their gene (GetNg,",
          "TrunThre = %d): %s."), tt, paste(sprintf("Ng=%s: %s", names(table(ng_all)),
          format(as.integer(table(ng_all)), big.mark = ",")), collapse = "; ")))
      } else ng_all <- NULL
      if (!is.null(ng_all)) {
        ng <- as.integer(ng_all[rownames(mat)])
        if (anyNA(ng)) { notes <- c(notes, sprintf("%d transcript(s) had no Ng value and were put in group 1.", sum(is.na(ng)))); ng[is.na(ng)] <- 1L }
        names(ng) <- rownames(mat)
      }
    }
  }

  lens <- NULL
  if (isTRUE(options$length_offset)) {
    lo <- length_offset_matrix(lengths, rownames(mat), colnames(mat))
    notes <- c(notes, lo$note)
    if (engine %in% c("deseq2", "edger")) lens <- lo$matrix else if (!is.null(lo$matrix))
      notes <- c(notes, sprintf("Gene-length offsets are not used by %s; they apply to DESeq2 and edgeR only.",
                                names(DE_ENGINES)[DE_ENGINES == engine]))
  }

  out <- switch(engine,
    deseq2 = .run_deseq2(mat, cd, main, covariates, types, contrasts, options, tick, lens),
    limma  = .run_limma(mat, cd, main, covariates, types, contrasts, options, tick),
    edger  = .run_edger(mat, cd, main, covariates, types, contrasts, options, tick, lens),
    ebseq  = .run_ebseq(mat, cd, main, contrasts, options, tick, ng),
    stop("Unknown engine: ", engine))

  # Transcript level: gene columns, and stage-wise testing of each comparison.
  if (level == "transcript") {
    tk <- t2g[match(rownames(mat), t2g$transcript_id), , drop = FALSE]
    st_note <- NULL
    for (id in names(out$results)) {
      df <- out$results[[id]]
      a <- tk[match(df$gene_id, tk$transcript_id), , drop = FALSE]
      df <- data.frame(gene_id = df$gene_id, parent_gene_id = a$gene_id, gene_name = a$gene_name,
                       transcript_name = a$transcript_name, transcript_biotype = a$transcript_biotype,
                       df[, -1, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
      if (isTRUE(options$stager) && engine != "ebseq" && !all(is.na(df$parent_gene_id))) {
        st <- stage_wise_dte(df, options$stager_alpha %||% 0.05)
        df <- st$df; st_note <- st$note
      }
      out$results[[id]] <- df
    }
    if (isTRUE(options$stager) && engine != "ebseq" && requireNamespace("stageR", quietly = TRUE))
      out$versions <- c(out$versions, stageR = as.character(utils::packageVersion("stageR")))
    if (isTRUE(options$stager) && engine == "ebseq")
      st_note <- "Stage-wise testing (stageR) needs p-values; EBSeq gives posterior probabilities, so it was skipped."
    notes <- c(notes, st_note)
  }

  tick(1, "Done")
  # Fields the engine returns that are merged with the shared ones below.
  out_qc <- out$qc; out_notes <- out$notes; out$qc <- NULL; out$notes <- NULL
  c(out, list(
    kept = rownames(mat), n_input = nrow(counts), counts = mat,
    formula = paste(deparse(design_formula(main, covariates)), collapse = ""),
    filter = flt$description, notes = c(notes, out_notes), coldata = cd, main = main,
    qc = c(qc, out_qc),
    covariates = covariates, contrasts = contrasts, settings = options,
    length_offset = !is.null(lens), lengths_used = lens,
    level = level, divided = divided, tx2gene = t2g, ngvector = ng,
    overdispersion = if (divided) transcripts$overdispersion else NULL,
    # what the run record, the R script and the bundle need to reproduce the run
    input_counts = counts[, cd$sample, drop = FALSE], filter_n = flt$n_required,
    sheet = sheet[match(cd$sample, sheet$sample), , drop = FALSE], types = types,
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

.run_deseq2 <- function(mat, cd, main, covariates, types, contrasts, options, tick, lens = NULL) {
  notes <- character(0)
  for (v in covariates) if (identical(types[[v]], "numeric")) {
    cd[[v]] <- as.numeric(scale(cd[[v]]))
    notes <- c(notes, sprintf("Numeric covariate '%s' was centred and scaled, as DESeq2 recommends.", v))
  }
  lrt <- identical(options$deseq2_test, "lrt")
  reduced <- if (length(covariates)) stats::reformulate(covariates) else ~ 1
  tick(0.15, "Building the DESeq2 data set")
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(
    countData = mat, colData = cd, design = design_formula(main, covariates)))
  if (!is.null(lens)) {
    # DESeq2 turns an 'avgTxLength' assay into gene- and sample-specific normalisation
    # factors (as for tximport input) when it estimates size factors.
    SummarizedExperiment::assay(dds, "avgTxLength", withDimnames = FALSE) <- lens
    notes <- c(notes, "Gene-length offsets: DESeq2 normalisation factors combine library size and each sample's gene lengths (as with tximport).")
  }
  tick(0.25, "Size factors, dispersions and model fit")
  dds <- if (lrt) {
    suppressMessages(DESeq2::DESeq(dds, test = "LRT", reduced = reduced,
                                   fitType = options$fit_type, quiet = TRUE))
  } else {
    suppressMessages(DESeq2::DESeq(dds, test = "Wald", fitType = options$fit_type, quiet = TRUE))
  }
  lfc_thr <- options$lfc_threshold
  if (lrt) {
    notes <- c(notes, sprintf(paste(
      "Likelihood-ratio test of %s against %s: the p-value asks whether '%s' matters at all, so it is the",
      "same for every comparison; each comparison's fold change comes from the full model."),
      paste(deparse(design_formula(main, covariates)), collapse = ""),
      paste(deparse(reduced), collapse = ""), main))
    if (lfc_thr > 0) {
      notes <- c(notes, "The |log2 FC| test threshold applies to Wald tests only and was ignored for the LRT.")
      lfc_thr <- 0
    }
  }

  shrink <- options$shrink %||% "none"
  if (!shrink %in% c("normal", "apeglm", "ashr")) shrink <- "none"
  relevelled <- list()      # apeglm needs a coefficient: refit once per reference level
  shrink_one <- function(r, cr) {
    if (shrink == "apeglm") {
      coef <- .deseq2_coef(dds, main, cr$target, cr$reference)
      d <- dds
      if (is.null(coef)) {
        key <- cr$reference
        if (is.null(relevelled[[key]])) {
          d2 <- dds
          d2[[main]] <- stats::relevel(d2[[main]], ref = cr$reference)
          relevelled[[key]] <<- suppressMessages(DESeq2::nbinomWaldTest(d2, quiet = TRUE))
        }
        d <- relevelled[[key]]
        coef <- .deseq2_coef(d, main, cr$target, cr$reference)
      }
      if (is.null(coef)) stop("Could not find the coefficient for ", cr$id, " for apeglm shrinkage.")
      suppressMessages(DESeq2::lfcShrink(d, coef = coef, res = r, type = "apeglm", quiet = TRUE))
    } else {
      suppressMessages(DESeq2::lfcShrink(dds, contrast = c(main, cr$target, cr$reference), res = r,
                                         type = shrink, quiet = TRUE))
    }
  }

  res <- list()
  n <- nrow(contrasts)
  for (i in seq_len(n)) {
    cr <- contrasts[i, ]
    tick(0.55 + 0.3 * (i - 1) / n, paste("Testing", cr$id))
    r <- DESeq2::results(dds, contrast = c(main, cr$target, cr$reference),
                         alpha = options$alpha, lfcThreshold = lfc_thr,
                         independentFiltering = isTRUE(options$independent_filtering),
                         cooksCutoff = isTRUE(options$cooks_cutoff))
    df <- data.frame(gene_id = rownames(r), as.data.frame(r), stringsAsFactors = FALSE,
                     check.names = FALSE)
    if (shrink != "none") {
      tick(0.55 + 0.3 * (i - 0.5) / n, paste("Shrinking fold changes:", cr$id))
      sr <- shrink_one(r, cr)
      df$log2FoldChange_mle <- df$log2FoldChange
      df$lfcSE_mle <- df$lfcSE
      df$log2FoldChange <- as.numeric(sr$log2FoldChange)
      df$lfcSE <- if ("lfcSE" %in% colnames(sr)) as.numeric(sr$lfcSE) else NA_real_
    }
    df <- df[order(df$pvalue, na.last = TRUE), , drop = FALSE]
    rownames(df) <- NULL
    res[[cr$id]] <- df
  }
  if (shrink != "none") notes <- c(notes, sprintf(paste(
    "log2FoldChange and lfcSE are %s-shrunken (lfcShrink); p-values are unchanged.",
    "The unshrunken estimates are kept as log2FoldChange_mle and lfcSE_mle."), shrink),
    if (shrink == "apeglm") paste(
      "apeglm is not symmetric: for low-count genes (especially with zeros in one group) the shrunken",
      "fold change depends on which group is the reference, so A vs B is not exactly minus B vs A.",
      "ashr and normal shrinkage are symmetric."))

  tick(0.88, sprintf("Computing the %s transform", if (identical(options$transform, "rlog")) "rlog" else "VST"))
  tr <- if (identical(options$transform, "rlog")) {
    DESeq2::rlog(dds, blind = isTRUE(options$blind))
  } else if (nrow(dds) >= 1000) {
    DESeq2::vst(dds, blind = isTRUE(options$blind))
  } else {
    notes <- c(notes, "Fewer than 1,000 genes: used varianceStabilizingTransformation() instead of vst().")
    DESeq2::varianceStabilizingTransformation(dds, blind = isTRUE(options$blind))
  }
  mc <- S4Vectors::mcols(dds)
  disp <- data.frame(baseMean = mc$baseMean, gene_est = mc$dispGeneEst, fit = mc$dispFit,
                     final = DESeq2::dispersions(dds),
                     outlier = if ("dispOutlier" %in% colnames(mc)) mc$dispOutlier else FALSE,
                     row.names = rownames(dds))
  sf <- DESeq2::sizeFactors(dds)
  if (is.null(sf)) {   # normalisation factors (length offsets): summarise per sample
    nf <- DESeq2::normalizationFactors(dds)
    sf <- exp(colMeans(log(nf)))
    sf <- sf / exp(mean(log(sf)))
  }
  vers <- c(DESeq2 = as.character(utils::packageVersion("DESeq2")))
  if (shrink %in% c("apeglm", "ashr")) vers[shrink] <- as.character(utils::packageVersion(shrink))
  list(engine = "DESeq2", results = res,
       qc = list(size_factors = sf, dispersion = disp),
       transformed = SummarizedExperiment::assay(tr),
       transform_label = if (identical(options$transform, "rlog")) "rlog" else "VST",
       normalised = DESeq2::counts(dds, normalized = TRUE),
       size_factors = sf,
       notes = notes,
       versions = vers)
}

#' Name of the DESeq2 coefficient for target vs reference, or NULL.
.deseq2_coef <- function(dds, main, target, reference) {
  rn <- DESeq2::resultsNames(dds)
  want <- make.names(paste0(main, "_", target, "_vs_", reference))
  hit <- rn[make.names(rn) == want]
  if (length(hit)) hit[1] else NULL
}

#' Gene-length matrix for offsets, aligned to `genes` x `samples`.
#'
#' Offsets only change the result when lengths differ between samples (RSEM or Salmon
#' effective lengths); one length per gene cancels out of every comparison.
#' Missing or non-positive lengths are replaced by the gene's geometric mean over the
#' other samples, or by 1 (no correction) when a gene has none.
#' @return list(matrix = NULL or genes x samples, note)
length_offset_matrix <- function(lengths, genes, samples) {
  if (is.null(lengths)) return(list(matrix = NULL,
    note = "Gene-length offsets were requested but no gene lengths are loaded; they were not used."))
  if (is.null(dim(lengths))) return(list(matrix = NULL, note = paste(
    "Gene-length offsets were not used: the loaded lengths are one value per gene (e.g. from a GTF),",
    "which cancels out of every comparison. Offsets need sample-specific lengths such as",
    "rsem.merged.gene_lengths.tsv or Salmon's effective lengths.")))
  miss_s <- setdiff(samples, colnames(lengths))
  if (length(miss_s)) return(list(matrix = NULL, note = sprintf(
    "Gene-length offsets were not used: the lengths file has no column for %s.",
    paste(utils::head(miss_s, 3), collapse = ", "))))
  m <- matrix(NA_real_, length(genes), length(samples), dimnames = list(genes, samples))
  hit <- intersect(genes, rownames(lengths))
  m[hit, ] <- suppressWarnings(as.numeric(as.matrix(lengths[hit, samples, drop = FALSE])))
  m[!is.finite(m) | m <= 0] <- NA
  if (all(apply(m, 1, function(x) all(is.na(x)) || stats::sd(x, na.rm = TRUE) %in% c(0, NA)))) {
    return(list(matrix = NULL, note = paste(
      "Gene-length offsets were not used: every gene has the same length in every sample,",
      "so they would cancel out of each comparison.")))
  }
  n_none <- sum(rowSums(!is.na(m)) == 0)
  gm <- exp(rowMeans(log(m), na.rm = TRUE)); gm[!is.finite(gm)] <- 1
  idx <- which(is.na(m), arr.ind = TRUE)
  if (nrow(idx)) m[idx] <- gm[idx[, 1]]
  m[rowSums(!is.na(m)) == 0, ] <- 1
  list(matrix = m, note = sprintf(
    "Gene-length offsets from sample-specific lengths%s.",
    if (n_none) sprintf("; %s gene(s) without a length were left uncorrected",
                        format(n_none, big.mark = ",")) else ""))
}

#' Results file name: <comparison>.<engine>.results.tsv
de_result_filename <- function(id, engine) {
  sprintf("%s.%s.results.tsv", id, tolower(engine))
}

#' Up / down / tested counts for one results table at the app's thresholds.
de_result_counts <- function(df, padj_cutoff = 0.05, lfc_cutoff = 1) {
  ok <- !is.na(df$padj) & df$padj <= padj_cutoff
  c(tested = sum(!is.na(df$pvalue)),
    up   = sum(ok & df$log2FoldChange >=  lfc_cutoff, na.rm = TRUE),
    down = sum(ok & df$log2FoldChange <= -lfc_cutoff, na.rm = TRUE))
}

# ---- limma-voom and edgeR ---------------------------------------------------------

#' Cell-means design (~ 0 + main + covariates): one column per level of the
#' main variable, so any two levels form a simple +1 / -1 contrast.
.cell_means_design <- function(cd, main, covariates) {
  f <- stats::as.formula(paste("~ 0 +", paste(c(main, covariates), collapse = " + ")))
  mm <- stats::model.matrix(f, cd)
  lv <- levels(cd[[main]])
  cols <- paste0(main, lv)
  if (!all(cols %in% colnames(mm))) stop("Could not build the design for '", main, "'.")
  list(design = mm, col = stats::setNames(cols, lv))
}

.contrast_vector <- function(dm, reference, target) {
  v <- stats::setNames(numeric(ncol(dm$design)), colnames(dm$design))
  v[dm$col[[target]]] <- 1; v[dm$col[[reference]]] <- -1
  v
}

#' TMM-style normalisation shared by limma-voom and edgeR.
.edger_normalise <- function(mat, method, lens = NULL) {
  # edgeR >= 4.x renamed calcNormFactors() to normLibSizes() and prints a notice for the old name
  norm_fun <- if (exists("normLibSizes", envir = asNamespace("edgeR"), inherits = FALSE))
    getExportedValue("edgeR", "normLibSizes") else edgeR::calcNormFactors
  y <- edgeR::DGEList(counts = mat)
  if (is.null(lens)) {
    y <- norm_fun(y, method = method)
    eff <- y$samples$lib.size * y$samples$norm.factors
    sf <- stats::setNames(eff / exp(mean(log(eff))), colnames(mat))
    return(list(y = y, size_factors = sf, normalised = t(t(mat) / sf)))
  }
  # tximport's recipe: length-correct, normalise, then use log(length x library) as the offset
  norm_mat <- lens / exp(rowMeans(log(lens)))
  norm_cts <- mat / norm_mat
  eff_lib <- norm_fun(edgeR::DGEList(norm_cts), method = method)
  eff_lib <- eff_lib$samples$norm.factors * colSums(norm_cts)
  off <- log(sweep(norm_mat, 2, eff_lib, "*"))
  y <- edgeR::scaleOffset(y, off)
  sf <- stats::setNames(eff_lib / exp(mean(log(eff_lib))), colnames(mat))
  list(y = y, size_factors = sf, normalised = norm_cts / rep(sf, each = nrow(norm_cts)))
}

.run_limma <- function(mat, cd, main, covariates, types, contrasts, options, tick) {
  notes <- character(0)
  block <- options$limma_block
  if (is.null(block) || !nzchar(block) || block == "(none)") block <- NULL
  fixed <- setdiff(covariates, block)
  tick(0.12, sprintf("Normalising (%s)", options$norm_method))
  nz <- .edger_normalise(mat, options$norm_method)
  dm <- .cell_means_design(cd, main, fixed)
  quality <- identical(options$limma_voom, "quality")
  vfun <- function(...) if (quality) limma::voomWithQualityWeights(nz$y, dm$design, ...) else
    limma::voom(nz$y, dm$design, ...)
  tick(0.25, "voom: mean-variance trend and precision weights")
  vqc <- limma::voom(nz$y, dm$design, save.plot = TRUE)          # for the QC plot
  if (is.null(block)) {
    v <- vfun()
    fit <- limma::lmFit(v, dm$design)
  } else {
    tick(0.35, sprintf("duplicateCorrelation on '%s'", block))
    b <- as.character(cd[[block]])
    v <- vfun()
    dc <- limma::duplicateCorrelation(v, dm$design, block = b)
    v <- vfun(block = b, correlation = dc$consensus.correlation)
    dc <- limma::duplicateCorrelation(v, dm$design, block = b)
    fit <- limma::lmFit(v, dm$design, block = b, correlation = dc$consensus.correlation)
    notes <- c(notes, sprintf("'%s' was modelled as a random block (duplicateCorrelation; consensus correlation %.2f).",
                              block, dc$consensus.correlation))
  }
  if (quality) notes <- c(notes, "Sample quality weights were estimated (voomWithQualityWeights).")
  base <- rowMeans(nz$normalised)
  res <- list(); n <- nrow(contrasts)
  for (i in seq_len(n)) {
    cr <- contrasts[i, ]
    tick(0.5 + 0.35 * (i - 1) / n, paste("Testing", cr$id))
    f2 <- limma::contrasts.fit(fit, .contrast_vector(dm, cr$reference, cr$target))
    if (options$limma_lfc > 0) {
      f2 <- limma::treat(f2, lfc = options$limma_lfc, trend = isTRUE(options$limma_trend),
                         robust = isTRUE(options$limma_robust))
      tt <- limma::topTreat(f2, coef = 1, number = Inf, sort.by = "none")
    } else {
      f2 <- limma::eBayes(f2, trend = isTRUE(options$limma_trend), robust = isTRUE(options$limma_robust))
      tt <- limma::topTable(f2, coef = 1, number = Inf, sort.by = "none")
    }
    se <- as.numeric(f2$stdev.unscaled[, 1] * sqrt(f2$s2.post))
    df <- data.frame(gene_id = rownames(tt), baseMean = unname(base[rownames(tt)]),
                     log2FoldChange = tt$logFC, lfcSE = se, stat = tt$t, pvalue = tt$P.Value,
                     padj = tt$adj.P.Val, AveExpr = tt$AveExpr, stringsAsFactors = FALSE)
    if (!is.null(tt$B)) df$B <- tt$B
    df <- df[order(df$pvalue, na.last = TRUE), , drop = FALSE]; rownames(df) <- NULL
    res[[cr$id]] <- df
  }
  if (options$limma_lfc > 0) notes <- c(notes, sprintf(
    "Tested against |log2 FC| > %g with limma::treat (p-values are for exceeding the threshold).", options$limma_lfc))
  list(engine = "limma", engine_label = "limma-voom", results = res,
       transformed = v$E, transform_label = "voom log2-CPM", normalised = nz$normalised,
       qc = list(size_factors = nz$size_factors,
                 voom = list(x = vqc$voom.xy$x, y = vqc$voom.xy$y,
                             line_x = vqc$voom.line$x, line_y = vqc$voom.line$y),
                 sample_weights = if (quality) { w <- v$targets$sample.weights; if (is.null(w)) w <- v$sample.weights; if (!is.null(w)) setNames(as.numeric(w), colnames(v$E)) else NULL } else NULL),
       notes = notes,
       versions = c(limma = as.character(utils::packageVersion("limma")),
                    edgeR = as.character(utils::packageVersion("edgeR"))))
}

.run_edger <- function(mat, cd, main, covariates, types, contrasts, options, tick, lens = NULL) {
  notes <- character(0)
  tick(0.12, sprintf("Normalising (%s)", options$norm_method))
  nz <- .edger_normalise(mat, options$norm_method, lens)
  if (!is.null(lens)) notes <- c(notes,
    "Gene-length offsets: edgeR offsets combine the normalised library sizes with each sample's gene lengths (tximport's recipe, scaleOffset).")
  dm <- .cell_means_design(cd, main, covariates)
  tick(0.25, "Estimating dispersions")
  y <- edgeR::estimateDisp(nz$y, dm$design, robust = isTRUE(options$edger_robust))
  ql <- !identical(options$edger_test, "lrt")
  tick(0.45, if (ql) "Quasi-likelihood fit" else "GLM fit")
  fit <- if (ql) edgeR::glmQLFit(y, dm$design, robust = isTRUE(options$edger_robust)) else
    edgeR::glmFit(y, dm$design)
  base <- rowMeans(nz$normalised)
  res <- list(); n <- nrow(contrasts)
  for (i in seq_len(n)) {
    cr <- contrasts[i, ]
    tick(0.6 + 0.25 * (i - 1) / n, paste("Testing", cr$id))
    cv <- .contrast_vector(dm, cr$reference, cr$target)
    tr <- if (options$edger_lfc > 0) edgeR::glmTreat(fit, contrast = cv, lfc = options$edger_lfc) else
      if (ql) edgeR::glmQLFTest(fit, contrast = cv) else edgeR::glmLRT(fit, contrast = cv)
    tt <- edgeR::topTags(tr, n = Inf, sort.by = "none")$table
    st <- if (!is.null(tt$F)) tt$F else if (!is.null(tt$LR)) tt$LR else NA_real_
    df <- data.frame(gene_id = rownames(tt), baseMean = unname(base[rownames(tt)]),
                     log2FoldChange = tt$logFC, stat = st, pvalue = tt$PValue, padj = tt$FDR,
                     logCPM = tt$logCPM, stringsAsFactors = FALSE)
    df <- df[order(df$pvalue, na.last = TRUE), , drop = FALSE]; rownames(df) <- NULL
    res[[cr$id]] <- df
  }
  notes <- c(notes, sprintf("edgeR %s test%s; common BCV %.3f.",
                            if (ql) "quasi-likelihood F" else "likelihood-ratio",
                            if (options$edger_lfc > 0) sprintf(" against |log2 FC| > %g (glmTreat)", options$edger_lfc) else "",
                            sqrt(y$common.dispersion)))
  list(engine = "edgeR", engine_label = "edgeR", results = res,
       transformed = edgeR::cpm(y, log = TRUE, prior.count = 2), transform_label = "log2-CPM (TMM)",
       normalised = nz$normalised,
       qc = list(size_factors = nz$size_factors,
                 bcv = data.frame(AveLogCPM = y$AveLogCPM, tagwise = sqrt(y$tagwise.dispersion),
                                  trended = sqrt(y$trended.dispersion), common = sqrt(y$common.dispersion),
                                  row.names = rownames(y))),
       notes = notes,
       versions = c(edgeR = as.character(utils::packageVersion("edgeR"))))
}

# ---- EBSeq ------------------------------------------------------------------------

#' EBSeq gives each gene a posterior probability of differential expression
#' (PPDE) instead of a p-value. For the Explorer, the probability of equal
#' expression (PPEE = 1 - PPDE) fills both p-value columns, so an adjusted-p
#' cut-off of 0.05 means PPDE >= 0.95 -- EBSeq's usual "FDR 5%" call.
.run_ebseq <- function(mat, cd, main, contrasts, options, tick, ng = NULL) {
  sf_all <- EBSeq::MedianNorm(mat)
  norm_all <- EBSeq::GetNormalizedMat(mat, sf_all)
  res <- list(); ppde <- list(); drop_notes <- character(0); n <- nrow(contrasts)
  for (i in seq_len(n)) {
    cr <- contrasts[i, ]
    tick(0.15 + 0.7 * (i - 1) / n, paste("EBSeq:", cr$id))
    keep <- as.character(cd[[main]]) %in% c(cr$reference, cr$target)
    sub <- mat[, keep, drop = FALSE]
    # Target first: EBSeq's fold change is the first condition over the second.
    cond <- factor(as.character(cd[[main]])[keep], levels = c(cr$target, cr$reference))
    sf <- EBSeq::MedianNorm(sub)
    eb <- NULL
    utils::capture.output(eb <- suppressMessages(EBSeq::EBTest(
      Data = sub, NgVector = if (is.null(ng)) NULL else as.numeric(ng[rownames(sub)]),
      Conditions = cond, sizeFactors = sf, maxround = options$ebseq_maxround)))
    pp <- EBSeq::GetPPMat(eb)
    fc <- EBSeq::PostFC(eb)
    g <- rownames(sub)
    # EBSeq drops genes it cannot test in this pair of groups (e.g. all-zero there),
    # so index by name and leave NA for the rest.
    mi <- match(g, rownames(pp))
    ppde_g <- as.numeric(pp[mi, "PPDE"]); ppee_g <- as.numeric(pp[mi, "PPEE"])
    lfc <- log2(as.numeric(fc$PostFC[match(g, names(fc$PostFC))]))
    real <- log2(as.numeric(fc$RealFC[match(g, names(fc$RealFC))]))
    lfc[!is.finite(lfc)] <- real[!is.finite(lfc)]
    if (anyNA(mi)) drop_notes <- c(drop_notes, sprintf(
      "%s gene(s) passed the overall filter but could not be tested by EBSeq for %s (for example, no counts in these two groups); their columns are NA there.",
      format(sum(is.na(mi)), big.mark = ","), cr$id))
    df <- data.frame(gene_id = g, baseMean = rowMeans(EBSeq::GetNormalizedMat(sub, sf)),
                     log2FoldChange = lfc, PPDE = ppde_g, PPEE = ppee_g,
                     pvalue = ppee_g, padj = ppee_g, log2RealFC = real, stringsAsFactors = FALSE)
    df <- df[order(df$PPEE, na.last = TRUE), , drop = FALSE]; rownames(df) <- NULL
    res[[cr$id]] <- df
    ppde[[cr$id]] <- ppde_g
  }
  list(engine = "EBSeq", engine_label = "EBSeq", results = res,
       transformed = log2(norm_all + 1), transform_label = "log2 median-normalised",
       normalised = norm_all,
       qc = list(size_factors = sf_all, ppde = ppde),
       notes = c(paste("EBSeq reports posterior probabilities, not p-values: the p-value and adjusted p",
                       "columns hold PPEE = 1 - PPDE, so adjusted p <= 0.05 means PPDE >= 0.95."),
                 "log2 fold changes are EBSeq's posterior fold changes (target / reference).",
                 drop_notes),
       versions = c(EBSeq = as.character(utils::packageVersion("EBSeq"))))
}
