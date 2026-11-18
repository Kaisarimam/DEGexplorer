# =============================================================================
# de_report.R -- run record, reproducible R script, design files and the export
# bundle for an Analysis-tab run (the `res` list from run_de_analysis()).
# =============================================================================

.dq <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")

# ---- result tables for export ------------------------------------------------------

#' One comparison's table as written to disk: gene symbols added for gene-level
#' runs (when an annotation is loaded); for transcript-level runs the id column is
#' named transcript_id (parent_gene_id, gene_name ... follow it).
#' Engine tag used in file names (deseq2, limma, edger, ebseq, edger_dtu, dexseq ...).
result_file_tag <- function(res) {
  if (is.null(res)) return("analysis")
  if (!is.null(res$engine_key)) return(res$engine_key)
  tolower(switch(res$engine, limma = "limma", edgeR = "edger", res$engine))
}

de_result_for_export <- function(res, id, annotation = NULL, df = NULL) {
  if (is.null(df)) df <- res$results[[id]]
  if (isTRUE(res$level %in% c("transcript", "dtu"))) {
    colnames(df)[1] <- "transcript_id"
    return(df)
  }
  if (!is.null(annotation) && all(c("gene_id", "gene_name") %in% colnames(annotation)) &&
      !"gene_name" %in% colnames(df)) {
    i <- annotation_match(df$gene_id, annotation$gene_id)
    df <- cbind(df[, 1, drop = FALSE], gene_name = annotation$gene_name[i], df[, -1, drop = FALSE])
  }
  df
}

# ---- design files (CSV) ---------------------------------------------------------

#' Sample sheet of the included samples: `sample` plus every variable.
de_samplesheet <- function(res) {
  sh <- res$sheet
  if (is.null(sh)) sh <- data.frame(sample = res$coldata$sample, res$coldata[, -1, drop = FALSE],
                                    check.names = FALSE, stringsAsFactors = FALSE)
  sh[] <- lapply(sh, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x })
  rownames(sh) <- NULL
  sh
}

#' Comparisons table: id, variable, reference, target, blocking.
de_contrasts_table <- function(contrasts, covariates = character(0)) {
  data.frame(id = contrasts$id, variable = contrasts$variable,
             reference = contrasts$reference, target = contrasts$target,
             blocking = rep(paste(covariates, collapse = ";"), nrow(contrasts)),
             stringsAsFactors = FALSE)
}

#' Read a sample sheet (CSV/TSV) for the editor: returns the sheet and variable types.
read_design_samplesheet <- function(path, samples) {
  sep <- if (grepl("\\.(tsv|txt)(\\.gz)?$", path, ignore.case = TRUE)) "\t" else ","
  df <- utils::read.table(path, sep = sep, header = TRUE, colClasses = "character",
                          check.names = FALSE, quote = "\"", comment.char = "", na.strings = character(0))
  colnames(df) <- trimws(colnames(df))
  id_col <- intersect(c("sample", "Sample", "sample_id", "id", "ID"), colnames(df))[1]
  if (is.na(id_col)) id_col <- colnames(df)[1]
  hit <- match(samples, df[[id_col]])
  if (!any(!is.na(hit))) stop("None of the sample names in the file match the count matrix columns (looked in '", id_col, "').")
  vars <- setdiff(colnames(df), c(id_col, "fastq_1", "fastq_2", "strandedness", "sample"))
  out <- data.frame(sample = samples, stringsAsFactors = FALSE)
  types <- list()
  for (v in vars) {
    x <- trimws(df[[v]][hit]); x[is.na(x)] <- ""
    nm <- make.names(v)
    out[[nm]] <- x
    num <- suppressWarnings(as.numeric(x[nzchar(x)]))
    types[[nm]] <- if (length(num) && !anyNA(num) && length(unique(num)) > 3) "numeric" else "factor"
  }
  list(sheet = out, types = types, matched = sum(!is.na(hit)), unmatched = samples[is.na(hit)],
       id_column = id_col)
}

#' Read a contrasts file: id (optional), variable, reference, target, blocking (optional).
read_design_contrasts <- function(path) {
  sep <- if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) "\t" else ","
  df <- utils::read.table(path, sep = sep, header = TRUE, colClasses = "character",
                          check.names = FALSE, quote = "\"", comment.char = "")
  colnames(df) <- tolower(trimws(colnames(df)))
  need <- c("variable", "reference", "target")
  if (!all(need %in% colnames(df))) stop("A contrasts file needs the columns variable, reference and target.")
  df[] <- lapply(df, trimws)
  id <- if ("id" %in% colnames(df)) df$id else rep("", nrow(df))
  id[!nzchar(id)] <- contrast_id(df$variable, df$target, df$reference)[!nzchar(id)]
  blocking <- if ("blocking" %in% colnames(df)) df$blocking else rep("", nrow(df))
  blk <- unique(unlist(strsplit(blocking[nzchar(blocking)], "[;,]")))
  list(contrasts = data.frame(id = make.names(id), variable = make.names(df$variable),
                              reference = df$reference, target = df$target, stringsAsFactors = FALSE),
       blocking = make.names(trimws(blk[nzchar(trimws(blk))])))
}

# ---- run record --------------------------------------------------------------

#' The engine settings that applied to this run, as label = value.
de_engine_settings <- function(res) {
  o <- res$settings; e <- res$engine
  yn <- function(x) if (isTRUE(x)) "yes" else "no"
  s <- switch(e,
    DESeq2 = c(
      "Test" = if (identical(o$deseq2_test, "lrt")) "likelihood-ratio (reduced model drops the main variable)" else "Wald",
      "Dispersion fit" = o$fit_type, "alpha (independent filtering)" = o$alpha,
      "|log2 FC| test threshold" = o$lfc_threshold, "Independent filtering" = yn(o$independent_filtering),
      "Cook's cutoff" = yn(o$cooks_cutoff), "Fold-change shrinkage" = o$shrink %||% "none",
      "Transform" = paste0(res$transform_label, if (isTRUE(o$blind)) " (blind)" else " (design-aware)")),
    limma = c(
      "Normalisation" = o$norm_method,
      "voom" = if (identical(o$limma_voom, "quality")) "voomWithQualityWeights" else "voom",
      "eBayes trend" = yn(o$limma_trend), "eBayes robust" = yn(o$limma_robust),
      "Random block" = o$limma_block %||% "(none)", "|log2 FC| threshold (treat)" = o$limma_lfc),
    edgeR = c(
      "Normalisation" = o$norm_method,
      "Test" = if (identical(o$edger_test, "lrt")) "likelihood-ratio (glmLRT)" else "quasi-likelihood F (glmQLFTest)",
      "Robust" = yn(o$edger_robust), "|log2 FC| threshold (glmTreat)" = o$edger_lfc),
    EBSeq = c("EM iterations (maxround)" = o$ebseq_maxround),
    character(0))
  if (identical(res$level, "dtu")) {
    d <- res$dtu_opts; k <- res$engine_key
    s <- c("Analysis" = "differential transcript usage (DTU)",
           "Method" = res$engine_label,
           "Counts" = DTU_COUNT_LABELS[[res$dtu_counts]],
           "Filter" = res$filter,
           if (k %in% c("edger_dtu", "limma_dtu")) c(
             "Gene-level test" = if (identical(d$gene_test, "F")) "diffSplice F-test" else "Simes",
             "Robust" = yn(d$robust)),
           if (k == "drimseq") c("Post-hoc proportion-SD filter" = if (d$sd_filter > 0) format(d$sd_filter) else "off",
                                 "Random seed (dmPrecision)" = d$seed),
           if (k == "saturn") c("satuRn p-values" = d$saturn_pvalues),
           if (k %in% c("dexseq", "drimseq")) c("Comparisons" = "each fitted on its own two groups"),
           "Stage-wise testing (stageR, dtu)" = if (isTRUE(d$stager)) sprintf("yes, alpha %g", d$stager_alpha) else "no")
    return(vapply(s, as.character, ""))
  }
  tx <- identical(res$level, "transcript")
  if (e %in% c("DESeq2", "edgeR")) s <- c(s, setNames(yn(res$length_offset),
                                             if (tx) "Effective-length offsets" else "Gene-length offsets"))
  if (tx) {
    s <- c("Feature level" = "transcripts (isoforms)", s)
    if (e %in% c("limma", "edgeR"))
      s <- c(s, "Counts divided by quantification overdispersion" = yn(res$divided))
    if (e == "EBSeq")
      s <- c(s, "Ng vector" = if (is.null(res$ngvector)) "none" else
               if (identical(o$ebseq_ng, "file")) "from file" else sprintf("isoforms per gene (GetNg, TrunThre %s)", o$ebseq_trunthre))
    if (e != "EBSeq") s <- c(s, "Stage-wise testing (stageR)" =
      if (any(vapply(res$results, function(d) "stageR_tx_padj" %in% colnames(d), logical(1))))
        sprintf("yes, alpha %g", o$stager_alpha) else "no")
  }
  vapply(s, as.character, "")
}

#' A Markdown record of everything that defines the run.
de_run_record <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1, source = NULL) {
  sh <- de_samplesheet(res)
  md_table <- function(df) {
    df[] <- lapply(df, function(x) gsub("|", "\\|", as.character(x), fixed = TRUE))
    c(paste0("| ", paste(colnames(df), collapse = " | "), " |"),
      paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
      apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  }
  cmp <- do.call(rbind, lapply(names(res$results), function(id) {
    n <- de_result_counts(res$results[[id]], padj_cutoff, lfc_cutoff)
    cr <- res$contrasts[res$contrasts$id == id, ]
    out <- data.frame(Comparison = id, Target = cr$target, Reference = cr$reference,
                      Tested = n[["tested"]], Up = n[["up"]], Down = n[["down"]])
    if (identical(res$level, "dtu")) out$`Genes with DTU` <- dtu_gene_count(res$results[[id]], padj_cutoff)
    out
  }))
  set <- de_engine_settings(res)
  vers <- c(res$versions, R = paste(R.version$major, R.version$minor, sep = "."))
  c("# DEGexplorer analysis record", "",
    sprintf("- **Run:** %s", res$run_at),
    sprintf("- **Engine:** %s", res$engine_label %||% res$engine),
    if (!is.null(source)) sprintf("- **Counts:** %s", source),
    sprintf("- **Model:** `%s` (main variable: %s%s)", res$formula, res$main,
            if (length(res$covariates)) paste0("; covariates: ", paste(res$covariates, collapse = ", ")) else ""),
    sprintf("- **Filter:** %s; %s of %s %s tested", res$filter,
            format(length(res$kept), big.mark = ","), format(res$n_input, big.mark = ","),
            if (isTRUE(res$level %in% c("transcript", "dtu"))) "transcripts" else "genes"),
    sprintf("- **Matrix for PCA and heatmaps:** %s", res$transform_label),
    "", "## Samples", "", md_table(sh),
    "", "## Comparisons", "",
    sprintf("Counts at adjusted p \u2264 %g and |log2 FC| \u2265 %g.", padj_cutoff, lfc_cutoff), "",
    md_table(cmp),
    "", "## Settings", "",
    md_table(data.frame(Setting = names(set), Value = unname(set))),
    if (length(res$notes)) c("", "## Notes", "", paste("-", res$notes)),
    "", "## Software", "",
    md_table(data.frame(Package = names(vers), Version = unname(vers))),
    "", "Reproduce with `Rscript reproduce.R` inside the export bundle.")
}

# ---- reproducible R script -------------------------------------------------------

#' A standalone R script that repeats the run from the bundle's files.
de_reproduce_script <- function(res) {
  if (identical(res$level, "dtu")) return(dtu_reproduce_script(res))
  o <- res$settings; e <- res$engine; main <- res$main; cov <- res$covariates
  cd <- res$coldata
  eng_file <- tolower(switch(e, limma = "limma", edgeR = "edger", e))
  pkgs <- switch(e, DESeq2 = "DESeq2", limma = c("limma", "edgeR"), edgeR = "edgeR", EBSeq = "EBSeq")
  sh <- o$shrink %||% "none"
  if (e == "DESeq2" && sh %in% c("apeglm", "ashr")) pkgs <- c(pkgs, sh)
  tx <- identical(res$level, "transcript")
  use_stager <- tx && any(vapply(res$results, function(d) "stageR_tx_padj" %in% colnames(d), logical(1)))
  if (use_stager) pkgs <- c(pkgs, "stageR")
  ng_mode <- if (tx && e == "EBSeq" && !is.null(res$ngvector)) (if (identical(o$ebseq_ng, "file")) "file" else "isoforms") else "none"
  vers <- paste(names(res$versions), res$versions, collapse = ", ")
  L <- character(0); add <- function(...) L <<- c(L, ...)

  add("# =============================================================================",
      sprintf("# Reproduces the DEGexplorer run of %s", res$run_at),
      sprintf("# Engine: %s (%s, R %s.%s)", res$engine_label %||% e, vers, R.version$major, R.version$minor),
      sprintf("# Model:  %s", res$formula),
      "#",
      "# Reads the files exported with it, from the same folder:",
      if (tx) "#   counts.tsv       raw transcript counts of the included samples, before filtering"
      else "#   counts.tsv       raw counts of the included samples, before filtering",
      "#   samplesheet.csv  sample sheet",
      "#   contrasts.csv    comparisons (id, variable, reference, target, blocking)",
      if (isTRUE(res$length_offset)) "#   lengths.tsv      sample-specific (effective) lengths for the offsets",
      if (isTRUE(res$divided)) "#   overdispersion.tsv  each transcript's quantification overdispersion (edgeR catchSalmon / catchKallisto)",
      if (tx) "#   tx2gene.tsv      transcript -> gene map",
      if (ng_mode == "file") "#   ngvector.tsv     EBSeq Ng vector",
      "# and writes reproduced/<comparison>.<engine>.results.tsv.",
      "#",
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
    if (is.numeric(cd[[v]])) {
      add(sprintf("coldata[[%s]] <- as.numeric(sheet[[%s]])", .dq(v), .dq(v)))
    } else {
      add(sprintf("coldata[[%s]] <- factor(trimws(sheet[[%s]]), levels = %s)", .dq(v), .dq(v),
                  .dq(levels(cd[[v]]))))
    }
  }
  add("counts <- counts[, rownames(coldata), drop = FALSE]",
      "counts[is.na(counts)] <- 0")
  if (tx) add("tx2gene <- read.delim(\"tx2gene.tsv\", colClasses = \"character\", check.names = FALSE)")
  if (isTRUE(res$divided)) {
    add("# Divide out read-to-transcript ambiguity (Baldoni et al. 2024); divided counts stay fractional",
        "od <- read.delim(\"overdispersion.tsv\", colClasses = c(\"character\", \"numeric\"))",
        "o <- od$overdispersion[match(rownames(counts), od$transcript_id)]",
        "o[is.na(o)] <- 1",
        "counts <- counts / pmax(o, 1)")
  } else {
    add("counts <- round(counts)   # RSEM / Salmon expected counts are rounded for the count models")
  }
  if (ng_mode == "isoforms") {
    g <- ifelse(is.na(res$tx2gene$gene_id), res$tx2gene$transcript_id, res$tx2gene$gene_id)
    tt <- min(as.integer(o$ebseq_trunthre %||% 3), max(table(g)))
    add("# EBSeq isoform groups: number of isoforms of each transcript's gene, truncated (GetNg)",
        "g <- ifelse(is.na(tx2gene$gene_id) | tx2gene$gene_id == \"\", tx2gene$transcript_id, tx2gene$gene_id)",
        sprintf("ng_all <- GetNg(tx2gene$transcript_id, g, TrunThre = %d)$IsoformNgTrun", tt))
  }
  add("",
      sprintf("# Filter: %s", res$filter))
  if (isTRUE(o$use_cpm)) {
    add(sprintf("abundance <- t(t(counts) / pmax(colSums(counts), 1) * 1e6)"))
  } else add("abundance <- counts")
  add(sprintf("keep <- rowSums(abundance >= %s) >= %d", format(o$min_abundance), as.integer(res$filter_n)),
      "counts <- counts[keep, , drop = FALSE]",
      sprintf("cat(sum(keep), \"%s pass the filter\\n\")", if (tx) "transcripts" else "genes"))
  if (ng_mode == "isoforms") add("ng <- as.integer(ng_all[rownames(counts)]); ng[is.na(ng)] <- 1L")
  if (ng_mode == "file") add("ngv <- read.delim(\"ngvector.tsv\", colClasses = c(\"character\", \"integer\"))",
                             "ng <- ngv$Ng[match(rownames(counts), ngv$transcript_id)]; ng[is.na(ng)] <- 1L")
  if (isTRUE(res$length_offset)) {
    add("lengths <- as.matrix(read.delim(\"lengths.tsv\", row.names = 1, check.names = FALSE))",
        "lengths <- lengths[rownames(counts), colnames(counts)]")
  }
  add("", "dir.create(\"reproduced\", showWarnings = FALSE)")
  if (use_stager) {
    add("# Stage-wise testing (stageR, method \"dte\"): screen genes (Simes of their transcripts,",
        "# BH across genes), then confirm transcripts within significant genes (Holm).",
        "stage_wise_dte <- function(df, alpha) {",
        "  ok <- !is.na(df$pvalue) & !is.na(df$parent_gene_id)",
        "  p <- setNames(df$pvalue[ok], df$gene_id[ok]); g <- df$parent_gene_id[ok]",
        "  sc <- tapply(p, g, function(x) { x <- sort(x); min(1, min(length(x) * x / seq_along(x))) })",
        "  screen <- setNames(as.numeric(sc), names(sc))   # stageR needs a plain numeric vector",
        "  conf <- matrix(p, ncol = 1, dimnames = list(names(p), \"transcript\"))",
        "  obj <- stageRTx(pScreen = screen, pConfirmation = conf, pScreenAdjusted = FALSE,",
        "                  tx2gene = data.frame(transcript = names(p), gene = g))",
        "  obj <- stageWiseAdjustment(obj, method = \"dte\", alpha = alpha)",
        "  adj <- as.data.frame(getAdjustedPValues(obj, order = FALSE, onlySignificantGenes = FALSE))",
        "  i <- match(df$gene_id, adj$txID)",
        "  df$stageR_gene_padj <- as.numeric(adj$gene[i]); df$stageR_tx_padj <- as.numeric(adj$transcript[i])",
        "  df",
        "}")
  }
  add("save_result <- function(df, id) {")
  if (tx) add(
      "  a <- tx2gene[match(df$gene_id, tx2gene$transcript_id), , drop = FALSE]",
      "  df <- data.frame(gene_id = df$gene_id, parent_gene_id = a$gene_id, gene_name = a$gene_name,",
      "                   transcript_name = a$transcript_name, transcript_biotype = a$transcript_biotype,",
      "                   df[, -1, drop = FALSE], check.names = FALSE)",
      if (use_stager) sprintf("  df <- stage_wise_dte(df, %s)", format(o$stager_alpha %||% 0.05)))
  add("  df <- df[order(df$pvalue, na.last = TRUE), , drop = FALSE]",
      if (tx) "  colnames(df)[1] <- \"transcript_id\"",
      sprintf("  write.table(df, file.path(\"reproduced\", paste0(id, \".%s.results.tsv\")),", eng_file),
      "              sep = \"\\t\", quote = FALSE, row.names = FALSE, na = \"NA\")",
      "}", "")

  if (e %in% c("limma", "edgeR")) {
    block <- o$limma_block
    if (e != "limma" || is.null(block) || block == "(none)") block <- NULL
    fixed <- setdiff(cov, block)
    add("# Cell-means design: one coefficient per level of the main variable",
        sprintf("design <- model.matrix(%s, coldata)",
                paste("~ 0 +", paste(c(main, fixed), collapse = " + "))),
        "contrast_vector <- function(target, reference) {",
        "  v <- setNames(numeric(ncol(design)), colnames(design))",
        "  v[paste0(main, target)] <- 1; v[paste0(main, reference)] <- -1",
        "  v",
        "}",
        "# edgeR >= 4 renamed calcNormFactors() to normLibSizes()",
        "norm_lib <- if (exists(\"normLibSizes\", envir = asNamespace(\"edgeR\"))) edgeR::normLibSizes else edgeR::calcNormFactors",
        "y <- DGEList(counts = counts)")
    if (isTRUE(res$length_offset) && e == "edgeR") {
      add("# Gene-length offsets (tximport's recipe)",
          "norm_mat <- lengths / exp(rowMeans(log(lengths)))",
          "norm_cts <- counts / norm_mat",
          sprintf("eff_lib <- norm_lib(DGEList(norm_cts), method = %s)", .dq(o$norm_method)),
          "eff_lib <- eff_lib$samples$norm.factors * colSums(norm_cts)",
          "y <- scaleOffset(y, log(sweep(norm_mat, 2, eff_lib, \"*\")))",
          "size_factors <- eff_lib / exp(mean(log(eff_lib)))",
          "base_mean <- rowMeans(norm_cts / rep(size_factors, each = nrow(norm_cts)))")
    } else {
      add(sprintf("y <- norm_lib(y, method = %s)", .dq(o$norm_method)),
          "eff <- y$samples$lib.size * y$samples$norm.factors",
          "size_factors <- eff / exp(mean(log(eff)))",
          "base_mean <- rowMeans(t(t(counts) / size_factors))")
    }
  }

  if (e == "DESeq2") {
    lrt <- identical(o$deseq2_test, "lrt")
    for (v in cov) if (is.numeric(cd[[v]]))
      add(sprintf("coldata[[%s]] <- as.numeric(scale(coldata[[%s]]))   # centred and scaled, as DESeq2 recommends", .dq(v), .dq(v)))
    add(sprintf("dds <- DESeqDataSetFromMatrix(counts, coldata, design = %s)", res$formula))
    if (isTRUE(res$length_offset))
      add("assay(dds, \"avgTxLength\", withDimnames = FALSE) <- lengths   # gene-length offsets, as with tximport")
    reduced <- if (length(cov)) paste("~", paste(cov, collapse = " + ")) else "~ 1"
    add(if (lrt) sprintf("dds <- DESeq(dds, test = \"LRT\", reduced = %s, fitType = %s)", reduced, .dq(o$fit_type))
        else sprintf("dds <- DESeq(dds, test = \"Wald\", fitType = %s)", .dq(o$fit_type)))
    if (sh == "apeglm") {
      add("", "# apeglm needs a model coefficient; refit with the reference as base level when needed",
          "shrink_apeglm <- function(r, target, reference) {",
          "  want <- make.names(paste0(main, \"_\", target, \"_vs_\", reference))",
          "  d <- dds",
          "  if (!want %in% make.names(resultsNames(d))) {",
          "    d[[main]] <- relevel(d[[main]], ref = reference)",
          "    d <- nbinomWaldTest(d)",
          "  }",
          "  coef <- resultsNames(d)[make.names(resultsNames(d)) == want][1]",
          "  lfcShrink(d, coef = coef, res = r, type = \"apeglm\")",
          "}")
    }
    add("", "for (i in seq_len(nrow(contrasts))) {",
        "  cr <- contrasts[i, ]",
        sprintf("  r <- results(dds, contrast = c(main, cr$target, cr$reference), alpha = %s,", format(o$alpha)),
        sprintf("               lfcThreshold = %s, independentFiltering = %s, cooksCutoff = %s)",
                format(if (lrt) 0 else o$lfc_threshold), isTRUE(o$independent_filtering), isTRUE(o$cooks_cutoff)),
        "  out <- data.frame(gene_id = rownames(r), as.data.frame(r), check.names = FALSE)")
    if (sh != "none") {
      add(if (sh == "apeglm") "  sr <- shrink_apeglm(r, cr$target, cr$reference)"
          else sprintf("  sr <- lfcShrink(dds, contrast = c(main, cr$target, cr$reference), res = r, type = %s)", .dq(sh)),
          "  out$log2FoldChange_mle <- out$log2FoldChange",
          "  out$lfcSE_mle <- out$lfcSE",
          "  out$log2FoldChange <- as.numeric(sr$log2FoldChange)",
          "  out$lfcSE <- if (\"lfcSE\" %in% colnames(sr)) as.numeric(sr$lfcSE) else NA_real_")
    }
    add("  save_result(out, cr$id)", "}")
    tr <- if (identical(o$transform, "rlog")) sprintf("rlog(dds, blind = %s)", isTRUE(o$blind)) else
      sprintf("vst(dds, blind = %s)", isTRUE(o$blind))
    add("", "write.table(data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE),",
        "            file.path(\"reproduced\", \"normalised_counts.tsv\"), sep = \"\\t\", quote = FALSE, row.names = FALSE)",
        sprintf("write.table(data.frame(gene_id = rownames(dds), assay(%s), check.names = FALSE),", tr),
        sprintf("            file.path(\"reproduced\", \"%s.tsv\"), sep = \"\\t\", quote = FALSE, row.names = FALSE)",
                tolower(res$transform_label)))
  }

  if (e == "limma") {
    quality <- identical(o$limma_voom, "quality")
    vf <- if (quality) "voomWithQualityWeights" else "voom"
    if (is.null(block)) {
      add(sprintf("v <- %s(y, design)", vf), "fit <- lmFit(v, design)")
    } else {
      add(sprintf("# '%s' is a random block: estimate the within-block correlation twice, as limma's guide does", block),
          sprintf("block <- as.character(coldata[[%s]])", .dq(block)),
          sprintf("v <- %s(y, design)", vf),
          "dc <- duplicateCorrelation(v, design, block = block)",
          sprintf("v <- %s(y, design, block = block, correlation = dc$consensus.correlation)", vf),
          "dc <- duplicateCorrelation(v, design, block = block)",
          "fit <- lmFit(v, design, block = block, correlation = dc$consensus.correlation)",
          "cat(\"consensus correlation:\", round(dc$consensus.correlation, 3), \"\\n\")")
    }
    add("", "for (i in seq_len(nrow(contrasts))) {",
        "  cr <- contrasts[i, ]",
        "  f2 <- contrasts.fit(fit, contrast_vector(cr$target, cr$reference))")
    if (o$limma_lfc > 0) {
      add(sprintf("  f2 <- treat(f2, lfc = %s, trend = %s, robust = %s)", format(o$limma_lfc), isTRUE(o$limma_trend), isTRUE(o$limma_robust)),
          "  tt <- topTreat(f2, coef = 1, number = Inf, sort.by = \"none\")")
    } else {
      add(sprintf("  f2 <- eBayes(f2, trend = %s, robust = %s)", isTRUE(o$limma_trend), isTRUE(o$limma_robust)),
          "  tt <- topTable(f2, coef = 1, number = Inf, sort.by = \"none\")")
    }
    add("  out <- data.frame(gene_id = rownames(tt), baseMean = unname(base_mean[rownames(tt)]),",
        "                    log2FoldChange = tt$logFC,",
        "                    lfcSE = as.numeric(f2$stdev.unscaled[, 1] * sqrt(f2$s2.post)),",
        "                    stat = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val, AveExpr = tt$AveExpr)",
        "  if (!is.null(tt$B)) out$B <- tt$B",
        "  save_result(out, cr$id)", "}",
        "write.table(data.frame(gene_id = rownames(v$E), v$E, check.names = FALSE),",
        "            file.path(\"reproduced\", \"voom_log2cpm.tsv\"), sep = \"\\t\", quote = FALSE, row.names = FALSE)")
  }

  if (e == "edgeR") {
    ql <- !identical(o$edger_test, "lrt")
    add(sprintf("y <- estimateDisp(y, design, robust = %s)", isTRUE(o$edger_robust)),
        if (ql) sprintf("fit <- glmQLFit(y, design, robust = %s)", isTRUE(o$edger_robust)) else "fit <- glmFit(y, design)",
        "cat(\"common BCV:\", round(sqrt(y$common.dispersion), 3), \"\\n\")",
        "", "for (i in seq_len(nrow(contrasts))) {",
        "  cr <- contrasts[i, ]",
        "  cv <- contrast_vector(cr$target, cr$reference)",
        if (o$edger_lfc > 0) sprintf("  tr <- glmTreat(fit, contrast = cv, lfc = %s)", format(o$edger_lfc))
        else if (ql) "  tr <- glmQLFTest(fit, contrast = cv)" else "  tr <- glmLRT(fit, contrast = cv)",
        "  tt <- topTags(tr, n = Inf, sort.by = \"none\")$table",
        "  stat <- if (!is.null(tt$F)) tt$F else if (!is.null(tt$LR)) tt$LR else NA_real_",
        "  out <- data.frame(gene_id = rownames(tt), baseMean = unname(base_mean[rownames(tt)]),",
        "                    log2FoldChange = tt$logFC, stat = stat, pvalue = tt$PValue, padj = tt$FDR,",
        "                    logCPM = tt$logCPM)",
        "  save_result(out, cr$id)", "}")
  }

  if (e == "EBSeq") {
    add("for (i in seq_len(nrow(contrasts))) {",
        "  cr <- contrasts[i, ]",
        "  keep_s <- as.character(coldata[[main]]) %in% c(cr$reference, cr$target)",
        "  sub <- counts[, keep_s, drop = FALSE]",
        "  # target first: EBSeq's fold change is the first condition over the second",
        "  cond <- factor(as.character(coldata[[main]])[keep_s], levels = c(cr$target, cr$reference))",
        "  sf <- MedianNorm(sub)",
        if (ng_mode != "none") sprintf("  eb <- EBTest(Data = sub, NgVector = as.numeric(ng[match(rownames(sub), rownames(counts))]), Conditions = cond, sizeFactors = sf, maxround = %d)", as.integer(o$ebseq_maxround))
        else sprintf("  eb <- EBTest(Data = sub, Conditions = cond, sizeFactors = sf, maxround = %d)", as.integer(o$ebseq_maxround)),
        "  pp <- GetPPMat(eb); fc <- PostFC(eb); g <- rownames(sub)",
        "  mi <- match(g, rownames(pp))   # EBSeq drops genes it cannot test in these two groups",
        "  lfc <- log2(as.numeric(fc$PostFC[match(g, names(fc$PostFC))]))",
        "  real <- log2(as.numeric(fc$RealFC[match(g, names(fc$RealFC))]))",
        "  lfc[!is.finite(lfc)] <- real[!is.finite(lfc)]",
        "  ppee <- as.numeric(pp[mi, \"PPEE\"])",
        "  out <- data.frame(gene_id = g, baseMean = rowMeans(GetNormalizedMat(sub, sf)),",
        "                    log2FoldChange = lfc, PPDE = as.numeric(pp[mi, \"PPDE\"]), PPEE = ppee,",
        "                    pvalue = ppee, padj = ppee, log2RealFC = real)   # PPEE = 1 - PPDE",
        "  save_result(out, cr$id)", "}")
  }
  add("", "writeLines(capture.output(sessionInfo()), file.path(\"reproduced\", \"session_info.txt\"))",
      "cat(\"Results written to\", normalizePath(\"reproduced\"), \"\\n\")")
  L
}

# ---- export bundle ---------------------------------------------------------------

.write_tsv <- function(df, file) utils::write.table(df, file, sep = "\t", quote = FALSE,
                                                     row.names = FALSE, na = "NA")
#' Like .write_tsv but numbers keep every digit (17 significant), so the inputs the
#' reproduce script reads back are bit-identical to what the app used.
.write_tsv_exact <- function(df, file) {
  for (cl in names(df)) if (is.double(df[[cl]]))
    df[[cl]] <- ifelse(is.na(df[[cl]]), NA_character_, sprintf("%.17g", df[[cl]]))
  .write_tsv(df, file)
}
.matrix_df <- function(m) data.frame(gene_id = rownames(m), m, check.names = FALSE,
                                     stringsAsFactors = FALSE)

#' Zip everything about a run: results, matrices, design, QC figures, record and script.
#'
#' @param annotation optional data frame (gene_id, gene_name) to add symbols to results
#' @param qc logical: render the QC figures (PNG) into qc/
write_de_bundle <- function(res, file, padj_cutoff = 0.05, lfc_cutoff = 1,
                            annotation = NULL, source = NULL, qc = TRUE) {
  dir <- tempfile("de_bundle_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  eng_file <- result_file_tag(res)
  dir.create(file.path(dir, "results"))
  for (id in names(res$results)) {
    .write_tsv(de_result_for_export(res, id, annotation),
               file.path(dir, "results", paste0(id, ".", eng_file, ".results.tsv")))
  }
  .write_tsv_exact(.matrix_df(res$input_counts), file.path(dir, "counts.tsv"))
  utils::write.csv(de_samplesheet(res), file.path(dir, "samplesheet.csv"), row.names = FALSE)
  utils::write.csv(de_contrasts_table(res$contrasts, res$covariates), file.path(dir, "contrasts.csv"),
                   row.names = FALSE)
  if (isTRUE(res$length_offset) && !is.null(res$lengths_used))
    .write_tsv_exact(.matrix_df(res$lengths_used), file.path(dir, "lengths.tsv"))
  if (isTRUE(res$divided) && !is.null(res$overdispersion))
    .write_tsv_exact(data.frame(transcript_id = names(res$overdispersion), overdispersion = unname(res$overdispersion)),
               file.path(dir, "overdispersion.tsv"))
  if (isTRUE(res$level %in% c("transcript", "dtu")) && !is.null(res$tx2gene))
    .write_tsv(res$tx2gene, file.path(dir, "tx2gene.tsv"))
  if (!is.null(res$dtu_lengths)) .write_tsv_exact(.matrix_df(res$dtu_lengths), file.path(dir, "lengths.tsv"))
  if (!is.null(res$ngvector) && identical(res$settings$ebseq_ng, "file"))
    .write_tsv(data.frame(transcript_id = names(res$ngvector), Ng = unname(res$ngvector)),
               file.path(dir, "ngvector.tsv"))
  dir.create(file.path(dir, "matrices"))
  .write_tsv(.matrix_df(round(res$normalised, 4)), file.path(dir, "matrices", "normalised_counts.tsv"))
  .write_tsv(.matrix_df(round(res$transformed, 5)),
             file.path(dir, "matrices", paste0(gsub("[^A-Za-z0-9]+", "_", tolower(res$transform_label)), ".tsv")))
  writeLines(de_run_record(res, padj_cutoff, lfc_cutoff, source), file.path(dir, "run_record.md"))
  writeLines(de_reproduce_script(res), file.path(dir, "reproduce.R"))
  if (qc) {
    dir.create(file.path(dir, "qc"))
    for (w in QC_PLOTS) {
      o <- tryCatch(build_qc_figure(w, res), error = function(e) NULL)
      if (is.null(o)) next
      f <- file.path(dir, "qc", paste0(w, ".png"))
      tryCatch(export_plot(o, f, "png", 7, 5.5, 150), error = function(e) NULL)
    }
    .write_tsv(qc_sample_summary(res), file.path(dir, "qc", "per_sample.tsv"))
  }
  writeLines(capture.output(utils::sessionInfo()), file.path(dir, "session_info.txt"))
  zip_directory(dir, file)
}
