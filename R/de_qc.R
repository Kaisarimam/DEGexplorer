# =============================================================================
# de_qc.R -- quality-control figures for a run from the Analysis tab
#
# Each function takes the list returned by run_de_analysis() (plus the raw
# counts where needed) and returns a ggplot, or a pheatmap object for the
# sample-distance heatmap.
# =============================================================================

QC_PLOTS <- c(
  "Library sizes"              = "library",
  "Count distributions"        = "distribution",
  "PCA"                        = "pca",
  "Sample distances"           = "distance",
  "Dispersion / mean-variance" = "dispersion",
  "P-value histograms"         = "pvalues"
)

#' Sample variables of a run as character columns (for colouring / annotation).
qc_sample_table <- function(res) {
  cd <- res$coldata
  out <- data.frame(sample = cd$sample, stringsAsFactors = FALSE)
  for (v in setdiff(colnames(cd), "sample")) out[[v]] <- as.character(cd[[v]])
  out
}

#' Per-sample numbers: library size, detected genes, size factor.
qc_sample_summary <- function(res) {
  q <- res$qc
  df <- data.frame(Sample = names(q$lib_size),
                   `Library size (M)` = round(q$lib_size / 1e6, 2),
                   `Genes detected` = as.integer(q$detected),
                   check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(q$size_factors)) df$`Size factor` <- round(unname(q$size_factors[df$Sample]), 3)
  if (!is.null(q$sample_weights)) df$`voom quality weight` <- round(as.numeric(q$sample_weights[match(df$Sample, names(q$sample_weights))]), 3)
  vt <- qc_sample_table(res)
  df[[res$main]] <- vt[[res$main]][match(df$Sample, vt$sample)]
  df
}

.qc_fill_scale <- function(values, name) {
  lv <- unique(stats::na.omit(values))
  ggplot2::scale_fill_manual(values = stats::setNames(group_palette(length(lv)), lv), name = name)
}

#' Library size per sample, coloured by the main variable.
qc_library_plot <- function(res, style = default_style()) {
  st <- resolve_style(style)
  vt <- qc_sample_table(res); q <- res$qc
  df <- data.frame(sample = names(q$lib_size), size = q$lib_size / 1e6,
                   detected = q$detected, stringsAsFactors = FALSE)
  df$group <- vt[[res$main]][match(df$sample, vt$sample)]
  df$sample <- factor(df$sample, levels = df$sample[order(df$group, df$sample)])
  ggplot2::ggplot(df, ggplot2::aes(x = sample, y = size, fill = group)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1fM\n%s genes", size,
                                                    format(detected, big.mark = ","))),
                       vjust = -0.25, size = st$label_size * 0.85, lineheight = 0.9) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    .qc_fill_scale(df$group, res$main) +
    deg_theme(st) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "Library sizes", x = NULL, y = "Total counts (millions)",
                  subtitle = "labels: total counts, and genes with at least one count")
}

#' Distribution of log2(count + 1) per sample, raw and normalised, for the
#' genes that passed filtering.
qc_distribution_plot <- function(res, counts = NULL, style = default_style()) {
  st <- resolve_style(style)
  vt <- qc_sample_table(res)
  smp <- vt$sample
  if (is.null(counts)) counts <- res$counts
  raw <- log2(counts[res$kept, smp, drop = FALSE] + 1)
  nrm <- log2(res$normalised[res$kept, smp, drop = FALSE] + 1)
  long <- function(m, what) data.frame(sample = rep(colnames(m), each = nrow(m)),
                                       value = as.vector(m), what = what, stringsAsFactors = FALSE)
  df <- rbind(long(raw, "Raw counts"), long(nrm, paste(res$engine_label %||% res$engine, "normalised")))
  df$what <- factor(df$what, levels = unique(df$what))
  df$group <- vt[[res$main]][match(df$sample, vt$sample)]
  df$sample <- factor(df$sample, levels = smp[order(vt[[res$main]], smp)])
  ggplot2::ggplot(df, ggplot2::aes(x = sample, y = value, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = st$line_width) +
    ggplot2::facet_wrap(~what, ncol = 1) +
    .qc_fill_scale(df$group, res$main) +
    deg_theme(st) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "Count distributions", x = NULL, y = "log2(count + 1)",
                  subtitle = sprintf("%s genes that passed filtering; medians should line up after normalisation",
                                     format(length(res$kept), big.mark = ",")))
}

#' PCA of the transformed matrix (VST / rlog / logCPM).
qc_pca_plot <- function(res, colour_by = res$main, shape_by = NULL, ntop = 500,
                        style = default_style()) {
  pr <- run_pca(res$transformed, ntop = ntop)
  make_pca_plot(pr, qc_sample_table(res), colour_by = colour_by,
                shape_by = if (identical(shape_by, colour_by)) NULL else shape_by,
                style = style,
                title = "PCA of the analysed samples",
                subtitle = sprintf("%s matrix, top %s most variable genes", res$transform_label,
                                   format(pr$ngenes, big.mark = ",")))
}

#' Euclidean sample-to-sample distances on the transformed matrix.
qc_distance_heatmap <- function(res, annotate_by = res$main, style = default_style()) {
  st <- resolve_style(style)
  m <- res$transformed
  d <- as.matrix(stats::dist(t(m)))
  vt <- qc_sample_table(res)
  ab <- intersect(annotate_by, setdiff(colnames(vt), "sample"))
  ann <- NULL; ann_col <- NULL
  if (length(ab)) {
    ann <- as.data.frame(vt[match(colnames(d), vt$sample), ab, drop = FALSE])
    rownames(ann) <- colnames(d)
    ann_col <- lapply(ab, function(a) {
      lv <- unique(stats::na.omit(ann[[a]])); stats::setNames(group_palette(length(lv), "Set2"), lv)
    })
    names(ann_col) <- ab
  }
  hc <- stats::hclust(stats::as.dist(d))
  pheatmap::pheatmap(
    d, clustering_distance_rows = stats::as.dist(d), clustering_distance_cols = stats::as.dist(d),
    cluster_rows = hc, cluster_cols = hc,
    annotation_col = ann, annotation_row = ann, annotation_colors = ann_col,
    color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100),
    main = sprintf("Sample distances (%s, Euclidean)", res$transform_label),
    fontsize = st$base_size, fontsize_row = st$axis_text_size, fontsize_col = st$axis_text_size,
    border_color = NA, display_numbers = ncol(d) <= 12, number_format = "%.0f",
    number_color = "grey20", silent = TRUE)
}

#' DESeq2 dispersion estimates against mean expression.
qc_dispersion_plot <- function(res, style = default_style()) {
  if (!is.null(res$qc$voom)) return(qc_voom_plot(res, style))
  if (!is.null(res$qc$bcv)) return(qc_bcv_plot(res, style))
  if (!is.null(res$qc$ql)) return(qc_ql_plot(res, style))
  if (!is.null(res$qc$precision)) return(qc_precision_plot(res, style))
  st <- resolve_style(style)
  dp <- res$qc$dispersion
  if (is.null(dp)) stop(if (identical(res$engine, "EBSeq"))
    "EBSeq fits no dispersion trend; see 'P-value histograms' for its PPDE distribution." else
    if (identical(res$engine, "satuRn")) "satuRn shrinks its dispersions internally and stores no plot data." else
    "No dispersion estimates are stored for this engine.")
  dp <- dp[is.finite(dp$baseMean) & dp$baseMean > 0, , drop = FALSE]
  pts <- rbind(
    data.frame(baseMean = dp$baseMean, value = dp$gene_est, what = "gene-wise estimate"),
    data.frame(baseMean = dp$baseMean, value = dp$final, what = "final (shrunk)"))
  # DESeq2 floors gene-wise estimates it cannot make at 1e-8; drawn, they
  # squash the informative range into the top of the plot.
  floor_hit <- sum(dp$gene_est <= 1e-7, na.rm = TRUE)
  pts <- pts[is.finite(pts$value) & pts$value > 1e-7, , drop = FALSE]
  pts$what <- factor(pts$what, levels = c("gene-wise estimate", "final (shrunk)"))
  fit <- dp[is.finite(dp$fit) & dp$fit > 0, , drop = FALSE]
  fit <- fit[order(fit$baseMean), , drop = FALSE]
  out <- dp[dp$outlier %in% TRUE & is.finite(dp$gene_est) & dp$gene_est > 0, , drop = FALSE]
  p <- ggplot2::ggplot(pts, ggplot2::aes(x = baseMean, y = value, colour = what)) +
    ggplot2::geom_point(size = st$point_size * 0.35, alpha = 0.45) +
    ggplot2::geom_line(data = fit, ggplot2::aes(x = baseMean, y = fit), inherit.aes = FALSE,
                       colour = "#D6604D", linewidth = st$line_width * 3) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = c("gene-wise estimate" = "grey35", "final (shrunk)" = "#4393C3"),
                                 name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    deg_theme(st) +
    ggplot2::labs(title = "Dispersion estimates", x = "Mean of normalised counts", y = "Dispersion",
                  subtitle = paste0("red line: fitted trend; outliers (circled) keep their gene-wise value",
                                    if (floor_hit) sprintf("; %s genes at DESeq2's minimum not drawn",
                                                           format(floor_hit, big.mark = ",")) else ""))
  if (nrow(out)) p <- p + ggplot2::geom_point(data = out, ggplot2::aes(x = baseMean, y = gene_est),
                                              inherit.aes = FALSE, shape = 1, colour = "#4393C3",
                                              size = st$point_size * 1.2)
  p
}

#' Histogram of raw p-values, one panel per comparison. A healthy test shows a
#' flat background with a peak near zero; a hump in the middle or a peak near
#' one points to a model problem.
qc_pvalue_plot <- function(res, ids = names(res$results), style = default_style()) {
  st <- resolve_style(style)
  ids <- intersect(ids, names(res$results))
  if (!length(ids)) stop("Choose at least one comparison.")
  if (identical(res$engine, "EBSeq")) return(qc_ppde_plot(res, ids, style))
  df <- do.call(rbind, lapply(ids, function(i) {
    p <- res$results[[i]]$pvalue; p <- p[is.finite(p)]
    data.frame(comparison = i, pvalue = p, stringsAsFactors = FALSE)
  }))
  df$comparison <- factor(df$comparison, levels = ids)
  ggplot2::ggplot(df, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.025), fill = "#4393C3", colour = "white",
                            linewidth = 0.2) +
    ggplot2::facet_wrap(~comparison, scales = "free_y") +
    deg_theme(st) +
    ggplot2::labs(title = "P-value histograms", x = "p-value",
                  y = if (isTRUE(res$level %in% c("transcript", "dtu"))) "Transcripts" else "Genes",
                  subtitle = "expected: flat background with a peak at 0")
}

#' edgeR quasi-likelihood dispersions (as plotQLDisp draws them).
qc_ql_plot <- function(res, style = default_style()) {
  st <- resolve_style(style); q <- res$qc$ql
  q <- q[is.finite(q$AveLogCPM), , drop = FALSE]
  pts <- rbind(data.frame(x = q$AveLogCPM, y = q$raw^0.25, what = "raw"),
               data.frame(x = q$AveLogCPM, y = q$squeezed^0.25, what = "squeezed"))
  pts$what <- factor(pts$what, levels = c("raw", "squeezed"))
  tr <- q[order(q$AveLogCPM), , drop = FALSE]
  ggplot2::ggplot(pts[is.finite(pts$y), ], ggplot2::aes(x = x, y = y, colour = what)) +
    ggplot2::geom_point(size = st$point_size * 0.35, alpha = 0.45) +
    ggplot2::geom_line(data = tr, ggplot2::aes(x = AveLogCPM, y = trend^0.25), inherit.aes = FALSE,
                       colour = "#D6604D", linewidth = st$line_width * 3) +
    ggplot2::scale_colour_manual(values = c(raw = "grey35", squeezed = "#4393C3"), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    deg_theme(st) +
    ggplot2::labs(title = "Quasi-likelihood dispersions", x = "Average log2 CPM",
                  y = "Quarter-root mean deviance", subtitle = "red line: trend (prior)")
}

#' DRIMSeq gene-wise precision against mean gene expression.
qc_precision_plot <- function(res, style = default_style()) {
  st <- resolve_style(style); p <- res$qc$precision
  p <- p[is.finite(p$precision) & p$mean_expression > 0, , drop = FALSE]
  ggplot2::ggplot(p, ggplot2::aes(x = mean_expression, y = precision)) +
    ggplot2::geom_point(size = st$point_size * 0.4, alpha = 0.45, colour = "grey30") +
    ggplot2::geom_hline(yintercept = p$common[1], colour = "#D6604D", linewidth = st$line_width * 3) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    deg_theme(st) +
    ggplot2::labs(title = "DRIMSeq precision", x = "Mean gene expression", y = "Gene-wise precision",
                  subtitle = "red line: common precision")
}

#' limma-voom mean-variance trend (the plot voom(plot = TRUE) draws).
qc_voom_plot <- function(res, style = default_style()) {
  st <- resolve_style(style); q <- res$qc$voom
  pts <- data.frame(x = q$x, y = q$y)
  ln <- data.frame(x = q$line_x, y = q$line_y)
  ggplot2::ggplot(pts, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(size = st$point_size * 0.35, alpha = 0.35, colour = "grey30") +
    ggplot2::geom_line(data = ln, colour = "#D6604D", linewidth = st$line_width * 3) +
    deg_theme(st) +
    ggplot2::labs(title = "voom: mean-variance trend", x = "log2(count size + 0.5)",
                  y = "Sqrt(standard deviation)",
                  subtitle = "each point a gene; red: fitted trend used for the precision weights")
}

#' edgeR biological coefficient of variation (as plotBCV draws it).
qc_bcv_plot <- function(res, style = default_style()) {
  st <- resolve_style(style); b <- res$qc$bcv
  b <- b[order(b$AveLogCPM), , drop = FALSE]
  ggplot2::ggplot(b, ggplot2::aes(x = AveLogCPM)) +
    ggplot2::geom_point(ggplot2::aes(y = tagwise), size = st$point_size * 0.35, alpha = 0.35, colour = "grey30") +
    ggplot2::geom_line(ggplot2::aes(y = trended), colour = "#D6604D", linewidth = st$line_width * 3) +
    ggplot2::geom_hline(yintercept = b$common[1], colour = "#4393C3", linewidth = st$line_width * 2) +
    deg_theme(st) +
    ggplot2::labs(title = "edgeR: biological coefficient of variation", x = "Average log2 CPM",
                  y = "BCV", subtitle = sprintf("grey: gene-wise; red: trend; blue: common BCV = %.3f", b$common[1]))
}

#' EBSeq posterior probability of differential expression, per comparison.
qc_ppde_plot <- function(res, ids = names(res$results), style = default_style()) {
  st <- resolve_style(style)
  df <- do.call(rbind, lapply(ids, function(i) {
    p <- res$results[[i]]$PPDE; p <- p[is.finite(p)]
    data.frame(comparison = i, ppde = p, stringsAsFactors = FALSE)
  }))
  df$comparison <- factor(df$comparison, levels = ids)
  ggplot2::ggplot(df, ggplot2::aes(x = ppde)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.025), fill = "#4393C3", colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = 0.95, linetype = "dashed", colour = "#D6604D") +
    ggplot2::facet_wrap(~comparison, scales = "free_y") +
    ggplot2::scale_y_sqrt() +
    deg_theme(st) +
    ggplot2::labs(title = "Posterior probability of differential expression (PPDE)", x = "PPDE",
                  y = "Genes (square-root scale)", subtitle = "dashed line: PPDE = 0.95 (FDR 5%)")
}

#' Dispatcher used by the Analysis tab.
build_qc_figure <- function(which, res, counts = NULL, colour_by = NULL, shape_by = NULL,
                            ntop = 500, annotate_by = NULL, ids = NULL, style = default_style()) {
  switch(which,
    library      = qc_library_plot(res, style),
    distribution = qc_distribution_plot(res, counts, style),
    pca          = qc_pca_plot(res, colour_by %||% res$main, shape_by, ntop, style),
    distance     = qc_distance_heatmap(res, annotate_by %||% res$main, style),
    dispersion   = qc_dispersion_plot(res, style),
    pvalues      = qc_pvalue_plot(res, ids %||% names(res$results), style),
    stop("Unknown QC plot: ", which))
}
