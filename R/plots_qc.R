# =============================================================================
# plots_qc.R -- PCA, scree plot, sample correlation and heatmaps
# =============================================================================

#' Run PCA on the most variable genes of an expression matrix.
#'
#' @param expr   genes x samples, already on a sensible (usually log) scale
#' @param ntop   number of most-variable genes to use (Inf = all)
#' @param scale_genes  scale each gene to unit variance before PCA
run_pca <- function(expr, ntop = 500, scale_genes = FALSE) {
  expr <- expr[stats::complete.cases(expr), , drop = FALSE]
  if (nrow(expr) < 2 || ncol(expr) < 2) {
    stop("PCA needs at least 2 genes and 2 samples after filtering.")
  }
  vars <- matrixStats::rowVars(as.matrix(expr))
  keep <- vars > 0 & is.finite(vars)
  expr <- expr[keep, , drop = FALSE]; vars <- vars[keep]
  if (nrow(expr) < 2) stop("No genes vary across samples, so PCA cannot be run.")

  ntop <- min(ntop, nrow(expr))
  sel  <- order(vars, decreasing = TRUE)[seq_len(ntop)]
  mat  <- t(expr[sel, , drop = FALSE])

  pca <- stats::prcomp(mat, center = TRUE, scale. = scale_genes)
  pct <- pca$sdev^2 / sum(pca$sdev^2) * 100
  list(pca = pca,
       scores = as.data.frame(pca$x),
       percent = pct,
       ngenes = ntop)
}

#' PCA score plot, coloured/shaped by sample metadata columns.
make_pca_plot <- function(pca_res, meta, colour_by = "condition", shape_by = NULL,
                          pc_x = 1, pc_y = 2, style = default_style(),
                          show_labels = TRUE, show_ellipse = FALSE,
                          palette = "Set1",
                          title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                          max_overlaps = 30) {
  st <- resolve_style(style)
  scores <- pca_res$scores
  npc <- ncol(scores)
  if (pc_x > npc || pc_y > npc) {
    stop(sprintf("Only %d principal component(s) are available.", npc))
  }
  df <- data.frame(sample = rownames(scores),
                   x = scores[[paste0("PC", pc_x)]],
                   y = scores[[paste0("PC", pc_y)]],
                   stringsAsFactors = FALSE)
  idx <- match(df$sample, meta$sample)
  for (cl in setdiff(colnames(meta), "sample")) df[[cl]] <- meta[[cl]][idx]

  has_colour <- !is.null(colour_by) && colour_by %in% colnames(df)
  has_shape  <- !is.null(shape_by) && shape_by %in% colnames(df) && !identical(shape_by, "(none)")

  aes_args <- list(x = quote(x), y = quote(y))
  if (has_colour) aes_args$colour <- as.name(colour_by)
  if (has_shape)  aes_args$shape  <- as.name(shape_by)

  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = st$point_size * 3, alpha = min(st$point_alpha + 0.25, 1)) +
    deg_theme(st)

  if (has_colour) {
    lv <- unique(stats::na.omit(df[[colour_by]]))
    p <- p + ggplot2::scale_colour_manual(values = stats::setNames(group_palette(length(lv), palette), lv),
                                          name = colour_by)
  }
  if (show_ellipse && has_colour) {
    tab <- table(df[[colour_by]])
    if (any(tab >= 3)) {
      p <- p + ggplot2::stat_ellipse(type = "norm", linetype = "dashed",
                                     linewidth = st$line_width, show.legend = FALSE)
    }
  }
  if (show_labels) {
    p <- p + gene_label_layer(df, st, max_overlaps,
                              mapping = ggplot2::aes(label = sample),
                              show.legend = FALSE)
  }
  apply_labels(p, title, subtitle, xlab, ylab,
               default_title = "Principal component analysis",
               default_subtitle = sprintf("top %s most variable genes",
                                          format(pca_res$ngenes, big.mark = ",")),
               default_x = sprintf("PC%d (%.1f%% variance)", pc_x, pca_res$percent[pc_x]),
               default_y = sprintf("PC%d (%.1f%% variance)", pc_y, pca_res$percent[pc_y]))
}

#' Scree plot: variance explained per PC, with a cumulative overlay.
make_scree_plot <- function(pca_res, n_pc = 10, style = default_style(),
                            show_cumulative = TRUE, bar_colour = "#4393C3",
                            line_colour = "#D6604D",
                            title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  pct <- pca_res$percent
  n_pc <- min(n_pc, length(pct))
  df <- data.frame(PC = factor(paste0("PC", seq_len(n_pc)),
                               levels = paste0("PC", seq_len(n_pc))),
                   variance = pct[seq_len(n_pc)],
                   cumulative = cumsum(pct)[seq_len(n_pc)],
                   stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC, y = variance)) +
    ggplot2::geom_col(fill = bar_colour, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", variance)),
                       vjust = -0.5, size = st$label_size, family = st$base_family) +
    deg_theme(st)

  if (show_cumulative) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = cumulative, group = 1),
                         colour = line_colour, linewidth = st$line_width * 2) +
      ggplot2::geom_point(ggplot2::aes(y = cumulative),
                          colour = line_colour, size = st$point_size * 2) +
      ggplot2::scale_y_continuous(
        limits = c(0, max(105, max(df$cumulative) * 1.05)),
        sec.axis = ggplot2::sec_axis(~., name = "Cumulative variance explained (%)"))
  }
  apply_labels(p, title, subtitle, xlab, ylab,
               default_title = "Scree plot", default_subtitle = NULL,
               default_x = "Principal component",
               default_y = "Variance explained (%)")
}

#' Sample-to-sample correlation heatmap (quick QC of clustering behaviour).
make_sample_correlation <- function(expr, meta, method = "pearson",
                                    annotate_by = "condition",
                                    palette = "RdBu", style = default_style(),
                                    title = "Sample correlation") {
  expr <- expr[stats::complete.cases(expr), , drop = FALSE]
  if (ncol(expr) < 2) stop("Need at least 2 samples for a correlation heatmap.")
  cm <- stats::cor(expr, method = method, use = "pairwise.complete.obs")
  ann <- NULL
  ab <- intersect(annotate_by, colnames(meta))
  if (length(ab)) {
    ann <- as.data.frame(meta[match(colnames(cm), meta$sample), ab, drop = FALSE])
    rownames(ann) <- colnames(cm)
  }
  st <- resolve_style(style)
  pheatmap::pheatmap(
    cm, annotation_col = ann, annotation_row = ann,
    color = heatmap_palette(palette), main = title,
    fontsize = st$base_size, fontsize_row = st$axis_text_size,
    fontsize_col = st$axis_text_size, border_color = NA, silent = TRUE
  )
}

# ---- heatmap ----------------------------------------------------------------

heatmap_palette <- function(name = "RdBu", n = 100) {
  switch(name,
    "RdBu"     = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(n),
    "RdYlBu"   = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(n),
    "Viridis"  = if (requireNamespace("viridisLite", quietly = TRUE))
                   viridisLite::viridis(n) else
                   grDevices::colorRampPalette(c("#440154", "#21908C", "#FDE725"))(n),
    "Magma"    = if (requireNamespace("viridisLite", quietly = TRUE))
                   viridisLite::magma(n) else
                   grDevices::colorRampPalette(c("#000004", "#B63679", "#FCFDBF"))(n),
    "navy"     = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(n),
    "GreenRed" = grDevices::colorRampPalette(c("#1B7837", "white", "#762A83"))(n),
    grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(n)
  )
}

#' Order samples so that the reference (control) group is on the left.
#'
#' pheatmap draws columns in the order given when `cluster_cols = FALSE`;
#' when clustering is on, the dendrogram decides, which is why control samples
#' could previously appear on the right. Passing the ordered sample vector here
#' and leaving clustering off keeps control-left / treated-right.
order_samples_by_group <- function(samples, meta, reference = NULL,
                                   group_col = "condition") {
  idx <- match(samples, meta$sample)
  grp <- as.character(meta[[group_col]][idx])
  grp[is.na(grp)] <- "ungrouped"
  levels_present <- unique(grp)
  if (!is.null(reference) && reference %in% levels_present) {
    lev <- c(reference, sort(setdiff(levels_present, reference)))
  } else {
    lev <- sort(levels_present)
  }
  ord <- order(factor(grp, levels = lev), samples)
  samples[ord]
}

#' Build a heatmap of selected genes.
#'
#' @param gene_ids  genes to show; use `all_genes = TRUE` to plot everything
#' @param expr      expression matrix (VST / log-normalised recommended)
#' @param meta      sample sheet
#' @param de        canonical DE table (for gene symbols)
#' @param reference reference/control condition -> drawn on the left
make_heatmap <- function(gene_ids, expr, meta, de = NULL,
                         all_genes = FALSE,
                         scale_rows = TRUE,
                         cluster_rows = TRUE, cluster_cols = FALSE,
                         show_rownames = TRUE, show_colnames = TRUE,
                         annotate_by = "condition",
                         reference = NULL,
                         palette = "RdBu",
                         style = default_style(),
                         max_genes = 5000,
                         title = "Heatmap",
                         gap_between_groups = TRUE,
                         group_col = "condition") {
  st <- resolve_style(style)
  if (all_genes) {
    gene_ids <- rownames(expr)
  }
  gene_ids <- unique(gene_ids)
  present <- intersect(gene_ids, rownames(expr))
  if (length(present) < 2) {
    stop(sprintf("Only %d of the requested genes are in the expression matrix; at least 2 are needed.",
                 length(present)))
  }
  truncated <- FALSE
  if (length(present) > max_genes) {
    v <- matrixStats::rowVars(as.matrix(expr[present, , drop = FALSE]))
    present <- present[order(v, decreasing = TRUE)[seq_len(max_genes)]]
    truncated <- TRUE
  }

  mat <- expr[present, , drop = FALSE]
  samples <- intersect(colnames(mat), meta$sample)
  if (length(samples) < 2) stop("Fewer than 2 samples of the sample sheet are present in the matrix.")

  # control-left ordering; only meaningful when column clustering is off
  ordered_samples <- order_samples_by_group(samples, meta, reference, group_col)
  mat <- mat[, ordered_samples, drop = FALSE]

  # drop zero-variance rows, which break row scaling
  if (scale_rows) {
    rv <- matrixStats::rowVars(as.matrix(mat))
    keep <- is.finite(rv) & rv > 0
    if (sum(keep) < 2) stop("Fewer than 2 genes vary across samples, so row scaling is not possible.")
    mat <- mat[keep, , drop = FALSE]
  }

  # readable row labels
  if (!is.null(de)) {
    sym <- de$display_gene[match(rownames(mat), de$gene_id)]
    sym[is.na(sym)] <- rownames(mat)[is.na(sym)]
    rownames(mat) <- make.unique(sym)
  }

  ab <- intersect(annotate_by, colnames(meta))
  ann_col <- NULL
  if (length(ab)) {
    ann_col <- as.data.frame(meta[match(ordered_samples, meta$sample), ab, drop = FALSE])
    rownames(ann_col) <- ordered_samples
  }

  gaps <- NULL
  if (gap_between_groups && !cluster_cols && group_col %in% colnames(meta)) {
    g <- as.character(meta[[group_col]][match(ordered_samples, meta$sample)])
    gaps <- cumsum(rle(g)$lengths)
    gaps <- gaps[-length(gaps)]
  }

  n_show <- nrow(mat)
  auto_rownames <- show_rownames && n_show <= 120

  hm <- pheatmap::pheatmap(
    mat,
    scale          = if (scale_rows) "row" else "none",
    cluster_rows   = cluster_rows && nrow(mat) > 1,
    cluster_cols   = cluster_cols,
    show_rownames  = auto_rownames,
    show_colnames  = show_colnames,
    annotation_col = ann_col,
    gaps_col       = gaps,
    color          = heatmap_palette(palette),
    fontsize       = st$base_size,
    fontsize_row   = st$axis_text_size,
    fontsize_col   = st$axis_text_size,
    main           = title,
    border_color   = NA,
    silent         = TRUE
  )
  list(plot = hm, n_genes = nrow(mat), n_requested = length(gene_ids),
       truncated = truncated, rownames_hidden = show_rownames && !auto_rownames,
       samples = ordered_samples)
}
