# =============================================================================
# plots_de.R -- volcano, MA and expression-vs-expression scatter plots
# =============================================================================

#' Label each gene Up / Down / NS against the current thresholds.
classify_degs <- function(df, padj_cutoff = 0.05, lfc_cutoff = 1,
                          basemean_cutoff = 0) {
  status <- rep("NS", nrow(df))
  bm <- df$baseMean
  bm_ok <- if (all(is.na(bm))) rep(TRUE, nrow(df)) else
    (!is.na(bm) & bm >= basemean_cutoff)
  sig <- !is.na(df$padj) & df$padj <= padj_cutoff & bm_ok
  status[sig & !is.na(df$log2FoldChange) & df$log2FoldChange >=  lfc_cutoff] <- "Up"
  status[sig & !is.na(df$log2FoldChange) & df$log2FoldChange <= -lfc_cutoff] <- "Down"
  df$status <- factor(status, levels = c("Up", "Down", "NS"))
  df
}

#' Significant subset, most significant first.
significant_degs <- function(df) {
  s <- df[df$status %in% c("Up", "Down"), , drop = FALSE]
  s[order(s$padj, -abs(s$log2FoldChange)), , drop = FALSE]
}

.status_scale <- function(plot_df, up_colour, down_colour, ns_colour, style) {
  n_up   <- sum(plot_df$status == "Up",   na.rm = TRUE)
  n_down <- sum(plot_df$status == "Down", na.rm = TRUE)
  list(
    colour = ggplot2::scale_colour_manual(
      values = c(Up = up_colour, Down = down_colour, NS = ns_colour),
      labels = c(Up = sprintf("Up (%d)", n_up),
                 Down = sprintf("Down (%d)", n_down),
                 NS = "Not significant"),
      drop = FALSE, name = NULL
    ),
    size = ggplot2::scale_size_manual(
      values = c(Up = style$point_size + 0.3, Down = style$point_size + 0.3,
                 NS = style$point_size), guide = "none", drop = FALSE),
    alpha = ggplot2::scale_alpha_manual(
      values = c(Up = min(style$point_alpha + 0.15, 1),
                 Down = min(style$point_alpha + 0.15, 1),
                 NS = style$point_alpha * 0.6), guide = "none", drop = FALSE)
  )
}

# ---- volcano ----------------------------------------------------------------

make_volcano <- function(df, padj_cutoff = 0.05, lfc_cutoff = 1,
                         basemean_cutoff = 0, label_ids = character(0),
                         up_colour = DEG_COLOURS$up,
                         down_colour = DEG_COLOURS$down,
                         ns_colour = DEG_COLOURS$ns,
                         style = default_style(),
                         show_cutoff_lines = TRUE,
                         title = NULL, subtitle = NULL,
                         xlab = NULL, ylab = NULL,
                         max_overlaps = 30) {
  st <- resolve_style(style)
  plot_df <- classify_degs(df[!is.na(df$padj) & !is.na(df$log2FoldChange), , drop = FALSE],
                           padj_cutoff, lfc_cutoff, basemean_cutoff)
  if (nrow(plot_df) == 0) stop("No genes have both a fold change and an adjusted p-value.")
  plot_df$neglog10padj <- -log10(pmax(plot_df$padj, 1e-300))

  sc <- .status_scale(plot_df, up_colour, down_colour, ns_colour, st)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FoldChange, y = neglog10padj,
                                              colour = status)) +
    ggplot2::geom_point(ggplot2::aes(size = status, alpha = status)) +
    sc$colour + sc$size + sc$alpha + deg_theme(st)

  if (show_cutoff_lines) {
    p <- p +
      ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                          linetype = "dashed", colour = "grey40", linewidth = st$line_width) +
      ggplot2::geom_hline(yintercept = -log10(padj_cutoff),
                          linetype = "dashed", colour = "grey40", linewidth = st$line_width)
  }
  if (length(label_ids)) {
    lab <- plot_df[plot_df$gene_id %in% label_ids, , drop = FALSE]
    if (nrow(lab)) {
      p <- p +
        ggplot2::geom_point(data = lab, size = st$highlight_size, shape = 21,
                            fill = NA, colour = "black", stroke = 0.5) +
        gene_label_layer(lab, st, max_overlaps,
                         mapping = ggplot2::aes(label = display_gene))
    }
  }
  apply_labels(p, title, subtitle, xlab, ylab,
               default_title = "Volcano plot", default_subtitle = NULL,
               default_x = expression(bold(log[2]~"fold change")),
               default_y = expression(bold(-log[10]~"adjusted "*italic(P))))
}

# ---- MA ---------------------------------------------------------------------

make_ma <- function(df, padj_cutoff = 0.05, lfc_cutoff = 1,
                    basemean_cutoff = 0, label_ids = character(0),
                    up_colour = DEG_COLOURS$up,
                    down_colour = DEG_COLOURS$down,
                    ns_colour = DEG_COLOURS$ns,
                    style = default_style(),
                    show_cutoff_lines = TRUE,
                    title = NULL, subtitle = NULL,
                    xlab = NULL, ylab = NULL,
                    max_overlaps = 30) {
  st <- resolve_style(style)
  if (all(is.na(df$baseMean))) {
    stop("This results table has no mean-expression column, so an MA plot cannot be drawn. ",
         "Load a count matrix on the Data tab, or map a mean-expression column.")
  }
  plot_df <- classify_degs(df[!is.na(df$baseMean) & !is.na(df$log2FoldChange), , drop = FALSE],
                           padj_cutoff, lfc_cutoff, basemean_cutoff)
  plot_df$log10baseMean <- log10(plot_df$baseMean + 1)

  sc <- .status_scale(plot_df, up_colour, down_colour, ns_colour, st)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log10baseMean, y = log2FoldChange,
                                              colour = status)) +
    ggplot2::geom_point(ggplot2::aes(size = status, alpha = status)) +
    sc$colour + sc$size + sc$alpha + deg_theme(st)

  if (show_cutoff_lines) {
    p <- p +
      ggplot2::geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff),
                          linetype = "dashed", colour = "grey40", linewidth = st$line_width) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey20", linewidth = st$line_width * 0.75)
  }
  if (length(label_ids)) {
    lab <- plot_df[plot_df$gene_id %in% label_ids, , drop = FALSE]
    if (nrow(lab)) {
      p <- p +
        ggplot2::geom_point(data = lab, size = st$highlight_size, shape = 21,
                            fill = NA, colour = "black", stroke = 0.5) +
        gene_label_layer(lab, st, max_overlaps,
                         mapping = ggplot2::aes(label = display_gene))
    }
  }
  apply_labels(p, title, subtitle, xlab, ylab,
               default_title = "MA plot", default_subtitle = NULL,
               default_x = expression(bold(log[10]~"(mean expression + 1)")),
               default_y = expression(bold(log[2]~"fold change")))
}

# ---- group-mean expression scatter -----------------------------------------

#' Mean expression per condition group, on the requested scale.
#'
#' @param expr   normalised expression matrix (genes x samples)
#' @param meta   data frame with `sample` and `condition`
#' @param group_x,group_y condition levels to place on each axis
#' @param log_transform  apply log2(x + prior) (set FALSE if already logged)
group_mean_expression <- function(expr, meta, group_x, group_y,
                                  log_transform = TRUE, prior = 1) {
  meta <- meta[meta$sample %in% colnames(expr), , drop = FALSE]
  sx <- meta$sample[meta$condition == group_x]
  sy <- meta$sample[meta$condition == group_y]
  if (length(sx) == 0 || length(sy) == 0) {
    stop("Both groups need at least one sample present in the expression matrix.")
  }
  mx <- rowMeans(expr[, sx, drop = FALSE], na.rm = TRUE)
  my <- rowMeans(expr[, sy, drop = FALSE], na.rm = TRUE)
  if (log_transform) { mx <- log2(mx + prior); my <- log2(my + prior) }
  data.frame(gene_id = rownames(expr), x = mx, y = my,
             stringsAsFactors = FALSE)
}

#' Expression-vs-expression scatter (control on x, treated on y) with the
#' identity line, in the style of published transcriptome-wide off-target figures.
#'
#' @param scatter_df  from `group_mean_expression()`
#' @param de          canonical DE table supplying status + gene names
#' @param offtarget_ids  gene ids to mark with open diamonds
#' @param focus_ids      gene ids to mark with a filled highlight dot
make_expression_scatter <- function(scatter_df, de,
                                    padj_cutoff = 0.05, lfc_cutoff = 1,
                                    basemean_cutoff = 0,
                                    label_ids = character(0),
                                    offtarget_ids = character(0),
                                    focus_ids = character(0),
                                    up_colour = "#E8A33D",
                                    down_colour = "#2C5AA0",
                                    ns_colour = "grey80",
                                    style = default_style(),
                                    show_identity = TRUE,
                                    identity_colour = "#2E7D32",
                                    show_counts = TRUE,
                                    group_x_name = "control",
                                    group_y_name = "treated",
                                    value_label = "log2 RPKM",
                                    title = NULL, subtitle = NULL,
                                    xlab = NULL, ylab = NULL,
                                    max_overlaps = 30) {
  st <- resolve_style(style)
  de_cl <- classify_degs(de, padj_cutoff, lfc_cutoff, basemean_cutoff)
  idx <- match(scatter_df$gene_id, de_cl$gene_id)
  scatter_df$status <- de_cl$status[idx]
  scatter_df$status[is.na(scatter_df$status)] <- "NS"
  scatter_df$status <- factor(scatter_df$status, levels = c("Up", "Down", "NS"))
  scatter_df$display_gene <- de_cl$display_gene[idx]
  scatter_df$display_gene[is.na(scatter_df$display_gene)] <- scatter_df$gene_id[is.na(scatter_df$display_gene)]
  scatter_df <- scatter_df[is.finite(scatter_df$x) & is.finite(scatter_df$y), , drop = FALSE]
  if (nrow(scatter_df) == 0) stop("No genes have finite expression values in both groups.")

  n_up   <- sum(scatter_df$status == "Up",   na.rm = TRUE)
  n_down <- sum(scatter_df$status == "Down", na.rm = TRUE)

  # draw NS first so coloured points sit on top
  ord <- order(factor(scatter_df$status, levels = c("NS", "Down", "Up")))
  scatter_df <- scatter_df[ord, , drop = FALSE]

  p <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = x, y = y, colour = status))
  if (show_identity) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
                                  colour = identity_colour,
                                  linewidth = st$line_width * 1.5)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(size = status, alpha = status)) +
    ggplot2::scale_colour_manual(
      values = c(Up = up_colour, Down = down_colour, NS = ns_colour),
      labels = c(
        Up   = bquote(Up~(P[adj] < .(padj_cutoff)~and~log[2]*FC > .(lfc_cutoff))),
        Down = bquote(Down~(P[adj] < .(padj_cutoff)~and~log[2]*FC < .(-lfc_cutoff))),
        NS   = "Not significant"),
      drop = FALSE, name = NULL) +
    ggplot2::scale_size_manual(
      values = c(Up = st$point_size, Down = st$point_size, NS = st$point_size * 0.8),
      guide = "none", drop = FALSE) +
    ggplot2::scale_alpha_manual(
      values = c(Up = st$point_alpha, Down = st$point_alpha, NS = st$point_alpha * 0.55),
      guide = "none", drop = FALSE) +
    deg_theme(st)

  if (length(offtarget_ids)) {
    ot <- scatter_df[scatter_df$gene_id %in% offtarget_ids, , drop = FALSE]
    if (nrow(ot)) {
      p <- p + ggplot2::geom_point(data = ot, shape = 23, fill = NA,
                                   colour = DEG_COLOURS$offtarget,
                                   size = st$highlight_size,
                                   stroke = 0.7, inherit.aes = TRUE,
                                   show.legend = FALSE)
    }
  }
  if (length(focus_ids)) {
    fc <- scatter_df[scatter_df$gene_id %in% focus_ids, , drop = FALSE]
    if (nrow(fc)) {
      p <- p + ggplot2::geom_point(data = fc, shape = 21,
                                   fill = DEG_COLOURS$highlight,
                                   colour = DEG_COLOURS$highlight,
                                   size = st$highlight_size + 0.6,
                                   inherit.aes = FALSE,
                                   ggplot2::aes(x = x, y = y),
                                   show.legend = FALSE)
    }
  }
  lab_ids <- unique(c(label_ids, offtarget_ids, focus_ids))
  if (length(lab_ids)) {
    lab <- scatter_df[scatter_df$gene_id %in% lab_ids, , drop = FALSE]
    if (nrow(lab)) {
      p <- p + gene_label_layer(lab, st, max_overlaps,
                                mapping = ggplot2::aes(label = display_gene),
                                fontface = if (isTRUE(st$label_bold)) "bold.italic" else "italic",
                                show.legend = FALSE)
    }
  }
  if (show_counts) {
    # Placed with npc coordinates so the labels can never be clipped by the
    # panel edge, and anchored on the side each group occupies: up-regulated
    # genes sit above the identity line (top-left), down-regulated below it
    # (bottom-right).
    up_lab   <- sprintf("%s genes up",   format(n_up,   big.mark = ","))
    down_lab <- sprintf("%s genes down", format(n_down, big.mark = ","))
    p <- p +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = up_lab,
                        colour = up_colour, hjust = -0.12, vjust = 1.8,
                        fontface = "bold", size = st$label_size + 0.8,
                        family = st$base_family) +
      ggplot2::annotate("text", x = Inf, y = -Inf, label = down_lab,
                        colour = down_colour, hjust = 1.12, vjust = -1.4,
                        fontface = "bold", size = st$label_size + 0.8,
                        family = st$base_family)
  }
  apply_labels(p, title, subtitle, xlab, ylab,
               default_title = sprintf("%s versus %s", group_y_name, group_x_name),
               default_subtitle = NULL,
               default_x = sprintf("%s, %s", value_label, group_x_name),
               default_y = sprintf("%s, %s", value_label, group_y_name))
}
