# =============================================================================
# plots_enrichment.R -- the Enrichment gallery figures
#
#   enrichment-score bar / dot plot        (one database)
#   GO results of three ontologies         (horizontal bars)
#   GO BP / CC / MF three-in-one           (vertical bars, framed slanted labels)
#   Sankey + dot plot                      (genes -> terms, with a dot panel)
#   KEGG pathway category summary          (pathways grouped by KEGG class)
#   gene-concept network (circular)        enrichment map
#   GO graph (enrichplot::goplot)          KEGG pathway map (pathview)
#
# Everything except the last two is plain ggplot2 + grid working on the tables
# from enrichment_tables.R, so it needs no Bioconductor package and draws
# imported results exactly as it draws freshly computed ones.
# =============================================================================

# ---- palettes ---------------------------------------------------------------

ENRICH_PALETTES <- list(
  classic   = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
                "#A6761D", "#666666"),
  npg       = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
                "#91D1C2", "#DC0000", "#7E6148", "#B09C85"),
  lancet    = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91",
                "#AD002A", "#ADB6B6", "#1B1919"),
  aaas      = c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021",
                "#5F559B", "#A20056", "#808180", "#1B1919"),
  jama      = c("#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599",
                "#80796B"),
  nejm      = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD",
                "#FFDC91", "#EE4C97"),
  jco       = c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC", "#003C67",
                "#8F7700", "#3B3B3B", "#A73030", "#4A6990"),
  frontiers = c("#D51317", "#F39200", "#EFD500", "#95C11F", "#007B3D", "#31B7BC",
                "#0094CD", "#164194", "#6F286A", "#706F6F")
)

ENRICH_PALETTE_CHOICES <- c(
  "Classic (green / orange / purple)" = "classic", "NPG" = "npg", "Lancet" = "lancet",
  "AAAS" = "aaas", "JAMA" = "jama", "NEJM" = "nejm", "JCO" = "jco",
  "Frontiers" = "frontiers")

#' n colours from a named palette; long requests are interpolated rather than
#' recycled so neighbouring terms never share a colour.
gallery_palette <- function(name = "classic", n = 3) {
  pal <- ENRICH_PALETTES[[name]]
  if (is.null(pal)) pal <- ENRICH_PALETTES[["classic"]]
  if (n <= length(pal)) return(pal[seq_len(n)])
  grDevices::colorRampPalette(pal)(n)
}

# ---- small helpers ----------------------------------------------------------

.trunc_label <- function(x, max_chars = 50) {
  x <- as.character(x)
  max_chars <- max(8L, as.integer(max_chars))
  long <- nchar(x) > max_chars
  x[long] <- paste0(trimws(substr(x[long], 1, max_chars - 3)), "...")
  x
}

.wrap_label <- function(x, width = 45) {
  vapply(as.character(x), function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

#' Duplicate term names (the same description in two tables) would collapse
#' into one factor level; make them unique without changing what is displayed
#' more than necessary.
.unique_labels <- function(x) {
  d <- duplicated(x)
  while (any(d)) { x[d] <- paste0(x[d], " "); d <- duplicated(x) }
  x
}

#' Rendered width of strings, in inches. A throw-away null device is opened
#' when none is active so this never leaves an Rplots.pdf behind.
.text_inches <- function(labels, fontsize = 10, family = "", face = "plain") {
  if (!length(labels)) return(numeric(0))
  opened <- FALSE
  if (grDevices::dev.cur() == 1L) { grDevices::pdf(NULL); opened <- TRUE }
  on.exit(if (opened) grDevices::dev.off(), add = TRUE)
  gp <- grid::gpar(fontsize = fontsize, fontfamily = family, fontface = face)
  vapply(as.character(labels), function(l)
    grid::convertWidth(grid::grobWidth(grid::textGrob(l, gp = gp)), "inches", valueOnly = TRUE),
    numeric(1), USE.NAMES = FALSE)
}

.need_terms <- function(df, what = "this figure") {
  if (is.null(df) || !nrow(df))
    stop("No terms pass the current cutoff for ", what,
         ". Raise the cutoff, or judge significance by p-value instead of adjusted p.")
  invisible(TRUE)
}

.term_genes <- function(df) {
  gl <- strsplit(as.character(df$geneID), "/", fixed = TRUE)
  lapply(gl, function(g) unique(g[nzchar(g)]))
}

#' Theme for the bar-style gallery figures: axis lines, no panel box.
.gallery_bar_theme <- function(st) {
  st$panel_border <- FALSE
  deg_theme(st) +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
                   axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
                   axis.text = ggplot2::element_text(colour = "black", size = st$axis_text_size),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}

.score_label <- function(rank_by) "Enrichment score (-log10(p-value))"

# ---- 1. enrichment-score bar plot -------------------------------------------

#' Horizontal bars of -log10(p) for the top terms of one database.
make_score_barplot <- function(df, top_n = 10, rank_by = "p.adjust", cutoff = 0.05,
                               fill = "#ADD8E6", outline = "black",
                               show_values = FALSE, digits = 2, label_chars = 60,
                               style = default_style(), title = NULL,
                               xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  d <- gallery_top(df, top_n, rank_by, cutoff); .need_terms(d, "the bar plot")
  d$label <- .unique_labels(.wrap_label(.trunc_label(d$Description, label_chars * 2), label_chars))
  d$label <- factor(d$label, levels = rev(d$label))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Score, y = label)) +
    ggplot2::geom_col(fill = fill, colour = outline, linewidth = 0.6, width = 0.8) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(
      mult = c(0, if (isTRUE(show_values)) 0.12 else 0.06))) +
    .gallery_bar_theme(st)
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = formatC(Score, format = "f", digits = digits)),
                                hjust = -0.15, size = st$axis_text_size * 0.8 / ggplot2::.pt,
                                family = st$base_family)
  }
  apply_labels(p, title, NULL, xlab, ylab, default_title = "Enrichment results",
               default_x = .score_label(rank_by), default_y = NULL)
}

# ---- 2. enrichment-score dot plot -------------------------------------------

#' Dot plot: x = -log10(p) (or gene ratio), size = gene count, colour = p.
make_score_dotplot <- function(df, top_n = 10, rank_by = "p.adjust", cutoff = 0.05,
                               x_axis = c("score", "ratio"), label_chars = 60,
                               style = default_style(), title = NULL,
                               xlab = NULL, ylab = NULL) {
  st <- resolve_style(style); x_axis <- match.arg(x_axis)
  d <- gallery_top(df, top_n, rank_by, cutoff); .need_terms(d, "the dot plot")
  if (x_axis == "ratio" && all(is.na(d$GeneRatioNum))) x_axis <- "score"
  d$x <- if (x_axis == "ratio") d$GeneRatioNum else d$Score
  d$label <- .unique_labels(.wrap_label(.trunc_label(d$Description, label_chars * 2), label_chars))
  d$label <- factor(d$label, levels = rev(d$label))
  d$colour_value <- d[[if (identical(rank_by, "p.adjust")) "p.adjust" else "pvalue"]]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = label, size = Count, colour = colour_value)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "red", high = "blue",
                                   name = if (identical(rank_by, "p.adjust")) "adj. P" else "pvalue") +
    ggplot2::scale_size_continuous(range = c(st$point_size * 1.2, st$point_size * 4.5),
                                   breaks = function(l) unique(round(pretty(l, 4)))) +
    deg_theme(st) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
                   panel.grid.minor = ggplot2::element_blank())
  if (isTRUE(st$grid)) p <- p + ggplot2::theme(
    panel.grid.major = ggplot2::element_line(colour = "grey35", linetype = "dotted", linewidth = 0.35))
  apply_labels(p, title, NULL, xlab, ylab, default_title = "Enrichment results",
               default_x = if (x_axis == "ratio") "Gene ratio" else .score_label(rank_by),
               default_y = NULL)
}

# ---- 3. GO: three ontologies, horizontal ------------------------------------

#' Stack the top terms of BP, CC and MF into one frame.
go_three_df <- function(tables, top_n = 10, rank_by = "p.adjust", cutoff = 0.05) {
  parts <- lapply(intersect(c("BP", "CC", "MF"), names(tables)), function(o) {
    d <- gallery_top(tables[[o]], top_n, rank_by, cutoff)
    if (is.null(d) || !nrow(d)) return(NULL)
    d$Ontology <- o
    d[, c("ID", "Description", "Score", "Count", "pvalue", "p.adjust", "Ontology")]
  })
  out <- do.call(rbind, parts)
  if (is.null(out) || !nrow(out))
    stop("No GO terms pass the current cutoff in any ontology. Raise the cutoff, ",
         "or judge significance by p-value instead of adjusted p.")
  out$Ontology <- factor(out$Ontology, levels = intersect(c("BP", "CC", "MF"), unique(out$Ontology)))
  out[order(out$Ontology, -out$Score), , drop = FALSE]
}

.ontology_colours <- function(palette, levels) {
  stats::setNames(gallery_palette(palette, 3), c("BP", "CC", "MF"))[levels]
}

make_go_three_horizontal <- function(tables, top_n = 10, rank_by = "p.adjust", cutoff = 0.05,
                                     palette = "classic", show_values = FALSE, digits = 2,
                                     label_chars = 70, style = default_style(),
                                     title = NULL, xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  d <- go_three_df(tables, top_n, rank_by, cutoff)
  d$key <- factor(seq_len(nrow(d)), levels = rev(seq_len(nrow(d))))
  labs <- stats::setNames(.trunc_label(d$Description, label_chars), as.character(d$key))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Score, y = key, fill = Ontology)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_discrete(labels = labs) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(
      mult = c(0, if (isTRUE(show_values)) 0.12 else 0.06))) +
    ggplot2::scale_fill_manual(values = .ontology_colours(palette, levels(d$Ontology)), name = NULL) +
    .gallery_bar_theme(st)
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = formatC(Score, format = "f", digits = digits)),
                                hjust = -0.15, size = st$axis_text_size * 0.75 / ggplot2::.pt,
                                family = st$base_family)
  }
  apply_labels(p, title, NULL, xlab, ylab, default_title = "GO Results of Three Ontologies",
               default_x = "Enrichment Score", default_y = NULL)
}

# ---- 4. GO: BP / CC / MF three-in-one with framed labels --------------------

#' The slanted term labels and their frames, drawn below the panel.
#'
#' Positions mix npc (where each bar sits across the panel) with inches (how
#' far the label run extends), so the frames stay parallel to the text at any
#' figure size.
.framed_labels_grob <- function(x_npc, half_npc, labels, group, colours, group_names,
                                angle, run_in, fontsize, family, name_fontsize) {
  th <- angle * pi / 180
  dx <- run_in * cos(th); dy <- run_in * sin(th)
  top <- 0.05                                   # inches below the axis line
  u_in <- function(v) grid::unit(v, "inches")
  u_x  <- function(v) grid::unit(v, "npc")
  y0   <- grid::unit(0, "npc")
  kids <- list()
  for (g in levels(group)) {
    idx <- which(group == g)
    col <- colours[[g]]
    x1 <- min(x_npc[idx]) - half_npc; x2 <- max(x_npc[idx]) + half_npc
    kids[[length(kids) + 1]] <- grid::polygonGrob(
      x = grid::unit.c(u_x(x1), u_x(x2), u_x(x2) - u_in(dx), u_x(x1) - u_in(dx)),
      y = grid::unit.c(y0 - u_in(top), y0 - u_in(top),
                       y0 - u_in(top + dy), y0 - u_in(top + dy)),
      gp = grid::gpar(col = col, fill = NA, lwd = 1))
    kids[[length(kids) + 1]] <- grid::textGrob(
      labels[idx], x = u_x(x_npc[idx]) - u_in(0.05 * cos(th)),
      y = y0 - u_in(top + 0.05 * sin(th)),
      rot = angle, just = c("right", "centre"),
      gp = grid::gpar(col = col, fontsize = fontsize, fontfamily = family))
    kids[[length(kids) + 1]] <- grid::textGrob(
      group_names[[g]], x = u_x((x1 + x2) / 2) - u_in(dx),
      y = y0 - u_in(top + dy + 0.06), just = c("centre", "top"),
      gp = grid::gpar(col = col, fontsize = name_fontsize, fontfamily = family))
  }
  grid::gTree(children = do.call(grid::gList, kids))
}

#' Vertical bars for BP, CC and MF side by side; each ontology's term labels
#' run diagonally beneath its bars inside a frame of the ontology's colour.
make_go_three_vertical <- function(tables, top_n = 10, rank_by = "p.adjust", cutoff = 0.05,
                                   palette = "classic", angle = 70, label_chars = 45,
                                   show_values = FALSE, digits = 2, y_breaks = NULL,
                                   term_size = NULL, style = default_style(),
                                   title = NULL, xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  d <- go_three_df(tables, top_n, rank_by, cutoff)
  angle <- min(max(suppressWarnings(as.numeric(angle)), 30), 90)
  if (is.na(angle)) angle <- 70
  d$x <- seq_len(nrow(d)) + (as.integer(d$Ontology) - 1) * 1.2   # gap between ontologies
  x_lo <- 0.2; x_hi <- max(d$x) + 0.8
  labs <- .trunc_label(d$Description, label_chars)
  fs <- if (is.null(term_size)) st$axis_text_size * 0.72 else term_size
  name_fs <- fs * 1.05
  run_in <- max(.text_inches(labs, fs, st$base_family)) + 0.18
  th <- angle * pi / 180
  cols <- .ontology_colours(palette, levels(d$Ontology))

  grob <- .framed_labels_grob(
    x_npc = (d$x - x_lo) / (x_hi - x_lo), half_npc = 0.5 / (x_hi - x_lo),
    labels = labs, group = d$Ontology, colours = as.list(cols),
    group_names = as.list(GALLERY_DB_LABELS[levels(d$Ontology)]),
    angle = angle, run_in = run_in, fontsize = fs, family = st$base_family,
    name_fontsize = name_fs)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = Score, fill = Ontology)) +
    ggplot2::geom_col(width = 0.82) +
    ggplot2::annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    ggplot2::scale_x_continuous(limits = c(x_lo, x_hi), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    .gallery_bar_theme(st) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.8, "lines"),
      plot.margin = grid::unit(c(0.15, 0.2, run_in * sin(th) + 0.42,
                                 max(0.1, run_in * cos(th) - 0.3)), "inches"))
  brk <- suppressWarnings(as.numeric(y_breaks)); brk <- brk[!is.na(brk)]
  p <- p + if (length(brk)) {
    ggplot2::scale_y_continuous(breaks = brk, expand = ggplot2::expansion(mult = c(0, 0.06)),
                                limits = c(0, max(c(d$Score, brk)) * 1.02))
  } else {
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.06)))
  }
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = formatC(Score, format = "f", digits = digits)),
                                vjust = -0.4, size = fs * 0.85 / ggplot2::.pt,
                                family = st$base_family)
  }
  apply_labels(p, title, NULL, NULL, ylab, default_title = NULL,
               default_x = NULL, default_y = "Enrichment score")
}

# ---- 5. Sankey + dot plot ---------------------------------------------------

.sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot, id, n = 40) {
  t <- seq(0, 1, length.out = n)
  s <- t * t * (3 - 2 * t)
  x <- x0 + (x1 - x0) * t
  data.frame(x = c(x, rev(x)),
             y = c(y0_top + (y1_top - y0_top) * s, rev(y0_bot + (y1_bot - y0_bot) * s)),
             id = id, stringsAsFactors = FALSE)
}

#' Left: a Sankey diagram linking each gene to the terms it belongs to.
#' Right: a dot plot for the same terms, row-aligned with the Sankey nodes
#' (x = gene ratio, size = gene count, colour = -log10 p).
#'
#' Both halves are drawn in one coordinate system, so the alignment is exact
#' and the result is an ordinary ggplot.
make_sankey_dot <- function(df, top_n = 10, rank_by = "p.adjust", cutoff = 0.05,
                            max_genes = 100, palette = "npg",
                            sankey_ratio = 7, dot_ratio = 5,
                            x_label = "Gene ratio", p_label = NULL,
                            label_chars = 38, term_size = NULL,
                            style = default_style(), title = NULL) {
  st <- resolve_style(style)
  d <- gallery_top(df, top_n, rank_by, cutoff); .need_terms(d, "the Sankey + dot plot")
  genes <- .term_genes(d)
  keep_term <- lengths(genes) > 0
  if (!any(keep_term)) stop("The selected terms carry no gene members, so there is nothing to link.")
  d <- d[keep_term, , drop = FALSE]; genes <- genes[keep_term]
  nT <- nrow(d)
  edges <- data.frame(term = rep(seq_len(nT), lengths(genes)), gene = unlist(genes),
                      stringsAsFactors = FALSE)

  # Cap the gene count, keeping the genes shared by most terms.
  max_genes <- max(5L, as.integer(max_genes))
  freq <- sort(table(edges$gene), decreasing = TRUE)
  dropped <- 0L
  if (length(freq) > max_genes) {
    dropped <- length(freq) - max_genes
    edges <- edges[edges$gene %in% names(freq)[seq_len(max_genes)], , drop = FALSE]
  }
  # Order genes by the mean row of their terms: fewer crossing ribbons.
  bary <- tapply(edges$term, edges$gene, mean)
  g_levels <- names(bary)[order(bary, names(bary))]
  edges$g <- match(edges$gene, g_levels)
  nG <- length(g_levels)
  edges <- edges[order(edges$g, edges$term), , drop = FALSE]

  n_link <- nrow(edges)
  H <- n_link * 1.3
  stack_nodes <- function(sizes) {                 # top-down stacking with gaps
    n <- length(sizes)
    gap <- if (n > 1) (H - sum(sizes)) / (n - 1) else 0
    off <- if (n > 1) 0 else (H - sum(sizes)) / 2
    top <- H - off - cumsum(c(0, sizes[-n])) - gap * (seq_len(n) - 1)
    data.frame(top = top, bottom = top - sizes)
  }
  g_size <- as.numeric(table(factor(edges$g, levels = seq_len(nG))))
  t_size <- as.numeric(table(factor(edges$term, levels = seq_len(nT))))
  g_pos <- stack_nodes(g_size); t_pos <- stack_nodes(t_size)

  # Slot each link into its gene node (ordered by term) and term node (by gene).
  edges$g_top <- g_pos$top[edges$g] - (stats::ave(edges$g, edges$g, FUN = seq_along) - 1)
  eo <- order(edges$term, edges$g)
  edges$t_top <- NA_real_
  edges$t_top[eo] <- t_pos$top[edges$term[eo]] -
    (stats::ave(edges$term[eo], edges$term[eo], FUN = seq_along) - 1)

  S <- max(1, as.numeric(sankey_ratio)); D <- max(1, as.numeric(dot_ratio))
  node_w <- 0.035 * S
  xg1 <- 0.20 * S; xg0 <- xg1 - node_w             # gene nodes
  xt0 <- S - node_w; xt1 <- S                      # term nodes
  ribbons <- do.call(rbind, lapply(seq_len(n_link), function(i)
    cbind(.sigmoid_ribbon(xg1, xt0, edges$g_top[i], edges$g_top[i] - 1,
                          edges$t_top[i], edges$t_top[i] - 1, id = i),
          term = edges$term[i])))
  t_cols <- gallery_palette(palette, nT)
  t_lab <- .trunc_label(d$Description, label_chars)
  term_key <- factor(seq_len(nT))
  first_term <- tapply(edges$term, edges$g, min)
  g_nodes <- data.frame(xmin = xg0, xmax = xg1, ymin = g_pos$bottom, ymax = g_pos$top,
                        term = factor(as.integer(first_term), levels = levels(term_key)),
                        label = g_levels, y = (g_pos$top + g_pos$bottom) / 2)
  t_nodes <- data.frame(xmin = xt0, xmax = xt1, ymin = t_pos$bottom, ymax = t_pos$top,
                        term = term_key, label = t_lab, y = (t_pos$top + t_pos$bottom) / 2)
  ribbons$term <- factor(ribbons$term, levels = levels(term_key))

  # ---- dot panel geometry
  use_ratio <- !all(is.na(d$GeneRatioNum))
  xv <- if (use_ratio) d$GeneRatioNum else d$Count
  rng <- range(xv, na.rm = TRUE)
  pad <- if (diff(rng) > 0) diff(rng) * 0.12 else max(abs(rng[1]) * 0.1, 0.01)
  lo <- rng[1] - pad; hi <- rng[2] + pad
  px0 <- S * 1.05; px1 <- px0 + D
  to_x <- function(v) px0 + (v - lo) / (hi - lo) * (px1 - px0)
  y_lo <- -0.02 * H; y_hi <- 1.02 * H
  ticks <- pretty(c(lo, hi), n = 3); ticks <- ticks[ticks >= lo & ticks <= hi]
  if (!length(ticks)) ticks <- mean(c(lo, hi))
  tick_len <- 0.012 * H
  dots <- data.frame(x = to_x(xv), y = t_nodes$y, Count = d$Count, score = d$Score)

  fs <- if (is.null(term_size)) st$axis_text_size * 0.8 else term_size
  gene_fs <- min(fs, max(3.5, 380 / nG * 0.55))   # shrink when the gene column is crowded
  mm <- function(pt) pt / ggplot2::.pt
  if (is.null(p_label) || !nzchar(trimws(p_label))) p_label <- "-log10(pvalue)"
  if (!use_ratio && identical(x_label, "Gene ratio")) x_label <- "Gene count"

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = ribbons, ggplot2::aes(x = x, y = y, group = id, fill = term),
                          alpha = 0.42, colour = NA) +
    ggplot2::geom_rect(data = g_nodes, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                                    ymax = ymax, fill = term)) +
    ggplot2::geom_rect(data = t_nodes, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                                    ymax = ymax, fill = term)) +
    ggplot2::geom_text(data = g_nodes, ggplot2::aes(x = xmin - 0.012 * S, y = y, label = label),
                       hjust = 1, size = mm(gene_fs), family = st$base_family) +
    ggplot2::geom_text(data = t_nodes, ggplot2::aes(x = xmin - 0.012 * S, y = y, label = label),
                       hjust = 1, size = mm(fs), family = st$base_family) +
    ggplot2::scale_fill_manual(values = stats::setNames(t_cols, levels(term_key)), guide = "none")

  if (isTRUE(st$grid)) {
    p <- p + ggplot2::annotate("segment", x = to_x(ticks), xend = to_x(ticks),
                               y = y_lo, yend = y_hi, colour = "grey90", linewidth = 0.3) +
      ggplot2::annotate("segment", x = px0, xend = px1, y = t_nodes$y, yend = t_nodes$y,
                        colour = "grey93", linewidth = 0.3)
  }
  p <- p +
    ggplot2::annotate("rect", xmin = px0, xmax = px1, ymin = y_lo, ymax = y_hi,
                      fill = NA, colour = "black", linewidth = 0.5) +
    ggplot2::annotate("segment", x = to_x(ticks), xend = to_x(ticks),
                      y = y_lo, yend = y_lo - tick_len, colour = "black", linewidth = 0.4) +
    ggplot2::annotate("text", x = to_x(ticks), y = y_lo - tick_len * 1.6,
                      label = format(ticks, trim = TRUE), vjust = 1,
                      size = mm(st$axis_text_size * 0.85), family = st$base_family) +
    ggplot2::annotate("text", x = (px0 + px1) / 2, y = y_lo - tick_len * 1.6 - 0.045 * H,
                      label = x_label, vjust = 1, fontface = "bold",
                      size = mm(st$axis_title_size * 0.9), family = st$base_family) +
    ggplot2::geom_point(data = dots, ggplot2::aes(x = x, y = y, size = Count, colour = score)) +
    ggplot2::scale_colour_gradient(low = "blue", high = "red", name = p_label) +
    ggplot2::scale_size_continuous(range = c(st$point_size * 1.1, st$point_size * 4), name = "Count",
                                   breaks = function(l) unique(round(pretty(l, 4)))) +
    ggplot2::scale_x_continuous(limits = c(0, px1 + 0.01 * D), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(y_lo - tick_len * 1.6 - 0.11 * H, y_hi),
                                expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = st$base_size, base_family = st$base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = st$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = st$subtitle_size, hjust = 0.5, colour = "grey35"),
      legend.title = ggplot2::element_text(size = st$legend_title_size),
      legend.text = ggplot2::element_text(size = st$legend_text_size),
      plot.margin = grid::unit(c(0.15, 0.15, 0.1, 0.15), "inches"))
  sub <- if (dropped > 0) sprintf(
    "%d of %d genes shown (those shared by the most terms)", nG, nG + dropped) else NULL
  p <- p + ggplot2::labs(title = if (is.null(title) || !nzchar(title)) NULL else title,
                         subtitle = sub)
  attr(p, "n_genes") <- nG; attr(p, "n_dropped") <- dropped
  p
}

# ---- 6. KEGG pathway category summary ---------------------------------------

#' Pathways as horizontal bars, grouped and coloured by top-level KEGG class,
#' with a class key running down the right-hand side.
make_kegg_category_plot <- function(df, top_n = 30, rank_by = "p.adjust", cutoff = 0.05,
                                    value = c("count", "score"), palette = "npg",
                                    show_values = TRUE, digits = 2, label_chars = 55,
                                    style = default_style(), title = NULL,
                                    xlab = NULL, ylab = NULL) {
  st <- resolve_style(style); value <- match.arg(value)
  d <- gallery_top(df, top_n, rank_by, cutoff); .need_terms(d, "the pathway category plot")
  d$Class <- kegg_pathway_class(d$ID)
  present <- intersect(KEGG_CLASSES, unique(d$Class))
  d$Class <- factor(d$Class, levels = present)
  d$value <- if (value == "count") d$Count else d$Score
  d <- d[order(d$Class, -d$value), , drop = FALSE]
  n <- nrow(d)
  d$y <- rev(seq_len(n))
  cols <- stats::setNames(gallery_palette(palette, length(KEGG_CLASSES)), KEGG_CLASSES)[present]
  x_max <- max(d$value) * 1.02
  key <- do.call(rbind, lapply(present, function(cl) {
    yy <- d$y[d$Class == cl]
    data.frame(Class = cl, ymin = min(yy) - 0.42, ymax = max(yy) + 0.42)
  }))
  key$Class <- factor(key$Class, levels = present)
  key_fs <- st$legend_text_size * 0.85
  right_in <- max(.text_inches(present, key_fs, st$base_family)) + 0.25

  p <- ggplot2::ggplot(d) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = value, ymin = y - 0.38, ymax = y + 0.38,
                                    fill = Class)) +
    ggplot2::geom_segment(data = key, ggplot2::aes(x = x_max * 1.035, xend = x_max * 1.035,
                                                   y = ymin, yend = ymax, colour = Class),
                          linewidth = 1.6, lineend = "butt") +
    ggplot2::geom_text(data = key, ggplot2::aes(x = x_max * 1.055, y = ymax, label = Class),
                       hjust = 0, vjust = 1, size = key_fs / ggplot2::.pt,
                       family = st$base_family) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_y_continuous(breaks = d$y, labels = .trunc_label(d$Description, label_chars),
                                limits = c(0.4, n + 0.6)) +
    ggplot2::scale_x_continuous(breaks = { b <- pretty(c(0, x_max), 6); b[b <= x_max] }) +
    ggplot2::coord_cartesian(xlim = c(0, x_max * 1.035), expand = FALSE, clip = "off") +
    .gallery_bar_theme(st) +
    ggplot2::theme(plot.margin = grid::unit(c(0.15, right_in, 0.1, 0.1), "inches"))
  if (isTRUE(show_values)) {
    lab <- if (value == "count") format(d$value, trim = TRUE) else
      formatC(d$value, format = "f", digits = digits)
    p <- p + ggplot2::geom_text(ggplot2::aes(x = value, y = y, label = lab), hjust = 1.25,
                                size = st$axis_text_size * 0.72 / ggplot2::.pt,
                                family = st$base_family)
  }
  apply_labels(p, title, NULL, xlab, ylab, default_title = NULL,
               default_x = if (value == "count") "Gene count" else .score_label(rank_by),
               default_y = NULL)
}

# ---- 7. gene-concept network (circular) -------------------------------------

#' Terms and their genes on one circle, each gene joined to its terms by a
#' curve through the middle. Genes are coloured by log2 fold change when one
#' is known for them.
#'
#' @param fc named numeric vector, names = gene symbols
make_circular_cnet <- function(df, top_n = 10, rank_by = "p.adjust", cutoff = 0.05,
                               fc = NULL, max_genes = 150, palette = "npg",
                               label_chars = 40, style = default_style(), title = NULL) {
  st <- resolve_style(style)
  d <- gallery_top(df, top_n, rank_by, cutoff); .need_terms(d, "the gene-concept network")
  genes <- .term_genes(d)
  keep <- lengths(genes) > 0
  if (!any(keep)) stop("The selected terms carry no gene members, so there is nothing to draw.")
  d <- d[keep, , drop = FALSE]; genes <- genes[keep]
  nT <- nrow(d)
  edges <- data.frame(term = rep(seq_len(nT), lengths(genes)), gene = unlist(genes),
                      stringsAsFactors = FALSE)
  freq <- sort(table(edges$gene), decreasing = TRUE)
  dropped <- 0L
  if (length(freq) > max_genes) {
    dropped <- length(freq) - max_genes
    edges <- edges[edges$gene %in% names(freq)[seq_len(max_genes)], , drop = FALSE]
  }
  bary <- tapply(edges$term, edges$gene, mean)
  g_levels <- names(bary)[order(bary, names(bary))]
  nG <- length(g_levels)

  # Terms take a wider slot than genes so their labels have room.
  slot_t <- 2.2
  total <- nG + nT * slot_t
  a_gene <- (seq_len(nG) - 0.5) / total
  a_term <- (nG + (seq_len(nT) - 0.5) * slot_t) / total
  # Clockwise, rotated so the term arc is centred on the left-hand side: term
  # labels then run out into free space and the legends sit on the right.
  start <- pi + (nG + nT * slot_t / 2) / total * 2 * pi
  ang <- function(f) start - f * 2 * pi
  gx <- cos(ang(a_gene)); gy <- sin(ang(a_gene))
  tx <- cos(ang(a_term)); ty <- sin(ang(a_term))

  tt <- seq(0, 1, length.out = 30)
  curves <- do.call(rbind, lapply(seq_len(nrow(edges)), function(i) {
    g <- match(edges$gene[i], g_levels); t <- edges$term[i]
    # quadratic Bezier with the control point at the centre
    data.frame(x = (1 - tt)^2 * gx[g] + tt^2 * tx[t], y = (1 - tt)^2 * gy[g] + tt^2 * ty[t],
               id = i, term = t)
  }))
  t_lab <- .unique_labels(.trunc_label(d$Description, label_chars))
  curves$term <- factor(t_lab[curves$term], levels = t_lab)

  radial <- function(x, y, label, r = 1.06) {
    a <- atan2(y, x) * 180 / pi
    left <- x < 0
    data.frame(x = x * r, y = y * r, label = label,
               angle = ifelse(left, a + 180, a), hjust = ifelse(left, 1, 0))
  }
  g_lab <- radial(gx, gy, g_levels)
  # Term labels stay horizontal (wrapped), anchored on the side they sit on.
  t_lbl <- data.frame(x = tx * 1.09, y = ty * 1.09, label = .wrap_label(t_lab, 24), angle = 0,
                      hjust = ifelse(tx < -0.35, 1, ifelse(tx > 0.35, 0, 0.5)))
  t_lbl$y <- t_lbl$y + ifelse(abs(tx) <= 0.35, sign(ty) * 0.1, 0)
  g_fc <- if (is.null(fc)) rep(NA_real_, nG) else unname(fc[g_levels])
  g_pts <- data.frame(x = gx, y = gy, fc = g_fc)
  t_pts <- data.frame(x = tx, y = ty, Count = d$Count)
  has_fc <- any(!is.na(g_fc))
  gene_fs <- min(st$label_size, max(1.6, 150 / nG))
  room_x <- 1.15 + 0.6 * st$label_size / 3.2
  room_y <- 1.12 + max(nchar(g_levels)) * 0.022 * gene_fs / 3

  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = curves, ggplot2::aes(x = x, y = y, group = id, colour = term),
                       linewidth = 0.45, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = stats::setNames(gallery_palette(palette, nT), t_lab),
                                 name = "category") +
    ggplot2::geom_point(data = t_pts, ggplot2::aes(x = x, y = y, size = Count),
                        shape = 21, fill = "#E5C494", colour = "grey35", stroke = 0.3) +
    ggplot2::scale_size_continuous(range = c(st$point_size * 1.6, st$point_size * 4.5),
                                   name = "size",
                                   breaks = function(l) unique(round(pretty(l, 4))))
  p <- p + if (has_fc) {
    ggplot2::geom_point(data = g_pts, ggplot2::aes(x = x, y = y, fill = fc), shape = 21,
                        size = st$point_size * 1.5, colour = "grey30", stroke = 0.2)
  } else {
    ggplot2::geom_point(data = g_pts, ggplot2::aes(x = x, y = y), shape = 21,
                        size = st$point_size * 1.5, fill = "#8DA0CB", colour = "grey30", stroke = 0.2)
  }
  if (has_fc) {
    p <- p + if (min(g_fc, na.rm = TRUE) < 0 && max(g_fc, na.rm = TRUE) > 0) {
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                    midpoint = 0, name = "log2 fold change", na.value = "grey80")
    } else {
      ggplot2::scale_fill_gradient(low = "blue", high = "red", name = "log2 fold change",
                                   na.value = "grey80")
    }
  }
  txt <- function(dat, size, face) ggplot2::geom_text(
    data = dat, ggplot2::aes(x = x, y = y, label = label, angle = angle, hjust = hjust),
    size = size, fontface = face, family = st$base_family, colour = st$label_colour,
    lineheight = 0.85)
  p <- p + txt(g_lab, gene_fs, "plain") +
    txt(t_lbl, st$label_size, if (isTRUE(st$label_bold)) "bold" else "plain") +
    ggplot2::coord_fixed(xlim = c(-room_x, room_x), ylim = c(-room_y, room_y), clip = "off") +
    ggplot2::theme_void(base_size = st$base_size, base_family = st$base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = st$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = st$subtitle_size, hjust = 0.5, colour = "grey35"),
      legend.title = ggplot2::element_text(size = st$legend_title_size),
      legend.text = ggplot2::element_text(size = st$legend_text_size * 0.9)) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 3, override.aes = list(linewidth = 1)),
                    size = ggplot2::guide_legend(order = 2))
  sub <- if (dropped > 0) sprintf(
    "%d of %d genes shown (those shared by the most terms)", nG, nG + dropped) else NULL
  p + ggplot2::labs(title = if (is.null(title) || !nzchar(title)) NULL else title, subtitle = sub)
}

# ---- 8. enrichment map ------------------------------------------------------

#' Force-directed layout (Fruchterman-Reingold) in base R, enough for the few
#' dozen nodes an enrichment map holds. Deterministic for a given seed.
.fr_layout <- function(n, from, to, weight, iterations = 400, seed = 42) {
  if (n == 1) return(cbind(0, 0))
  a0 <- 2 * pi * (seq_len(n) - 1) / n
  if (exists(".Random.seed", envir = globalenv())) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  pos <- cbind(cos(a0), sin(a0)) + matrix(stats::runif(2 * n, -0.05, 0.05), n)
  k <- 2.4 / sqrt(n)
  temp <- 0.25
  for (it in seq_len(iterations)) {
    dx <- outer(pos[, 1], pos[, 1], "-"); dy <- outer(pos[, 2], pos[, 2], "-")
    dist <- pmax(sqrt(dx^2 + dy^2), 1e-3)
    rep_f <- k^2 / dist^2                                   # repulsion / distance
    disp <- cbind(rowSums(dx * rep_f), rowSums(dy * rep_f))
    if (length(from)) {
      ex <- pos[from, 1] - pos[to, 1]; ey <- pos[from, 2] - pos[to, 2]
      ed <- pmax(sqrt(ex^2 + ey^2), 1e-3)
      f <- ed / k * (0.15 + 0.5 * weight)                        # attraction / distance
      for (j in seq_along(from)) {
        disp[from[j], ] <- disp[from[j], ] - c(ex[j], ey[j]) * f[j]
        disp[to[j], ]   <- disp[to[j], ]   + c(ex[j], ey[j]) * f[j]
      }
    }
    disp <- disp - pos * 1.4                                # gentle pull to the centre
    len <- pmax(sqrt(rowSums(disp^2)), 1e-9)
    pos <- pos + disp / len * pmin(len, temp)
    temp <- temp * 0.992
  }
  pos <- sweep(pos, 2, colMeans(pos))
  pos / max(abs(pos))
}

#' Terms as nodes (size = gene count, colour = p), joined when their gene sets
#' overlap (Jaccard similarity >= `min_similarity`).
make_native_emap <- function(df, top_n = 20, rank_by = "p.adjust", cutoff = 0.05,
                             min_similarity = 0.2, label_chars = 45,
                             style = default_style(), title = NULL) {
  st <- resolve_style(style)
  d <- gallery_top(df, top_n, rank_by, cutoff); .need_terms(d, "the enrichment map")
  genes <- .term_genes(d)
  n <- nrow(d)
  from <- integer(0); to <- integer(0); w <- numeric(0)
  if (n > 1) for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    u <- length(union(genes[[i]], genes[[j]]))
    s <- if (u > 0) length(intersect(genes[[i]], genes[[j]])) / u else 0
    if (s >= min_similarity) { from <- c(from, i); to <- c(to, j); w <- c(w, s) }
  }
  pos <- .fr_layout(n, from, to, w)
  nodes <- data.frame(x = pos[, 1], y = pos[, 2], Count = d$Count,
                      label = .trunc_label(d$Description, label_chars),
                      p = d[[if (identical(rank_by, "p.adjust")) "p.adjust" else "pvalue"]])
  p <- ggplot2::ggplot()
  if (length(from)) {
    ed <- data.frame(x = pos[from, 1], y = pos[from, 2], xend = pos[to, 1], yend = pos[to, 2],
                     similarity = w)
    p <- p + ggplot2::geom_segment(data = ed, ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                                                           linewidth = similarity),
                                   colour = "grey70", lineend = "round") +
      ggplot2::scale_linewidth_continuous(range = c(0.5, 3), name = "similarity")
  }
  p <- p +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y, size = Count, colour = p)) +
    ggplot2::scale_colour_gradient(low = "red", high = "blue",
                                   name = if (identical(rank_by, "p.adjust")) "adj. P" else "pvalue") +
    ggplot2::scale_size_continuous(range = c(st$point_size * 1.6, st$point_size * 5), name = "size",
                                   breaks = function(l) unique(round(pretty(l, 4)))) +
    ggrepel::geom_text_repel(data = nodes, ggplot2::aes(x = x, y = y, label = label),
                             size = st$label_size, family = st$base_family,
                             fontface = if (isTRUE(st$label_bold)) "bold" else "plain",
                             colour = st$label_colour, max.overlaps = Inf, seed = 42,
                             point.padding = 0.4, box.padding = 0.5,
                             segment.colour = "grey60", segment.size = 0.25,
                             min.segment.length = 0.4) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.22)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.1)) +
    ggplot2::theme_void(base_size = st$base_size, base_family = st$base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = st$title_size, face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(size = st$legend_title_size),
      legend.text = ggplot2::element_text(size = st$legend_text_size))
  p + ggplot2::labs(title = if (is.null(title) || !nzchar(title)) NULL else title)
}

# ---- 9. GO graph (needs enrichplot) -----------------------------------------

#' Rebuild a clusterProfiler `enrichResult` from a gallery table, for the one
#' figure (the GO graph) that cannot be drawn from a plain table.
gallery_table_to_enrich_result <- function(df, ontology) {
  if (!requireNamespace("DOSE", quietly = TRUE))
    stop('The GO graph needs clusterProfiler / enrichplot: BiocManager::install(c("clusterProfiler", "enrichplot"))')
  res <- df[, GALLERY_STD_COLS, drop = FALSE]
  res <- res[!duplicated(res$ID), , drop = FALSE]
  rownames(res) <- res$ID
  gs <- stats::setNames(strsplit(res$geneID, "/", fixed = TRUE), res$ID)
  methods::new("enrichResult", result = res, pvalueCutoff = 1, pAdjustMethod = "BH",
               qvalueCutoff = 1, organism = "Homo sapiens", ontology = ontology,
               gene = unique(unlist(gs)), keytype = "SYMBOL", universe = character(0),
               gene2Symbol = character(0), geneSets = gs, readable = TRUE)
}

#' The induced GO graph of the top terms (enrichplot::goplot).
#'
#' @param obj an `enrichResult` from enrichGO, or NULL to rebuild one from `df`
make_go_graph <- function(obj = NULL, df = NULL, ontology = "BP", top_n = 10,
                          rank_by = "p.adjust", cutoff = 0.05, title = NULL) {
  if (!ontology %in% c("BP", "CC", "MF")) stop("The GO graph is only defined for GO ontologies.")
  if (!HAS_ENRICHPLOT)
    stop('The GO graph needs the enrichplot package: BiocManager::install("enrichplot")')
  if (!HAS_GGARCHERY)
    stop('The GO graph needs the ggarchery package: install.packages("ggarchery")')
  if (is.null(obj)) {
    .need_terms(gallery_sig(df, rank_by, cutoff), "the GO graph")
    obj <- gallery_table_to_enrich_result(gallery_sig(df, rank_by, cutoff), ontology)
  }
  if (!nrow(as.data.frame(obj))) stop("No significant terms to plot.")
  p <- tryCatch(enrichplot::goplot(obj, showCategory = top_n),
                error = function(e) stop("enrichplot::goplot failed: ", conditionMessage(e),
                                         " (it needs GO.db and GOSemSim, and GO ids as term ids)."))
  if (!is.null(title) && nzchar(title)) p <- p + ggplot2::ggtitle(title)
  p
}

# ---- 10. KEGG pathway map (needs pathview + internet) -----------------------

HAS_PATHVIEW <- has_pkg("pathview")

#' Render a KEGG pathway map coloured by log2 fold change, via pathview.
#' pathview downloads the map from KEGG, so this needs an internet connection.
#'
#' @param fc named numeric vector of log2 fold changes
#' @param id_type what the names of `fc` are: "SYMBOL", "ENTREZ" or "ENSEMBL"
#' @return path to the rendered PNG
render_pathview_map <- function(pathway_id, fc, id_type = "SYMBOL", species = NULL,
                                limit = NULL, out_dir = tempfile("pathview_")) {
  if (!HAS_PATHVIEW)
    stop('KEGG pathway maps need the pathview package: BiocManager::install("pathview")')
  pid <- trimws(as.character(pathway_id))
  num <- sub("^[A-Za-z]*", "", pid)
  if (!grepl("^[0-9]{5}$", num)) stop("'", pid, "' is not a KEGG pathway id (expected e.g. hsa04110).")
  if (is.null(species) || !nzchar(species)) {
    species <- tolower(sub("[0-9]+$", "", pid))
    if (!nzchar(species) || species %in% c("map", "ko")) species <- "hsa"
  }
  fc <- fc[!is.na(fc) & !is.na(names(fc)) & nzchar(names(fc))]
  fc <- fc[!duplicated(names(fc))]
  if (!length(fc)) stop("No fold changes are available to colour the map with.")
  if (is.null(limit) || is.na(limit) || limit <= 0) limit <- max(1, ceiling(max(abs(fc))))
  limit <- min(limit, 10)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  owd <- setwd(out_dir); on.exit(setwd(owd), add = TRUE)
  # pathview reads its bundled data sets by name, which only works once the
  # package is attached.
  if (!"package:pathview" %in% search())
    suppressPackageStartupMessages(attachNamespace("pathview"))
  res <- tryCatch(suppressWarnings(suppressMessages(
    pathview::pathview(gene.data = fc, pathway.id = num, species = species,
                       gene.idtype = id_type, kegg.dir = out_dir,
                       out.suffix = "pathview", kegg.native = TRUE,
                       limit = list(gene = limit, cpd = 1)))),
    error = function(e) e)
  png <- file.path(out_dir, sprintf("%s%s.pathview.png", species, num))
  if (!file.exists(png)) {
    why <- if (inherits(res, "error")) conditionMessage(res) else "KEGG did not return the map"
    stop("pathview could not render ", species, num, " (", why,
         "). It downloads the map from KEGG, so check the internet connection.")
  }
  png
}

#' A drawing function for a PNG on disk, so a pathway map goes through the
#' same TIFF / PDF / PNG export path as every other figure.
png_draw_function <- function(png_path) {
  force(png_path)
  function() {
    if (!requireNamespace("png", quietly = TRUE))
      stop('Re-exporting a pathway map needs the png package: install.packages("png")')
    img <- png::readPNG(png_path)
    grid::grid.newpage()
    grid::grid.raster(img, interpolate = TRUE)
  }
}

# ---- dispatcher -------------------------------------------------------------

GALLERY_FIGURES <- c(
  "Enrichment-score bar plot"                    = "score_bar",
  "Enrichment-score dot plot"                    = "score_dot",
  "GO: three ontologies (horizontal bars)"       = "go3_h",
  "GO: BP / CC / MF three-in-one (framed labels)" = "go3_v",
  "Sankey + dot plot"                            = "sankey",
  "KEGG pathway category summary"                = "kegg_class",
  "Gene-concept network (circular)"              = "cnet",
  "Enrichment map"                               = "emap",
  "GO graph (goplot)"                            = "goplot",
  "KEGG pathway map (pathview)"                  = "pathview"
)

#' Which database(s) a figure draws from.
gallery_figure_scope <- function(figure) {
  switch(figure,
    go3_h = , go3_v = "go3",
    kegg_class = , pathview = "kegg",
    goplot = "go1",
    "any")
}

#' Suggested export size (inches) for each figure.
gallery_figure_size <- function(figure) {
  switch(figure,
    score_bar = c(8, 6), score_dot = c(8, 6), go3_h = c(10, 7.5), go3_v = c(10, 7.5),
    sankey = c(10, 8), kegg_class = c(10, 6.5), cnet = c(11, 8.5), emap = c(10, 8),
    goplot = c(12, 10), c(8, 6))
}

#' Build one gallery figure.
#'
#' @param tables named list of gallery tables (BP / CC / MF / KEGG)
#' @param db     the database for single-database figures
#' @param opts   list of figure options; anything missing takes its default
build_gallery_figure <- function(figure, tables, db = NULL, opts = list(),
                                 style = default_style(), labels = list(),
                                 fc = NULL, objects = list()) {
  o <- function(name, default) { v <- opts[[name]]; if (is.null(v)) default else v }
  top_n <- o("top_n", 10); rank_by <- o("rank_by", "p.adjust"); cutoff <- o("cutoff", 0.05)
  ttl <- function(default) axis_label(labels$title, default)
  one <- function() {
    if (is.null(db) || is.null(tables[[db]]))
      stop("No ", if (is.null(db)) "enrichment" else GALLERY_DB_LABELS[[db]],
           " results are loaded for this figure.")
    tables[[db]]
  }
  switch(figure,
    score_bar = make_score_barplot(one(), top_n, rank_by, cutoff,
      fill = o("bar_colour", "#ADD8E6"), show_values = o("show_values", FALSE),
      digits = o("digits", 2), label_chars = o("label_chars", 60), style = style,
      title = ttl(GALLERY_DB_TITLES[[db]]), xlab = labels$xlab, ylab = labels$ylab),
    score_dot = make_score_dotplot(one(), top_n, rank_by, cutoff,
      x_axis = o("dot_x", "score"), label_chars = o("label_chars", 60), style = style,
      title = ttl(GALLERY_DB_TITLES[[db]]), xlab = labels$xlab, ylab = labels$ylab),
    go3_h = make_go_three_horizontal(tables, top_n, rank_by, cutoff,
      palette = o("palette", "classic"), show_values = o("show_values", FALSE),
      digits = o("digits", 2), label_chars = o("label_chars", 60) + 10, style = style,
      title = labels$title, xlab = labels$xlab, ylab = labels$ylab),
    go3_v = make_go_three_vertical(tables, top_n, rank_by, cutoff,
      palette = o("palette", "classic"), angle = o("angle", 70),
      label_chars = min(o("label_chars", 60), o("slant_chars", 45)),
      show_values = o("show_values", FALSE), digits = o("digits", 2),
      y_breaks = o("y_breaks", NULL), style = style,
      title = labels$title, ylab = labels$ylab),
    sankey = make_sankey_dot(one(), top_n, rank_by, cutoff,
      max_genes = o("max_genes", 100), palette = o("palette", "npg"),
      sankey_ratio = o("sankey_ratio", 7), dot_ratio = o("dot_ratio", 5),
      x_label = axis_label(labels$xlab, "Gene ratio"), p_label = labels$ylab,
      label_chars = min(o("label_chars", 60), 42), style = style, title = labels$title),
    kegg_class = make_kegg_category_plot(tables[["KEGG"]], top_n, rank_by, cutoff,
      value = o("kegg_value", "count"), palette = o("palette", "npg"),
      show_values = o("show_values", TRUE), digits = o("digits", 2),
      label_chars = o("label_chars", 60), style = style,
      title = labels$title, xlab = labels$xlab, ylab = labels$ylab),
    cnet = make_circular_cnet(one(), top_n, rank_by, cutoff, fc = fc,
      max_genes = o("max_genes", 100), palette = o("palette", "npg"),
      label_chars = min(o("label_chars", 60), 45), style = style, title = labels$title),
    emap = make_native_emap(one(), top_n, rank_by, cutoff,
      min_similarity = o("min_similarity", 0.2),
      label_chars = min(o("label_chars", 60), 50), style = style, title = labels$title),
    goplot = make_go_graph(objects[[db]], one(), db, top_n, rank_by, cutoff, labels$title),
    stop("Unknown gallery figure: ", figure))
}

# ---- one-click results bundle -----------------------------------------------

#' Write the full set of tables and figures into one folder:
#'
#'   go/       go.BP.all.xlsx  go.BP.sig.xlsx  BP_Enrichment_Score_barplot.*
#'             BP_Enrichment_Score_dotplot.*  BP_cnetplot.*  BP_emapplot.*
#'             [BP_goplot.*]   (same for CC, MF)   GO_Three_Ontologies.*
#'             GO_BP_CC_MF_three_in_one.*
#'   pathway/  pathway.all.xlsx  pathway.sig.xlsx  Pathway_Enrichment_Score_*.*
#'             Pathway_cnetplot.*  Pathway_emapplot.*  Pathway_sankey_dot.*
#'             Pathway_category_summary.*  [hsaNNNNN.pathview.png]
#'
#' A figure that cannot be drawn is skipped and recorded in bundle_log.txt;
#' one failure never costs the rest of the bundle.
#'
#' @return character vector: the log lines
write_gallery_bundle <- function(tables, out_dir, opts = list(), style = default_style(),
                                 fc = NULL, objects = list(),
                                 formats = c("png", "pdf", "tiff"), dpi = 300,
                                 pathview = FALSE, pathview_fc = fc,
                                 pathview_id_type = "SYMBOL", progress = NULL) {
  log <- character(0)
  note <- function(...) log <<- c(log, sprintf(...))
  tick <- function(msg) if (is.function(progress)) progress(msg)
  rank_by <- if (is.null(opts$rank_by)) "p.adjust" else opts$rank_by
  cutoff  <- if (is.null(opts$cutoff)) 0.05 else opts$cutoff

  save_fig <- function(dir, stem, figure, db = NULL, size = gallery_figure_size(figure)) {
    obj <- tryCatch(build_gallery_figure(figure, tables, db, opts, style, list(), fc, objects),
                    error = function(e) e)
    if (inherits(obj, "error")) { note("SKIPPED %s: %s", stem, conditionMessage(obj)); return() }
    for (fmt in formats) {
      f <- file.path(dir, paste0(stem, ".", fmt))
      r <- tryCatch({ export_plot(obj, f, fmt, size[1], size[2], dpi); NULL },
                    error = function(e) conditionMessage(e))
      if (!is.null(r)) note("FAILED %s.%s: %s", stem, fmt, r)
    }
    note("ok      %s", stem)
  }
  save_tables <- function(dir, stem, df) {
    write_table_xlsx_or_csv(gallery_export_table(df), file.path(dir, paste0(stem, ".all")))
    write_table_xlsx_or_csv(gallery_export_table(gallery_sig(df, rank_by, cutoff)),
                            file.path(dir, paste0(stem, ".sig")))
    note("ok      %s tables (%d terms, %d significant)", stem, nrow(df),
         nrow(gallery_sig(df, rank_by, cutoff)))
  }

  go <- intersect(c("BP", "CC", "MF"), names(tables))
  if (length(go)) {
    gd <- file.path(out_dir, "go"); dir.create(gd, recursive = TRUE, showWarnings = FALSE)
    for (o in go) {
      tick(paste("GO", o))
      save_tables(gd, paste0("go.", o), tables[[o]])
      save_fig(gd, paste0(o, "_Enrichment_Score_barplot"), "score_bar", o)
      save_fig(gd, paste0(o, "_Enrichment_Score_dotplot"), "score_dot", o)
      save_fig(gd, paste0(o, "_cnetplot"), "cnet", o)
      save_fig(gd, paste0(o, "_emapplot"), "emap", o)
      if (HAS_ENRICHPLOT && HAS_GGARCHERY) save_fig(gd, paste0(o, "_goplot"), "goplot", o)
    }
    tick("GO, three ontologies")
    save_fig(gd, "GO_Three_Ontologies", "go3_h")
    save_fig(gd, "GO_BP_CC_MF_three_in_one", "go3_v")
    if (!HAS_ENRICHPLOT) note("note    GO graphs (goplot) skipped: enrichplot is not installed")
    else if (!HAS_GGARCHERY) note("note    GO graphs (goplot) skipped: ggarchery is not installed")
  }
  if (!is.null(tables[["KEGG"]])) {
    tick("KEGG pathways")
    pd <- file.path(out_dir, "pathway"); dir.create(pd, recursive = TRUE, showWarnings = FALSE)
    save_tables(pd, "pathway", tables[["KEGG"]])
    save_fig(pd, "Pathway_Enrichment_Score_barplot", "score_bar", "KEGG")
    save_fig(pd, "Pathway_Enrichment_Score_dotplot", "score_dot", "KEGG")
    save_fig(pd, "Pathway_cnetplot", "cnet", "KEGG")
    save_fig(pd, "Pathway_emapplot", "emap", "KEGG")
    save_fig(pd, "Pathway_sankey_dot", "sankey", "KEGG")
    save_fig(pd, "Pathway_category_summary", "kegg_class")
    if (isTRUE(pathview)) {
      if (!HAS_PATHVIEW) {
        note("note    pathway maps skipped: pathview is not installed")
      } else {
        ids <- gallery_sig(tables[["KEGG"]], rank_by, cutoff)$ID
        ids <- ids[.is_kegg_id(ids)]
        for (id in ids) {
          tick(paste("pathway map", id))
          r <- tryCatch({
            png <- render_pathview_map(id, pathview_fc, pathview_id_type)
            file.copy(png, file.path(pd, basename(png)), overwrite = TRUE); NULL
          }, error = function(e) conditionMessage(e))
          if (is.null(r)) note("ok      %s.pathview", id) else note("SKIPPED %s.pathview: %s", id, r)
        }
      }
    }
  }
  if (!length(go) && is.null(tables[["KEGG"]])) stop("No enrichment results are loaded.")
  writeLines(c(sprintf("DEGexplorer enrichment bundle, written %s", format(Sys.time())),
               sprintf("Significance: %s <= %s; top %s terms per figure.", rank_by, cutoff,
                       if (is.null(opts$top_n)) 10 else opts$top_n), "", log),
             file.path(out_dir, "bundle_log.txt"))
  log
}
