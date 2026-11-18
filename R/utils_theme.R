# =============================================================================
# utils_theme.R -- shared plot styling
#
# Every plot in the app draws its typography and point geometry from a single
# `style` list, so the sidebar's font/point-size controls behave identically
# across the volcano, MA, scatter, PCA, scree and enrichment panels.
# =============================================================================

#' Default styling. The UI overrides individual entries; anything omitted
#' falls back to these values.
default_style <- function() {
  list(
    base_size      = 13,    # global text size (pt)
    title_size     = 15,
    subtitle_size  = 12,
    axis_title_size= 13,
    axis_text_size = 11,
    legend_title_size = 12,
    legend_text_size  = 11,
    label_size     = 3.2,   # gene tag size (ggplot mm units)
    label_bold     = TRUE,
    label_colour   = "#000000",
    point_size     = 1.6,
    point_alpha    = 0.7,
    highlight_size = 2.8,   # points that carry a gene tag
    line_width     = 0.4,
    base_family    = "",    # "" = device default
    legend_position= "right",
    grid           = TRUE,
    panel_border   = TRUE
  )
}

#' Merge user overrides onto the defaults.
resolve_style <- function(...) {
  st <- default_style()
  ov <- list(...)
  if (length(ov) == 1 && is.list(ov[[1]]) && is.null(names(ov))) ov <- ov[[1]]
  for (nm in names(ov)) if (!is.null(ov[[nm]])) st[[nm]] <- ov[[nm]]
  st
}

#' The ggplot2 theme used by every panel.
deg_theme <- function(style = default_style()) {
  st <- resolve_style(style)
  th <- ggplot2::theme_bw(base_size = st$base_size, base_family = st$base_family) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = st$title_size, face = "bold", hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(size = st$subtitle_size, hjust = 0.5),
      axis.title      = ggplot2::element_text(size = st$axis_title_size, face = "bold"),
      axis.text       = ggplot2::element_text(size = st$axis_text_size),
      legend.title    = ggplot2::element_text(size = st$legend_title_size),
      legend.text     = ggplot2::element_text(size = st$legend_text_size),
      legend.position = st$legend_position,
      strip.text      = ggplot2::element_text(size = st$axis_title_size, face = "bold"),
      strip.background= ggplot2::element_rect(fill = "grey92", colour = NA)
    )
  th <- th + if (st$grid) {
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92"))
  } else {
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  if (!st$panel_border) th <- th + ggplot2::theme(panel.border = ggplot2::element_blank())
  th
}

#' Repelled gene tags with styling taken from `style`.
#' Any argument passed through `...` (e.g. `fontface`) overrides the styled
#' default rather than being supplied twice.
gene_label_layer <- function(data, style = default_style(),
                             max_overlaps = 30, seed = 42, ...) {
  st <- resolve_style(style)
  dots <- list(...)
  args <- list(
    data = data,
    colour   = st$label_colour,
    size     = st$label_size,
    fontface = if (isTRUE(st$label_bold)) "bold" else "plain",
    family   = st$base_family,
    max.overlaps = max_overlaps,
    min.segment.length = 0,
    segment.size = 0.3,
    segment.colour = "grey40",
    box.padding = 0.35,
    seed = seed
  )
  args[names(dots)] <- dots
  do.call(ggrepel::geom_text_repel, args)
}

#' Resolve an axis label: user text wins, otherwise the computed default.
#' Accepts NULL/"" as "use the default".
axis_label <- function(user_value, default_value) {
  if (is.null(user_value)) return(default_value)
  uv <- trimws(as.character(user_value)[1])
  if (!nzchar(uv)) return(default_value)
  uv
}

#' Apply user-supplied title / axis labels to a plot in one step.
apply_labels <- function(p, title = NULL, subtitle = NULL,
                         x = NULL, y = NULL,
                         default_title = NULL, default_subtitle = NULL,
                         default_x = NULL, default_y = NULL) {
  ttl <- axis_label(title,    default_title)
  sub <- axis_label(subtitle, default_subtitle)
  xl  <- axis_label(x,        default_x)
  yl  <- axis_label(y,        default_y)
  p + ggplot2::labs(
    title = if (is.null(ttl) || identical(ttl, "")) NULL else ttl,
    subtitle = if (is.null(sub) || identical(sub, "")) NULL else sub,
    x = xl, y = yl
  )
}

# Colour-blind-safe defaults used across panels
DEG_COLOURS <- list(
  up   = "#D6604D",
  down = "#4393C3",
  ns   = "grey75",
  highlight = "#E4211C",
  offtarget = "#B2182B"
)

#' Discrete palette for grouping variables (PCA, heatmap annotations).
group_palette <- function(n, name = "Set1") {
  if (n <= 0) return(character(0))
  maxn <- RColorBrewer::brewer.pal.info[name, "maxcolors"]
  if (n <= maxn) {
    cols <- RColorBrewer::brewer.pal(max(3, n), name)[seq_len(n)]
  } else {
    cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(maxn, name))(n)
  }
  cols
}
