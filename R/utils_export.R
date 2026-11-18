# =============================================================================
# utils_export.R -- figure export and gene-set selection shared by all tabs
# =============================================================================

# --- device helpers ----------------------------------------------------------

#' Sanity-check a user-supplied figure dimension.
.safe_dim <- function(x, default) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1 || is.na(x) || !is.finite(x) || x <= 0) return(default)
  min(max(x, 1), 60)
}

.safe_dpi <- function(x, default = 300) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1 || is.na(x) || !is.finite(x) || x < 36) return(default)
  min(x, 1200)
}

#' Open a graphics device, trying the best option first and degrading
#' gracefully. Returns the name of the device that actually opened.
#'
#' PDF notes: `cairo_pdf()` defaults to `onefile = FALSE`, which means that if a
#' figure ever draws more than one page (grid::grid.newpage(), a faceted grob
#' tree, enrichplot's network layouts) every page is written to the *same*
#' path and the earlier pages are clobbered -- the usual symptom being a PDF
#' that is empty, truncated or simply never finishes downloading. Forcing
#' `onefile = TRUE` fixes that, and the pdf() fallback keeps the button alive on
#' R builds without cairo.
.open_export_device <- function(file, format, width, height, dpi) {
  attempts <- if (format == "tiff") {
    list(
      list(name = "tiff", fun = function()
        grDevices::tiff(filename = file, width = width, height = height,
                        units = "in", res = dpi, compression = "lzw", bg = "white")),
      list(name = "tiff (no compression)", fun = function()
        grDevices::tiff(filename = file, width = width, height = height,
                        units = "in", res = dpi, bg = "white")),
      list(name = "png", fun = function()
        grDevices::png(filename = file, width = width, height = height,
                       units = "in", res = dpi, bg = "white"))
    )
  } else if (format == "png") {
    list(
      list(name = "png", fun = function()
        grDevices::png(filename = file, width = width, height = height,
                       units = "in", res = dpi, bg = "white"))
    )
  } else {
    c(
      if (isTRUE(capabilities("cairo"))) list(
        list(name = "cairo_pdf", fun = function()
          grDevices::cairo_pdf(filename = file, width = width, height = height,
                               onefile = TRUE, bg = "white"))
      ),
      list(
        list(name = "pdf", fun = function()
          grDevices::pdf(file = file, width = width, height = height,
                         onefile = TRUE, bg = "white", useDingbats = FALSE))
      )
    )
  }

  problems <- character(0)
  for (a in attempts) {
    ok <- tryCatch({ a$fun(); TRUE },
                   error   = function(e) { problems <<- c(problems, paste0(a$name, ": ", conditionMessage(e))); FALSE },
                   warning = function(w) { problems <<- c(problems, paste0(a$name, ": ", conditionMessage(w))); FALSE })
    if (isTRUE(ok)) return(a$name)
  }
  stop("No ", toupper(format), " graphics device could be opened. ",
       paste(problems, collapse = "; "))
}

#' Draw any figure object this app can produce onto the current device.
#' Everything goes through one code path -- open device, draw, close -- so the
#' TIFF and PDF buttons behave identically and no ggsave()/device-argument
#' version differences can creep in.
.draw_export_object <- function(obj) {
  if (is.null(obj)) stop("there is no figure to export yet")
  if (inherits(obj, "pheatmap")) {
    grid::grid.newpage(); grid::grid.draw(obj$gtable); return(invisible(TRUE))
  }
  if (inherits(obj, "recordedplot")) {
    grDevices::replayPlot(obj); return(invisible(TRUE))
  }
  if (is.function(obj)) { obj(); return(invisible(TRUE)) }
  # ggplot, ggraph (cnetplot / emapplot), patchwork -- all print()able
  if (inherits(obj, c("ggplot", "ggraph", "patchwork", "gg"))) {
    print(obj); return(invisible(TRUE))
  }
  if (inherits(obj, c("gtable", "gTree", "grob"))) {
    grid::grid.newpage(); grid::grid.draw(obj); return(invisible(TRUE))
  }
  print(obj)
  invisible(TRUE)
}

#' Save a ggplot, ggraph, patchwork, pheatmap, gtable or raw grob to TIFF,
#' PDF or PNG at publication quality.
#'
#' The figure is rendered to a private temp file carrying the correct extension
#' and only copied to the download target once it has been verified non-empty,
#' so a half-written file is never served to the browser.
export_plot <- function(plot_obj, file, format = c("tiff", "pdf", "png"),
                        width = 7, height = 6, dpi = 300) {
  format <- match.arg(format)
  width  <- .safe_dim(width, 7)
  height <- .safe_dim(height, 6)
  dpi    <- .safe_dpi(dpi)

  tmp <- tempfile(fileext = paste0(".", format))
  on.exit(unlink(tmp), add = TRUE)

  dev_before <- grDevices::dev.cur()
  dev_name   <- .open_export_device(tmp, format, width, height, dpi)

  drawn <- tryCatch({ .draw_export_object(plot_obj); TRUE },
                    error = function(e) e)

  # Always close, even if drawing blew up half way, or the device stays open
  # and the file on disk is left truncated.
  if (!identical(grDevices::dev.cur(), dev_before)) {
    try(grDevices::dev.off(), silent = TRUE)
  }

  if (inherits(drawn, "error")) {
    stop("The figure could not be rendered for export (", conditionMessage(drawn), ").")
  }
  if (!file.exists(tmp) || file.size(tmp) < 500) {
    stop("The ", toupper(format), " device (", dev_name, ") produced an empty file. ",
         "Try a smaller figure, or switch plot type.")
  }
  if (!file.copy(tmp, file, overwrite = TRUE)) {
    stop("The exported ", toupper(format), " could not be written to the download location.")
  }
  invisible(file)
}

#' Export wrapper for download handlers: on failure it still writes a valid
#' one-page file carrying the error text, so the browser always receives
#' something and the user learns why the figure did not come out, instead of
#' the button appearing to do nothing.
#'
#' @return NULL on success, or the error message (so the caller can also raise
#'   a Shiny notification).
export_plot_safe <- function(plot_obj, file, format = c("tiff", "pdf", "png"),
                             width = 7, height = 6, dpi = 300) {
  format <- match.arg(format)
  res <- tryCatch({ export_plot(plot_obj, file, format, width, height, dpi); NULL },
                  error = function(e) conditionMessage(e))
  if (is.null(res)) return(invisible(NULL))

  msg <- paste(strwrap(paste("Export failed:", res), width = 70), collapse = "\n")
  note <- tryCatch({
    dev_before <- grDevices::dev.cur()
    .open_export_device(file, format, .safe_dim(width, 7), .safe_dim(height, 6), .safe_dpi(dpi))
    grid::grid.newpage()
    grid::grid.text(msg, x = 0.5, y = 0.5,
                    gp = grid::gpar(fontsize = 11, col = "#b02a37", fontfamily = "sans"))
    if (!identical(grDevices::dev.cur(), dev_before)) try(grDevices::dev.off(), silent = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(note)) writeLines(msg, file)
  invisible(res)
}

#' Resolve pasted free text into gene ids, matching on symbol first then id.
#' Returns both the hits and whatever could not be matched, so the UI can tell
#' the user which of their genes were not found.
match_user_genes <- function(text, de) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(list(ids = character(0), missing = character(0)))
  }
  raw <- trimws(unlist(strsplit(text, "[,;\n\r\t ]+")))
  raw <- unique(raw[nzchar(raw)])
  if (!length(raw)) return(list(ids = character(0), missing = character(0)))

  by_name <- de$gene_id[match(toupper(raw), toupper(de$gene_name))]
  by_disp <- de$gene_id[match(toupper(raw), toupper(de$display_gene))]
  by_id   <- de$gene_id[match(toupper(raw), toupper(de$gene_id))]
  hit <- ifelse(!is.na(by_name), by_name, ifelse(!is.na(by_disp), by_disp, by_id))
  list(ids = unique(hit[!is.na(hit)]), missing = raw[is.na(hit)])
}

#' The gene-selection modes offered everywhere a gene set is needed.
GENE_SELECT_MODES <- c(
  "None"                                  = "none",
  "All significant DEGs"                  = "all_sig",
  "Top N by significance"                 = "top_n",
  "Top N by |log2 fold change|"           = "top_fc",
  "Up-regulated only (top N)"             = "top_up",
  "Down-regulated only (top N)"           = "top_down",
  "Select genes (searchable)"             = "dropdown",
  "Paste my own gene list"                = "user_list",
  "Genes from an enriched pathway"        = "pathway"
)

#' Central resolver: turn a selection mode plus its parameters into gene ids.
#'
#' @param mode     one of `GENE_SELECT_MODES`
#' @param de_all   canonical DE table (already classified)
#' @param top_n    N for the top-N modes
#' @param dropdown ids chosen in a selectize input
#' @param pasted   free text from a textarea
#' @param pathway  ids belonging to a chosen enrichment term
resolve_gene_selection <- function(mode, de_all, top_n = 15,
                                   dropdown = character(0), pasted = "",
                                   pathway = character(0)) {
  sig <- de_all[de_all$status %in% c("Up", "Down"), , drop = FALSE]
  switch(mode,
    none      = character(0),
    all_sig   = sig$gene_id,
    top_n     = head(sig$gene_id[order(sig$padj, -abs(sig$log2FoldChange))], top_n),
    top_fc    = head(sig$gene_id[order(-abs(sig$log2FoldChange))], top_n),
    top_up    = {
      u <- sig[sig$status == "Up", , drop = FALSE]
      head(u$gene_id[order(u$padj, -u$log2FoldChange)], top_n)
    },
    top_down  = {
      d <- sig[sig$status == "Down", , drop = FALSE]
      head(d$gene_id[order(d$padj, d$log2FoldChange)], top_n)
    },
    dropdown  = as.character(dropdown),
    user_list = match_user_genes(pasted, de_all)$ids,
    pathway   = as.character(pathway),
    character(0)
  )
}

#' Tidy DE table for on-screen display / CSV export.
deg_display_table <- function(df, digits = 3) {
  out <- data.frame(
    Gene            = df$display_gene,
    `Ensembl ID`    = df$gene_id,
    `Mean expr.`    = round(df$baseMean, 1),
    `log2 FC`       = round(df$log2FoldChange, digits),
    `P value`       = ifelse(is.na(df$pvalue), NA_character_,
                             formatC(df$pvalue, format = "e", digits = 2)),
    `adj. P`        = ifelse(is.na(df$padj), NA_character_,
                             formatC(df$padj, format = "e", digits = 2)),
    Status          = as.character(df$status),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  out <- out[order(df$padj, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Timestamped filename helper.
stamp_name <- function(prefix, ext) {
  sprintf("%s_%s.%s", prefix, format(Sys.time(), "%Y%m%d_%H%M%S"), ext)
}
