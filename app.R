# =============================================================================
# DEGexplorer -- an interactive workbench for RNA-seq differential expression
#
#   Volcano | MA | expression scatter | PCA + scree | heatmap |
#   GO/KEGG/MSigDB enrichment | enrichment gallery | gene intelligence
#
# Run with:   shiny::runApp()      (from this directory)
#        or:  R -e 'shiny::runApp("path/to/DEGexplorer", launch.browser = TRUE)'
#
# Licence: MIT. See LICENSE.
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(DT)
  library(readr)
  library(RColorBrewer)
  library(matrixStats)
  library(grid)
})

# Optional at load time; the app reports precisely what is missing where used.
for (p in c("shinyjs", "colourpicker", "shinythemes", "openxlsx", "gridExtra")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    warning(sprintf("Package '%s' is not installed; some controls will be unavailable. Run install.packages('%s').", p, p))
  }
}

# ---- locate files relative to the app directory ----------------------------
# shiny::runApp() evaluates app.R in its own environment while the modules
# below are sourced into the global environment, so these two helpers are
# assigned globally to keep them visible to every module.
local({
  d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
  if (is.null(d) || !nzchar(d)) d <- getwd()
  assign("APP_DIR", normalizePath(d, mustWork = FALSE), envir = globalenv())
})

#' Resolve a path inside the app directory, falling back to the working
#' directory so the app runs from either location.
assign("app_file", function(...) {
  p <- file.path(get("APP_DIR", envir = globalenv()), ...)
  if (file.exists(p)) return(p)
  p2 <- file.path(getwd(), ...)
  if (file.exists(p2)) return(p2)
  p
}, envir = globalenv())

# ---- load modules -----------------------------------------------------------
for (f in c("utils_io.R", "utils_normalise.R", "utils_theme.R", "utils_export.R",
            "plots_de.R", "plots_qc.R", "enrichment.R", "enrichment_tables.R",
            "plots_enrichment.R", "gene_report.R",
            "sheet_editor.R", "organisms.R", "annotation_sources.R", "transcripts.R",
            "de_analysis.R", "dtu.R", "de_qc.R", "de_report.R", "analysis_ui.R", "analysis_server.R",
            "gene_explorer.R",
            "ui.R", "server.R")) {
  source(app_file("R", f), local = FALSE, encoding = "UTF-8")
}

shinyApp(ui = app_ui(), server = app_server)
