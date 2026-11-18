# ---------------------------------------------------------------------------
# install_dependencies.R -- one-shot setup for DEGexplorer
#
#   Rscript install_dependencies.R            # required + recommended
#   Rscript install_dependencies.R --all      # also enrichment + live lookups
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
full <- "--all" %in% args

required <- c("shiny", "ggplot2", "ggrepel", "pheatmap", "DT", "readr",
              "RColorBrewer", "matrixStats")
recommended <- c("shinyjs", "colourpicker", "shinythemes", "openxlsx",
                 "gridExtra", "jsonlite", "curl", "viridisLite", "zip", "png")
# DESeq2, limma, edgeR and EBSeq are the Analysis-tab engines (limma-voom also uses edgeR)
bioc <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "pathview",
          "DESeq2", "limma", "edgeR", "EBSeq", "apeglm",
          "stageR",   # stage-wise testing for transcript-level analysis
          "DEXSeq", "DRIMSeq", "satuRn")   # transcript usage (DTU); edgeR / limma diffSplice need nothing extra
# Other organisms: install that organism's OrgDb for GO / KEGG enrichment, e.g.
# BiocManager::install("org.Mm.eg.db")   # mouse; org.Rn.eg.db rat, org.Dr.eg.db zebrafish, ...
cran_extra <- c("msigdbr", "ggarchery", "ashr")   # ashr: DESeq2 fold-change shrinkage

install_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(miss)) { message("  all present"); return(invisible()) }
  message("  installing: ", paste(miss, collapse = ", "))
  install.packages(miss, repos = "https://cloud.r-project.org")
}

message("Required packages:");    install_missing(required)
message("Recommended packages:"); install_missing(recommended)

if (full) {
  message("Optional CRAN packages:"); install_missing(cran_extra)
  message("Bioconductor packages:")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  miss <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) BiocManager::install(miss, update = FALSE, ask = FALSE)
  else message("  all present")
} else {
  message("\nEnrichment and live-lookup packages were skipped.")
  message("Run  Rscript install_dependencies.R --all  to add them.")
}

message("\nDone. Start the app with:  shiny::runApp()")
