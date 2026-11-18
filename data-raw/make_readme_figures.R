# =============================================================================
# make_readme_figures.R -- redraw the gallery figures shown in the README
# (man/figures/gallery_*.png) from the bundled example enrichment tables, with the
# app's own gallery code at its default sizes (80 dpi). Run from the repository
# root after data-raw/make_example_enrichment.R:
#
#   Rscript data-raw/make_readme_figures.R
# =============================================================================
src <- readLines("app.R"); src <- src[!grepl("^shinyApp\\(", src)]
suppressMessages(eval(parse(text = src), envir = globalenv()))
de <- standardise_de_table(read_table_auto("data/example_deseq2_results.tsv.gz"))$data
tabs <- import_enrichment_files(example_enrichment_files("data/example_enrichment"), de = de)$tables
fc <- stats::setNames(de$log2FoldChange, de$display_gene)
fc <- fc[!is.na(fc) & !duplicated(names(fc))]
figs <- list(score_bar = "BP", score_dot = "BP", go3_h = NULL, go3_v = NULL, sankey = "BP",
             kegg_class = NULL, cnet = "BP", emap = "BP")
# The KEGG demo table is built locally only (KEGG data are not redistributed);
# without it the KEGG figure is left as it is.
if (!"KEGG" %in% names(tabs)) figs$kegg_class <- NULL
for (f in names(figs)) {
  obj <- build_gallery_figure(f, tabs, figs[[f]], fc = fc)
  sz <- gallery_figure_size(f)
  out <- file.path("man", "figures", paste0("gallery_", f, ".png"))
  export_plot(obj, out, "png", sz[1], sz[2], 80)
  cat("wrote", out, "\n")
}
