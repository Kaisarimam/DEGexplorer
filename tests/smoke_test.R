# =============================================================================
# smoke_test.R -- end-to-end checks against the bundled example data.
#
#   Rscript tests/smoke_test.R        (from the repository root)
#
# Exits non-zero on the first failure so CI catches regressions.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(pheatmap)
  library(readr); library(RColorBrewer); library(matrixStats); library(grid)
})

root <- if (file.exists("app.R")) "." else ".."
for (f in c("utils_io.R", "utils_normalise.R", "utils_theme.R", "utils_export.R", "organisms.R",
            "plots_de.R", "plots_qc.R", "enrichment.R", "enrichment_tables.R",
            "plots_enrichment.R", "gene_report.R", "de_analysis.R",
            "annotation_sources.R", "transcripts.R", "dtu.R", "de_qc.R", "de_report.R", "gene_explorer.R")) {
  source(file.path(root, "R", f))
}

ok <- 0L
check <- function(label, expr) {
  res <- tryCatch(isTRUE(expr), error = function(e) {
    cat("FAIL ", label, " -- ", conditionMessage(e), "\n", sep = ""); FALSE
  })
  if (res) { cat("ok   ", label, "\n", sep = ""); ok <<- ok + 1L }
  else { cat("FAIL ", label, "\n", sep = ""); quit(status = 1L) }
}

d <- function(f) file.path(root, "data", f)

# ---- 1. readers and format detection ---------------------------------------

deseq_raw <- read_table_auto(d("example_deseq2_results.tsv.gz"))
limma_raw <- read_table_auto(d("example_limma_results.tsv.gz"))

check("DESeq2 file is recognised", guess_de_pipeline(deseq_raw) == "DESeq2")
check("limma file is recognised", grepl("limma", guess_de_pipeline(limma_raw)))

# result-table layouts (annotated, plain, DEGs only)
plain_raw <- read_table_auto(d("example_limma_results_plain.tsv.gz"))
filt_raw  <- read_table_auto(d("example_limma_results_filtered.tsv.gz"))

i_ann  <- describe_de_table(limma_raw, "x.limma.annotated.tsv")
i_res  <- describe_de_table(plain_raw, "x.limma.results.tsv")
i_filt <- describe_de_table(filt_raw,  "x.limma.results_filtered.tsv")
i_ds2  <- describe_de_table(deseq_raw, "x.deseq2.annotated.tsv")

check("limma annotated layout recognised",
      i_ann$nfcore && i_ann$pipeline == "limma" && i_ann$flavour == "annotated")
check("limma plain results layout recognised",
      i_res$nfcore && i_res$flavour == "results" && !i_res$filtered)
check("limma DEGs-only layout recognised",
      i_filt$nfcore && i_filt$flavour == "results_filtered" && i_filt$filtered)
check("DESeq2 layout recognised", i_ds2$nfcore && i_ds2$pipeline == "DESeq2")
check("filtered table raises a warning note", length(i_filt$notes) > 0)
check("symbol-less table is flagged",
      any(grepl("No gene symbols", i_res$notes)))
check("all three limma flavours share the statistics columns",
      identical(detect_de_columns(plain_raw)$log2FoldChange, "logFC") &&
      identical(detect_de_columns(filt_raw)$padj, "adj.P.Val"))

# annotation back-fill
ann_map <- read_gene_annotation(d("example_limma_results.tsv.gz"))
check("annotation map built", nrow(ann_map) > 20000 &&
      all(c("gene_id", "gene_name") %in% colnames(ann_map)))
plain_std <- standardise_de_table(plain_raw)
check("plain results have no symbols",
      all(plain_std$data$display_gene == plain_std$data$gene_id))
filled <- apply_gene_annotation(plain_std$data, ann_map)
check("symbols filled in from the annotation file", filled$n_filled > 10000)
check("filled table now shows symbols",
      sum(filled$data$display_gene != filled$data$gene_id) == filled$n_filled)
check("annotation never overwrites existing symbols",
      apply_gene_annotation(standardise_de_table(limma_raw)$data, ann_map)$n_filled == 0)

# enrichment background
filt_std <- standardise_de_table(filt_raw)$data
check("filtered table alone gives a tiny universe",
      length(enrichment_universe(filt_std)) < 500)
check("unfiltered file restores the universe",
      length(enrichment_universe(filt_std, plain_raw$gene_id)) > 20000)
check("count matrix is used as a fallback background",
      length(enrichment_universe(filt_std, NULL, paste0("G", 1:9000))) == 9000)

de_d <- standardise_de_table(deseq_raw)
de_l <- standardise_de_table(limma_raw)

check("DESeq2 maps to the canonical schema",
      all(c("gene_id", "log2FoldChange", "pvalue", "padj", "baseMean", "display_gene")
          %in% colnames(de_d$data)))
check("limma maps logFC -> log2FoldChange", de_l$mapping$log2FoldChange == "logFC")
check("limma maps adj.P.Val -> padj", de_l$mapping$padj == "adj.P.Val")
check("limma AveExpr is linearised", any(grepl("log-scaled", de_l$notes)))
check("both tables have the same gene count", nrow(de_d$data) == nrow(de_l$data))
check("gene ids are unique", !anyDuplicated(de_d$data$gene_id))
check("display_gene never empty", all(nzchar(de_d$data$display_gene)))

# ---- 2. matrices ------------------------------------------------------------

cnt <- read_matrix_auto(d("example_rsem_gene_counts.tsv.gz"))
len <- read_matrix_auto(d("example_rsem_gene_lengths.tsv.gz"))
meta <- read_sample_metadata(d("example_samplesheet.csv"))

check("count matrix parsed", is.matrix(cnt$matrix) && ncol(cnt$matrix) == 6)
check("gene_name column captured separately", !is.null(cnt$annotation))
check("length matrix aligns to counts",
      length(intersect(rownames(cnt$matrix), rownames(len$matrix))) > 20000)
check("sample sheet has a condition column", "condition" %in% colnames(meta))
check("sample ids match matrix columns", all(meta$sample %in% colnames(cnt$matrix)))

# ---- 3. normalisation -------------------------------------------------------

tpm  <- normalise_counts(cnt$matrix, len$matrix, "tpm")
cpm  <- normalise_counts(cnt$matrix, NULL, "cpm")
fpkm <- normalise_counts(cnt$matrix, len$matrix, "fpkm")
ds2  <- normalise_counts(cnt$matrix, NULL, "deseq2")
tmm  <- normalise_counts(cnt$matrix, NULL, "tmm")

check("TPM sums to 1e6 per sample", all(abs(colSums(tpm) - 1e6) < 1))
check("CPM sums to 1e6 per sample", all(abs(colSums(cpm) - 1e6) < 1))
check("FPKM is finite and non-negative", all(is.finite(fpkm)) && min(fpkm) >= 0)
check("DESeq2 size factors are positive", all(deseq2_size_factors(cnt$matrix) > 0))
check("TMM factors multiply to 1", abs(prod(tmm_norm_factors(cnt$matrix)) - 1) < 1e-8)
check("DESeq2 normalisation preserves dimensions", all(dim(ds2) == dim(cnt$matrix)))
check("TPM without lengths is refused",
      inherits(tryCatch(normalise_counts(cnt$matrix, NULL, "tpm"), error = function(e) e), "error"))

# ---- 4. DE classification ---------------------------------------------------

cls <- classify_degs(de_d$data, 0.05, 1)
check("classification yields Up/Down/NS", all(levels(cls$status) == c("Up", "Down", "NS")))
check("significant genes found", sum(cls$status != "NS") > 0)
check("stricter cutoff yields no more genes",
      sum(classify_degs(de_d$data, 0.01, 2)$status != "NS") <= sum(cls$status != "NS"))

sig <- significant_degs(cls)
check("significant genes sorted by adjusted p", !is.unsorted(sig$padj, na.rm = TRUE))

# ---- 5. gene selection ------------------------------------------------------

check("top-N selection honours N", length(resolve_gene_selection("top_n", cls, 10)) == 10)
check("up-only selection is all Up",
      all(cls$status[match(resolve_gene_selection("top_up", cls, 5), cls$gene_id)] == "Up"))
check("down-only selection is all Down",
      all(cls$status[match(resolve_gene_selection("top_down", cls, 5), cls$gene_id)] == "Down"))
m <- match_user_genes("IFIT1, MX1, NOT_A_REAL_GENE", de_d$data)
check("pasted symbols matched", length(m$ids) == 2)
check("unmatched genes reported", identical(m$missing, "NOT_A_REAL_GENE"))

# ---- 6. plots ---------------------------------------------------------------

labs <- head(sig$gene_id, 8)
check("volcano builds", inherits(make_volcano(cls, label_ids = labs), "ggplot"))
check("MA builds", inherits(make_ma(cls, label_ids = labs), "ggplot"))

sdf <- group_mean_expression(fpkm, meta, "mock", "stim")
check("group means computed for every gene", nrow(sdf) == nrow(fpkm))
check("scatter builds",
      inherits(make_expression_scatter(sdf, cls, label_ids = labs), "ggplot"))

pca <- run_pca(log2(cpm + 1), ntop = 500)
check("PCA returns scores for every sample", nrow(pca$scores) == ncol(cpm))
check("PCA variance sums to ~100%", abs(sum(pca$percent) - 100) < 1e-6)
check("PCA plot builds", inherits(make_pca_plot(pca, meta), "ggplot"))
check("scree plot builds", inherits(make_scree_plot(pca), "ggplot"))

# ---- 7. heatmap, including the reference-left ordering -----------------------

vst <- log2(cpm + 1)
hm <- make_heatmap(head(sig$gene_id, 30), vst, meta, de_d$data,
                   reference = "mock", cluster_cols = FALSE)
cond <- meta$condition[match(hm$samples, meta$sample)]
n_ctrl <- sum(cond == "mock")
check("heatmap plots the requested genes", hm$n_genes == 30)
check("reference (mock) samples occupy the leftmost columns",
      all(cond[seq_len(n_ctrl)] == "mock"))
check("stim samples occupy the rightmost columns",
      all(cond[(n_ctrl + 1):length(cond)] == "stim"))

hm_all <- make_heatmap(NULL, vst, meta, de_d$data, all_genes = TRUE,
                       reference = "mock", max_genes = 500)
check("all-genes heatmap respects the row cap", hm_all$n_genes == 500)
check("all-genes heatmap reports truncation", isTRUE(hm_all$truncated))
check("all-genes mode includes non-significant genes",
      hm_all$n_requested > sum(cls$status != "NS"))

# ---- 8. export --------------------------------------------------------------

tmp_t <- tempfile(fileext = ".tiff"); tmp_p <- tempfile(fileext = ".pdf")
export_plot(make_volcano(cls), tmp_t, "tiff", 6, 5, 150)
export_plot(hm$plot, tmp_p, "pdf", 6, 5)
check("TIFF written", file.exists(tmp_t) && file.size(tmp_t) > 5000)
check("PDF written from a pheatmap", file.exists(tmp_p) && file.size(tmp_p) > 2000)

# Shiny hands downloadHandler a temp path with NO extension, and the enrichment
# tab can produce a stacked gtable rather than a ggplot -- both must export to
# PDF as readily as to TIFF.
is_pdf <- function(f) identical(readBin(f, "raw", 4L), as.raw(c(0x25, 0x50, 0x44, 0x46)))
stacked <- gridExtra::arrangeGrob(make_volcano(cls), make_ma(cls), ncol = 1)
for (nm in c("ggplot", "gtable", "pheatmap")) {
  obj <- switch(nm, ggplot = make_volcano(cls), gtable = stacked, pheatmap = hm$plot)
  fp <- tempfile(); ft <- tempfile()            # deliberately extension-less
  export_plot(obj, fp, "pdf", 8, 7, 300)
  export_plot(obj, ft, "tiff", 8, 7, 150)
  check(paste0("PDF export from a ", nm, " is a real PDF"),
        file.size(fp) > 2000 && is_pdf(fp))
  check(paste0("TIFF export from a ", nm), file.size(ft) > 5000)
  unlink(c(fp, ft))
}

fe <- tempfile()
err <- export_plot_safe(NULL, fe, "pdf", 6, 5, 300)
check("a failed export still returns a readable PDF",
      is.character(err) && file.size(fe) > 500 && is_pdf(fe))
unlink(fe)

tmp_csv <- tempfile(fileext = ".csv")
write_matrix_file(head(tpm, 100), tmp_csv, "csv", annotation = cnt$annotation)
back <- read_table_auto(tmp_csv)
check("matrix round-trips through CSV", nrow(back) == 100)
unlink(c(tmp_t, tmp_p, tmp_csv))

# ---- 8b. GTF / GFF3 annotation ---------------------------------------------

# Small GENCODE-style fixtures: versioned ids, gene + transcript + exon rows,
# header comments; plus an Ensembl-style GFF3 and a transcript-only GTF.
gtf_src <- utils::head(limma_raw[!is.na(limma_raw$gene_name) & nzchar(limma_raw$gene_name), ], 120)
gtf_lines <- c("#!genome-build GRCh38")
for (i in seq_len(nrow(gtf_src))) {
  at <- sprintf(paste0('gene_id "%s.%d"; gene_version "%d"; gene_name "%s"; ',
                       'gene_source "havana"; gene_biotype "%s"; level 2;'),
                gtf_src$gene_id[i], i %% 9 + 1, i %% 9 + 1, gtf_src$gene_name[i],
                if (i %% 3 == 0) "lncRNA" else "protein_coding")
  row <- function(ft, a) paste(gtf_src$chromosome[i], "HAVANA", ft, 1000 * i,
                               1000 * i + 800, ".", "+", ".", a, sep = "\t")
  gtf_lines <- c(gtf_lines, row("gene", at),
                 row("transcript", paste0(at, sprintf(' transcript_id "ENST%08d.1";', i))),
                 row("exon", paste0(at, sprintf(' transcript_id "ENST%08d.1"; exon_number 1;', i))))
}
f_gtf <- tempfile(fileext = ".gtf"); writeLines(gtf_lines, f_gtf)
f_gz  <- tempfile(fileext = ".gtf.gz")
local({ cn <- gzfile(f_gz, "wt"); writeLines(gtf_lines, cn); close(cn) })
f_tx  <- tempfile(fileext = ".gtf")
writeLines(grep("\ttranscript\t|^#", gtf_lines, value = TRUE), f_tx)
f_gff <- tempfile(fileext = ".gff3")
writeLines(c("##gff-version 3", vapply(seq_len(60), function(i)
  paste(gtf_src$chromosome[i], "ensembl", "gene", 1000 * i, 1000 * i + 800, ".", "-", ".",
        sprintf("ID=gene:%s;Name=%s;biotype=protein_coding;description=test%%20gene",
                gtf_src$gene_id[i], gtf_src$gene_name[i]), sep = "\t"),
  character(1))), f_gff)

check("plain GTF recognised by content", is_gtf_file(f_gtf))
check("gzipped GTF recognised", is_gtf_file(f_gz))
check("GFF3 recognised", is_gtf_file(f_gff))
check("a TSV is not mistaken for a GTF", !is_gtf_file(d("example_limma_results.tsv.gz")))

g <- read_gtf_annotation(f_gtf)
check("GTF keeps one row per gene", nrow(g) == nrow(gtf_src))
check("GTF yields gene ids and symbols",
      all(c("gene_id", "gene_name") %in% colnames(g)) && all(!is.na(g$gene_name)))
check("GTF yields the locus", all(c("chromosome", "start", "end", "strand") %in% colnames(g)))
check("GTF yields extra attributes",
      all(c("gene_biotype", "gene_version") %in% colnames(g)))
check("transcript and exon rows did not become genes", !any(duplicated(g$gene_id)))
check("GTF source recorded", identical(attr(g, "source"), "GTF"))
check("gzipped GTF parses identically", identical(read_gtf_annotation(f_gz), g))
check("transcript-only GTF still yields genes",
      nrow(read_gtf_annotation(f_tx)) == nrow(gtf_src))

gf <- read_gtf_annotation(f_gff)
check("GFF3 strips the 'gene:' id prefix", all(grepl("^ENSG[0-9]+$", gf$gene_id)))
check("GFF3 reads Name as the symbol", all(!is.na(gf$gene_name)))
check("GFF3 percent-encoding decoded", any(grepl(" ", gf$description, fixed = TRUE)))
check("GFF3 source recorded", identical(attr(gf, "source"), "GFF3"))

check("read_gene_annotation dispatches on content",
      identical(read_gene_annotation(f_gtf), g))

# version-suffix tolerance: GENCODE ids vs unversioned pipeline ids
check("GTF ids carry versions, DE ids do not",
      all(grepl("\\.[0-9]+$", g$gene_id)) && !any(grepl("\\.", plain_std$data$gene_id)))
gfill <- apply_gene_annotation(plain_std$data, g)
check("versioned GTF ids still match unversioned DE ids", gfill$n_filled == nrow(gtf_src))
check("GTF attributes land on the DE table",
      all(c("chromosome", "strand", "gene_biotype") %in% colnames(gfill$data)))
check("a GTF alone removes the need for an id/symbol table",
      sum(gfill$data$display_gene != gfill$data$gene_id) == nrow(gtf_src))

# ---- 9. enrichment results table -------------------------------------------

set.seed(42)
fake_enrich <- local({
  n <- 6
  gl <- replicate(n, paste(sample(de_d$data$gene_id, sample(c(2, 5, 200), 1)), collapse = "/"))
  data.frame(ID = paste0("GO:", 1:n), Description = paste("term", 1:n),
             GeneRatio = paste0(6:(5 + n), "/1420"), RichFactor = runif(n),
             pvalue = runif(n, 1e-7, 1e-3), p.adjust = runif(n, 1e-5, 0.01),
             geneID = gl, Count = sample(2:200, n, TRUE), stringsAsFactors = FALSE)
})

etb <- enrich_result_table(fake_enrich, de_d$data, "Up")
check("enrichment table gains a gene symbol column", "geneName" %in% colnames(etb))
check("enrichment table keeps the Ensembl id column", "geneID" %in% colnames(etb))
check("gene columns sit at the end next to Count",
      tail(colnames(etb), 4)[1] == "Count" && tail(colnames(etb), 1) == "Direction")
check("symbols resolved where the DE table knows them",
      any(!grepl("^ENSG", unlist(strsplit(etb$geneName, "/")))))
check("unknown ids fall back to the Ensembl accession",
      all(lengths(strsplit(etb$geneName, "/")) == lengths(strsplit(etb$geneID, "/"))))
check("export table carries no HTML",
      !any(grepl("<", c(etb$geneName, etb$geneID), fixed = TRUE)))

fmt <- format_enrich_numbers(etb)
check("wide numerics rounded for display", all(nchar(format(fmt$RichFactor)) <= 6))
check("p-values shown in scientific notation", all(grepl("e[+-]", fmt$p.adjust)))

cells <- collapsible_gene_cell(etb$geneName)
long  <- lengths(strsplit(etb$geneName, "/")) > 3
check("long gene lists collapse into <details>", all(grepl("^<details", cells[long])))
check("short gene lists stay plain", all(grepl("^<span", cells[!long])))
check("collapsed cells still contain every gene",
      all(mapply(function(cell, src) {
        all(vapply(strsplit(src, "/")[[1]], function(g) grepl(g, cell, fixed = TRUE), logical(1)))
      }, cells, etb$geneName)))
check("the toggle reports how many are hidden", all(grepl("more</span>", cells[long])))
check("friendly headers applied",
      identical(enrich_table_colnames(etb)[match(c("geneName", "geneID"), colnames(etb))],
                c("Genes", "Ensembl IDs")))

# KEGG is the one database clusterProfiler runs on Entrez ids: those native
# accessions get their own column, and no other database gains an empty one.
kegg_fake <- local({
  ens <- fake_enrich$geneID
  ent <- vapply(strsplit(ens, "/"), function(x)
    paste(seq_along(x) + 1000L, collapse = "/"), character(1))
  cbind(fake_enrich[, setdiff(colnames(fake_enrich), "geneID")],
        geneID = ent, geneID_ensembl = ens, stringsAsFactors = FALSE)
})
ktb <- enrich_result_table(kegg_fake, de_d$data, "Up")
check("KEGG keeps its Entrez ids in their own column", "entrezID" %in% colnames(ktb))
check("Entrez column header is labelled",
      enrich_table_colnames(ktb)[match("entrezID", colnames(ktb))] == "Entrez IDs")
check("KEGG geneID column still holds Ensembl accessions",
      all(grepl("^ENS", vapply(strsplit(ktb$geneID, "/"), `[`, character(1), 1))))
check("Entrez ids align one-to-one with the symbols",
      all(lengths(strsplit(ktb$entrezID, "/")) == lengths(strsplit(ktb$geneName, "/"))))
check("Entrez ids survive untouched",
      identical(ktb$entrezID, kegg_fake$geneID))
check("non-KEGG results gain no empty Entrez column", !"entrezID" %in% colnames(etb))
check("every gene column is collapsible",
      all(vapply(intersect(ENRICH_GENE_COLS, colnames(ktb)),
                 function(cl) all(grepl("^<(details|span)", collapsible_gene_cell(ktb[[cl]]))),
                 logical(1))))
check("Entrez column exports as plain text",
      !any(grepl("<", ktb$entrezID, fixed = TRUE)))

# ---- 10. gene report --------------------------------------------------------

L <- gene_resource_links("IFIT1", "ENSG00000185745")
check("resource links generated", length(L) >= 12 && all(grepl("^https://", unlist(L))))
check("literature query includes context",
      grepl("T cell", literature_query("IFIT1", "T cell, interferon")))
tbl <- gene_link_table(head(sig$gene_id, 5), de_d$data)
check("link table built", nrow(tbl) == 5 && "Resources" %in% colnames(tbl))
html <- build_gene_report(head(sig$gene_id, 5), de_d$data)
check("HTML report built", nchar(html) > 5000 && grepl("</html>", html))
check("report contains gene symbols", grepl(sig$display_gene[1], html, fixed = TRUE))

# ---- 10b. bundled GO release ------------------------------------------------

check("GO release databases are offered", all(GO_RELEASE_DATABASES %in% ENRICH_DATABASES))
check("GO release version is read", grepl("^\\d{4}-\\d{2}-\\d{2}$", go_release_version()))
gr <- go_release_term2gene("GOREL_BP", stats::setNames(de_d$data$display_gene, de_d$data$gene_id))
check("GO release gene sets are keyed by the table's Ensembl ids",
      nrow(gr$t2g) > 10000 && all(grepl("^ENSG", gr$t2g$gene)))
check("GO release sets carry term names", "GO:0060337" %in% gr$t2n$term)
check("IFIT1 is not in GO:0060337 in this release",
      !any(de_d$data$display_gene[match(gr$t2g$gene[gr$t2g$term == "GO:0060337"],
                                        de_d$data$gene_id)] %in% "IFIT1"))
check("GO release needs no org.Hs.eg.db",
      !length(enrichment_missing_packages("GOREL_BP")) || !HAS_CLUSTERPROFILER)

# ---- 11. enrichment gallery -------------------------------------------------

demo_f <- example_enrichment_files(d("example_enrichment"))
check("three GO demo tables are shipped and no KEGG data", length(demo_f) == 3 &&
        !any(grepl("pathway", basename(demo_f))))
gimp <- import_enrichment_files(demo_f, de = de_d$data)
gt   <- gimp$tables
check("demo tables import as GO BP/CC/MF from their file names",
      identical(names(gt), c("BP", "CC", "MF")) && all(vapply(gt, nrow, 1L) > 50))
check("demo tables are in the standard nine-column layout",
      all(vapply(demo_f, function(f) identical(colnames(utils::read.csv(f, check.names = FALSE)),
                                               GALLERY_STD_COLS), logical(1))))
demo_genes <- unique(unlist(strsplit(unlist(lapply(gt, `[[`, "geneID")), "/", fixed = TRUE)))
demo_sig <- significant_degs(cls)
check("every gene in the demo tables is a DEG of the bundled demo data",
      length(demo_genes) > 50 && all(demo_genes %in% demo_sig$display_gene))
check("significance defaults to adjusted p, and p-value is the looser option",
      all(gallery_sig(gt$BP)$p.adjust <= 0.05) &&
      nrow(gallery_sig(gt$BP, "pvalue")) > nrow(gallery_sig(gt$BP)))
check("enrichment score is -log10(p)",
      isTRUE(all.equal(gt$BP$Score[1], -log10(gt$BP$pvalue[1]))))
check("gene ratios are evaluated without eval()",
      isTRUE(all.equal(ratio_to_numeric(c("8/90", "0.25", "x")), c(8 / 90, 0.25, NA))))
gal_fc <- stats::setNames(de_d$data$log2FoldChange, de_d$data$display_gene)

# importing tables written by other tools
alien <- data.frame(Category = c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "KEGG_PATHWAY"),
                    Term = c("GO:0006955~immune response", "GO:0005886~plasma membrane",
                             "hsa04110:Cell cycle"),
                    PValue = c(1e-4, 2e-3, 0.03), Genes = c("A, B, C", "A, D", "E, F"),
                    FDR = c(0.01, 0.04, 0.2), stringsAsFactors = FALSE)
al <- split_enrichment_by_database(standardise_enrich_table(alien))
check("a DAVID-style table is split into its databases",
      setequal(names(al), c("BP", "CC", "KEGG")) && al$BP$ID == "GO:0006955" &&
      al$BP$Description == "immune response" && al$BP$geneID == "A/B/C" && al$BP$Count == 3 &&
      al$KEGG$ID == "hsa04110")
tmpx <- file.path(tempdir(), "my.go.MF.all.csv")
utils::write.csv(gallery_export_table(head(gt$MF, 40)), tmpx, row.names = FALSE)
imp <- import_enrichment_files(tmpx)
check("database is recognised from the file name on import",
      identical(names(imp$tables), "MF") && nrow(imp$tables$MF) == 40)
check("an unidentifiable table is reported, not guessed",
      { f2 <- file.path(tempdir(), "results.csv"); file.copy(tmpx, f2, overwrite = TRUE)
        r <- import_enrichment_files(f2); !length(r$tables) && any(grepl("could not tell", r$notes)) })
check("the sidebar hint settles an unidentifiable table",
      identical(names(import_enrichment_files(file.path(tempdir(), "results.csv"), hint = "CC")$tables), "CC"))
unlink(c(tmpx, file.path(tempdir(), "results.csv")))

check("KEGG classes are derived from pathway ids",
      identical(kegg_pathway_class(c("hsa00010", "hsa03010", "hsa04514", "hsa04110", "hsa04713",
                                     "hsa03320", "hsa05217", "hsa04141", "Cell cycle")),
                c("Metabolism", "Genetic Information Processing",
                  "Environmental Information Processing", "Cellular Processes",
                  "Organismal Systems", "Organismal Systems", "Human Diseases",
                  "Genetic Information Processing", "Unclassified")))

# every table-driven figure builds and renders. No KEGG table is shipped (KEGG data
# may not be redistributed), so the KEGG-class figure gets a stand-in table: real
# pathway-id formats (which carry the class), placeholder names, genes from GO BP.
if (is.null(gt$KEGG)) {
  kk <- head(gt$BP, 12)
  kk$ID <- c("hsa00010", "hsa03010", "hsa04514", "hsa04110", "hsa04713", "hsa05217",
             "hsa04141", "hsa03320", "hsa04060", "hsa04620", "hsa05164", "hsa04630")[seq_len(nrow(kk))]
  kk$Description <- paste("Test pathway", seq_len(nrow(kk)))
  gt$KEGG <- kk
}
gal_figs <- c(score_bar = "BP", score_dot = "CC", go3_h = NA, go3_v = NA, sankey = "BP",
              kegg_class = NA, cnet = "MF", emap = "BP")
for (fg in names(gal_figs)) {
  db <- if (is.na(gal_figs[[fg]])) NULL else gal_figs[[fg]]
  fp <- tempfile()
  r <- tryCatch({
    p <- build_gallery_figure(fg, gt, db, list(show_values = TRUE), fc = gal_fc)
    export_plot(p, fp, "png", 8, 6, 72); inherits(p, "ggplot") && file.size(fp) > 5000
  }, error = function(e) { cat("     ", conditionMessage(e), "\n"); FALSE })
  check(paste0("gallery figure renders: ", fg), r)
  unlink(fp)
}
check("PNG export is a real PNG",
      { fp <- tempfile(); export_plot(make_volcano(cls), fp, "png", 5, 4, 72)
        identical(readBin(fp, "raw", 4L), as.raw(c(0x89, 0x50, 0x4e, 0x47))) })
check("three-ontology frame keeps BP, CC, MF order with n terms each",
      { g3 <- go_three_df(gt, 5); identical(as.integer(table(g3$Ontology)), c(5L, 5L, 5L)) &&
        identical(levels(g3$Ontology), c("BP", "CC", "MF")) })
check("Sankey caps the gene count and says so",
      { p <- make_sankey_dot(gt$BP, top_n = 30, max_genes = 12)
        attr(p, "n_genes") == 12 && attr(p, "n_dropped") > 0 })
check("a cutoff nothing passes gives a clear message",
      grepl("No terms pass", tryCatch(make_score_barplot(gt$KEGG, cutoff = 1e-30),
                                      error = function(e) conditionMessage(e))))
check("figures cope with a single term",
      inherits(make_native_emap(gt$KEGG, top_n = 1), "ggplot") &&
      inherits(make_sankey_dot(gt$KEGG, top_n = 1), "ggplot") &&
      inherits(make_circular_cnet(gt$KEGG, top_n = 1), "ggplot"))
check("missing optional packages are reported by name",
      HAS_ENRICHPLOT || grepl("enrichplot", tryCatch(make_go_graph(NULL, gt$BP, "BP"),
                                                     error = function(e) conditionMessage(e))))
check("pathway maps report what they need",
      HAS_PATHVIEW || grepl("pathview", tryCatch(render_pathview_map("hsa05160", gal_fc),
                                                 error = function(e) conditionMessage(e))))

# the one-click bundle
bdir <- tempfile("bundle_"); dir.create(bdir)
blog <- write_gallery_bundle(gt[c("CC", "KEGG")], bdir, list(top_n = 5), fc = gal_fc,
                             formats = "png", dpi = 50)
bfiles <- list.files(bdir, recursive = TRUE)
check("bundle is laid out as go/ and pathway/ with tables and figures",
      all(c("go/CC_Enrichment_Score_barplot.png", "go/CC_cnetplot.png", "go/GO_Three_Ontologies.png",
            "pathway/Pathway_sankey_dot.png", "pathway/Pathway_category_summary.png",
            "bundle_log.txt") %in% bfiles) &&
      sum(grepl("^(go/go\\.CC|pathway/pathway)\\.(all|sig)\\.(xlsx|csv)$", bfiles)) == 4)
check("bundle log records no failures", !any(grepl("^(FAILED|SKIPPED)", blog)))
bzip <- tempfile(fileext = ".zip")
check("bundle zips", { zip_directory(bdir, bzip); file.size(bzip) > 10000 })
unlink(c(bdir, bzip), recursive = TRUE)

# ---- 12. differential expression in the app (Analysis tab) ------------------

cnt_r <- read_matrix_auto(d("example_rsem_gene_counts.tsv.gz"))
meta_ex <- read_sample_metadata(d("example_samplesheet.csv"))
ss <- design_scaffold(colnames(cnt_r$matrix))
for (v in c("condition", "donor")) ss[[v]] <- as.character(meta_ex[[v]][match(ss$sample, meta_ex$sample)])
ty  <- list(condition = "factor", donor = "factor")
inc <- ss$sample
cons <- data.frame(id = contrast_id("condition", "stim", "mock"), variable = "condition",
                   reference = "mock", target = "stim", stringsAsFactors = FALSE)

check("sample names split into the requested piece",
      identical(guess_from_names(c("A_x_1", "B_y_2"), "_", 2), c("x", "y")) &&
      identical(guess_from_names("A_x_1", "_", -1), "1") &&
      is.na(guess_from_names("A", "_", 3)))
check("variable names are cleaned and must be new",
      clean_variable_name("cell type") == "cell.type" &&
      inherits(tryCatch(clean_variable_name("donor", "donor"), error = function(e) e), "error"))
check("the model puts the main variable last",
      identical(paste(deparse(design_formula("condition", "donor")), collapse = ""), "~donor + condition"))
check("the example's paired design passes the checks",
      !length(validate_design(ss, inc, "condition", "donor", ty, cons)$errors))
check("a covariate confounded with the main variable is caught", {
  s2 <- ss; s2$batch <- s2$condition
  any(grepl("confounded", validate_design(s2, inc, "condition", "batch",
                                          c(ty, list(batch = "factor")), cons)$errors)) })
check("a sample with no value is caught", {
  s3 <- ss; s3$condition[1] <- ""
  any(grepl("no value", validate_design(s3, inc, "condition", character(0), ty, cons)$errors)) })
check("a design with no comparisons is caught",
      any(grepl("comparison", validate_design(ss, inc, "condition", character(0), ty, cons[0, ])$errors)))
check("a level missing from the included samples is caught", {
  keep <- ss$sample[ss$condition == "stim"]
  any(grepl("only one level", validate_design(ss, keep, "condition", character(0), ty, cons)$errors)) })
flt <- filter_counts(round(cnt_r$matrix[, inc]), 10, "group", group = ss$condition)
check("smallest-group filtering requires 3 of the 6 samples (3 per group)",
      flt$n_required == 3 && sum(flt$keep) > 5000 && sum(flt$keep) < nrow(cnt_r$matrix))
check("proportion filtering rounds up",
      filter_counts(matrix(1:8, 2), 1, "proportion", min_proportion = 0.6)$n_required == 3)

if (HAS_DESEQ2) {
  dres <- run_de_analysis(cnt_r$matrix, ss, inc, "condition", "donor", ty, cons, "deseq2")
  r1 <- dres$results[[cons$id]]
  check("DESeq2 returns one table per comparison, one row per kept gene",
        identical(names(dres$results), cons$id) && nrow(r1) == length(dres$kept))
  check("DESeq2 output loads as a standard DE table", {
    st <- standardise_de_table(r1)
    all(c("gene_id", "log2FoldChange", "pvalue", "padj") %in% colnames(st$data)) })
  check("the transformed matrix covers kept genes and included samples",
        identical(dim(dres$transformed), c(length(dres$kept), length(inc))))
  check("non-integer counts are reported in the run notes",
        any(grepl("rounded", dres$notes)))
  check("library sizes and detected genes cover every included sample",
        identical(names(dres$qc$lib_size), inc) && all(dres$qc$detected > 0))
  check("DESeq2 dispersions are kept for the QC plot",
        nrow(dres$qc$dispersion) == length(dres$kept) && !is.null(dres$qc$size_factors))
  for (w in QC_PLOTS) check(paste("QC figure builds:", w), {
    o <- build_qc_figure(w, dres)
    f <- tempfile(fileext = ".png")
    export_plot(o, f, "png", 6, 5, 60); file.size(f) > 1000 })
  check("the per-sample QC table has one row per sample",
        nrow(qc_sample_summary(dres)) == length(inc))
  check("the run finds differentially expressed genes",
        de_result_counts(r1, 0.05, 1)[["up"]] + de_result_counts(r1, 0.05, 1)[["down"]] > 0)
} else {
  cat("skip DESeq2 is not installed: the engine run is not tested\n")
}

# ---- 13. annotation download and gene lengths from a GTF -----------------------

gtf <- tempfile(fileext = ".gtf")
writeLines(c(
  "#!genome-build GRCh38",
  paste(sep = "\t", "1\tensembl\tgene\t100\t400\t.\t+\t.", 'gene_id "ENSG00000000001"; gene_version "3"; gene_name "AAA1"; gene_biotype "protein_coding";'),
  paste(sep = "\t", "1\tensembl\ttranscript\t100\t400\t.\t+\t.", 'gene_id "ENSG00000000001"; transcript_id "T1"; gene_name "AAA1";'),
  # transcript 1: exons 100-199 and 300-400; transcript 2: 150-249 (overlaps the first)
  paste(sep = "\t", "1\tensembl\texon\t100\t199\t.\t+\t.", 'gene_id "ENSG00000000001"; transcript_id "T1";'),
  paste(sep = "\t", "1\tensembl\texon\t300\t400\t.\t+\t.", 'gene_id "ENSG00000000001"; transcript_id "T1";'),
  paste(sep = "\t", "1\tensembl\texon\t150\t249\t.\t+\t.", 'gene_id "ENSG00000000001"; transcript_id "T2";'),
  paste(sep = "\t", "2\tensembl\tgene\t10\t59\t.\t-\t.", 'gene_id "ENSG00000000002"; gene_name "BBB2"; gene_biotype "lncRNA";'),
  paste(sep = "\t", "2\tensembl\texon\t10\t59\t.\t-\t.", 'gene_id "ENSG00000000002"; transcript_id "T3";')), gtf)
gl <- gtf_gene_lengths(gtf)
check("GTF gene length is the union of overlapping exons",
      identical(unname(gl[c("ENSG00000000001", "ENSG00000000002")]), c(251, 50)))
check("GTF lengths match version-suffixed matrix ids",
      identical(unname(match_gene_lengths(gl, c("ENSG00000000002.4", "ENSG00000000009.1"))), 50) )
check("a GTF gives symbols and biotypes too",
      identical(read_gtf_annotation(gtf)$gene_name, c("AAA1", "BBB2")))
check("GENCODE ids are recognised as GENCODE human",
      identical(guess_annotation_source(c("ENSG00000000003.15", "ENSG00000000005.6"))[c("species", "provider")],
                list(species = "human", provider = "gencode")))
check("Ensembl mouse ids are recognised",
      identical(guess_annotation_source(c("ENSMUSG00000000001", "ENSMUSG00000000003"))$provider, "ensembl"))
fake_listing <- function(u) paste('<a href="Homo_sapiens.GRCh38.115.gtf.gz">', '<a href="Homo_sapiens.GRCh38.115.chr.gtf.gz">',
                                  '<a href="gencode.v48.annotation.gtf.gz">', '<a href="gencode.v48.primary_assembly.annotation.gtf.gz">')
check("current Ensembl GTF is found in the listing",
      identical(annotation_url("ensembl", "human", "current", fake_listing)$url,
                "https://ftp.ensembl.org/pub/current/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"))
check("latest GENCODE GTF is found in the listing",
      identical(annotation_url("gencode", "human", "latest", fake_listing)$file, "gencode.v48.annotation.gtf.gz"))
check("pinned releases build the archive address",
      identical(annotation_url("ensembl", "mouse", "110")$url,
                "https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz") &&
      identical(annotation_url("gencode", "mouse", "M37")$url,
                "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.annotation.gtf.gz"))
check("lengths from a GTF give TPM", {
  cm <- matrix(c(10, 20, 30, 40), 2, dimnames = list(c("ENSG00000000001", "ENSG00000000002"), c("s1", "s2")))
  tp <- normalise_counts(cm, match_gene_lengths(gl, rownames(cm)), "tpm")
  isTRUE(all.equal(unname(colSums(tp)), c(1e6, 1e6))) })

# ---- 14. full gene lists and the Gene explorer ------------------------------------

de_c <- classify_degs(de_d$data, 0.05, 1)
full <- de_full_table(de_c)
check("the full gene list keeps every gene, most significant first",
      nrow(full) == nrow(de_c) && identical(which.min(full$padj), 1L) &&
      all(c("gene_id", "gene_name", "log2FoldChange", "pvalue", "padj", "status") %in% colnames(full)))
check("the full gene list carries annotation columns", "gene_biotype" %in% colnames(full))
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- tempfile(fileext = ".xlsx"); write_de_workbook(de_c, wb, 0.05, 1)
  check("the Excel workbook has all, significant, up, down and thresholds sheets",
        identical(openxlsx::getSheetNames(wb), c("All genes", "Significant", "Up", "Down", "Thresholds")) &&
        nrow(openxlsx::read.xlsx(wb, "Up")) == sum(de_c$status == "Up"))
}
cnt_all <- read_matrix_auto(d("example_rsem_gene_counts.tsv.gz"))
catg <- gene_catalog(rownames(cnt_all$matrix), de_c, NULL, cnt_all$annotation)
check("the gene catalogue covers every matrix gene and every tested gene",
      all(rownames(cnt_all$matrix) %in% catg$gene_id) && all(de_c$gene_id %in% catg$gene_id))
check("matrix genes get symbols from the matrix's gene_name column when the DE table lacks them",
      sum(catg$symbol != catg$gene_id) > sum(de_c$display_gene != de_c$gene_id))
ex_n <- normalise_counts(cnt_all$matrix, NULL, "deseq2")
gtab <- gene_expression_table(catg, ex_n, meta_ex, de_c)
check("the all-genes table has ids, symbols, stats, group means and every sample",
      nrow(gtab) == nrow(catg) &&
      all(c("Gene ID", "Symbol", "log2 FC", "P value", "adj. P", "Status", "mean: mock", "mean: stim") %in% colnames(gtab)) &&
      all(colnames(ex_n) %in% colnames(gtab)))
check("genes outside the DE results are marked 'Not tested'",
      all(gtab$Status[!gtab$`Gene ID` %in% de_c$gene_id] == "Not tested"))
ifit1 <- catg$gene_id[match("IFIT1", catg$symbol)]
check("group means match the samples", {
  ctl <- meta_ex$sample[meta_ex$condition == "mock"]
  isTRUE(all.equal(gtab$`mean: mock`[gtab$`Gene ID` == ifit1], round(mean(ex_n[ifit1, ctl]), 3)))
})
for (lay in c("samples", "groups")) check(paste("gene expression plot builds:", lay), {
  p <- gene_bar_plot(ex_n, meta_ex, c(ifit1, catg$gene_id[match("MX1", catg$symbol)]), catg, de_c,
                     layout = lay, log_scale = lay == "groups")
  f <- tempfile(fileext = ".png"); export_plot(p, f, "png", 6, 4, 60); file.size(f) > 1000 })
check("an unknown gene gives a clear message",
      grepl("None of the selected", tryCatch(gene_bar_plot(ex_n, meta_ex, "NOPE", catg),
                                              error = function(e) conditionMessage(e))))

# ---- 15. other DE engines: limma-voom, edgeR, EBSeq ---------------------------

check("EBSeq refuses covariates with a clear message",
      any(grepl("EBSeq", engine_design_errors("ebseq", "donor", ty, de_default_options()))))
check("a limma blocking factor must be a ticked factor covariate", {
  o <- de_default_options(); o$limma_block <- "batch"
  length(engine_design_errors("limma", "donor", ty, o)) > 0 })
check("a ticked blocking factor is accepted", {
  o <- de_default_options(); o$limma_block <- "donor"
  !length(engine_design_errors("limma", "donor", ty, o)) })
check("the checks list engine-specific errors",
      any(grepl("EBSeq", c(validate_design(ss, inc, "condition", "donor", ty, cons)$errors,
                           engine_design_errors("ebseq", "donor", ty, de_default_options())))))

engine_runs <- list(
  limma = list(engine = "limma", covars = "donor", opts = list()),
  limma_quality = list(engine = "limma", covars = "donor", opts = list(limma_voom = "quality")),
  edger_ql = list(engine = "edger", covars = "donor", opts = list(edger_test = "ql")),
  edger_lrt = list(engine = "edger", covars = "donor", opts = list(edger_test = "lrt")),
  ebseq = list(engine = "ebseq", covars = character(0), opts = list(ebseq_maxround = 3)))
for (nm in names(engine_runs)) {
  er <- engine_runs[[nm]]
  if (!is.null(de_engine_missing(er$engine))) { cat("skip", nm, "- not installed:", de_engine_missing(er$engine), "\n"); next }
  o <- modifyList(de_default_options(), er$opts)
  # a gene subset keeps the smoke test quick; the engines are validated on the full matrix elsewhere
  sub <- cnt_r$matrix[seq_len(min(4000, nrow(cnt_r$matrix))), ]
  rr <- tryCatch(run_de_analysis(sub, ss, inc, "condition", er$covars, ty, cons, er$engine, o),
                 error = function(e) e)
  check(paste(nm, "runs"), !inherits(rr, "error"))
  if (inherits(rr, "error")) { cat("  ", conditionMessage(rr), "\n"); next }
  r1 <- rr$results[[cons$id]]
  check(paste(nm, "returns a standard table for every kept gene"),
        nrow(r1) == length(rr$kept) &&
        all(c("gene_id", "log2FoldChange", "pvalue", "padj") %in% colnames(r1)) &&
        all(!is.na(r1$log2FoldChange)))
  check(paste(nm, "output loads into the explorer"),
        all(c("gene_id", "log2FoldChange", "padj") %in% colnames(standardise_de_table(r1)$data)))
  check(paste(nm, "transformed matrix and size factors cover every sample"),
        identical(dim(rr$transformed), c(length(rr$kept), length(inc))) &&
        length(rr$qc$size_factors) == length(inc))
  check(paste(nm, "finds DE genes in the stim-vs-mock demo"),
        sum(unlist(de_result_counts(r1, 0.05, 1)[c("up", "down")])) > 0)
  check(paste(nm, "agrees with the direction of IFIT1 (up after treatment)"), {
    i <- match(catg$gene_id[match("IFIT1", catg$symbol)], r1$gene_id)
    is.na(i) || r1$log2FoldChange[i] > 2 })
  check(paste(nm, "records engine label and package versions"),
        nzchar(rr$engine_label) && length(rr$versions) > 0)
  for (w in QC_PLOTS) {
    fig <- tryCatch(build_qc_figure(w, rr), error = function(e) e)
    if (er$engine == "ebseq" && w == "dispersion") {
      check("EBSeq dispersion panel explains why it is empty",
            inherits(fig, "error") && grepl("EBSeq", conditionMessage(fig)))
      next
    }
    check(paste(nm, "QC figure builds:", w), !inherits(fig, "error") && {
      f <- tempfile(fileext = ".png"); export_plot(fig, f, "png", 6, 5, 60); file.size(f) > 1000 })
  }
  qs <- qc_sample_summary(rr)
  check(paste(nm, "per-sample QC table has one row per sample"), nrow(qs) == length(inc))
  if (nm == "limma_quality")
    check("voom quality weights appear in the per-sample table",
          "voom quality weight" %in% colnames(qs) && all(!is.na(qs$`voom quality weight`)))
}
if (is.null(de_engine_missing("limma"))) {
  o <- modifyList(de_default_options(), list(limma_block = "donor"))
  sub <- cnt_r$matrix[seq_len(min(3000, nrow(cnt_r$matrix))), ]
  rb <- run_de_analysis(sub, ss, inc, "condition", "donor", ty, cons, "limma", o)
  check("limma with a blocking factor reports the consensus correlation",
        any(grepl("correlation", rb$notes, ignore.case = TRUE)))
  o$limma_lfc <- 1
  rt <- run_de_analysis(sub, ss, inc, "condition", "donor", ty, cons, "limma", modifyList(o, list(limma_block = "(none)")))
  check("limma TREAT (lfc threshold) runs and is no more liberal than the plain test",
        sum(rt$results[[1]]$padj < 0.05, na.rm = TRUE) <=
        sum(run_de_analysis(sub, ss, inc, "condition", "donor", ty, cons, "limma")$results[[1]]$padj < 0.05, na.rm = TRUE))
}

# ---- 16. shrinkage, LRT, length offsets, run record, script, bundle, design files ----

len_ex <- read_matrix_auto(d("example_rsem_gene_lengths.tsv.gz"))$matrix
sub16 <- cnt_r$matrix[seq_len(min(4000, nrow(cnt_r$matrix))), ]
check("length offsets need sample-specific lengths", {
  g <- rownames(sub16)[1:50]
  is.null(length_offset_matrix(NULL, g, inc)$matrix) &&
    is.null(length_offset_matrix(setNames(len_ex[, 1], rownames(len_ex)), g, inc)$matrix) &&
    identical(dim(length_offset_matrix(len_ex, g, inc)$matrix), c(50L, length(inc)))
})
check("missing lengths are filled and never zero", {
  l2 <- len_ex; l2[rownames(sub16)[1], 1] <- NA; l2[rownames(sub16)[2], ] <- 0
  m <- length_offset_matrix(l2, rownames(sub16)[1:5], inc)$matrix
  all(is.finite(m)) && all(m > 0) })
check("the recommended preset filters 10 counts in 1 sample and shrinks with ashr",
      with(recommended_options(), min_abundance == 10 && filter_mode == "samples" &&
             min_samples == 1 && shrink == "ashr" && isTRUE(length_offset)))
check("a missing shrinkage package is reported before running", {
  o <- de_default_options(); o$shrink <- "apeglm"
  requireNamespace("apeglm", quietly = TRUE) || any(grepl("apeglm", engine_design_errors("deseq2", "donor", ty, o))) })

cons16 <- rbind(cons, data.frame(id = "condition_mock_vs_stim", variable = "condition",
                                 reference = "stim", target = "mock"))
if (HAS_DESEQ2) {
  base16 <- run_de_analysis(sub16, ss, inc, "condition", "donor", ty, cons, "deseq2")
  for (sh in c("normal", "ashr", "apeglm")) {
    if (sh != "normal" && !requireNamespace(sh, quietly = TRUE)) { cat("skip", sh, "is not installed\n"); next }
    r16 <- run_de_analysis(sub16, ss, inc, "condition", "donor", ty, cons16, "deseq2",
                           list(shrink = sh))
    a <- r16$results[[1]]; b <- base16$results[[1]]; m <- match(a$gene_id, b$gene_id)
    check(paste(sh, "shrinkage keeps p-values and the unshrunken fold change"),
          isTRUE(all.equal(a$pvalue, b$pvalue[m])) && isTRUE(all.equal(a$log2FoldChange_mle, b$log2FoldChange[m])))
    check(paste(sh, "shrinkage pulls fold changes towards zero"),
          median(abs(a$log2FoldChange), na.rm = TRUE) < median(abs(a$log2FoldChange_mle), na.rm = TRUE))
    if (sh != "apeglm") {
      check(paste(sh, "shrinkage gives mirror-image fold changes for the reversed comparison"), {
        x <- r16$results[[1]]; y <- r16$results[[2]]
        cor(x$log2FoldChange, -y$log2FoldChange[match(x$gene_id, y$gene_id)], use = "c") > 0.99 })
    } else {
      # apeglm is not symmetric (its estimate for low-count genes depends on the base level),
      # so check the reversed comparison against a fresh fit with 'stim' as the base level.
      check("apeglm on a reversed comparison matches a fresh DESeq2 fit with that base level", {
        y <- r16$results[[2]]
        cdf <- data.frame(row.names = inc, donor = factor(ss$donor), condition = relevel(factor(ss$condition), "stim"))
        kk <- y$gene_id
        f2 <- suppressMessages(DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(round(sub16[kk, inc]), cdf, ~ donor + condition), quiet = TRUE))
        s2 <- suppressMessages(DESeq2::lfcShrink(f2, coef = "condition_mock_vs_stim", type = "apeglm", quiet = TRUE))
        cor(y$log2FoldChange, s2[kk, "log2FoldChange"], use = "c") > 0.9999 })
      check("apeglm runs carry a note about asymmetry", any(grepl("not symmetric", r16$notes)))
    }
  }
  lrt16 <- run_de_analysis(sub16, ss, inc, "condition", "donor", ty, cons16, "deseq2", list(deseq2_test = "lrt"))
  check("the LRT gives one p-value per gene, shared by every comparison", {
    x <- lrt16$results[[1]]; y <- lrt16$results[[2]]
    isTRUE(all.equal(x$pvalue, y$pvalue[match(x$gene_id, y$gene_id)])) && any(grepl("Likelihood-ratio", lrt16$notes)) })
  off16 <- run_de_analysis(sub16, ss, inc, "condition", "donor", ty, cons, "deseq2",
                           list(length_offset = TRUE), lengths = len_ex)
  check("DESeq2 length offsets change the fit and are recorded",
        isTRUE(off16$length_offset) && !isTRUE(all.equal(off16$results[[1]]$pvalue, base16$results[[1]]$pvalue)) &&
        length(off16$qc$size_factors) == length(inc))
  check("the run record lists samples, comparisons, settings and versions", {
    rec <- de_run_record(off16)
    any(grepl("^## Samples", rec)) && any(grepl("condition_stim_vs_mock", rec)) &&
      any(grepl("Gene-length offsets \\| yes", rec)) && any(grepl("DESeq2", rec)) })
}
if (is.null(de_engine_missing("edger"))) {
  e0 <- run_de_analysis(sub16, ss, inc, "condition", "donor", ty, cons, "edger")
  e1 <- run_de_analysis(sub16, ss, inc, "condition", "donor", ty, cons, "edger", list(length_offset = TRUE), lengths = len_ex)
  check("edgeR length offsets change the fit",
        isTRUE(e1$length_offset) && !isTRUE(all.equal(e0$results[[1]]$pvalue, e1$results[[1]]$pvalue)))
}
check("limma reports that it does not use length offsets", is.null(de_engine_missing("limma")) && {
  l1 <- run_de_analysis(sub16, ss, inc, "condition", "donor", ty, cons, "limma", list(length_offset = TRUE), lengths = len_ex)
  !isTRUE(l1$length_offset) && any(grepl("not used by limma", l1$notes)) } || !is.null(de_engine_missing("limma")))

# The generated script must reproduce the app's numbers from the bundle alone.
repro_cfgs <- list(deseq2 = list(shrink = if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal", length_offset = TRUE),
                   limma = list(limma_voom = "quality"), edger = list(length_offset = TRUE), ebseq = list(ebseq_maxround = 3))
for (eng in names(repro_cfgs)) {
  if (!is.null(de_engine_missing(eng))) next
  rr <- run_de_analysis(sub16, ss, inc, "condition", if (eng == "ebseq") character(0) else "donor", ty, cons16,
                        eng, repro_cfgs[[eng]], lengths = len_ex)
  z <- tempfile(fileext = ".zip")
  write_de_bundle(rr, z, qc = eng == "deseq2")
  ex <- tempfile(); dir.create(ex); utils::unzip(z, exdir = ex)
  check(paste(eng, "bundle has results, design, script and record"),
        all(file.exists(file.path(ex, c("reproduce.R", "run_record.md", "samplesheet.csv", "contrasts.csv", "counts.tsv")))) &&
        length(list.files(file.path(ex, "results"))) == nrow(cons16) &&
        (eng != "deseq2" || (file.exists(file.path(ex, "lengths.tsv")) && file.exists(file.path(ex, "qc", "pca.png")))))
  owd <- setwd(ex)
  st <- system2(file.path(R.home("bin"), "Rscript"), "reproduce.R", stdout = FALSE, stderr = FALSE)
  setwd(owd)
  check(paste(eng, "generated R script runs and reproduces every comparison"), st == 0 && all(vapply(names(rr$results), function(id) {
    f <- list.files(file.path(ex, "reproduced"), paste0("^", id, "\\."), full.names = TRUE)
    if (length(f) != 1) return(FALSE)
    b <- utils::read.delim(f, check.names = FALSE); a <- rr$results[[id]]; m <- match(a$gene_id, b$gene_id)
    # ashr / apeglm shrinkage is an optimisation that converges to ~1e-8, so shrunken
    # fold changes are compared at 1e-6; the unshrunken estimates must match exactly.
    shr <- "log2FoldChange_mle" %in% colnames(a)
    !anyNA(m) && isTRUE(all.equal(a$log2FoldChange, b$log2FoldChange[m], tolerance = if (shr) 1e-6 else 1e-8)) &&
      (!shr || isTRUE(all.equal(a$log2FoldChange_mle, b$log2FoldChange_mle[m], tolerance = 1e-8))) &&
      isTRUE(all.equal(a$padj, b$padj[m], tolerance = 1e-8))
  }, logical(1))))
}

# A main variable with a third level: each EBSeq comparison tests its own two groups,
# and genes it cannot test there (e.g. all-zero in those samples) come back as NA.
if (is.null(de_engine_missing("ebseq"))) {
  ss3 <- ss; ss3$condition <- c("mock", "stim", "mock", "stim", "high", "high")
  cons3 <- rbind(cons, data.frame(id = "condition_high_vs_mock", variable = "condition",
                                  reference = "mock", target = "high"))
  r3 <- tryCatch(run_de_analysis(sub16, ss3, inc, "condition", character(0),
                                 list(condition = "factor"), cons3, "ebseq", list(ebseq_maxround = 3)),
                 error = function(e) e)
  check("EBSeq runs on a 3-level main variable (pairwise)", !inherits(r3, "error") &&
        length(r3$results) == 2 && nrow(r3$results[[1]]) == length(r3$kept))
  if (!inherits(r3, "error"))
    check("EBSeq keeps every filtered gene in each comparison, NA where untestable",
          all(vapply(r3$results, function(df) all(!is.na(df$gene_id)) &&
                       all(is.na(df$padj) == is.na(df$log2FoldChange)), logical(1))))
}

check("a saved design reads back: sample sheet and contrasts", {
  f1 <- tempfile(fileext = ".csv"); f2 <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(sample = rev(inc), condition = rev(ss$condition), age = c(30, 41, 52, 63, 70, 81)), f1, row.names = FALSE)
  utils::write.csv(de_contrasts_table(cons16, "donor"), f2, row.names = FALSE)
  s1 <- read_design_samplesheet(f1, inc); c1 <- read_design_contrasts(f2)
  identical(s1$sheet$condition, ss$condition) && identical(s1$types$age, "numeric") && s1$matched == length(inc) &&
    identical(c1$contrasts$id, cons16$id) && identical(c1$blocking, "donor") })
check("a contrasts file without ids gets automatic names", {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(variable = "condition", reference = "mock", target = "stim"), f, row.names = FALSE)
  identical(read_design_contrasts(f)$contrasts$id, "condition_stim_vs_mock") })
check("a sample sheet that matches no sample is refused", {
  f <- tempfile(fileext = ".csv"); utils::write.csv(data.frame(sample = "x", condition = "a"), f, row.names = FALSE)
  inherits(tryCatch(read_design_samplesheet(f, inc), error = function(e) e), "error") })

# ---- 17. other organisms -------------------------------------------------------------

check("organisms are recognised from gene and transcript ids", identical(
  vapply(list(c("ENSG00000000003", "ENSG00000000005.6"), c("ENSMUSG00000000001"), c("ENSRNOG00000000001"),
              c("ENSDARG00000000001"), c("FBgn0000003"), c("WBGene00000001"), c("AT1G01010"),
              c("YAL001C", "YAL002W"), c("ENST00000456328.2"), c("foo")), guess_organism, ""),
  c("human", "mouse", "rat", "zebrafish", "fly", "worm", "arabidopsis", "yeast", "human", NA)))
check("transcript ids are told apart from gene ids",
      looks_like_transcripts(c("ENSMUST00000000001", "ENSMUST00000000010")) &&
      !looks_like_transcripts(c("ENSMUSG00000000001")) && looks_like_transcripts(c("AT1G01010.1", "AT1G01020.2")))
fake_rat <- function(u) paste('<a href="Rattus_norvegicus.GRCr8.114.gtf.gz">', '<a href="Rattus_norvegicus.GRCr8.115.gtf.gz">',
                             '<a href="Rattus_norvegicus.GRCr8.115.chr.gtf.gz">', '<a href="Rattus_norvegicus.GRCr8.115.abinitio.gtf.gz">')
check("Ensembl GTF for any species is found from the directory listing", {
  a <- annotation_url("ensembl", "rat", "current", fake_rat)
  identical(a$file, "Rattus_norvegicus.GRCr8.115.gtf.gz") && grepl("/pub/current/gtf/rattus_norvegicus/", a$url) })
check("Ensembl Genomes divisions (plants) use their own server",
      grepl("ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.61.gtf.gz",
            annotation_url("ensembl", "arabidopsis", "61", function(u) '<a href="Arabidopsis_thaliana.TAIR10.61.gtf.gz">')$url, fixed = TRUE))
check("'Other' takes any Ensembl species name", {
  a <- annotation_url("ensembl", "other", "current", function(u) '<a href="Ovis_aries.ARS-UI_Ramb_v2.0.115.gtf.gz">',
                      other = list(ensembl_dir = "ovis_aries"))
  identical(a$file, "Ovis_aries.ARS-UI_Ramb_v2.0.115.gtf.gz") })
check("GENCODE is refused for organisms it does not cover",
      inherits(tryCatch(annotation_url("gencode", "rat", "current", fake_rat), error = function(e) e), "error"))
check("enrichment asks for the organism's OrgDb and flags human-only gene sets", {
  m <- enrichment_missing_packages("BP", organism = organism_info("mouse"))
  g <- enrichment_missing_packages("GOREL_BP", organism = organism_info("mouse"))
  (requireNamespace("org.Mm.eg.db", quietly = TRUE) || any(grepl("org.Mm.eg.db", m))) && any(grepl("human GAF", g)) })
check("gene links follow the organism (KEGG code; no human-only resources for mouse)", {
  l <- gene_resource_links("Ifit1", "ENSMUSG00000034459", "Mus musculus")
  grepl("dbkey=mmu", l$KEGG) && is.null(l$GTEx) && is.null(l[["Human Protein Atlas"]]) && !is.null(l$Ensembl) })

# ---- 18. transcript (isoform)-level analysis ------------------------------------------

check("GENCODE pipe ids split into transcript, gene and names", {
  x <- split_tx_ids(c("ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|", "ENST1"))
  identical(x$transcript_id, c("ENST00000456328.2", "ENST1")) && identical(x$gene_name[1], "DDX11L2") &&
    identical(x$transcript_name[1], "DDX11L2-202") && is.na(x$gene_id[2]) })

gtf_tx <- tempfile(fileext = ".gtf")
writeLines(c(paste("1", "t", "exon", "100", "199", ".", "+", ".", 'gene_id "G1"; transcript_id "T1"; gene_name "AAA"; transcript_name "AAA-201"; transcript_biotype "protein_coding";', sep = "\t"),
             paste("1", "t", "exon", "300", "349", ".", "+", ".", 'gene_id "G1"; transcript_id "T1"; gene_name "AAA"; transcript_name "AAA-201"; transcript_biotype "protein_coding";', sep = "\t"),
             paste("1", "t", "exon", "100", "149", ".", "+", ".", 'gene_id "G1"; transcript_id "T2"; gene_name "AAA"; transcript_name "AAA-202"; transcript_biotype "retained_intron";', sep = "\t")), gtf_tx)
check("a GTF gives the transcript table (gene, names, biotype, exon-sum length)", {
  t <- gtf_transcript_table(gtf_tx)
  identical(t$transcript_id, c("T1", "T2")) && identical(t$tx_length, c(150, 50)) &&
    identical(t$transcript_biotype[2], "retained_intron") && identical(t$gene_name[1], "AAA") })

# A small Salmon project with Gibbs resamples: genes split into 1-4 isoforms,
# resample variance growing with the number of isoforms (read-to-transcript ambiguity).
set.seed(7)
tx_root <- file.path(tempdir(), "salmon_fixture"); unlink(tx_root, recursive = TRUE); dir.create(tx_root)
g_m <- round(cnt_r$matrix[order(-rowMeans(cnt_r$matrix))[1:500], ])
niso <- sample(1:4, nrow(g_m), replace = TRUE)
tx_rows <- unlist(lapply(seq_len(nrow(g_m)), function(i) rep(i, niso[i])))
tx_id <- unlist(lapply(seq_len(nrow(g_m)), function(i) sprintf("%sT%02d", rownames(g_m)[i], seq_len(niso[i]))))
prop <- unlist(lapply(niso, function(n) { p <- rgamma(n, 1); p / sum(p) }))
tx_cnt <- matrix(rpois(length(tx_rows) * ncol(g_m), g_m[tx_rows, ] * prop), length(tx_rows), dimnames = list(tx_id, colnames(g_m)))
tx_len <- round(runif(nrow(tx_cnt), 500, 4000))
for (smp in colnames(tx_cnt)) {
  d <- file.path(tx_root, smp); dir.create(file.path(d, "aux_info", "bootstrap"), recursive = TRUE)
  write.table(data.frame(Name = tx_id, Length = tx_len, EffectiveLength = tx_len - 150 + sample(-3:3, length(tx_len), TRUE),
                         TPM = 1, NumReads = tx_cnt[, smp]), file.path(d, "quant.sf"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(sprintf('{"num_targets": %d, "num_bootstraps": 20, "samp_type": "gibbs"}', nrow(tx_cnt)), file.path(d, "aux_info", "meta_info.json"))
  phi <- 1 + 2 * (niso[tx_rows] - 1)
  b <- sapply(1:20, function(k) pmax(0, tx_cnt[, smp] + rnorm(nrow(tx_cnt), 0, sqrt(phi * pmax(tx_cnt[, smp], 0.5)))))
  con <- gzcon(file(file.path(d, "aux_info", "bootstrap", "bootstraps.gz"), "wb")); writeBin(as.vector(b), con); close(con)
}
tx_map <- data.frame(transcript_id = tx_id, gene_id = rownames(g_m)[tx_rows],
                     gene_name = catg$symbol[match(rownames(g_m)[tx_rows], catg$gene_id)],
                     transcript_name = paste0(catg$symbol[match(rownames(g_m)[tx_rows], catg$gene_id)], "-", 200 + sequence(niso)))
q_tx <- read_quant_dirs(tx_root)
# Overdispersion comes from edgeR's catchSalmon; without edgeR the counts and
# lengths are still read and the overdispersion is left out.
check("Salmon folders give counts, effective lengths and overdispersion", {
  identical(dim(q_tx$counts), dim(tx_cnt)) && identical(q_tx$tool, "salmon") &&
    if (requireNamespace("edgeR", quietly = TRUE))
      !is.null(q_tx$overdispersion) && all(q_tx$overdispersion >= 1) &&
        # more isoforms -> more ambiguity
        median(q_tx$overdispersion[niso[tx_rows] == 4]) > 2 * median(q_tx$overdispersion[niso[tx_rows] == 1])
    else is.null(q_tx$overdispersion) })
check("without resamples there is no overdispersion and the note says how to get it", {
  r2 <- file.path(tempdir(), "salmon_noboot"); unlink(r2, recursive = TRUE); dir.create(r2)
  for (smp in colnames(tx_cnt)[1:2]) {
    d <- file.path(r2, smp); dir.create(file.path(d, "aux_info"), recursive = TRUE)
    file.copy(file.path(tx_root, smp, "quant.sf"), d)
    writeLines(sprintf('{"num_targets": %d, "num_bootstraps": 0}', nrow(tx_cnt)), file.path(d, "aux_info", "meta_info.json"))
  }
  q0 <- read_quant_dirs(r2); is.null(q0$overdispersion) && grepl("numGibbsSamples", q0$note) })
check("a merged transcript table carries its transcript-to-gene map", {
  f <- tempfile(fileext = ".tsv")
  write.table(data.frame(tx = tx_id, gene_id = tx_map$gene_id, tx_cnt, check.names = FALSE), f, sep = "\t", quote = FALSE, row.names = FALSE)
  t <- read_transcript_matrix(f); identical(dim(t$matrix), dim(tx_cnt)) && identical(t$tx2gene$gene_id, tx_map$gene_id) })
check("Ng vectors read as one column (row order) or two columns (id, Ng)", {
  f1 <- tempfile(); writeLines(as.character(pmin(niso[tx_rows], 3)), f1)
  f2 <- tempfile(); write.table(data.frame(tx_id, pmin(niso[tx_rows], 3)), f2, quote = FALSE, row.names = FALSE, col.names = FALSE)
  f3 <- tempfile(); writeLines(c("1", "2"), f3)
  identical(unname(read_ngvector(f1, tx_id)), as.integer(pmin(niso[tx_rows], 3))) &&
    identical(read_ngvector(f2, rev(tx_id)), read_ngvector(f1, tx_id)[rev(tx_id)]) &&
    inherits(tryCatch(read_ngvector(f3, tx_id), error = function(e) e), "error") })

tx_in <- list(tx2gene = tx_map, overdispersion = q_tx$overdispersion)
ss_tx <- design_scaffold(colnames(q_tx$counts))
for (v in c("condition", "donor")) ss_tx[[v]] <- ss[[v]][match(ss_tx$sample, ss$sample)]
tx_runs <- list()
for (eng in c("deseq2", "limma", "edger", "ebseq")) {
  if (!is.null(de_engine_missing(eng))) next
  o <- list(feature_level = "transcript", length_offset = eng == "deseq2", ebseq_maxround = 3)
  r <- tryCatch(run_de_analysis(q_tx$counts, ss_tx, inc, "condition", if (eng == "ebseq") character(0) else "donor",
                                ty, cons16, eng, o, lengths = q_tx$lengths, transcripts = tx_in), error = function(e) e)
  check(paste(eng, "runs at transcript level"), !inherits(r, "error"))
  if (inherits(r, "error")) { cat("  ", conditionMessage(r), "\n"); next }
  tx_runs[[eng]] <- r
  d1 <- r$results[[1]]
  check(paste(eng, "transcript results carry the parent gene and transcript names"),
        all(c("parent_gene_id", "gene_name", "transcript_name") %in% colnames(d1)) &&
        identical(d1$parent_gene_id, tx_map$gene_id[match(d1$gene_id, tx_map$transcript_id)]))
  if (eng %in% c("limma", "edger"))
    check(paste(eng, "analyses counts divided by their overdispersion, unrounded"),
          isTRUE(r$divided) && any(r$counts != round(r$counts)) &&
          isTRUE(all.equal(r$counts[1, ], q_tx$counts[rownames(r$counts)[1], inc] / q_tx$overdispersion[rownames(r$counts)[1]])))
  if (eng == "deseq2") check("DESeq2 transcript run uses effective-length offsets", isTRUE(r$length_offset))
  if (eng == "ebseq") check("EBSeq groups isoforms by their gene's isoform count (GetNg)",
                            identical(unname(r$ngvector), as.integer(pmin(niso[tx_rows], 3))[match(r$kept, tx_id)]))
  check(paste(eng, "transcript results reload as a transcript table"), {
    std <- standardise_de_table(de_result_for_export(r, names(r$results)[1]))$data
    identical(std$gene_id, d1$gene_id) && "parent_gene_id" %in% colnames(std) })
}
if (!is.null(tx_runs$edger)) {
  check("dividing out ambiguity changes the edgeR result", {
    r0 <- run_de_analysis(q_tx$counts, ss_tx, inc, "condition", "donor", ty, cons16, "edger",
                          list(feature_level = "transcript", divide_overdispersion = FALSE), transcripts = tx_in)
    !isTRUE(r0$divided) && !isTRUE(all.equal(r0$results[[1]]$pvalue, tx_runs$edger$results[[1]]$pvalue)) })
  check("enrichment maps DE transcripts to their unique genes", {
    st <- standardise_de_table(de_result_for_export(tx_runs$edger, names(tx_runs$edger$results)[1]))$data
    g <- enrichment_ids(st, st$gene_id[1:50])
    !anyDuplicated(g) && all(g %in% tx_map$gene_id) &&
      identical(sort(enrichment_universe(st)), sort(unique(st$parent_gene_id))) })
}
if (requireNamespace("stageR", quietly = TRUE) && !is.null(tx_runs$edger)) {
  d1 <- tx_runs$edger$results[[1]]
  check("stageR adds gene-level screen and transcript confirmation columns",
        all(c("stageR_gene_padj", "stageR_tx_padj") %in% colnames(d1)) &&
        all(is.na(d1$stageR_tx_padj[!is.na(d1$stageR_gene_padj) & d1$stageR_gene_padj > 0.05])) &&
        (!any(d1$stageR_gene_padj <= 0.05, na.rm = TRUE) || any(!is.na(d1$stageR_tx_padj))) &&
        any(grepl("Stage-wise testing \\(stageR", tx_runs$edger$notes)))
} else cat("skip stageR is not installed\n")
for (eng in names(tx_runs)) {
  rr <- tx_runs[[eng]]
  z <- tempfile(fileext = ".zip"); write_de_bundle(rr, z, qc = FALSE)
  ex <- tempfile(); dir.create(ex); utils::unzip(z, exdir = ex)
  owd <- setwd(ex); st <- system2(file.path(R.home("bin"), "Rscript"), "reproduce.R", stdout = FALSE, stderr = FALSE); setwd(owd)
  check(paste(eng, "transcript-level script reproduces every comparison"), st == 0 && all(vapply(names(rr$results), function(id) {
    f <- list.files(file.path(ex, "reproduced"), paste0("^", id, "\\."), full.names = TRUE)
    if (length(f) != 1) return(FALSE)
    b <- utils::read.delim(f, check.names = FALSE); a <- rr$results[[id]]; m <- match(a$gene_id, b$transcript_id)
    !anyNA(m) && isTRUE(all.equal(a$log2FoldChange, b$log2FoldChange[m], tolerance = 1e-8)) &&
      isTRUE(all.equal(a$padj, b$padj[m], tolerance = 1e-8)) &&
      (!"stageR_tx_padj" %in% colnames(a) || isTRUE(all.equal(a$stageR_tx_padj, as.numeric(b$stageR_tx_padj[m]), tolerance = 1e-8)))
  }, logical(1))))
}

# ---- 19. differential transcript usage (DTU) -------------------------------------------

# Planted isoform switches: 40 multi-isoform genes change their transcripts' shares in stim.
set.seed(11)
ng19 <- 400; niso19 <- sample(1:4, ng19, TRUE); gid19 <- sprintf("G%03d", seq_len(ng19))
rows19 <- rep(seq_len(ng19), niso19)
tx19 <- unlist(lapply(seq_len(ng19), function(i) sprintf("%sT%d", gid19[i], seq_len(niso19[i]))))
grp19 <- rep(c("mock", "stim"), each = 3); don19 <- rep(c("d1", "d2", "d3"), 2)
mu19 <- rgamma(ng19, 1.2, 1 / 300); dtu_true <- sample(which(niso19 >= 2), 40)
cnt19 <- matrix(0, length(tx19), 6, dimnames = list(tx19, paste0(don19, "_", grp19)))
for (i in seq_len(ng19)) {
  r <- which(rows19 == i); p0 <- rgamma(niso19[i], 2); p0 <- p0 / sum(p0)
  for (j in 1:6) {
    p <- p0
    if (i %in% dtu_true && grp19[j] == "stim") { p <- rev(p0); if (all(abs(p - p0) < 0.1)) { p <- p0; p[1:2] <- pmax(p0[2:1] + c(0.3, -0.3), 0.01); p <- p / sum(p) } }
    cnt19[r, j] <- rnbinom(length(r), mu = mu19[i] * p * runif(1, 0.8, 1.2), size = 20)
  }
}
od19 <- setNames(1 + 1.5 * (niso19[rows19] - 1) * runif(length(tx19), 0.5, 1.5), tx19)
len19 <- matrix(round(runif(length(tx19), 400, 3000)), length(tx19), 6, dimnames = dimnames(cnt19)) +
  matrix(sample(-5:5, length(cnt19), TRUE), nrow(cnt19))
map19 <- data.frame(transcript_id = tx19, gene_id = gid19[rows19], gene_name = paste0("SYM", rows19),
                    transcript_name = paste0("SYM", rows19, "-", 200 + sequence(niso19)), stringsAsFactors = FALSE)
ss19 <- data.frame(sample = colnames(cnt19), condition = grp19, donor = don19, stringsAsFactors = FALSE)
ty19 <- list(condition = "factor", donor = "factor")
cons19 <- data.frame(id = c("condition_stim_vs_mock", "condition_mock_vs_stim"), variable = "condition",
                     reference = c("mock", "stim"), target = c("stim", "mock"), stringsAsFactors = FALSE)
tx19_in <- list(tx2gene = map19, overdispersion = od19)

check("scaledTPM keeps each sample's total and follows tximport's formula", {
  x <- dtu_counts(cnt19, "scaledTPM", len19, NULL, map19$gene_id)
  ab <- cnt19 / len19; man <- t(t(ab) / colSums(ab) * colSums(cnt19))
  isTRUE(all.equal(x, man)) && isTRUE(all.equal(colSums(x), colSums(cnt19))) })
check("dtuScaledTPM gives every transcript of a gene its gene's median length", {
  x <- dtu_counts(cnt19, "dtuScaledTPM", len19, NULL, map19$gene_id)
  g <- map19$gene_id; med <- tapply(rowMeans(len19), g, median)
  ab <- cnt19 / len19 * as.numeric(med[g]); man <- t(t(ab) / colSums(ab) * colSums(cnt19))
  isTRUE(all.equal(x, man)) })
check("divided counts are counts / overdispersion", isTRUE(all.equal(dtu_counts(cnt19, "divided", NULL, od19), cnt19 / od19)))
check("count type falls back when lengths or resamples are missing, and says so", {
  a <- dtu_resolve_counts("dexseq", "auto", FALSE, TRUE); b <- dtu_resolve_counts("edger_dtu", "auto", TRUE, FALSE)
  identical(a$mode, "estimated") && grepl("lengths", a$note) && identical(b$mode, "estimated") && grepl("resamples", b$note) &&
    identical(dtu_resolve_counts("dexseq", "auto", TRUE, FALSE)$mode, "scaledTPM") &&
    identical(dtu_resolve_counts("limma_dtu", "auto", FALSE, TRUE)$mode, "divided") })
check("DRIMSeq-style filter: gene, transcript count and share rules; genes need two transcripts", {
  m <- rbind(a1 = c(50, 50, 50, 50), a2 = c(20, 20, 20, 20), a3 = c(1, 1, 1, 1),   # a3 too small
             b1 = c(100, 100, 100, 100), b2 = c(5, 5, 5, 5),                        # b keeps one -> dropped
             c1 = c(3, 3, 3, 3), c2 = c(3, 3, 3, 3))                                # gene below 10
  colnames(m) <- paste0("s", 1:4)
  k <- dtu_filter_drimseq(m, c("a", "a", "a", "b", "b", "c", "c"), 4, 10, 2, 10, 2, 0.1)
  identical(unname(k), c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)) })
if (requireNamespace("DRIMSeq", quietly = TRUE)) check("the filter matches DRIMSeq::dmFilter", {
  x <- dtu_counts(cnt19, "scaledTPM", len19, NULL, map19$gene_id)
  k <- dtu_filter_drimseq(x, map19$gene_id, 6, 10, 3, 10, 3, 0.1)
  d <- DRIMSeq::dmDSdata(counts = data.frame(gene_id = map19$gene_id, feature_id = tx19, x, check.names = FALSE),
                         samples = data.frame(sample_id = colnames(x), condition = grp19))
  d <- DRIMSeq::dmFilter(d, min_samps_gene_expr = 6, min_samps_feature_expr = 3, min_samps_feature_prop = 3,
                         min_gene_expr = 10, min_feature_expr = 10, min_feature_prop = 0.1)
  setequal(names(k)[k], DRIMSeq::counts(d)$feature_id) })
check("per-gene q-values: one per gene in [0, 1], ranked like the genes' smallest p", {
  set.seed(3); p <- runif(300)^2; g <- rep(sprintf("g%02d", 1:100), each = 3)
  q <- dtu_per_gene_qvalue(p, g); pm <- tapply(p, g, min)[names(q)]
  length(q) == 100 && all(q >= 0 & q <= 1) && cor(q, pm, method = "spearman") > 0.95 })
if (requireNamespace("DEXSeq", quietly = TRUE)) check("per-gene q-values equal DEXSeq's perGeneQValueExact", {
  set.seed(3); p <- runif(300)^2; g <- rep(sprintf("g%02d", 1:100), each = 3); g[1:30] <- rep(sprintf("h%02d", 1:15), each = 2)
  q <- dtu_per_gene_qvalue(p, g)
  gs <- split(seq_along(p), factor(g)); pG <- sapply(gs, function(i) min(p[i])); th <- unique(sort(pG))
  q2 <- pmin(1, DEXSeq:::perGeneQValueExact(pG, th, gs)[match(pG, th)])
  isTRUE(all.equal(unname(q[names(gs)]), unname(q2))) })

dtu_runs <- list()
for (eng in DTU_ENGINES) {
  if (!is.null(de_engine_missing(eng))) { cat("skip ", engine_label(eng), " is not installed\n", sep = ""); next }
  r <- tryCatch(run_de_analysis(cnt19, ss19, ss19$sample, "condition", "donor", ty19, cons19, eng,
                                list(feature_level = "dtu"), lengths = len19, transcripts = tx19_in),
                error = function(e) e)
  check(paste(engine_label(eng), "runs (DTU)"), !inherits(r, "error"))
  if (inherits(r, "error")) { cat("  ", conditionMessage(r), "\n"); next }
  dtu_runs[[eng]] <- r
  d1 <- r$results[[1]]; d2 <- r$results[[2]]
  check(paste(engine_label(eng), "results carry transcript, gene, usage and gene-level columns"),
        all(c("gene_id", "parent_gene_id", "gene_name", "transcript_name", "log2FoldChange", "pvalue", "padj",
              "gene_padj", "prop_reference", "prop_target", "delta_prop") %in% colnames(d1)) &&
          identical(r$level, "dtu") && all(table(d1$parent_gene_id) >= 2))
  check(paste(engine_label(eng), "uses its workflow's counts"),
        identical(r$dtu_counts, if (eng %in% c("edger_dtu", "limma_dtu")) "divided" else "scaledTPM"))
  sig <- unique(d1$parent_gene_id[!is.na(d1$gene_padj) & d1$gene_padj <= 0.05])
  cat(sprintf("     %s: %d genes called, %d planted of %d\n", engine_label(eng), length(sig), sum(sig %in% gid19[dtu_true]), length(dtu_true)))
  check(paste(engine_label(eng), "finds planted switches and few others"),
        if (eng == "saturn") length(sig) == 0 || mean(sig %in% gid19[dtu_true]) >= 0.7
        else sum(sig %in% gid19[dtu_true]) >= 10 && mean(sig %in% gid19[dtu_true]) >= 0.8)
  check(paste(engine_label(eng), "reversed comparison: same p-values, opposite share change"), {
    m <- match(d1$gene_id, d2$gene_id)
    !anyNA(m) && isTRUE(all.equal(d1$pvalue, d2$pvalue[m], tolerance = 1e-5)) &&
      isTRUE(all.equal(d1$delta_prop, -d2$delta_prop[m], tolerance = 1e-8)) &&
      cor(d1$log2FoldChange, -d2$log2FoldChange[m], use = "complete.obs") > 0.99 })
  if (requireNamespace("stageR", quietly = TRUE))
    check(paste(engine_label(eng), "stage-wise DTU columns: transcripts confirmed only in genes that pass the screen"),
          all(c("stageR_gene_padj", "stageR_tx_padj") %in% colnames(d1)) &&
            all(is.na(d1$stageR_tx_padj[!is.na(d1$stageR_gene_padj) & d1$stageR_gene_padj > 0.05])) &&
            (!any(d1$stageR_gene_padj <= 0.05, na.rm = TRUE) || any(!is.na(d1$stageR_tx_padj))))
  check(paste(engine_label(eng), "usage plot, QC figures and run record build"), {
    p <- dtu_usage_plot(r, names(r$results)[1], dtu_gene_choices(r, names(r$results)[1])[[1]])
    q <- lapply(c("library", "pca", "pvalues"), function(w) build_qc_figure(w, r))
    rec <- de_run_record(r)
    inherits(p, "ggplot") && any(grepl("differential transcript usage", rec)) && any(grepl("Genes with DTU", rec)) })
}
if (!is.null(dtu_runs$edger_dtu)) check("edgeR diffSplice log2FoldChange is on the log2 scale", {
  d <- dtu_runs$edger_dtu$results[[1]]
  if (exists("diffSplice.DGEGLM", envir = asNamespace("edgeR"), inherits = FALSE)) {
    # edgeR >= 4.6: relative to the gene's other transcripts, so for two-transcript
    # genes it is the log2 odds ratio of the shares
    two <- names(which(table(d$parent_gene_id) == 2))
    s <- d$parent_gene_id %in% two & d$padj < 0.01 & pmin(d$prop_reference, d$prop_target) > 0.05 &
      pmax(d$prop_reference, d$prop_target) < 0.95
    x <- log2(d$prop_target[s] / (1 - d$prop_target[s])) - log2(d$prop_reference[s] / (1 - d$prop_reference[s]))
  } else {
    s <- d$padj < 0.01 & d$prop_reference > 0.05 & d$prop_target > 0.05
    x <- log2(d$prop_target[s] / d$prop_reference[s])
  }
  b <- coef(lm(d$log2FoldChange[s] ~ x))[2]
  cat(sprintf("     slope %.3f over %d transcripts\n", b, sum(s)))
  sum(s) >= 6 && b > 0.8 && b < 1.25 })
if (is.null(de_engine_missing("edger_dtu"))) check("a DTU run without a transcript-to-gene map stops with an explanation", {
  e <- tryCatch(run_de_analysis(cnt19, ss19, ss19$sample, "condition", character(0), ty19, cons19[1, ], "edger_dtu",
                                list(feature_level = "dtu"), transcripts = list(tx2gene = NULL)), error = function(e) conditionMessage(e))
  grepl("transcript-to-gene map", e) })
if (length(dtu_runs)) check("DTU results reload as a transcript table and map to unique genes for enrichment", {
  r <- dtu_runs[[1]]; st <- standardise_de_table(de_result_for_export(r, names(r$results)[1]))$data
  g <- enrichment_ids(st, st$gene_id[1:40])
  "parent_gene_id" %in% colnames(st) && !anyDuplicated(g) && all(g %in% gid19) })
for (eng in names(dtu_runs)) {
  rr <- dtu_runs[[eng]]
  z <- tempfile(fileext = ".zip"); write_de_bundle(rr, z, qc = FALSE)
  ex <- tempfile(); dir.create(ex); utils::unzip(z, exdir = ex)
  owd <- setwd(ex); st <- system2(file.path(R.home("bin"), "Rscript"), "reproduce.R", stdout = FALSE, stderr = FALSE); setwd(owd)
  check(paste(engine_label(eng), "DTU script reproduces every comparison"), st == 0 && all(vapply(names(rr$results), function(id) {
    f <- file.path(ex, "reproduced", paste0(id, ".dtu.results.tsv"))
    if (!file.exists(f)) return(FALSE)
    b <- utils::read.delim(f, check.names = FALSE); a <- rr$results[[id]]; m <- match(a$gene_id, b$transcript_id)
    num <- function(x) suppressWarnings(as.numeric(x))
    !anyNA(m) && nrow(a) == nrow(b) &&
      all(vapply(intersect(c("log2FoldChange", "pvalue", "padj", "gene_padj", "delta_prop", "stageR_tx_padj"), colnames(a)),
                 function(cl) isTRUE(all.equal(a[[cl]], num(b[[cl]][m]), tolerance = 1e-8)), logical(1)))
  }, logical(1))))
}

# runApp() would otherwise source R/*.R alphabetically before app.R, and
# annotation_sources.R needs organisms.R first.
check("Shiny's R/ autoloading is switched off (app.R sources the modules in order)",
      file.exists(file.path(root, "R", "_disable_autoload.R")))

cat(sprintf("\nAll %d checks passed.\n", ok))
