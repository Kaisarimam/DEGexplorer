# =============================================================================
# make_example_de_tables.R -- rebuild the bundled example result tables from the
# bundled example counts (data/example_rsem_gene_counts.tsv.gz, its gene lengths
# and data/example_samplesheet.csv). Run from the repository root:
#
#   Rscript data-raw/make_example_de_tables.R
#
# Writes, for the comparison stim vs mock with donor as a covariate:
#   data/example_deseq2_results.tsv.gz         DESeq2 (Wald, ashr-shrunken fold changes,
#                                              gene-length offsets), annotated
#   data/example_limma_results.tsv.gz          limma-voom, annotated
#   data/example_limma_results_plain.tsv.gz    limma-voom, results only (no symbols)
#   data/example_limma_results_filtered.tsv.gz limma-voom, DEGs only (adj. P <= 0.05, |logFC| >= 1)
#   data/example_vst.tsv.gz                    DESeq2 VST matrix (blind)
# Genes are kept with at least 10 counts in at least one sample.
# The annotation columns come from data/example_annotation.tsv.gz (Ensembl GTF fields).
# =============================================================================
suppressMessages({
  for (f in c("utils_io.R", "utils_normalise.R", "de_analysis.R")) source(file.path("R", f))
  library(DESeq2)
})
cnt  <- read_matrix_auto("data/example_rsem_gene_counts.tsv.gz")
len  <- read_matrix_auto("data/example_rsem_gene_lengths.tsv.gz")$matrix
meta <- read.csv("data/example_samplesheet.csv", colClasses = "character")
ann  <- read.delim(gzfile("data/example_annotation.tsv.gz"), check.names = FALSE, colClasses = "character")

ss <- design_scaffold(colnames(cnt$matrix))
for (v in c("condition", "donor")) ss[[v]] <- meta[[v]][match(ss$sample, meta$sample)]
types <- list(condition = "factor", donor = "factor")
cons <- data.frame(id = "condition_stim_vs_mock", variable = "condition",
                   reference = "mock", target = "stim", stringsAsFactors = FALSE)
shrink <- if (requireNamespace("ashr", quietly = TRUE)) "ashr" else {
  message("ashr is not installed: using normal shrinkage"); "normal" }
opts <- modifyList(recommended_options(), list(shrink = shrink))

ds <- run_de_analysis(cnt$matrix, ss, ss$sample, "condition", "donor", types, cons, "deseq2",
                      opts, lengths = len)
lm <- run_de_analysis(cnt$matrix, ss, ss$sample, "condition", "donor", types, cons, "limma",
                      modifyList(opts, list(length_offset = FALSE)))

annotate <- function(df) {
  a <- ann[match(df$gene_id, ann$gene_id), setdiff(colnames(ann), "gene_id"), drop = FALSE]
  # genes missing from the annotation file take the counts file's symbol, if it has a real one
  miss <- is.na(a$gene_name) | a$gene_name %in% c("", "NA")
  nm <- cnt$annotation$gene_name[match(df$gene_id[miss], cnt$annotation$gene_id)]
  nm[!is.na(nm) & nm == df$gene_id[miss]] <- NA
  a$gene_name[miss] <- nm
  cbind(df, a)
}
write_gz <- function(df, f) {
  con <- gzfile(f, "w"); on.exit(close(con))
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
d <- ds$results[[1]]; d <- d[order(d$gene_id), ]
write_gz(annotate(d[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]),
         "data/example_deseq2_results.tsv.gz")
l <- lm$results[[1]]; l <- l[order(l$gene_id), ]
lt <- data.frame(gene_id = l$gene_id, logFC = l$log2FoldChange, AveExpr = l$AveExpr, t = l$stat,
                 P.Value = l$pvalue, adj.P.Val = l$padj, B = l$B)
write_gz(annotate(lt), "data/example_limma_results.tsv.gz")
write_gz(lt, "data/example_limma_results_plain.tsv.gz")
write_gz(lt[lt$adj.P.Val <= 0.05 & abs(lt$logFC) >= 1, ], "data/example_limma_results_filtered.tsv.gz")
v <- ds$transformed[order(rownames(ds$transformed)), ]
write_gz(data.frame(gene_id = rownames(v), v, check.names = FALSE), "data/example_vst.tsv.gz")

deg <- function(lfc, p) c(up = sum(p <= .05 & lfc >= 1, na.rm = TRUE), down = sum(p <= .05 & lfc <= -1, na.rm = TRUE))
cat(sprintf("DESeq2 (%s): %d genes tested; DEGs up %d / down %d\n", shrink, nrow(d),
            deg(d$log2FoldChange, d$padj)[1], deg(d$log2FoldChange, d$padj)[2]))
cat(sprintf("limma-voom: DEGs up %d / down %d\n", deg(lt$logFC, lt$adj.P.Val)[1], deg(lt$logFC, lt$adj.P.Val)[2]))
