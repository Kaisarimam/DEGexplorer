# =============================================================================
# transcripts.R -- inputs for transcript (isoform)-level differential expression
#
# Transcript-level analysis needs, besides the counts:
#   - the transcript -> gene map   (symbols, EBSeq's Ng vector, stage-wise testing,
#                                   and gene-level enrichment of the results)
#   - effective transcript lengths (DESeq2 / edgeR offsets, as tximport txOut)
#   - quantification uncertainty   (edgeR / limma: counts divided by the
#                                   read-to-transcript ambiguity overdispersion,
#                                   Baldoni et al. 2024, estimated from Salmon /
#                                   kallisto bootstrap or Gibbs resamples)
#   - optionally an Ng vector      (EBSeq; rsem-generate-ngvector)
# They come from a merged transcript count table, from the quantification
# folders themselves (Salmon, kallisto, RSEM), from a GTF, or from small files.
# =============================================================================

#' Split GENCODE-style pipe ids ("ENST...|ENSG...|OTTHUMG|OTTHUMT|NAME-201|NAME|len|biotype|").
#' @return data.frame(transcript_id, gene_id, transcript_name, gene_name); NA when absent
split_tx_ids <- function(ids) {
  ids <- as.character(ids)
  piped <- grepl("|", ids, fixed = TRUE)
  out <- data.frame(transcript_id = ids, gene_id = NA_character_, transcript_name = NA_character_,
                    gene_name = NA_character_, stringsAsFactors = FALSE)
  if (any(piped)) {
    parts <- strsplit(ids[piped], "|", fixed = TRUE)
    get <- function(k) vapply(parts, function(p) if (length(p) >= k && nzchar(p[k])) p[k] else NA_character_, "")
    out$transcript_id[piped] <- get(1); out$gene_id[piped] <- get(2)
    out$transcript_name[piped] <- get(5); out$gene_name[piped] <- get(6)
  }
  out
}

TX_ID_COLUMNS   <- c("tx", "transcript_id", "transcript", "target_id", "name", "txid", "tx_id", "isoform_id")
GENE_ID_COLUMNS <- c("gene_id", "gene", "geneid", "gene_ids")

#' Read a merged transcript count table (Salmon / RSEM merged tables, or any
#' transcript x sample matrix with optional gene columns).
#' @return list(matrix, tx2gene) with tx2gene = data.frame(transcript_id, gene_id,
#'   gene_name, transcript_name) or NULL when the file carries no gene column
read_transcript_matrix <- function(path) {
  df <- read_table_auto(path)
  key <- tolower(gsub("[^A-Za-z0-9_]", "", colnames(df)))
  num <- vapply(df, function(x) is.numeric(x) || all(suppressWarnings(!is.na(as.numeric(x[!is.na(x)])))), logical(1))
  id_i <- match(TX_ID_COLUMNS, key); id_i <- id_i[!is.na(id_i)][1]
  if (is.na(id_i)) id_i <- which(!num)[1]
  if (is.na(id_i)) id_i <- 1L
  gene_i <- match(GENE_ID_COLUMNS, key); gene_i <- gene_i[!is.na(gene_i) & gene_i != id_i][1]
  name_i <- match(c("gene_name", "genename", "symbol"), key); name_i <- name_i[!is.na(name_i)][1]
  sample_cols <- setdiff(which(num), c(id_i, gene_i, name_i))
  if (!length(sample_cols)) stop("No numeric sample columns found in ", basename(path), ".")
  sp <- split_tx_ids(df[[id_i]])
  m <- as.matrix(as.data.frame(lapply(df[sample_cols], function(x) suppressWarnings(as.numeric(x)))))
  colnames(m) <- colnames(df)[sample_cols]
  rownames(m) <- sp$transcript_id
  if (anyDuplicated(rownames(m))) stop("Transcript ids are not unique in ", basename(path), ".")
  t2g <- sp
  if (!is.na(gene_i)) t2g$gene_id <- ifelse(is.na(t2g$gene_id), as.character(df[[gene_i]]), t2g$gene_id)
  if (!is.na(name_i)) t2g$gene_name <- ifelse(is.na(t2g$gene_name), as.character(df[[name_i]]), t2g$gene_name)
  list(matrix = m, tx2gene = if (all(is.na(t2g$gene_id))) NULL else t2g)
}

#' Read a transcript -> gene table (two or more columns; a GTF also works).
read_tx2gene <- function(path) {
  if (grepl("\\.(gtf|gff|gff3)(\\.gz)?$", path, ignore.case = TRUE)) {
    t <- gtf_transcript_table(path)
    return(t[, c("transcript_id", "gene_id", "gene_name", "transcript_name"), drop = FALSE])
  }
  df <- read_table_auto(path)
  key <- tolower(gsub("[^A-Za-z0-9_]", "", colnames(df)))
  ti <- match(TX_ID_COLUMNS, key); ti <- ti[!is.na(ti)][1]; if (is.na(ti)) ti <- 1L
  gi <- match(GENE_ID_COLUMNS, key); gi <- gi[!is.na(gi) & gi != ti][1]; if (is.na(gi)) gi <- setdiff(seq_along(df), ti)[1]
  ni <- match(c("gene_name", "genename", "symbol"), key); ni <- ni[!is.na(ni)][1]
  tn <- match(c("transcript_name", "txname", "tx_name"), key); tn <- tn[!is.na(tn)][1]
  if (is.na(gi)) stop("A transcript-to-gene table needs at least two columns (transcript, gene).")
  data.frame(transcript_id = sub("\\|.*$", "", as.character(df[[ti]])), gene_id = as.character(df[[gi]]),
             gene_name = if (!is.na(ni)) as.character(df[[ni]]) else NA_character_,
             transcript_name = if (!is.na(tn)) as.character(df[[tn]]) else NA_character_,
             stringsAsFactors = FALSE)
}

#' Combine transcript -> gene sources; earlier sources win, later ones fill gaps.
merge_tx2gene <- function(ids, ...) {
  out <- data.frame(transcript_id = ids, gene_id = NA_character_, gene_name = NA_character_,
                    transcript_name = NA_character_, transcript_biotype = NA_character_,
                    stringsAsFactors = FALSE)
  for (src in list(...)) {
    if (is.null(src) || !nrow(src)) next
    i <- match(sub("\\.[0-9]+$", "", ids), sub("\\.[0-9]+$", "", src$transcript_id))
    i2 <- match(ids, src$transcript_id); i[!is.na(i2)] <- i2[!is.na(i2)]
    for (cl in intersect(c("gene_id", "gene_name", "transcript_name", "transcript_biotype"), colnames(src))) {
      v <- as.character(src[[cl]])[i]
      fill <- is.na(out[[cl]]) & !is.na(v) & nzchar(v) & v != "NA"
      out[[cl]][fill] <- v[fill]
    }
  }
  out
}

# ---- quantification folders ---------------------------------------------------------

#' Find per-sample quantification folders under a directory (Salmon quant.sf,
#' kallisto abundance.tsv/.h5, RSEM *.isoforms.results).
find_quant_dirs <- function(root) {
  if (!dir.exists(root)) stop("Folder not found: ", root)
  f <- list.files(root, pattern = "^(quant\\.sf|abundance\\.tsv|abundance\\.h5)$|\\.isoforms\\.results$",
                  recursive = TRUE, full.names = TRUE)
  f <- f[!grepl("/(aux_info|libParams|logs)/", f)]
  if (!length(f)) stop("No Salmon (quant.sf), kallisto (abundance.tsv) or RSEM (*.isoforms.results) files under ", root, ".")
  tool <- ifelse(basename(f) == "quant.sf", "salmon",
          ifelse(grepl("^abundance", basename(f)), "kallisto", "rsem"))
  d <- data.frame(dir = dirname(f), file = f, tool = tool, stringsAsFactors = FALSE)
  # one entry per folder; kallisto: prefer abundance.tsv over .h5
  d <- d[order(d$dir, -grepl("\\.tsv$|quant\\.sf$|results$", d$file)), , drop = FALSE]
  d <- d[!duplicated(d$dir), , drop = FALSE]
  if (length(unique(d$tool)) > 1) stop("The folder mixes outputs of ", paste(unique(d$tool), collapse = " and "), ".")
  d$sample <- ifelse(d$tool == "rsem", sub("\\.isoforms\\.results$", "", basename(d$file)), basename(d$dir))
  if (anyDuplicated(d$sample)) d$sample <- make.unique(d$sample)
  d
}

#' Read per-sample quantifications into counts, effective lengths and (when
#' resamples exist) the edgeR overdispersion of each transcript.
#'
#' @return list(counts, lengths, overdispersion (named numeric or NULL), tx2gene,
#'   tool, note)
read_quant_dirs <- function(root) {
  d <- find_quant_dirs(root)
  tool <- d$tool[1]
  rd <- function(f) as.data.frame(readr::read_tsv(f, progress = FALSE, show_col_types = FALSE))
  tabs <- lapply(d$file[!grepl("\\.h5$", d$file)], rd)
  if (length(tabs) != nrow(d)) stop("kallisto .h5 files need abundance.tsv beside them (run kallisto with --plaintext or keep both).")
  col <- switch(tool,
    salmon   = list(id = "Name", count = "NumReads", len = "EffectiveLength"),
    kallisto = list(id = "target_id", count = "est_counts", len = "eff_length"),
    rsem     = list(id = "transcript_id", count = "expected_count", len = "effective_length"))
  ids <- as.character(tabs[[1]][[col$id]])
  for (t in tabs[-1]) if (!identical(as.character(t[[col$id]]), ids))
    stop("The samples were not quantified against the same transcripts.")
  counts <- sapply(tabs, function(t) as.numeric(t[[col$count]]))
  lens   <- sapply(tabs, function(t) as.numeric(t[[col$len]]))
  sp <- split_tx_ids(ids)
  dimnames(counts) <- dimnames(lens) <- list(sp$transcript_id, d$sample)
  t2g <- sp
  if (tool == "rsem") t2g$gene_id <- as.character(tabs[[1]]$gene_id)

  # Read-to-transcript ambiguity: edgeR's catch* functions estimate it from the resamples.
  od <- NULL; note <- NULL
  od_raw <- tryCatch(switch(tool,
      salmon   = edgeR::catchSalmon(d$dir, verbose = FALSE)$annotation$Overdispersion,
      kallisto = edgeR::catchKallisto(d$dir, verbose = FALSE)$annotation$Overdispersion,
      rsem     = if (exists("catchRSEM", envir = asNamespace("edgeR")))
                   getExportedValue("edgeR", "catchRSEM")(d$file, verbose = FALSE)$annotation$Overdispersion
                 else stop("this edgeR has no catchRSEM()")),
    error = function(e) structure(conditionMessage(e), class = "od_error"))
  if (inherits(od_raw, "od_error") || all(is.na(od_raw))) {
    note <- paste0("No quantification resamples found (", if (inherits(od_raw, "od_error")) unclass(od_raw) else "0 resamples",
                   "): counts cannot be divided by their overdispersion. Rerun ",
                   switch(tool, salmon = "Salmon with --numGibbsSamples 50 (or --numBootstraps)",
                          kallisto = "kallisto with -b 50", rsem = "RSEM with --calc-pme / Gibbs sampling"),
                   " to account for read-to-transcript ambiguity in edgeR / limma.")
  } else {
    od <- stats::setNames(pmax(as.numeric(od_raw), 1), sp$transcript_id)
  }
  list(counts = counts, lengths = lens, overdispersion = od, tx2gene = t2g, tool = tool,
       n_samples = nrow(d), note = note)
}

#' Read an Ng vector for EBSeq: two columns (transcript, Ng), or one column in the
#' order of the count matrix rows (as rsem-generate-ngvector writes it).
read_ngvector <- function(path, ids) {
  x <- utils::read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(x) >= 2 && !is.numeric(x[[1]])) {
    if (!is.numeric(x[[2]])) x <- x[-1, , drop = FALSE]          # header row
    v <- stats::setNames(as.integer(x[[2]]), sub("\\|.*$", "", x[[1]]))
    out <- v[ids]
    if (mean(!is.na(out)) < 0.9) stop("Fewer than 90% of the transcripts have an Ng value in this file.")
    return(out)
  }
  v <- suppressWarnings(as.integer(x[[1]]))
  v <- v[!is.na(v)]
  if (length(v) != length(ids)) stop(sprintf(
    "The Ng vector has %d values but the count matrix has %d transcripts; a one-column file must follow the matrix's row order (the RSEM reference order).",
    length(v), length(ids)))
  stats::setNames(v, ids)
}
