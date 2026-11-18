# =============================================================================
# utils_io.R -- file reading, format auto-detection, DE-result adapters
#
# The app accepts differential-expression tables from DESeq2, limma(-voom),
# edgeR and generic pipelines. Rather than demanding one fixed column naming
# scheme, every table is passed through `standardise_de_table()`, which maps
# whichever dialect it finds onto one canonical set of columns:
#
#   gene_id, gene_name, display_gene, baseMean, log2FoldChange, pvalue, padj
#
# Everything downstream (plots, heatmaps, enrichment) speaks only canonical.
# =============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- generic delimited reader ----------------------------------------------

#' Read a delimited text file, guessing the separator from the first line.
#' Handles .tsv/.txt/.csv and gzipped versions of each.
read_table_auto <- function(path, n_max = Inf) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  first <- tryCatch(readLines(con, n = 1L, warn = FALSE), finally = close(con))
  if (length(first) == 0 || !nzchar(first[1])) stop("File appears to be empty: ", basename(path))

  n_tab <- lengths(regmatches(first, gregexpr("\t", first)))
  n_com <- lengths(regmatches(first, gregexpr(",",  first)))
  n_sem <- lengths(regmatches(first, gregexpr(";",  first)))
  delim <- c("\t", ",", ";")[which.max(c(n_tab, n_com, n_sem))]
  if (max(n_tab, n_com, n_sem) == 0) delim <- "\t"

  df <- readr::read_delim(
    path, delim = delim, show_col_types = FALSE, progress = FALSE,
    n_max = n_max, na = c("", "NA", "NaN", "null", "None", "."),
    name_repair = "minimal"
  )
  df <- df[, colnames(df) != "" & !is.na(colnames(df)), drop = FALSE]
  if (ncol(df) < 2) {
    stop("Only one column was parsed from ", basename(path),
         " -- the file may use an unusual delimiter.")
  }
  as.data.frame(df, stringsAsFactors = FALSE)
}

# ---- DE-result dialect detection -------------------------------------------

# Candidate column names per canonical field, in priority order.
# Matching is case-insensitive and ignores dots/underscores/spaces, so
# "adj.P.Val", "adj_p_val" and "AdjPVal" all resolve to the same field.
DE_SYNONYMS <- list(
  gene_id = c("gene_id", "geneid", "ensembl", "ensemblid", "ensemblgeneid",
              "id", "row", "rowname", "rownames", "feature", "featureid",
              "gene", "probe", "probeid", "target_id",
              # transcript-level tables (only when no gene id column exists)
              "transcript_id", "transcriptid", "tx", "txid", "transcript", "isoform_id", "isoform"),
  gene_name = c("gene_name", "genename", "symbol", "genesymbol", "hgncsymbol",
                "externalgenename", "name", "gene_symbol"),
  log2FoldChange = c("log2foldchange", "logfc", "log2fc", "log2foldchanges",
                     "log2ratio", "logratio", "fold_change_log2", "lfc",
                     "log2fold_change", "coefficient", "estimate"),
  baseMean = c("basemean", "aveexpr", "logcpm", "averageexpression",
               "meanexpression", "baseexpr", "avgexpr", "abundance",
               "meanexpr", "expression", "amean"),
  pvalue = c("pvalue", "p.value", "pval", "p", "praw", "rawp",
             "pvalues", "p_value", "nominalp"),
  padj = c("padj", "adj.p.val", "fdr", "qvalue", "q", "adjp", "adjpval",
           "adjustedpvalue", "padjusted", "bh", "fdr_bh", "svalue",
           "adj_pvalue", "q_value")
)

# Column sets written by nf-core/differentialabundance, used to recognise which
# of its outputs has been handed to the app.
NFCORE_LIMMA_CORE  <- c("geneid", "logfc", "aveexpr", "t", "pvalue", "adjpval", "b")
NFCORE_DESEQ2_CORE <- c("geneid", "basemean", "log2foldchange", "pvalue", "padj")
# Columns the GTF-based annotation step appends to the *.annotated.tsv flavour.
NFCORE_ANNOTATION_COLS <- c("chromosome", "start", "end", "width", "source", "type",
                            "score", "phase", "geneversion", "transcriptid",
                            "transcriptversion", "genesource", "genebiotype",
                            "transcriptsource", "transcriptbiotype", "tag",
                            "transcriptsupportlevel", "genename", "transcriptname",
                            "ccdsid")

.norm_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))

#' Work out which column of `df` supplies each canonical field.
#' Returns a named list (canonical field -> actual column name or NA).
detect_de_columns <- function(df) {
  keys <- .norm_key(colnames(df))
  out <- list()
  used <- character(0)
  for (field in names(DE_SYNONYMS)) {
    hit <- NA_character_
    for (cand in .norm_key(DE_SYNONYMS[[field]])) {
      idx <- which(keys == cand & !(colnames(df) %in% used))
      if (length(idx) > 0) { hit <- colnames(df)[idx[1]]; break }
    }
    # gene_id fallback: first character-like column that looks like an accession
    if (is.na(hit) && field == "gene_id") {
      for (j in seq_along(df)) {
        v <- df[[j]]
        if ((is.character(v) || is.factor(v)) && !(colnames(df)[j] %in% used)) {
          if (any(grepl("^(ENS[A-Z]*G|ENS[A-Z]*T)[0-9]", as.character(v)[1:min(50, nrow(df))]))) {
            hit <- colnames(df)[j]; break
          }
        }
      }
      if (is.na(hit) && (is.character(df[[1]]) || is.factor(df[[1]]))) hit <- colnames(df)[1]
    }
    if (!is.na(hit)) used <- c(used, hit)
    out[[field]] <- hit
  }
  out
}

#' Name the pipeline a DE table most likely came from (for display only).
guess_de_pipeline <- function(df) {
  keys <- .norm_key(colnames(df))
  if (all(c("log2foldchange", "basemean") %in% keys)) return("DESeq2")
  if (all(c("logfc", "aveexpr") %in% keys) || "adjpval" %in% keys ||
      "adj.p.val" %in% colnames(df))                  return("limma / limma-voom")
  if (all(c("logfc", "logcpm") %in% keys))            return("edgeR")
  "generic / unknown"
}

#' Recognise which nf-core/differentialabundance output a table is, if any.
#'
#' The pipeline emits, per contrast and per method:
#'   <contrast>.limma.results.tsv            -- every gene tested
#'   <contrast>.limma.results_filtered.tsv   -- significant genes only
#'   <contrast>.limma.annotated.tsv          -- every gene plus GTF annotation
#' and the equivalents for deseq2. All three limma flavours share the same
#' statistics columns, so the flavour is told apart by annotation columns and by
#' the file name, and the filtered flavour additionally by its p-value profile.
#'
#' @return list(pipeline, flavour, nfcore, annotated, filtered, label, notes)
describe_de_table <- function(df, filename = NULL) {
  keys  <- .norm_key(colnames(df))
  fname <- tolower(basename(filename %||% ""))

  is_limma  <- all(NFCORE_LIMMA_CORE  %in% keys)
  is_deseq2 <- all(NFCORE_DESEQ2_CORE %in% keys)
  nfcore    <- is_limma || is_deseq2

  method <- if (is_limma) "limma" else if (is_deseq2) "DESeq2" else guess_de_pipeline(df)
  annotated <- sum(NFCORE_ANNOTATION_COLS %in% keys) >= 5 || grepl("annotated", fname, fixed = TRUE)

  # Filtered tables are short and carry no non-significant rows.
  padj_col <- detect_de_columns(df)$padj
  padj_vals <- if (!is.na(padj_col)) suppressWarnings(as.numeric(df[[padj_col]])) else NA_real_
  padj_vals <- padj_vals[is.finite(padj_vals)]
  looks_filtered <- length(padj_vals) > 0 &&
    nrow(df) < 5000 && max(padj_vals) <= 0.25
  filtered <- grepl("filtered", fname, fixed = TRUE) || looks_filtered

  flavour <- if (filtered) "results_filtered" else if (annotated) "annotated" else "results"

  label <- if (nfcore) {
    sprintf("%s results (%s)", method, flavour)
  } else method

  notes <- character(0)
  if (filtered) {
    notes <- c(notes, sprintf(
      paste("This looks like a pre-filtered results table (%s rows, all adjusted p <= %.3g).",
            "Volcano / MA plots will show only the genes that survived the pipeline's own cut-off,",
            "and enrichment run against it would use those genes as their own background.",
            "Upload the unfiltered '*.%s.results.tsv' as the enrichment background to fix that."),
      format(nrow(df), big.mark = ","), max(padj_vals),
      if (is_deseq2) "deseq2" else "limma"))
  }
  if (nfcore && !annotated && !any(c("genename", "symbol") %in% keys)) {
    notes <- c(notes, paste(
      "No gene symbols in this file -- plots and reports will fall back to Ensembl ids.",
      "Add the matching '*.annotated.tsv' as a gene-annotation file to get symbols."))
  }

  list(pipeline = method, flavour = flavour, nfcore = nfcore,
       annotated = annotated, filtered = filtered, label = label, notes = notes)
}

#' Convert any supported DE table into the canonical schema.
#'
#' @param df         data frame straight from `read_table_auto()`
#' @param overrides  optional named list to force particular column choices,
#'                   e.g. list(padj = "FDR"); supplied by the UI mapper.
#' @return list(data = <canonical data frame>, mapping = <list>,
#'              pipeline = <chr>, notes = <chr vector>)
standardise_de_table <- function(df, overrides = list()) {
  mapping <- detect_de_columns(df)
  for (nm in names(overrides)) {
    v <- overrides[[nm]]
    if (!is.null(v) && !identical(v, "") && !identical(v, "(auto)") && v %in% colnames(df)) {
      mapping[[nm]] <- v
    }
  }
  notes <- character(0)

  need <- c("gene_id", "log2FoldChange", "pvalue")
  miss <- need[vapply(need, function(f) is.na(mapping[[f]]), logical(1))]
  if (length(miss) > 0) {
    stop("Could not identify required column(s): ", paste(miss, collapse = ", "),
         ".\nColumns found in the file: ", paste(colnames(df), collapse = ", "),
         "\nUse the 'Column mapping' controls to point the app at the right columns.")
  }

  # Transcript-level table whose gene_id repeats: the transcript column is the row id
  # and the gene becomes parent_gene_id (used to map transcripts to genes for enrichment).
  if (anyDuplicated(df[[mapping$gene_id]]) && !"parent_gene_id" %in% colnames(df)) {
    k <- .norm_key(colnames(df))
    txc <- colnames(df)[match(.norm_key(c("transcript_id", "tx", "txid", "target_id", "isoform_id", "transcript")), k)]
    txc <- txc[!is.na(txc) & txc != mapping$gene_id][1]
    if (!is.na(txc) && !anyDuplicated(df[[txc]])) {
      df$parent_gene_id <- as.character(df[[mapping$gene_id]])
      mapping$gene_id <- txc
      notes <- c(notes, sprintf("Transcript-level table: rows are identified by '%s'; its gene column is kept as parent_gene_id.", txc))
    }
  }

  out <- data.frame(
    gene_id        = as.character(df[[mapping$gene_id]]),
    log2FoldChange = suppressWarnings(as.numeric(df[[mapping$log2FoldChange]])),
    pvalue         = suppressWarnings(as.numeric(df[[mapping$pvalue]])),
    stringsAsFactors = FALSE
  )

  out$padj <- if (!is.na(mapping$padj)) {
    suppressWarnings(as.numeric(df[[mapping$padj]]))
  } else {
    notes <- c(notes, "No adjusted p-value column found; Benjamini-Hochberg FDR was computed from the raw p-values.")
    p.adjust(out$pvalue, method = "BH")
  }

  out$baseMean <- if (!is.na(mapping$baseMean)) {
    bm <- suppressWarnings(as.numeric(df[[mapping$baseMean]]))
    # limma AveExpr / edgeR logCPM are already on a log scale; put them back on
    # a linear scale so the MA plot's log10(x+1) axis stays meaningful.
    key <- .norm_key(mapping$baseMean)
    if (key %in% c("aveexpr", "logcpm") ||
        (all(is.finite(bm)) && max(bm, na.rm = TRUE) < 40 && any(bm < 0, na.rm = TRUE))) {
      notes <- c(notes, sprintf("Column '%s' looks log-scaled; converted to a linear scale for expression axes.",
                                mapping$baseMean))
      2^bm
    } else bm
  } else {
    notes <- c(notes, "No mean-expression column found; MA-plot x-axis will be unavailable unless a count matrix is loaded.")
    rep(NA_real_, nrow(out))
  }

  out$gene_name <- if (!is.na(mapping$gene_name)) {
    as.character(df[[mapping$gene_name]])
  } else rep(NA_character_, nrow(out))

  # carry through any extra annotation columns that are not already mapped
  taken <- unlist(mapping[!vapply(mapping, is.na, logical(1))], use.names = FALSE)
  extra <- setdiff(colnames(df), taken)
  if (length(extra) > 0) out <- cbind(out, df[, extra, drop = FALSE])

  out <- out[!is.na(out$gene_id) & out$gene_id != "", , drop = FALSE]
  # collapse duplicated gene ids, keeping the most significant row
  if (anyDuplicated(out$gene_id)) {
    n_dup <- sum(duplicated(out$gene_id))
    ord <- order(out$pvalue, na.last = TRUE)
    out <- out[ord, , drop = FALSE]
    out <- out[!duplicated(out$gene_id), , drop = FALSE]
    notes <- c(notes, sprintf("%d duplicated gene id(s) collapsed, keeping the most significant row for each.", n_dup))
  }

  out$display_gene <- ifelse(
    is.na(out$gene_name) | trimws(out$gene_name) %in% c("", "NA", "-"),
    out$gene_id, out$gene_name
  )
  # transcripts: label by transcript name (IFIT1-201) where there is one
  if ("parent_gene_id" %in% colnames(out) && "transcript_name" %in% colnames(out)) {
    tn <- as.character(out$transcript_name)
    ok <- !is.na(tn) & nzchar(tn) & tn != "NA"
    out$display_gene[ok] <- tn[ok]
  }
  rownames(out) <- NULL

  list(data = out, mapping = mapping, pipeline = guess_de_pipeline(df), notes = notes)
}

#' The bundled example datasets offered on the Data tab.
EXAMPLE_DATASETS <- c(
  "DESeq2 \u2014 annotated results"                 = "deseq2",
  "limma \u2014 annotated results"                  = "limma_annotated",
  "limma \u2014 plain results.tsv (no symbols)"     = "limma_results",
  "limma \u2014 results_filtered.tsv (DEGs only)"   = "limma_filtered"
)

# ---- GTF / GFF annotation ---------------------------------------------------

#' Open a text file, transparently handling gzip.
.open_text_conn <- function(path) {
  magic <- tryCatch(readBin(path, "raw", 2L), error = function(e) raw(0))
  if (length(magic) == 2L && magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)) {
    gzfile(path, "rt")
  } else file(path, "rt")
}

#' First non-comment, non-empty line of a file.
.first_data_line <- function(path, max_scan = 200L) {
  con <- .open_text_conn(path); on.exit(close(con), add = TRUE)
  for (i in seq_len(max_scan)) {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) return(NA_character_)
    if (nzchar(trimws(ln)) && !startsWith(ln, "#")) return(ln)
  }
  NA_character_
}

#' Is this a GTF / GFF annotation rather than a delimited table?
#'
#' Decided from the content, because Shiny stores uploads under a temporary
#' name and the original extension cannot be relied on.
is_gtf_file <- function(path, filename = NULL) {
  nm <- tolower(basename(filename %||% path))
  if (grepl("\\.(gtf|gff|gff3)(\\.gz|\\.bgz)?$", nm)) return(TRUE)
  ln <- .first_data_line(path)
  if (is.na(ln)) return(FALSE)
  f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  length(f) == 9L && grepl('gene_id |gene_id=|ID=|gene_name |gene=', f[9])
}

#' GFF3 writes attributes as key=value, GTF as key "value".
.is_gff3 <- function(attrs) {
  s <- utils::head(attrs[nzchar(attrs)], 50)
  if (!length(s)) return(FALSE)
  mean(grepl("^[^;=\\s]+=", trimws(s))) > 0.5
}

#' Which attribute keys does this file actually carry?
.gtf_attr_keys <- function(attrs, gff3, n_sample = 500L) {
  toks <- unlist(strsplit(utils::head(attrs, n_sample), ";", fixed = TRUE))
  toks <- trimws(toks); toks <- toks[nzchar(toks)]
  keys <- if (gff3) sub("=.*$", "", toks) else sub("[ \t].*$", "", toks)
  keys <- keys[grepl("^[A-Za-z_][A-Za-z0-9_.:-]*$", keys)]
  names(sort(table(keys), decreasing = TRUE))
}

#' Pull one attribute out of a GTF / GFF attribute column.
#' Uses a single perl regexpr pass with a capture group, which stays fast
#' across the ~60,000 gene rows of a human annotation.
.gtf_attr <- function(attrs, key, gff3) {
  pat <- if (gff3) paste0("(?:^|;)[ \t]*", key, "=([^;]*)")
         else      paste0("(?:^|;)[ \t]*", key, "[ \t]+\"?([^\";]*)\"?")
  m <- regexpr(pat, attrs, perl = TRUE)
  out <- rep(NA_character_, length(attrs))
  hit <- m > 0L
  if (any(hit)) {
    cs <- attr(m, "capture.start")[hit, 1L]
    cl <- attr(m, "capture.length")[hit, 1L]
    out[hit] <- trimws(substring(attrs[hit], cs, cs + cl - 1L))
  }
  out[!is.na(out) & !nzchar(out)] <- NA_character_
  # GFF3 percent-encodes reserved characters, and "ID=gene:ENSG..." is common.
  if (gff3 && any(hit)) {
    out[hit] <- sub("^(gene|transcript):", "", out[hit])
    if (any(grepl("%", out[hit], fixed = TRUE)))
      out[hit] <- utils::URLdecode(out[hit])
  }
  out
}

#' Read gene-level annotation straight out of a GTF or GFF3.
#'
#' Removes the need for a separate id/symbol table: point the app at the same
#' annotation the aligner and nf-core/differentialabundance used and every
#' gene id, symbol, biotype and locus comes from one file.
#'
#' Only "gene" rows are kept (falling back to "transcript" rows, deduplicated,
#' for annotations that have none), and the file is streamed in chunks so a
#' full human GTF never has to sit in memory at once.
#'
#' @param max_extra how many additional attribute columns beyond the known set
#'   to carry through.
read_gtf_annotation <- function(path, max_extra = 10L) {
  cols <- readr::cols_only(
    X1 = readr::col_character(), X3 = readr::col_character(),
    X4 = readr::col_double(),    X5 = readr::col_double(),
    X7 = readr::col_character(), X9 = readr::col_character())

  grab <- function(types) {
    cb <- readr::DataFrameCallback$new(function(x, pos) {
      if (!nrow(x)) return(x[0, , drop = FALSE])
      x[!is.na(x$X3) & x$X3 %in% types, , drop = FALSE]
    })
    suppressWarnings(readr::read_tsv_chunked(
      path, cb, chunk_size = 2e5, col_names = FALSE, col_types = cols,
      comment = "#", progress = FALSE))
  }

  df <- grab(c("gene", "transcript"))
  if (is.null(df) || !nrow(df)) df <- grab("exon")
  if (is.null(df) || !nrow(df)) {
    stop("No gene, transcript or exon records found in ", basename(path),
         " -- is this really a GTF/GFF file?")
  }
  if (any(df$X3 == "gene")) df <- df[df$X3 == "gene", , drop = FALSE]

  attrs <- df$X9
  gff3  <- .is_gff3(attrs)
  keys  <- .gtf_attr_keys(attrs, gff3)

  id_key <- intersect(c("gene_id", "geneID", "ID", "gene"), keys)[1]
  if (is.na(id_key)) stop("No gene_id attribute found in ", basename(path), ".")
  gene_id <- .gtf_attr(attrs, id_key, gff3)

  name_keys <- intersect(c("gene_name", "gene", "Name", "gene_symbol", "symbol"), keys)
  name_keys <- setdiff(name_keys, id_key)
  gene_name <- rep(NA_character_, length(gene_id))
  for (k in name_keys) {
    todo <- is.na(gene_name)
    if (!any(todo)) break
    gene_name[todo] <- .gtf_attr(attrs[todo], k, gff3)
  }

  out <- data.frame(gene_id = gene_id, gene_name = gene_name,
                    chromosome = df$X1, start = df$X4, end = df$X5,
                    strand = df$X7, stringsAsFactors = FALSE)

  known <- c("gene_biotype", "gene_type", "biotype", "gene_version",
             "gene_source", "level", "hgnc_id", "havana_gene", "tag")
  extra <- c(intersect(known, keys),
             utils::head(setdiff(keys, c(known, id_key, name_keys,
                                         "transcript_id", "transcript_name",
                                         "transcript_version", "exon_id",
                                         "exon_number", "Parent", "ID", "Name")),
                         max_extra))
  for (k in unique(extra)) {
    v <- .gtf_attr(attrs, k, gff3)
    if (any(!is.na(v))) out[[k]] <- v
  }

  out <- out[!is.na(out$gene_id) & nzchar(out$gene_id), , drop = FALSE]
  out <- out[!duplicated(out$gene_id), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "source") <- if (gff3) "GFF3" else "GTF"
  attr(out, "n_symbols") <- sum(!is.na(out$gene_name))
  out
}

# ---- gene annotation --------------------------------------------------------

#' Read a gene annotation table and reduce it to the columns worth carrying.
#'
#' Accepts anything with a gene id column plus a symbol column: an
#' nf-core/differentialabundance '*.annotated.tsv', the pipeline's own
#' gene-info table, a two-column id/symbol map, or a GTF-derived table.
read_gene_annotation <- function(path, filename = NULL) {
  if (is_gtf_file(path, filename)) return(read_gtf_annotation(path))
  df  <- read_table_auto(path)
  map <- detect_de_columns(df)
  if (is.na(map$gene_id)) {
    stop("No gene id column found in ", basename(path), ".")
  }
  if (is.na(map$gene_name)) {
    stop("No gene symbol column (gene_name / symbol) found in ", basename(path),
         ". Columns present: ", paste(utils::head(colnames(df), 15), collapse = ", "))
  }
  keep <- c(map$gene_id, map$gene_name,
            intersect(colnames(df),
                      c("gene_biotype", "biotype", "chromosome", "description",
                        "gene_version", "entrezid", "entrez_id")))
  out <- df[, unique(keep), drop = FALSE]
  colnames(out)[1:2] <- c("gene_id", "gene_name")
  out$gene_id   <- as.character(out$gene_id)
  out$gene_name <- as.character(out$gene_name)
  out <- out[!is.na(out$gene_id) & out$gene_id != "", , drop = FALSE]
  out <- out[!duplicated(out$gene_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Match DE gene ids against annotation ids, tolerating version suffixes.
#'
#' GENCODE writes "ENSG00000223972.5" while Ensembl and most count pipelines
#' write "ENSG00000223972"; without this, a perfectly good GTF would silently
#' match nothing.
annotation_match <- function(ids, ann_ids) {
  idx  <- match(ids, ann_ids)
  miss <- is.na(idx)
  if (any(miss)) {
    strip <- function(x) sub("\\.[0-9]+$", "", x)
    idx[miss] <- match(strip(ids[miss]), strip(ann_ids))
  }
  idx
}

#' Fill in missing gene symbols on a canonical DE table from an annotation map.
#' Existing symbols win; only blanks are filled, so an annotated results file
#' is never overwritten by a sparser map.
#'
#' @return list(data = <de with symbols>, n_filled = <int>)
apply_gene_annotation <- function(de, ann) {
  if (is.null(de) || is.null(ann) || !nrow(ann)) return(list(data = de, n_filled = 0L))
  blank <- is.na(de$gene_name) | trimws(de$gene_name) %in% c("", "NA", "-")
  idx   <- annotation_match(de$gene_id, ann$gene_id)
  hit   <- ann$gene_name[idx]
  fill  <- blank & !is.na(hit) & trimws(hit) != ""
  de$gene_name[fill] <- hit[fill]

  extra <- setdiff(colnames(ann), c("gene_id", "gene_name"))
  for (cl in extra) {
    if (!cl %in% colnames(de)) de[[cl]] <- ann[[cl]][idx]
  }
  de$display_gene <- ifelse(
    is.na(de$gene_name) | trimws(de$gene_name) %in% c("", "NA", "-"),
    de$gene_id, de$gene_name)
  list(data = de, n_filled = sum(fill))
}

#' Gene ids to use as the statistical background for enrichment.
#' Prefers an explicitly supplied unfiltered table, then a count matrix, then
#' whatever the DE table itself holds.
enrichment_universe <- function(de, universe_ids = NULL, count_ids = NULL) {
  # Transcript-level results: the background is the genes of the tested transcripts.
  if ("parent_gene_id" %in% colnames(de)) {
    u <- de$parent_gene_id
    return(unique(as.character(u[!is.na(u) & u != ""])))
  }
  u <- universe_ids
  if (is.null(u) || !length(u)) u <- count_ids
  if (is.null(u) || !length(u)) u <- de$gene_id[!is.na(de$padj)]
  unique(as.character(u[!is.na(u) & u != ""]))
}

# ---- count / length matrices ------------------------------------------------

#' Read a gene x sample numeric matrix.
#' Tolerates an extra `gene_name` (or similar) second column, as produced by
#' nf-core/rnaseq's rsem.merged.gene_counts.tsv, and returns that separately.
read_matrix_auto <- function(path) {
  df <- read_table_auto(path)
  id_col <- colnames(df)[1]
  ids <- as.character(df[[1]])

  rest <- df[, -1, drop = FALSE]
  # drop leading non-numeric annotation columns (e.g. gene_name)
  ann <- NULL
  while (ncol(rest) > 0) {
    v <- rest[[1]]
    num <- suppressWarnings(as.numeric(as.character(v)))
    if (all(is.na(num)) && !all(is.na(v))) {
      if (is.null(ann)) ann <- data.frame(gene_id = ids,
                                          gene_name = as.character(v),
                                          stringsAsFactors = FALSE)
      rest <- rest[, -1, drop = FALSE]
    } else break
  }
  if (ncol(rest) < 1) stop("No numeric sample columns found in ", basename(path))

  mat <- suppressWarnings(as.matrix(vapply(rest, function(x) as.numeric(as.character(x)),
                                            numeric(nrow(rest)))))
  colnames(mat) <- colnames(rest)
  rownames(mat) <- ids

  keep <- !is.na(ids) & ids != ""
  mat <- mat[keep, , drop = FALSE]
  if (!is.null(ann)) ann <- ann[keep, , drop = FALSE]

  if (anyDuplicated(rownames(mat))) {
    o <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
    mat <- mat[o, , drop = FALSE]
    if (!is.null(ann)) ann <- ann[o, , drop = FALSE]
    dup <- duplicated(rownames(mat))
    mat <- mat[!dup, , drop = FALSE]
    if (!is.null(ann)) ann <- ann[!dup, , drop = FALSE]
  }
  mat[!is.finite(mat)] <- NA_real_
  list(matrix = mat, annotation = ann, id_column = id_col)
}

# ---- sample metadata --------------------------------------------------------

#' Read a sample sheet. First column = sample id; a `condition`-like column is
#' located case-insensitively and renamed to `condition`.
read_sample_metadata <- function(path) {
  df <- read_table_auto(path)
  colnames(df)[1] <- "sample"
  df$sample <- as.character(df$sample)

  keys <- .norm_key(colnames(df))
  cand <- c("condition", "group", "treatment", "sampletype", "genotype",
            "status", "class", "arm")
  hit <- NA_character_
  for (cd in cand) { i <- which(keys == cd); if (length(i)) { hit <- colnames(df)[i[1]]; break } }
  if (is.na(hit)) {
    stop("The sample sheet needs a condition/group column (tried: ",
         paste(cand, collapse = ", "), "). Columns found: ",
         paste(colnames(df), collapse = ", "))
  }
  if (hit != "condition") colnames(df)[colnames(df) == hit] <- "condition"
  df$condition <- as.character(df$condition)
  df <- df[!is.na(df$sample) & df$sample != "", , drop = FALSE]
  rownames(df) <- NULL
  df
}

#' Build a plausible sample sheet from matrix column names when the user has
#' not supplied one: everything starts in a single "ungrouped" group so the
#' app still runs, and the user can edit groups in the UI.
scaffold_metadata <- function(sample_ids) {
  data.frame(sample = sample_ids,
             condition = rep("ungrouped", length(sample_ids)),
             stringsAsFactors = FALSE)
}

# ---- writing ----------------------------------------------------------------

#' Write a matrix (with a gene id column) as CSV or TSV.
write_matrix_file <- function(mat, file, format = c("csv", "tsv"),
                              id_name = "gene_id", annotation = NULL) {
  format <- match.arg(format)
  out <- data.frame(x = rownames(mat), stringsAsFactors = FALSE, check.names = FALSE)
  colnames(out) <- id_name
  if (!is.null(annotation)) {
    gn <- annotation$gene_name[match(rownames(mat), annotation$gene_id)]
    out$gene_name <- gn
  }
  out <- cbind(out, as.data.frame(mat, check.names = FALSE))
  if (format == "csv") {
    utils::write.csv(out, file, row.names = FALSE)
  } else {
    utils::write.table(out, file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  invisible(file)
}
