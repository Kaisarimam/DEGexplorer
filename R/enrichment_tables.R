# =============================================================================
# enrichment_tables.R -- enrichment results as plain tables
#
# The Enrichment gallery works from ordinary data frames in clusterProfiler's
# column layout (ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
# geneID, Count), one per database: GO BP, GO CC, GO MF and KEGG. That keeps
# every gallery figure independent of Bioconductor, and means results computed
# elsewhere -- clusterProfiler, an online GO/pathway enrichment service, DAVID,
# g:Profiler -- can be imported and drawn without re-running anything.
# =============================================================================

GALLERY_DATABASES <- c("BP", "CC", "MF", "KEGG")

GALLERY_DB_LABELS <- c(BP = "Biological process", CC = "Cellular component",
                       MF = "Molecular function", KEGG = "KEGG pathway")

GALLERY_DB_TITLES <- c(BP = "Biological Process", CC = "Cellular Component",
                       MF = "Molecular Function", KEGG = "Pathway Analysis")

GALLERY_STD_COLS <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue",
                      "p.adjust", "qvalue", "geneID", "Count")

# ---- column detection -------------------------------------------------------

.norm_colname <- function(x) gsub("[^a-z0-9]", "", tolower(x))

.ENRICH_COL_ALIASES <- list(
  ID          = c("id", "termid", "goid", "go", "pathwayid", "keggid", "native"),
  Description = c("description", "term", "termname", "pathway", "pathwayname",
                  "name", "goterm"),
  GeneRatio   = c("generatio", "ratio", "richfactor", "precision"),
  BgRatio     = c("bgratio", "backgroundratio"),
  pvalue      = c("pvalue", "pval", "p", "rawp", "nominalp"),
  p.adjust    = c("padjust", "padj", "fdr", "adjpval", "adjpvalue", "adjustedpvalue",
                  "benjamini", "bh", "qvalue"),
  qvalue      = c("qvalue", "qval"),
  geneName    = c("genename", "genesymbol", "symbols", "genesymbols"),
  geneID      = c("geneid", "genes", "gene", "geneids", "coreenrichment",
                  "intersection", "intersections", "members"),
  Count       = c("count", "genecount", "genenumber", "intersectionsize", "hits",
                  "overlap", "numgenes"),
  ONTOLOGY    = c("ontology", "category", "source", "subgroup", "namespace",
                  "database", "db", "class")
)

.find_enrich_col <- function(cols, key) {
  hit <- match(.ENRICH_COL_ALIASES[[key]], .norm_colname(cols))
  hit <- hit[!is.na(hit)]
  if (length(hit)) cols[hit[1]] else NA_character_
}

#' Evaluate "8/90"-style ratios (or plain numbers) without eval(parse()).
ratio_to_numeric <- function(x) {
  x <- as.character(x)
  vapply(strsplit(x, "/", fixed = TRUE), function(p) {
    v <- suppressWarnings(as.numeric(p))
    if (length(v) == 2 && !anyNA(v) && v[2] != 0) v[1] / v[2]
    else if (length(v) == 1) v[1] else NA_real_
  }, numeric(1))
}

#' Translate the many spellings of an ontology into BP / CC / MF / KEGG.
normalise_ontology <- function(x) {
  k <- .norm_colname(x)
  out <- rep(NA_character_, length(k))
  out[grepl("^(bp|gobp|gotermbp.*|biologicalprocess)$", k)] <- "BP"
  out[grepl("^(cc|gocc|gotermcc.*|cellularcomponent)$", k)] <- "CC"
  out[grepl("^(mf|gomf|gotermmf.*|molecularfunction)$", k)] <- "MF"
  out[grepl("^(kegg|keggpathway|pathway|keggpathways)$", k)] <- "KEGG"
  out
}

.is_kegg_id <- function(id) grepl("^[a-z]{2,5}[0-9]{5}$", id)

# ---- standardising ----------------------------------------------------------

#' Coerce any enrichment table (or a clusterProfiler result) to the standard
#' layout, with numeric statistics, "/"-separated gene symbols and two derived
#' columns: `GeneRatioNum` and `Score` (= -log10 p-value).
#'
#' @param x   data frame, or anything `as.data.frame()` understands
#' @param de  optional canonical DE table; Ensembl ids in the gene column are
#'            translated to symbols through it
standardise_enrich_table <- function(x, de = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!nrow(df)) return(.empty_enrich_table())
  cols <- colnames(df)
  pick <- function(key) .find_enrich_col(cols, key)

  c_desc <- pick("Description"); c_id <- pick("ID")
  if (is.na(c_desc) && is.na(c_id))
    stop("No term / pathway name column found (expected e.g. 'Description', 'Term' or 'pathway').")
  c_p <- pick("pvalue"); c_padj <- pick("p.adjust")
  if (is.na(c_p) && is.na(c_padj))
    stop("No p-value column found (expected e.g. 'pvalue', 'p.adjust' or 'FDR').")

  desc <- as.character(df[[if (is.na(c_desc)) c_id else c_desc]])
  id   <- if (is.na(c_id)) rep(NA_character_, nrow(df)) else as.character(df[[c_id]])
  # DAVID writes "GO:0006955~immune response" / "hsa04110:Cell cycle" in one field
  tilde <- grepl("^(GO:[0-9]{7})~", desc)
  id[tilde] <- sub("~.*$", "", desc[tilde]); desc[tilde] <- sub("^[^~]*~", "", desc[tilde])
  colon <- grepl("^[a-z]{2,5}[0-9]{5}:", desc)
  id[colon] <- sub(":.*$", "", desc[colon]); desc[colon] <- sub("^[^:]*:", "", desc[colon])
  id[is.na(id) | !nzchar(id)] <- desc[is.na(id) | !nzchar(id)]

  num <- function(col) if (is.na(col)) rep(NA_real_, nrow(df)) else
    suppressWarnings(as.numeric(df[[col]]))
  pval <- num(c_p); padj <- num(c_padj)
  if (is.na(c_p))    pval <- padj
  if (is.na(c_padj)) padj <- stats::p.adjust(pval, "BH")
  c_q <- pick("qvalue")
  qval <- if (is.na(c_q)) padj else num(c_q)

  # Gene members: prefer an explicit symbol column (this app's own export has
  # both geneName and geneID), otherwise whatever gene column exists.
  c_gene <- pick("geneName"); if (is.na(c_gene)) c_gene <- pick("geneID")
  genes <- if (is.na(c_gene)) rep("", nrow(df)) else as.character(df[[c_gene]])
  genes[is.na(genes)] <- ""
  genes <- gsub("\\s*[,;|]\\s*|\\s+", "/", trimws(genes))
  genes <- gsub("/+", "/", genes)
  if (!is.null(de) && any(grepl("ENS[A-Z]*G[0-9]{6,}", genes))) {
    genes <- ensembl_to_symbols(genes, de)
  }
  n_genes <- lengths(strsplit(genes, "/", fixed = TRUE))
  n_genes[!nzchar(genes)] <- 0L

  c_cnt <- pick("Count")
  count <- if (is.na(c_cnt)) n_genes else num(c_cnt)
  count[is.na(count)] <- n_genes[is.na(count)]

  c_gr <- pick("GeneRatio"); c_bg <- pick("BgRatio")
  gr <- if (is.na(c_gr)) rep(NA_character_, nrow(df)) else as.character(df[[c_gr]])
  bg <- if (is.na(c_bg)) rep(NA_character_, nrow(df)) else as.character(df[[c_bg]])

  out <- data.frame(ID = id, Description = desc, GeneRatio = gr, BgRatio = bg,
                    pvalue = pval, p.adjust = padj, qvalue = qval,
                    geneID = genes, Count = count,
                    stringsAsFactors = FALSE)
  out$GeneRatioNum <- ratio_to_numeric(gr)
  out$Score <- -log10(pmax(out$pvalue, .Machine$double.xmin))
  c_ont <- pick("ONTOLOGY")
  if (!is.na(c_ont)) out$ONTOLOGY <- normalise_ontology(df[[c_ont]])
  out <- out[!is.na(out$pvalue) & !is.na(out$Description) & nzchar(out$Description), , drop = FALSE]
  out <- out[order(out$pvalue), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.empty_enrich_table <- function() {
  data.frame(ID = character(0), Description = character(0), GeneRatio = character(0),
             BgRatio = character(0), pvalue = numeric(0), p.adjust = numeric(0),
             qvalue = numeric(0), geneID = character(0), Count = numeric(0),
             GeneRatioNum = numeric(0), Score = numeric(0), stringsAsFactors = FALSE)
}

#' The significant terms of a table, most significant first.
#'
#' @param rank_by "p.adjust" (the default) or "pvalue"
gallery_sig <- function(df, rank_by = "p.adjust", cutoff = 0.05) {
  if (is.null(df) || !nrow(df)) return(df)
  rank_by <- if (identical(rank_by, "pvalue")) "pvalue" else "p.adjust"
  cutoff <- suppressWarnings(as.numeric(cutoff))
  if (length(cutoff) != 1 || is.na(cutoff)) cutoff <- 0.05
  df <- df[!is.na(df[[rank_by]]) & df[[rank_by]] <= cutoff, , drop = FALSE]
  df[order(df[[rank_by]], df$pvalue), , drop = FALSE]
}

gallery_top <- function(df, top_n = 10, rank_by = "p.adjust", cutoff = 0.05) {
  s <- gallery_sig(df, rank_by, cutoff)
  if (is.null(s)) return(s)
  utils::head(s, max(1L, as.integer(top_n)))
}

#' The table as written to disk: the nine standard columns only.
gallery_export_table <- function(df) {
  if (is.null(df)) return(NULL)
  df[, intersect(GALLERY_STD_COLS, colnames(df)), drop = FALSE]
}

# ---- reading files ----------------------------------------------------------

#' Read one enrichment results file (xlsx, csv, tsv, txt; optionally gzipped).
read_enrichment_file <- function(path, name = basename(path)) {
  if (grepl("\\.xlsx?$", name, ignore.case = TRUE)) {
    if (!requireNamespace("openxlsx", quietly = TRUE))
      stop('Reading Excel files needs the openxlsx package: install.packages("openxlsx")')
    return(openxlsx::read.xlsx(path, check.names = FALSE))
  }
  read_table_auto(path)
}

#' Decide which database a table belongs to.
#' Order of evidence: an ontology column, KEGG-style ids, the file name, and
#' finally the caller's hint. Returns a named list of tables (a single file
#' holding all three ontologies is split).
split_enrichment_by_database <- function(std, name = "", hint = "auto") {
  if (!nrow(std)) return(list())
  db <- rep(NA_character_, nrow(std))
  if ("ONTOLOGY" %in% colnames(std)) db <- std$ONTOLOGY
  db[is.na(db) & .is_kegg_id(std$ID)] <- "KEGG"
  if (anyNA(db)) {
    from_name <- NA_character_
    nm <- toupper(name)
    if (grepl("(^|[^A-Z])(KEGG|PATHWAY)([^A-Z]|$)", nm)) from_name <- "KEGG"
    for (o in c("BP", "CC", "MF"))
      if (grepl(paste0("(^|[^A-Z])", o, "([^A-Z]|$)"), nm)) from_name <- o
    if (is.na(from_name) && hint %in% GALLERY_DATABASES) from_name <- hint
    db[is.na(db)] <- from_name
  }
  if (all(is.na(db))) return(list())
  std$ONTOLOGY <- NULL
  out <- split(std[!is.na(db), , drop = FALSE], db[!is.na(db)])
  lapply(out, function(d) { rownames(d) <- NULL; d })
}

#' Import one or more enrichment result files into a gallery table set.
#'
#' @return list(tables = named list keyed by BP/CC/MF/KEGG, notes = character)
import_enrichment_files <- function(paths, names = basename(paths), hint = "auto", de = NULL) {
  tables <- list(); notes <- character(0)
  for (i in seq_along(paths)) {
    one <- tryCatch({
      std <- standardise_enrich_table(read_enrichment_file(paths[i], names[i]), de)
      split_enrichment_by_database(std, names[i], hint)
    }, error = function(e) e)
    if (inherits(one, "error")) {
      notes <- c(notes, sprintf("%s: %s", names[i], conditionMessage(one))); next
    }
    if (!length(one)) {
      notes <- c(notes, sprintf(
        "%s: could not tell which database this is. Put BP, CC, MF or KEGG in the file name, or choose it in the sidebar.",
        names[i])); next
    }
    for (db in names(one)) {
      if (!is.null(tables[[db]]))
        notes <- c(notes, sprintf("%s: replaces the %s table loaded from an earlier file.", names[i], db))
      tables[[db]] <- one[[db]]
      notes <- c(notes, sprintf("%s: %s terms read as %s.", names[i],
                                format(nrow(one[[db]]), big.mark = ","), GALLERY_DB_LABELS[[db]]))
    }
  }
  list(tables = tables[intersect(GALLERY_DATABASES, names(tables))], notes = notes)
}

#' The demo tables shipped for the import option (data/example_enrichment/).
#' They are computed from the app's own bundled demo data -- see
#' data-raw/make_example_enrichment.R -- and are ordinary files in the import
#' layout, so they double as a template for your own.
example_enrichment_files <- function(dir) {
  f <- file.path(dir, c("go.BP.all.csv", "go.CC.all.csv", "go.MF.all.csv", "pathway.all.csv"))
  f[file.exists(f)]
}

#' clusterProfiler results -> gallery tables (gene members as symbols).
enrich_results_to_gallery <- function(results, de = NULL) {
  out <- list()
  for (db in intersect(GALLERY_DATABASES, names(results))) {
    r <- results[[db]]
    if (is.null(r)) next
    df <- as.data.frame(r)
    if (!nrow(df)) next
    if ("geneID_ensembl" %in% colnames(df)) {      # KEGG ran on Entrez ids
      df$geneID <- df$geneID_ensembl; df$geneID_ensembl <- NULL
    }
    out[[db]] <- standardise_enrich_table(df, de)
  }
  out
}

# ---- KEGG pathway classes ---------------------------------------------------

KEGG_CLASSES <- c("Metabolism", "Genetic Information Processing",
                  "Environmental Information Processing", "Cellular Processes",
                  "Organismal Systems", "Human Diseases", "Drug Development",
                  "Unclassified")

# KEGG map numbers encode the top-level class for most of the range (00/01 =
# metabolism, 03 = genetic information processing, 05 = human diseases, ...).
# The 04xxx block mixes four classes, so those maps are listed explicitly;
# anything in 04xxx not listed here is an Organismal Systems map.
.KEGG_04_ENVIRONMENTAL <- c(
  "04010", "04011", "04012", "04013", "04014", "04015", "04016", "04020", "04022",
  "04024", "04060", "04061", "04064", "04066", "04068", "04070", "04071", "04072",
  "04075", "04080", "04081", "04082", "04150", "04151", "04152", "04310", "04330",
  "04340", "04341", "04350", "04370", "04371", "04390", "04391", "04392", "04512",
  "04514", "04630", "04668")
.KEGG_04_GENETIC <- c("04120", "04122", "04130", "04141")
.KEGG_04_CELLULAR <- c(
  "04110", "04111", "04112", "04113", "04114", "04115", "04136", "04137", "04138",
  "04139", "04140", "04142", "04144", "04145", "04146", "04148", "04210", "04214",
  "04215", "04216", "04217", "04218", "04510", "04520", "04530", "04540", "04550",
  "04810", "04814", "04820")
.KEGG_02_CELLULAR <- c("02024", "02025", "02026", "02030", "02040")

#' Top-level KEGG class of a pathway, from its id (hsa04110, map00010, 04110).
kegg_pathway_class <- function(id) {
  num <- sub("^[A-Za-z]*", "", trimws(as.character(id)))
  ok <- grepl("^[0-9]{5}$", num)
  pre <- substr(num, 1, 2)
  cls <- rep("Unclassified", length(num))
  cls[ok & pre %in% c("00", "01")] <- "Metabolism"
  cls[ok & pre == "02"] <- "Environmental Information Processing"
  cls[ok & num %in% .KEGG_02_CELLULAR] <- "Cellular Processes"
  cls[ok & pre == "03"] <- "Genetic Information Processing"
  cls[ok & num == "03320"] <- "Organismal Systems"          # PPAR signalling
  cls[ok & pre == "04"] <- "Organismal Systems"
  cls[ok & num %in% .KEGG_04_ENVIRONMENTAL] <- "Environmental Information Processing"
  cls[ok & num %in% .KEGG_04_GENETIC] <- "Genetic Information Processing"
  cls[ok & num %in% .KEGG_04_CELLULAR] <- "Cellular Processes"
  cls[ok & pre == "05"] <- "Human Diseases"
  cls[ok & pre == "07"] <- "Drug Development"
  cls
}

# ---- writing ----------------------------------------------------------------

#' Write a table as .xlsx when openxlsx is available, otherwise as .csv.
#' Returns the path actually written.
write_table_xlsx_or_csv <- function(df, path_no_ext) {
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    f <- paste0(path_no_ext, ".xlsx"); openxlsx::write.xlsx(df, f); return(f)
  }
  f <- paste0(path_no_ext, ".csv"); utils::write.csv(df, f, row.names = FALSE); f
}

#' Zip a directory's contents (zip package if present, else the system zip).
zip_directory <- function(dir, zipfile) {
  files <- list.files(dir, recursive = TRUE)
  if (!length(files)) stop("Nothing was produced to put in the archive.")
  target <- normalizePath(zipfile, mustWork = FALSE)
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zip(target, files = list.files(dir), root = dir)
  } else {
    owd <- setwd(dir); on.exit(setwd(owd), add = TRUE)
    st <- utils::zip(target, files = files, flags = "-r9Xq")
    if (!identical(st, 0L) || !file.exists(target))
      stop('No zip program is available. Run install.packages("zip").')
  }
  invisible(target)
}
