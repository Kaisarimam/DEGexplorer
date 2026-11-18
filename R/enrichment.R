# =============================================================================
# enrichment.R -- GO / KEGG / MSigDB over-representation analysis
#
# clusterProfiler, the organism's OrgDb (org.Hs.eg.db for human, org.Mm.eg.db
# for mouse, ...), enrichplot and msigdbr are all optional. The rest of the app
# runs without them; this module reports precisely which package is missing for
# the current selection rather than failing opaquely. The organism comes from
# the app-wide Organism setting (R/organisms.R); human is the default.
# =============================================================================

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

HAS_CLUSTERPROFILER <- has_pkg("clusterProfiler")
HAS_ENRICHPLOT      <- has_pkg("enrichplot")
# enrichplot::goplot() draws its arrows with ggarchery in recent enrichplot releases.
HAS_GGARCHERY       <- has_pkg("ggarchery")
HAS_ORGDB           <- has_pkg("org.Hs.eg.db")
HAS_MSIGDBR         <- has_pkg("msigdbr")
HAS_GRIDEXTRA       <- has_pkg("gridExtra")

MSIGDB_DATABASES <- c("HALLMARK", "REACTOME", "WIKIPATHWAYS", "IMMUNESIGDB", "GOBP_MSIG", "TFT")

ENRICH_DATABASES <- c(
  "GO (organism OrgDb): Biological Process" = "BP",
  "GO (organism OrgDb): Molecular Function" = "MF",
  "GO (organism OrgDb): Cellular Component" = "CC",
  "KEGG pathway"             = "KEGG",
  "MSigDB: Hallmark"         = "HALLMARK",
  "MSigDB: Reactome"         = "REACTOME",
  "MSigDB: WikiPathways"     = "WIKIPATHWAYS",
  "MSigDB: ImmuneSigDB"      = "IMMUNESIGDB",
  "MSigDB: TF targets"       = "TFT"
)

# ---- GO, official Gene Ontology Consortium release (bundled) ----------------
#
# data/go_release/ holds gene sets built from one dated GO release
# (go-basic.obo + the UniProt human GAF), propagated over is_a + part_of.
# They identify genes by symbol and are matched to the loaded DE table by
# symbol. Rebuild for a newer release with
#   python3 data-raw/build_go_table.py go-basic.obo HUMAN-uniprot.gaf.gz OUTDIR --app-dir data/go_release

GO_RELEASE_DATABASES <- c("GOREL_BP", "GOREL_MF", "GOREL_CC")

.go_release_file <- function(f) {
  if (exists("app_file", mode = "function")) {
    p <- app_file("data", "go_release", f)
    if (file.exists(p)) return(p)
  }
  for (p in c(file.path("data", "go_release", f), file.path("..", "data", "go_release", f))) {
    if (file.exists(p)) return(p)
  }
  file.path("data", "go_release", f)
}

#' The bundled release date ("2026-07-26"), or NA when the files are absent.
go_release_version <- function() {
  f <- .go_release_file("VERSION.txt")
  if (!file.exists(f)) return(NA_character_)
  x <- grep("^release:", readLines(f, warn = FALSE), value = TRUE)
  if (!length(x)) return(NA_character_)
  trimws(sub("^release:", "", x[1]))
}

local({
  v <- go_release_version()
  if (!is.na(v)) {
    lab <- paste0("GO release ", v, " (human): ",
                  c("Biological Process", "Molecular Function", "Cellular Component"))
    ENRICH_DATABASES <<- append(ENRICH_DATABASES,
                                stats::setNames(GO_RELEASE_DATABASES, lab), after = 3)
  }
})

.go_release_cache <- new.env(parent = emptyenv())

#' Read the bundled term -> gene and term -> name tables once per session.
go_release_sets <- function() {
  if (!is.null(.go_release_cache$t2g)) return(.go_release_cache)
  f1 <- .go_release_file("term2gene.tsv.gz"); f2 <- .go_release_file("term2name.tsv.gz")
  if (!file.exists(f1) || !file.exists(f2)) {
    stop("The bundled GO release files are missing from data/go_release/. Rebuild them with ",
         "data-raw/build_go_table.py --app-dir data/go_release.")
  }
  rd <- function(f) as.data.frame(readr::read_tsv(f, col_types = readr::cols(.default = "c"),
                                                  progress = FALSE))
  .go_release_cache$t2g <- rd(f1)
  .go_release_cache$t2n <- rd(f2)
  .go_release_cache
}

#' Gene sets for one GO release ontology, re-keyed from symbols to the
#' Ensembl ids of the loaded DE table.
#'
#' @param symbol_map named character vector: names = Ensembl ids, values = symbols.
go_release_term2gene <- function(database, symbol_map) {
  ont <- sub("^GOREL_", "", database)
  s <- go_release_sets()
  t2g <- s$t2g[s$t2g$ontology == ont, c("go_id", "gene_symbol"), drop = FALSE]
  keep <- !is.na(symbol_map) & nzchar(symbol_map) & !grepl("^ENS", symbol_map)
  if (!any(keep)) {
    stop("The GO release databases match genes by symbol, and the loaded table has no gene symbols.")
  }
  idx <- split(names(symbol_map)[keep], unname(symbol_map[keep]))   # symbol -> Ensembl id(s)
  t2g <- t2g[t2g$gene_symbol %in% names(idx), , drop = FALSE]
  ids <- idx[t2g$gene_symbol]
  out <- unique(data.frame(term = rep(t2g$go_id, lengths(ids)),
                           gene = unlist(ids, use.names = FALSE),
                           stringsAsFactors = FALSE))
  t2n <- s$t2n[s$t2n$ontology == ont, c("go_id", "name"), drop = FALSE]
  colnames(t2n) <- c("term", "name")
  list(t2g = out, t2n = t2n)
}

.default_organism <- function() organism_info("human")

#' Which optional packages are missing for a given database / plot choice?
#'
#' Also reports, as a line starting with "#", a database that cannot serve the
#' organism (the bundled GO release is human only; MSigDB covers the organisms
#' msigdbr knows).
enrichment_missing_packages <- function(database, plot_type = "dot", organism = NULL) {
  org <- organism %||% .default_organism()
  miss <- character(0)
  need_orgdb <- database %in% c("BP", "MF", "CC", "KEGG")
  orgdb_line <- if (need_orgdb) organism_orgdb_missing(org) else NULL
  if (!HAS_CLUSTERPROFILER && !is.null(orgdb_line) && startsWith(orgdb_line, "BiocManager")) {
    miss <- c(miss, sprintf('BiocManager::install(c("clusterProfiler", "%s"))', org$orgdb))
  } else {
    if (!HAS_CLUSTERPROFILER) miss <- c(miss, 'BiocManager::install("clusterProfiler")')
    if (!is.null(orgdb_line)) miss <- c(miss, if (startsWith(orgdb_line, "BiocManager")) orgdb_line
                                              else paste("#", orgdb_line, "Set its OrgDb package under Organism (Other)."))
  }
  if (database %in% GO_RELEASE_DATABASES && !identical(org$key, "human")) {
    miss <- c(miss, sprintf("# The bundled GO release is built from the human GAF; for %s choose a 'GO (organism OrgDb)' database.", org$label))
  }
  if (database %in% MSIGDB_DATABASES) {
    if (!HAS_MSIGDBR) miss <- c(miss, 'install.packages("msigdbr")')
    if (!nzchar(org$msigdbr %||% "")) miss <- c(miss, sprintf("# MSigDB has no gene sets for %s.", org$label))
  }
  if (plot_type %in% c("cnet", "emap", "heat") && !HAS_ENRICHPLOT) {
    miss <- c(miss, 'BiocManager::install("enrichplot")')
  }
  miss
}

# ---- MSigDB retrieval -------------------------------------------------------

.msigdbr_cache <- new.env(parent = emptyenv())

#' Fetch an MSigDB collection, tolerating the v10 rename of
#' category/subcategory to collection/subcollection.
fetch_msigdbr_raw <- function(collection, subcollection_options = NULL, species = "Homo sapiens") {
  try_one <- function(subcoll) {
    a1 <- list(species = species, collection = collection)
    if (!is.null(subcoll)) a1$subcollection <- subcoll
    r <- tryCatch(do.call(msigdbr::msigdbr, a1), error = function(e) NULL)
    if (is.null(r)) {
      a2 <- list(species = species, category = collection)
      if (!is.null(subcoll)) a2$subcategory <- subcoll
      r <- tryCatch(do.call(msigdbr::msigdbr, a2), error = function(e) NULL)
    }
    r
  }
  if (is.null(subcollection_options)) {
    r <- try_one(NULL)
    if (is.null(r) || nrow(r) == 0) stop("msigdbr returned no gene sets for collection '", collection, "'.")
    return(r)
  }
  for (sc in subcollection_options) {
    r <- try_one(sc)
    if (!is.null(r) && nrow(r) > 0) return(r)
  }
  stop("Could not retrieve MSigDB collection '", collection, "'. ",
       "Check msigdbr::msigdbr_collections() for the names your version uses.")
}

fetch_msigdbr_cached <- function(collection, subcollection_options = NULL, species = "Homo sapiens") {
  key <- paste(species, collection, paste(subcollection_options, collapse = "|"), sep = "::")
  if (!is.null(.msigdbr_cache[[key]])) return(.msigdbr_cache[[key]])
  r <- fetch_msigdbr_raw(collection, subcollection_options, species)
  .msigdbr_cache[[key]] <- r
  r
}

pretty_msigdb_name <- function(x) {
  x <- gsub("^(HALLMARK_|REACTOME_|WP_|GOBP_|GOMF_|GOCC_|BIOCARTA_|PID_)", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

msigdb_term2gene <- function(database, species = "Homo sapiens") {
  spec <- switch(database,
    "HALLMARK"     = list(c = "H",  s = NULL),
    "REACTOME"     = list(c = "C2", s = c("REACTOME", "CP:REACTOME")),
    "WIKIPATHWAYS" = list(c = "C2", s = c("WIKIPATHWAYS", "CP:WIKIPATHWAYS")),
    "IMMUNESIGDB"  = list(c = "C7", s = c("IMMUNESIGDB")),
    "TFT"          = list(c = "C3", s = c("TFT:GTRD", "GTRD", "TFT")),
    "GOBP_MSIG"    = list(c = "C5", s = c("GO:BP", "BP")),
    stop("Unknown MSigDB database: ", database)
  )
  raw <- as.data.frame(fetch_msigdbr_cached(spec$c, spec$s, species))
  gcol <- intersect(c("ensembl_gene", "human_ensembl_gene"), colnames(raw))[1]
  if (is.na(gcol)) stop("msigdbr output has no Ensembl gene id column.")
  t2g <- unique(raw[, c("gs_name", gcol)])
  colnames(t2g) <- c("term", "gene")
  t2g <- t2g[!is.na(t2g$gene) & t2g$gene != "", , drop = FALSE]
  scol <- intersect(c("gene_symbol", "human_gene_symbol"), colnames(raw))[1]
  t2s <- if (is.na(scol)) NULL else {
    x <- unique(raw[, c("gs_name", scol)]); colnames(x) <- c("term", "gene")
    x[!is.na(x$gene) & x$gene != "", , drop = FALSE]
  }
  t2n <- unique(data.frame(term = t2g$term, name = pretty_msigdb_name(t2g$term),
                           stringsAsFactors = FALSE))
  list(t2g = t2g, t2s = t2s, t2n = t2n)
}

# ---- run --------------------------------------------------------------------

#' Over-representation analysis for one gene set.
#'
#' `genes` and `universe` are the table's gene ids (Ensembl ids for most
#' organisms; version suffixes are fine). `symbol_map` (names = ids, values =
#' symbols) is used for the bundled GO release, and as a fallback when the
#' organism's databases do not know the ids. Results carry a `geneID_ensembl`
#' column with the table's own ids, so pathway genes map back to the plots.
run_enrichment <- function(genes, universe, database = "BP",
                           pvalue_cutoff = 0.05, qvalue_cutoff = 0.2,
                           symbol_map = NULL, organism = NULL) {
  org <- organism %||% .default_organism()
  miss <- enrichment_missing_packages(database, organism = org)
  if (length(miss)) {
    stop(if (all(startsWith(miss, "#"))) paste(sub("^# ", "", miss), collapse = " ")
         else paste0("Required package(s) not installed. Run:\n", paste(miss, collapse = "\n")))
  }
  genes <- unique(genes[!is.na(genes)])
  universe <- unique(universe[!is.na(universe)])
  if (length(genes) < 3) {
    stop("Need at least 3 genes to run enrichment (", length(genes), " supplied).")
  }
  back_map <- function(res, key_to_id) {
    if (is.null(res) || !nrow(as.data.frame(res))) return(res)
    res@result$geneID_ensembl <- vapply(strsplit(res@result$geneID, "/"), function(k)
      paste(unique(stats::na.omit(unlist(key_to_id[k], use.names = FALSE))), collapse = "/"), character(1))
    res
  }
  tag <- function(res, keytype) {
    if (!is.null(res)) { attr(res, "dex_orgdb") <- org$orgdb; attr(res, "dex_keytype") <- keytype }
    res
  }

  if (database %in% c("BP", "MF", "CC", "KEGG")) {
    kt <- orgdb_keytype(org$orgdb, unique(c(genes, universe)), symbol_map)
    key_to_id <- split(names(kt$map), unname(kt$map))
    g2 <- unique(stats::na.omit(unname(kt$map[genes])))
    u2 <- unique(stats::na.omit(unname(kt$map[universe])))
    if (database != "KEGG") {
      res <- clusterProfiler::enrichGO(
        gene = g2, universe = u2, OrgDb = org$orgdb,
        keyType = kt$keytype, ont = database, pAdjustMethod = "BH",
        pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff, readable = FALSE)
      return(tag(back_map(res, key_to_id), kt$keytype))
    }
    map <- suppressMessages(suppressWarnings(clusterProfiler::bitr(g2, kt$keytype, "ENTREZID", org$orgdb)))
    uni <- suppressMessages(suppressWarnings(clusterProfiler::bitr(u2, kt$keytype, "ENTREZID", org$orgdb)))
    if (nrow(map) < 3) stop("Too few genes mapped to Entrez IDs for KEGG enrichment.")
    if (!nzchar(org$kegg %||% "")) stop("No KEGG organism code is set for ", org$label, ".")
    # KEGG gene ids are NCBI Gene ids for most vertebrates; other organisms are converted.
    kegg_kt <- if (org$kegg %in% c("hsa", "mmu", "rno", "dre", "gga", "ssc", "bta", "cfa", "mcc")) "kegg" else "ncbi-geneid"
    res <- clusterProfiler::enrichKEGG(
      gene = unique(map$ENTREZID), universe = unique(uni$ENTREZID),
      organism = org$kegg, keyType = kegg_kt, pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff)
    ez_to_id <- lapply(split(map[[kt$keytype]], map$ENTREZID), function(k) unlist(key_to_id[k], use.names = FALSE))
    return(tag(back_map(res, ez_to_id), "ENTREZID"))
  }

  if (database %in% GO_RELEASE_DATABASES) {
    if (is.null(symbol_map)) stop("Gene symbols are needed for the GO release databases.")
    t2 <- go_release_term2gene(database, symbol_map)
    # Same gene-set size limits as enrichGO's defaults, so the two GO sources
    # are directly comparable.
    return(clusterProfiler::enricher(
      gene = genes, universe = universe,
      TERM2GENE = t2$t2g, TERM2NAME = t2$t2n, pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff,
      minGSSize = 10, maxGSSize = 500))
  }

  # MSigDB: match by Ensembl id, or by symbol when the ids are not Ensembl ids.
  t2 <- msigdb_term2gene(database, org$msigdbr)
  base <- function(x) sub("\\.[0-9]+$", "", x)
  if (mean(base(universe) %in% t2$t2g$gene) >= 0.3) {
    key <- stats::setNames(base(universe), universe); sets <- t2$t2g
  } else if (!is.null(symbol_map) && !is.null(t2$t2s)) {
    key <- stats::setNames(unname(symbol_map[universe]), universe); sets <- t2$t2s
  } else stop("Few of the gene ids match MSigDB's ", org$label, " gene sets.")
  key_to_id <- split(names(key), unname(key))
  res <- clusterProfiler::enricher(
    gene = unique(stats::na.omit(unname(key[genes]))), universe = unique(stats::na.omit(unname(key))),
    TERM2GENE = sets, TERM2NAME = t2$t2n, pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff,
    minGSSize = 5, maxGSSize = 500)
  res <- back_map(res, key_to_id)
  if (!is.null(res)) {
    attr(res, "dex_orgdb") <- org$orgdb
    attr(res, "dex_keytype") <- if (identical(sets, t2$t2s)) "SYMBOL" else "ENSEMBL"
  }
  res
}

#' Gene ids to test: transcript-level tables map each transcript to its gene
#' (`parent_gene_id`), so a gene counts once however many of its transcripts changed.
enrichment_ids <- function(de, ids) {
  if (!"parent_gene_id" %in% colnames(de)) return(ids)
  p <- de$parent_gene_id[match(ids, de$gene_id)]
  unique(p[!is.na(p) & nzchar(p)])
}

#' id -> symbol map for enrichment (gene ids, also for transcript-level tables).
enrichment_symbol_map <- function(de) {
  if ("parent_gene_id" %in% colnames(de)) {
    sym <- if ("gene_name" %in% colnames(de)) de$gene_name else de$display_gene
    m <- stats::setNames(as.character(sym), de$parent_gene_id)
    return(m[!duplicated(names(m)) & !is.na(names(m))])
  }
  stats::setNames(de$display_gene, de$gene_id)
}

#' A gene-level view of a transcript-level table (one row per gene: its most
#' significant transcript), for labelling enrichment members; gene tables pass through.
enrichment_gene_view <- function(de) {
  if (!"parent_gene_id" %in% colnames(de)) return(de)
  o <- de[order(de$padj, na.last = TRUE), , drop = FALSE]
  o <- o[!is.na(o$parent_gene_id) & !duplicated(o$parent_gene_id), , drop = FALSE]
  o$gene_id <- o$parent_gene_id
  o$display_gene <- ifelse(!is.na(o$gene_name) & nzchar(o$gene_name), o$gene_name, o$parent_gene_id)
  o
}

#' Pathway gene sets re-keyed to the table's rows: for transcript-level tables
#' each gene expands to all of its tested transcripts.
enrichment_sets_to_rows <- function(sets, de) {
  if (!"parent_gene_id" %in% colnames(de) || !length(sets)) return(sets)
  lapply(sets, function(g) de$gene_id[de$parent_gene_id %in% g])
}

#' term label -> Ensembl gene ids, for feeding pathway genes back into plots.
enrichment_gene_sets <- function(enrich_result, database, tag = NULL) {
  df <- as.data.frame(enrich_result)
  if (!nrow(df)) return(list())
  idc <- if ("geneID_ensembl" %in% colnames(df)) "geneID_ensembl" else "geneID"
  gl <- strsplit(df[[idc]], "/")
  lab <- paste0(df$Description, "  (", df$ID, ")")
  if (!is.null(tag)) lab <- paste0("[", tag, "] ", lab)
  names(gl) <- lab
  gl
}

# ---- results table ----------------------------------------------------------

#' Map a "/"-separated list of Ensembl ids onto gene symbols.
#' Ids with no symbol in the DE table are kept as-is, so nothing is lost.
ensembl_to_symbols <- function(ids_slashed, de) {
  if (!length(ids_slashed)) return(character(0))
  if (is.null(de) || !nrow(de)) return(as.character(ids_slashed))
  vapply(strsplit(as.character(ids_slashed), "/", fixed = TRUE), function(ids) {
    ids <- ids[nzchar(ids)]
    if (!length(ids)) return("")
    s <- de$display_gene[match(ids, de$gene_id)]
    paste(ifelse(is.na(s) | !nzchar(s), ids, s), collapse = "/")
  }, character(1), USE.NAMES = FALSE)
}

#' Turn a long "/"-separated gene list into a click-to-expand table cell.
#'
#' A GO term can carry hundreds of genes. Printed in full, one such cell makes
#' its table row hundreds of pixels tall and pushes everything below it -- the
#' enrichment figure included -- far off the bottom of the page. Rendering the
#' list inside <details> keeps every row one line high until the user opens it.
collapsible_gene_cell <- function(items, n_preview = 3) {
  esc <- function(x) htmltools::htmlEscape(x)
  vapply(items, function(x) {
    if (is.na(x) || !nzchar(x)) return("")
    g <- unlist(strsplit(x, "/", fixed = TRUE)); g <- g[nzchar(g)]
    if (!length(g)) return("")
    if (length(g) <= n_preview) {
      return(sprintf('<span class="gene-cell-plain">%s</span>', esc(paste(g, collapse = ", "))))
    }
    sprintf(paste0('<details class="gene-cell"><summary>%s',
                   '<span class="gene-more"> +%d more</span></summary>',
                   '<div class="gene-list">%s</div></details>'),
            esc(paste(utils::head(g, n_preview), collapse = ", ")),
            length(g) - n_preview,
            esc(paste(g, collapse = ", ")))
  }, character(1), USE.NAMES = FALSE)
}

#' Round the wide numeric columns clusterProfiler returns so the table stays
#' readable (raw values run to 17 significant figures).
format_enrich_numbers <- function(df) {
  for (cl in intersect(colnames(df),
                       c("RichFactor", "FoldEnrichment", "oddsRatio", "zScore"))) {
    df[[cl]] <- round(suppressWarnings(as.numeric(df[[cl]])), 3)
  }
  for (cl in intersect(colnames(df), c("pvalue", "p.adjust", "qvalue"))) {
    v <- suppressWarnings(as.numeric(df[[cl]]))
    df[[cl]] <- ifelse(is.na(v), NA_character_, formatC(v, format = "e", digits = 2))
  }
  df
}

#' One enrichment result -> a tidy data frame, ready for display or for
#' CSV / Excel export.
#'
#' Gene membership is reported in up to three columns. `geneName` and `geneID`
#' (symbols and Ensembl accessions) are present for every database, so the
#' column means the same thing wherever you look. KEGG is the one database
#' clusterProfiler runs on Entrez ids rather than Ensembl, and those native
#' accessions are kept alongside in `entrezID` rather than being discarded.
#' Databases that never produce Entrez ids get no empty column.
enrich_result_table <- function(res, de = NULL, direction = NULL) {
  df <- as.data.frame(res)
  if (!nrow(df)) return(NULL)

  has_ens <- "geneID_ensembl" %in% colnames(df) && any(nzchar(df$geneID_ensembl))
  native  <- df$geneID                       # as the database reported it
  ens     <- if (has_ens) df$geneID_ensembl else native
  df$geneID_ensembl <- NULL

  df$geneID   <- ens
  df$geneName <- ensembl_to_symbols(ens, de)
  # `native` differs from `ens` only when the run was done on Entrez ids (KEGG).
  if (has_ens && !identical(as.character(native), as.character(ens))) {
    df$entrezID <- native
  }
  if (!is.null(direction)) df$Direction <- direction

  tail_cols <- intersect(c("Count", "geneName", "geneID", "entrezID", "Direction"),
                         colnames(df))
  df[, c(setdiff(colnames(df), tail_cols), tail_cols), drop = FALSE]
}

#' The gene-membership columns present in an enrichment table, in display order.
ENRICH_GENE_COLS <- c("geneName", "geneID", "entrezID")

#' Friendlier column headers for the on-screen table.
ENRICH_TABLE_HEADERS <- c(geneName = "Genes", geneID = "Ensembl IDs",
                          entrezID = "Entrez IDs",
                          p.adjust = "adj. P", qvalue = "q value",
                          pvalue = "P value")

enrich_table_colnames <- function(df) {
  nm <- colnames(df)
  hit <- match(nm, names(ENRICH_TABLE_HEADERS))
  ifelse(is.na(hit), nm, ENRICH_TABLE_HEADERS[hit])
}

# ---- plots ------------------------------------------------------------------

enrich_top_df <- function(enrich_result, top_n) {
  df <- as.data.frame(enrich_result)
  if (!nrow(df)) return(df)
  df <- head(df[order(df$p.adjust), ], top_n)
  df$GeneRatioNum <- vapply(df$GeneRatio, function(x) eval(parse(text = x)), numeric(1))
  df
}

make_dotplot <- function(res, top_n = 20, style = default_style(),
                         title = NULL, xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  if (is.null(res) || !nrow(as.data.frame(res))) stop("No significant terms to plot.")
  df <- enrich_top_df(res, top_n)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = GeneRatioNum, y = Description,
                                        size = Count, colour = p.adjust)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "red", high = "blue", name = "adj. P") +
    deg_theme(st)
  apply_labels(p, title, NULL, xlab, ylab,
               default_title = "Enrichment results", default_x = "Gene ratio", default_y = NULL)
}

make_barplot <- function(res, top_n = 20, style = default_style(),
                         title = NULL, xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  if (is.null(res) || !nrow(as.data.frame(res))) stop("No significant terms to plot.")
  df <- enrich_top_df(res, top_n)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Count, y = Description, fill = p.adjust)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient(low = "red", high = "blue", name = "adj. P") +
    deg_theme(st)
  apply_labels(p, title, NULL, xlab, ylab,
               default_title = "Enrichment results", default_x = "Gene count", default_y = NULL)
}

#' Combine up- and down-regulated results into one tidy frame with a
#' collision-proof row key (the same term can appear in both directions).
split_enrich_df <- function(res_up, res_down, top_n) {
  prep <- function(r, d) {
    if (is.null(r)) return(NULL)
    df <- as.data.frame(r)
    if (!nrow(df)) return(NULL)
    df <- head(df[order(df$p.adjust), ], top_n)
    df$GeneRatioNum <- vapply(df$GeneRatio, function(x) eval(parse(text = x)), numeric(1))
    df$Direction <- d
    df$row_key <- paste(d, df$Description, sep = "\r")
    df
  }
  out <- rbind(prep(res_up, "Up-regulated"), prep(res_down, "Down-regulated"))
  if (is.null(out) || !nrow(out)) stop("No significant terms in either direction to plot.")
  out$Direction <- factor(out$Direction, levels = c("Up-regulated", "Down-regulated"))
  out <- out[order(out$Direction, out$p.adjust), ]
  out$row_key <- factor(out$row_key, levels = rev(out$row_key))
  out
}

.row_key_labels <- function(cmb) {
  stats::setNames(sub("^(Up-regulated|Down-regulated)\r", "", as.character(cmb$row_key)),
                  as.character(cmb$row_key))
}

make_split_dotplot <- function(res_up, res_down, top_n = 15, style = default_style(),
                               title = NULL, xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  cmb <- split_enrich_df(res_up, res_down, top_n)
  p <- ggplot2::ggplot(cmb, ggplot2::aes(x = GeneRatioNum, y = row_key,
                                         size = Count, colour = p.adjust)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~Direction, scales = "free_y", ncol = 1) +
    ggplot2::scale_y_discrete(labels = .row_key_labels(cmb)) +
    ggplot2::scale_colour_gradient(low = "red", high = "blue", name = "adj. P") +
    deg_theme(st)
  apply_labels(p, title, NULL, xlab, ylab,
               default_title = "Enrichment results", default_x = "Gene ratio", default_y = NULL)
}

make_split_barplot <- function(res_up, res_down, top_n = 15, style = default_style(),
                               up_colour = DEG_COLOURS$up, down_colour = DEG_COLOURS$down,
                               title = NULL, xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  cmb <- split_enrich_df(res_up, res_down, top_n)
  p <- ggplot2::ggplot(cmb, ggplot2::aes(x = Count, y = row_key, fill = Direction)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~Direction, scales = "free_y", ncol = 1) +
    ggplot2::scale_y_discrete(labels = .row_key_labels(cmb)) +
    ggplot2::scale_fill_manual(values = c("Up-regulated" = up_colour,
                                          "Down-regulated" = down_colour), guide = "none") +
    deg_theme(st)
  apply_labels(p, title, NULL, xlab, ylab,
               default_title = "Enrichment results", default_x = "Gene count", default_y = NULL)
}

make_readable <- function(res, database) {
  # Swap the ids in geneID for gene symbols (cnet/heat plots label genes with
  # them). The keytype slot clusterProfiler stores does not always describe the
  # ids in geneID, so the ids themselves decide: each candidate keytype is
  # tried against the OrgDb and the one that maps the most ids is used.
  db <- attr(res, "dex_orgdb") %||% "org.Hs.eg.db"
  if (identical(attr(res, "dex_keytype"), "SYMBOL")) return(res)
  if (is.null(res) || !nrow(as.data.frame(res))) return(res)
  if (isTRUE(tryCatch(res@readable, error = function(e) FALSE))) return(res)
  if (!requireNamespace(db, quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) return(res)
  orgdb <- getExportedValue(db, db)
  ids <- unique(unlist(strsplit(res@result$geneID, "/", fixed = TRUE)))
  ids <- ids[nzchar(ids)]
  slot_kt <- tryCatch(res@keytype, error = function(e) NULL)
  cands <- unique(c(attr(res, "dex_keytype"), slot_kt,
                    if (database == "KEGG") "ENTREZID" else c("ENSEMBL", "ENTREZID")))
  cands <- intersect(cands, setdiff(AnnotationDbi::keytypes(orgdb), "SYMBOL"))
  best <- NULL
  for (kt in cands) {
    sym <- tryCatch(suppressMessages(AnnotationDbi::mapIds(orgdb, keys = ids, column = "SYMBOL",
                                                           keytype = kt, multiVals = "first")),
                    error = function(e) NULL)
    if (!is.null(sym) && (is.null(best) || sum(!is.na(sym)) > sum(!is.na(best$sym))))
      best <- list(kt = kt, sym = sym)
  }
  if (is.null(best) || !any(!is.na(best$sym))) return(res)
  g2s <- stats::setNames(ifelse(is.na(best$sym), names(best$sym), unname(best$sym)), names(best$sym))
  res@result$geneID <- vapply(strsplit(res@result$geneID, "/", fixed = TRUE),
                              function(k) paste(g2s[k], collapse = "/"), character(1))
  res@gene2Symbol <- g2s
  res@keytype <- best$kt
  res@readable <- TRUE
  res
}

make_cnetplot <- function(res, database, fc_named = NULL, top_n = 8,
                          style = default_style(), title = "") {
  if (!HAS_ENRICHPLOT) stop("The enrichplot package is required for this plot type.")
  if (is.null(res) || !nrow(as.data.frame(res))) stop("No significant terms to plot.")
  r <- make_readable(res, database)
  p <- tryCatch(
    enrichplot::cnetplot(r, showCategory = top_n, color.params = list(foldChange = fc_named)),
    error = function(e) tryCatch(
      enrichplot::cnetplot(r, showCategory = top_n, foldChange = fc_named),
      error = function(e2) enrichplot::cnetplot(r, showCategory = top_n)))
  p + ggplot2::ggtitle(title)
}

make_emapplot <- function(res, top_n = 20, style = default_style(), title = "") {
  if (!HAS_ENRICHPLOT) stop("The enrichplot package is required for this plot type.")
  df <- as.data.frame(res)
  if (nrow(df) < 2) stop("At least 2 significant terms are needed for an enrichment map.")
  r2 <- tryCatch(enrichplot::pairwise_termsim(res),
                 error = function(e) stop("Could not compute term similarity: ", conditionMessage(e)))
  enrichplot::emapplot(r2, showCategory = min(top_n, nrow(df))) + ggplot2::ggtitle(title)
}

make_heatplot <- function(res, database, fc_named = NULL, top_n = 20,
                          style = default_style(), title = "") {
  if (!HAS_ENRICHPLOT) stop("The enrichplot package is required for this plot type.")
  if (is.null(res) || !nrow(as.data.frame(res))) stop("No significant terms to plot.")
  r <- make_readable(res, database)
  p <- tryCatch(enrichplot::heatplot(r, showCategory = top_n, foldChange = fc_named),
                error = function(e) enrichplot::heatplot(r, showCategory = top_n))
  st <- resolve_style(style)
  p + ggplot2::ggtitle(title) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1,
                                                       size = st$axis_text_size * 0.7))
}

combine_stacked <- function(p_top, p_bottom, title = NULL) {
  if (is.null(p_bottom)) return(p_top)
  if (!HAS_GRIDEXTRA) {
    warning("gridExtra is not installed; showing the up-regulated panel only.")
    return(p_top)
  }
  gridExtra::arrangeGrob(p_top, p_bottom, ncol = 1,
    top = if (!is.null(title)) grid::textGrob(title, gp = grid::gpar(fontface = "bold", fontsize = 14)) else NULL)
}

#' Dispatcher used by the Enrichment tab.
build_enrichment_figure <- function(plot_type, database, res_single = NULL,
                                    res_up = NULL, res_down = NULL, split = FALSE,
                                    fc_named = NULL, top_n = 15,
                                    style = default_style(),
                                    title = "Enrichment results",
                                    xlab = NULL, ylab = NULL) {
  if (split) {
    if (plot_type == "dot") return(make_split_dotplot(res_up, res_down, top_n, style, title, xlab, ylab))
    if (plot_type == "bar") return(make_split_barplot(res_up, res_down, top_n, style, title = title, xlab = xlab, ylab = ylab))
    one <- function(r, d) {
      if (is.null(r) || !nrow(as.data.frame(r))) return(NULL)
      sub <- paste0(title, " \u2014 ", d)
      switch(plot_type,
        cnet = make_cnetplot(r, database, fc_named, ceiling(top_n / 2), style, sub),
        emap = tryCatch(make_emapplot(r, top_n, style, sub), error = function(e) NULL),
        heat = make_heatplot(r, database, fc_named, top_n, style, sub),
        stop("Unknown plot type: ", plot_type))
    }
    pu <- one(res_up, "Up-regulated"); pd <- one(res_down, "Down-regulated")
    if (is.null(pu) && is.null(pd)) stop("No significant terms in either direction to plot.")
    if (is.null(pu)) return(pd)
    if (is.null(pd)) return(pu)
    return(combine_stacked(pu, pd))
  }
  switch(plot_type,
    dot  = make_dotplot(res_single, top_n, style, title, xlab, ylab),
    bar  = make_barplot(res_single, top_n, style, title, xlab, ylab),
    cnet = make_cnetplot(res_single, database, fc_named, top_n, style, title),
    emap = make_emapplot(res_single, top_n, style, title),
    heat = make_heatplot(res_single, database, fc_named, top_n, style, title),
    stop("Unknown plot type: ", plot_type))
}
