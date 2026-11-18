# =============================================================================
# gene_report.R -- external gene knowledge and literature
#
# Two layers, deliberately separated:
#
#   1. Link layer (always available, no network needed). Builds deep links
#      into the major gene resources for every gene of interest. This is the
#      dependable core of the feature.
#
#   2. Live layer (optional, needs internet). Pulls gene summaries from
#      MyGene.info and publications from Europe PMC. Every call is wrapped so
#      a blocked network, proxy or offline laptop degrades to the link layer
#      instead of erroring out.
#
# `build_gene_report()` combines whatever is available into a single
# self-contained HTML document the user can download.
# =============================================================================

HAS_JSONLITE <- requireNamespace("jsonlite", quietly = TRUE)
HAS_CURL     <- requireNamespace("curl", quietly = TRUE)

# ---- link layer -------------------------------------------------------------

.url_encode <- function(x) utils::URLencode(as.character(x), reserved = TRUE)

#' Deep links for one gene. `symbol` is preferred where a resource is
#' symbol-keyed; `ensembl` is used where an Ensembl id is the better key.
gene_resource_links <- function(symbol, ensembl = NA_character_,
                                species = "Homo sapiens") {
  sym <- if (is.na(symbol) || !nzchar(symbol)) ensembl else symbol
  s <- .url_encode(sym)
  ens_ok <- !is.na(ensembl) && grepl("^ENS", ensembl)
  e <- .url_encode(ensembl)
  sp <- .url_encode(species)

  human <- identical(species, "Homo sapiens")
  org_key <- names(ORGANISMS)[vapply(ORGANISMS, function(o) identical(o$species, species), logical(1))][1]
  kegg <- if (!is.na(org_key)) ORGANISMS[[org_key]]$kegg else ""
  ens_ok <- ens_ok || (!is.na(ensembl) && !is.na(org_key) &&
    nzchar(ORGANISMS[[org_key]]$gene_prefix) && startsWith(ensembl, ORGANISMS[[org_key]]$gene_prefix))
  links <- list(
    "NCBI Gene"    = sprintf("https://www.ncbi.nlm.nih.gov/gene/?term=%s+AND+%s[orgn]", s, .url_encode(gsub(" ", "+", species))),
    "GeneCards"    = if (human) sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", s),
    "UniProt"      = sprintf("https://www.uniprot.org/uniprotkb?query=%s+AND+organism_name:%%22%s%%22", s, sp),
    "Open Targets" = if (human) { if (ens_ok) sprintf("https://platform.opentargets.org/target/%s", e)
                     else sprintf("https://platform.opentargets.org/search?q=%s", s) },
    "STRING"       = sprintf("https://string-db.org/cgi-bin/network.pl?identifier=%s&species=%s", s, .url_encode(species)),
    "GTEx"         = if (human) sprintf("https://gtexportal.org/home/gene/%s", s),
    "Human Protein Atlas" = if (human) { if (ens_ok) sprintf("https://www.proteinatlas.org/%s", e)
                            else sprintf("https://www.proteinatlas.org/search/%s", s) },
    "PubMed"       = sprintf("https://pubmed.ncbi.nlm.nih.gov/?term=%s", s),
    "Europe PMC"   = sprintf("https://europepmc.org/search?query=%s", s),
    "Reactome"     = sprintf("https://reactome.org/content/query?q=%s", s),
    "KEGG"         = if (nzchar(kegg)) sprintf("https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=50&dbkey=%s&keywords=%s", kegg, s),
    "AlphaFold"    = sprintf("https://alphafold.ebi.ac.uk/search/text/%s", s)
  )
  if (ens_ok) {
    host <- if (!is.na(org_key) && identical(ORGANISMS[[org_key]]$division, "plants")) "plants.ensembl.org" else "www.ensembl.org"
    links[["Ensembl"]] <- sprintf("https://%s/%s/Gene/Summary?g=%s", host, gsub(" ", "_", species), e)
  }
  Filter(Negate(is.null), links)
}

#' A PubMed / Europe PMC query for a gene in the context of the study.
#' `context_terms` lets the user bias the search, e.g. "T cell, interferon".
literature_query <- function(symbol, context_terms = "") {
  ct <- trimws(context_terms)
  if (!nzchar(ct)) return(symbol)
  parts <- trimws(unlist(strsplit(ct, "[,;]")))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(symbol)
  sprintf('%s AND (%s)', symbol, paste(sprintf('"%s"', parts), collapse = " OR "))
}

#' Table of gene + links, suitable for rendering with DT.
gene_link_table <- function(genes, de = NULL, species = "Homo sapiens",
                            context_terms = "") {
  if (!length(genes)) return(data.frame())
  if (!is.null(de)) {
    idx <- match(genes, de$gene_id)
    sym <- de$display_gene[idx]
    ens <- de$gene_id[idx]
    lfc <- de$log2FoldChange[idx]
    padj <- de$padj[idx]
  } else {
    sym <- genes; ens <- genes
    lfc <- rep(NA_real_, length(genes)); padj <- rep(NA_real_, length(genes))
  }
  sym[is.na(sym)] <- genes[is.na(sym)]

  mk <- function(i) {
    L <- gene_resource_links(sym[i], ens[i], species)
    paste(sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', unlist(L), names(L)),
          collapse = " &middot; ")
  }
  pm <- vapply(seq_along(genes), function(i) {
    q <- .url_encode(literature_query(sym[i], context_terms))
    sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/?term=%s" target="_blank" rel="noopener">search</a>', q)
  }, character(1))

  fmt_p <- function(x) {
    ifelse(is.na(x), "-",
           ifelse(x < 1e-4, formatC(x, format = "e", digits = 2),
                  formatC(x, format = "g", digits = 3)))
  }
  data.frame(
    Gene = sym,
    `Ensembl ID` = ens,
    log2FC = round(lfc, 3),
    `adj. P` = fmt_p(padj),
    Literature = pm,
    Resources = vapply(seq_along(genes), mk, character(1)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# ---- live layer -------------------------------------------------------------

.http_get_json <- function(url, timeout = 15) {
  if (!HAS_JSONLITE) stop("The jsonlite package is required for live lookups.")
  txt <- NULL
  if (HAS_CURL) {
    h <- curl::new_handle(timeout = timeout, connecttimeout = timeout,
                          useragent = "DEGexplorer/1.0 (R Shiny)")
    res <- curl::curl_fetch_memory(url, handle = h)
    if (res$status_code >= 400) stop("HTTP ", res$status_code)
    txt <- rawToChar(res$content)
  } else {
    con <- url(url, "rt")
    txt <- paste(readLines(con, warn = FALSE), collapse = "")
    close(con)
  }
  jsonlite::fromJSON(txt, simplifyVector = TRUE)
}

#' Gene summaries from MyGene.info (https://mygene.info).
#' Returns a data frame; on any failure returns NULL and sets an attribute
#' describing what went wrong, so the caller can inform the user.
fetch_gene_summaries <- function(symbols, species = "human", timeout = 20) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) return(NULL)
  symbols <- head(symbols, 200)
  out <- tryCatch({
    q <- paste(symbols, collapse = ",")
    url <- sprintf("https://mygene.info/v3/query?q=%s&scopes=symbol,ensembl.gene&species=%s&fields=symbol,name,summary,entrezgene,type_of_gene,genomic_pos&size=1",
                   .url_encode(q), .url_encode(species))
    res <- .http_get_json(url, timeout)
    if (is.null(res) || !length(res)) return(NULL)
    df <- as.data.frame(res, stringsAsFactors = FALSE)
    keep <- intersect(c("query", "symbol", "name", "summary", "entrezgene", "type_of_gene"),
                      colnames(df))
    df <- df[, keep, drop = FALSE]
    df
  }, error = function(e) {
    structure(list(), class = "degx_fetch_error", failure = conditionMessage(e))
  })
  out
}

#' Recent publications for a gene from Europe PMC.
fetch_publications <- function(symbol, context_terms = "", n = 5, timeout = 20) {
  q <- literature_query(symbol, context_terms)
  tryCatch({
    url <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%s&format=json&pageSize=%d&sort=CITED%%20desc",
                   .url_encode(q), as.integer(n))
    res <- .http_get_json(url, timeout)
    hits <- res$resultList$result
    if (is.null(hits) || !length(hits)) return(NULL)
    df <- as.data.frame(hits, stringsAsFactors = FALSE)
    cols <- intersect(c("id", "source", "pmid", "doi", "title", "authorString",
                        "journalTitle", "pubYear", "citedByCount"), colnames(df))
    df <- df[, cols, drop = FALSE]
    df$url <- if ("pmid" %in% colnames(df)) {
      ifelse(is.na(df$pmid) | df$pmid == "",
             sprintf("https://europepmc.org/article/%s/%s", df$source, df$id),
             sprintf("https://pubmed.ncbi.nlm.nih.gov/%s/", df$pmid))
    } else sprintf("https://europepmc.org/article/%s/%s", df$source, df$id)
    df
  }, error = function(e) structure(list(), class = "degx_fetch_error",
                                   failure = conditionMessage(e)))
}

#' TRUE when a live-fetch helper returned an error sentinel.
is_fetch_error <- function(x) inherits(x, "degx_fetch_error")

#' The message attached to a fetch-error sentinel (or NULL).
fetch_error_message <- function(x) if (is_fetch_error(x)) attr(x, "failure") else NULL

# ---- combined report --------------------------------------------------------

.esc <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' Assemble a standalone HTML dossier for a set of genes.
#'
#' @param genes        gene ids (matching `de$gene_id`)
#' @param de           canonical DE table
#' @param summaries    optional data frame from `fetch_gene_summaries()`
#' @param publications optional named list of data frames from `fetch_publications()`
#' @param meta_info    named character vector of study details for the header
build_gene_report <- function(genes, de, summaries = NULL, publications = NULL,
                              species = "Homo sapiens", context_terms = "",
                              meta_info = character(0),
                              title = "Differentially expressed gene report") {
  idx <- match(genes, de$gene_id)
  sym <- de$display_gene[idx]; sym[is.na(sym)] <- genes[is.na(sym)]
  lfc <- de$log2FoldChange[idx]; padj <- de$padj[idx]; bm <- de$baseMean[idx]

  css <- "
  body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;
       line-height:1.55;margin:0;padding:32px;color:#1f2933;background:#fff;max-width:1100px}
  h1{font-size:1.7rem;margin:0 0 4px;color:#14202b}
  h2{font-size:1.15rem;margin:34px 0 10px;padding-bottom:6px;border-bottom:2px solid #e3e8ee;color:#14202b}
  h3{font-size:1rem;margin:22px 0 6px;color:#1f3d5a}
  .sub{color:#66737f;margin:0 0 20px;font-size:.92rem}
  table{border-collapse:collapse;width:100%;margin:10px 0 18px;font-size:.9rem}
  th,td{border:1px solid #dde3ea;padding:7px 9px;text-align:left;vertical-align:top}
  th{background:#f4f7fa;font-weight:600}
  tr:nth-child(even) td{background:#fafbfc}
  .up{color:#c0392b;font-weight:600}.down{color:#2471a3;font-weight:600}
  .meta{background:#f4f7fa;border-left:4px solid #5b8db8;padding:10px 14px;margin:14px 0;font-size:.9rem}
  .links a{margin-right:9px;white-space:nowrap}
  a{color:#1f6feb;text-decoration:none}a:hover{text-decoration:underline}
  .summary{background:#fcfdfe;border:1px solid #e3e8ee;border-radius:5px;padding:9px 12px;margin:7px 0;font-size:.9rem}
  .note{color:#8a6d3b;background:#fcf8e3;border-left:4px solid #dcc27a;padding:9px 13px;margin:14px 0;font-size:.88rem}
  footer{margin-top:40px;padding-top:14px;border-top:1px solid #e3e8ee;color:#8a959f;font-size:.82rem}
  "
  h <- c('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">',
         sprintf('<title>%s</title>', .esc(title)),
         sprintf('<style>%s</style></head><body>', css),
         sprintf('<h1>%s</h1>', .esc(title)),
         sprintf('<p class="sub">%d gene(s) &middot; %s &middot; generated %s</p>',
                 length(genes), .esc(species), format(Sys.time(), "%Y-%m-%d %H:%M")))

  if (length(meta_info)) {
    h <- c(h, '<div class="meta">',
           paste(sprintf("<strong>%s:</strong> %s", .esc(names(meta_info)), .esc(meta_info)),
                 collapse = "<br>"),
           '</div>')
  }

  # summary table
  h <- c(h, "<h2>Differential expression summary</h2>",
         "<table><tr><th>#</th><th>Gene</th><th>Ensembl ID</th><th>log2 FC</th><th>adj. P</th><th>Mean expr.</th><th>Direction</th></tr>")
  for (i in seq_along(genes)) {
    dir <- if (!is.na(lfc[i]) && lfc[i] > 0) '<span class="up">Up</span>'
           else if (!is.na(lfc[i])) '<span class="down">Down</span>' else "-"
    h <- c(h, sprintf("<tr><td>%d</td><td><strong>%s</strong></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
                      i, .esc(sym[i]), .esc(genes[i]),
                      ifelse(is.na(lfc[i]), "-", sprintf("%.3f", lfc[i])),
                      ifelse(is.na(padj[i]), "-", formatC(padj[i], format = "e", digits = 2)),
                      ifelse(is.na(bm[i]), "-", format(round(bm[i], 1), big.mark = ",")),
                      dir))
  }
  h <- c(h, "</table>")

  live_failed <- is_fetch_error(summaries)
  have_summaries <- !live_failed && is.data.frame(summaries) && nrow(summaries) > 0

  # per-gene detail
  h <- c(h, "<h2>Gene dossiers</h2>")
  for (i in seq_along(genes)) {
    h <- c(h, sprintf('<h3>%d. %s <span style="font-weight:400;color:#8a959f">(%s)</span></h3>',
                      i, .esc(sym[i]), .esc(genes[i])))
    if (have_summaries) {
      row <- summaries[summaries$symbol == sym[i] |
                       (!is.na(summaries$query) & summaries$query == genes[i]), , drop = FALSE]
      if (nrow(row)) {
        if ("name" %in% colnames(row) && !is.na(row$name[1]))
          h <- c(h, sprintf('<p><em>%s</em></p>', .esc(row$name[1])))
        if ("summary" %in% colnames(row) && !is.na(row$summary[1]) && nzchar(row$summary[1]))
          h <- c(h, sprintf('<div class="summary">%s</div>', .esc(row$summary[1])))
      }
    }
    L <- gene_resource_links(sym[i], genes[i], species)
    h <- c(h, '<p class="links">',
           paste(sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>',
                         unlist(L), .esc(names(L))), collapse = " "),
           '</p>')

    if (!is.null(publications) && !is.null(publications[[sym[i]]])) {
      pubs <- publications[[sym[i]]]
      if (is.data.frame(pubs) && nrow(pubs)) {
        h <- c(h, "<table><tr><th>Publication</th><th>Journal</th><th>Year</th><th>Cited</th></tr>")
        for (j in seq_len(nrow(pubs))) {
          h <- c(h, sprintf('<tr><td><a href="%s" target="_blank" rel="noopener">%s</a><br><span style="color:#8a959f;font-size:.85em">%s</span></td><td>%s</td><td>%s</td><td>%s</td></tr>',
                            .esc(pubs$url[j]), .esc(pubs$title[j]),
                            .esc(if ("authorString" %in% colnames(pubs)) pubs$authorString[j] else ""),
                            .esc(if ("journalTitle" %in% colnames(pubs)) pubs$journalTitle[j] else ""),
                            .esc(if ("pubYear" %in% colnames(pubs)) pubs$pubYear[j] else ""),
                            .esc(if ("citedByCount" %in% colnames(pubs)) pubs$citedByCount[j] else "")))
        }
        h <- c(h, "</table>")
      }
    } else {
      q <- .url_encode(literature_query(sym[i], context_terms))
      h <- c(h, sprintf('<p style="font-size:.9rem">Literature: <a href="https://pubmed.ncbi.nlm.nih.gov/?term=%s" target="_blank" rel="noopener">PubMed</a> &middot; <a href="https://europepmc.org/search?query=%s" target="_blank" rel="noopener">Europe PMC</a></p>', q, q))
    }
  }

  if (live_failed) {
    h <- c(h, sprintf('<div class="note"><strong>Live lookups unavailable.</strong> Gene descriptions and publication lists could not be retrieved (%s). All resource links above remain fully functional.</div>',
                      .esc(fetch_error_message(summaries))))
  }

  h <- c(h, '<footer>Generated by DEGexplorer. Gene descriptions courtesy of MyGene.info; publication metadata courtesy of Europe PMC. Always confirm findings against the primary source.</footer>',
         "</body></html>")
  paste(h, collapse = "\n")
}
