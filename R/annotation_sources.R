# =============================================================================
# annotation_sources.R -- fetch a GTF from Ensembl or GENCODE, and derive gene
# lengths from any GTF
#
# One GTF gives both things the app needs when the pipeline's own files are
# missing: gene symbols / biotypes (via read_gtf_annotation) and gene lengths
# for TPM / FPKM (union of each gene's exons).
# =============================================================================

# Species come from the shared ORGANISMS registry (R/organisms.R); "other" is any
# Ensembl or Ensembl Genomes species, named by its directory (e.g. ovis_aries).
ANNOTATION_SPECIES <- ORGANISMS

ANNOTATION_PROVIDERS <- c("Ensembl" = "ensembl", "GENCODE" = "gencode")

#' Suggest the annotation a count matrix was built with, from its gene ids.
#'
#' Human vs mouse comes from the Ensembl id prefix; GENCODE writes version
#' suffixes (ENSG00000223972.5) where Ensembl does not. The release cannot be
#' read from the ids, so the newest one is suggested.
guess_annotation_source <- function(ids) {
  ids <- utils::head(as.character(ids[!is.na(ids)]), 5000)
  sp <- guess_organism(ids)
  versioned <- mean(grepl("^ENS[A-Z]*[GT][0-9]+\\.[0-9]+$", ids)) > 0.5
  gencode <- versioned && !is.na(sp) && sp %in% c("human", "mouse")
  list(species = sp, provider = if (gencode) "gencode" else "ensembl",
       note = if (is.na(sp)) "The ids do not say which organism this is; choose it yourself (Other = any Ensembl species)." else
         sprintf("Ids look like %s %s%s.", if (gencode) "GENCODE" else "Ensembl",
                 ORGANISMS[[sp]]$label, if (gencode) " (version suffixes)" else ""))
}

#' Download address of a GTF.
#'
#' @param release "current" (newest) or a release number: Ensembl 103+ (e.g.
#'   "115"), GENCODE human (e.g. "48") or mouse (e.g. "M37" or "37").
#' @param listing for "current": a function(url) returning the directory
#'   listing as text, used to find the newest file name (replaceable in tests).
annotation_url <- function(provider, species, release = "current",
                           listing = function(u) paste(readLines(u, warn = FALSE), collapse = "\n"),
                           other = NULL) {
  sp <- organism_info(species, other)
  if (identical(species, "other") && !nzchar(sp$ensembl_dir))
    stop("Give the organism's Ensembl name, e.g. ovis_aries (see ensembl.org/info/about/species.html).")
  release <- trimws(as.character(release %||% "current"))
  if (!nzchar(release)) release <- "current"
  latest <- tolower(release) %in% c("current", "latest")

  if (identical(provider, "ensembl")) {
    div <- sp$division %||% ""
    div_label <- if (nzchar(div)) paste0("Ensembl ", tools::toTitleCase(div)) else "Ensembl"
    rel_dir <- if (latest) "current" else {
      rel <- suppressWarnings(as.integer(release))
      if (is.na(rel)) stop("An Ensembl release is a number, e.g. 115.")
      if (identical(species, "mouse") && !nzchar(div) && rel < 103) stop("Mouse GRCm39 starts at Ensembl release 103.")
      paste0("release-", rel)
    }
    dir <- if (nzchar(div)) sprintf("https://ftp.ensemblgenomes.ebi.ac.uk/pub/%s/%s/gtf/%s/", div, rel_dir, sp$ensembl_dir)
           else sprintf("https://ftp.ensembl.org/pub/%s/gtf/%s/", rel_dir, sp$ensembl_dir)
    # Known assembly and a pinned release: the address is fixed, no listing needed.
    if (!latest && !is.null(sp$ensembl_file) && !nzchar(div)) {
      f <- sprintf("%s.%s.gtf.gz", sp$ensembl_file, sub("^release-", "", rel_dir))
      return(list(url = paste0(dir, f), file = f,
                  label = sprintf("%s release %s, %s", div_label, sub("^release-", "", rel_dir), sp$label)))
    }
    # Otherwise read the directory: files are <Species>.<assembly>.<release>.gtf.gz
    # (the .chr / .abinitio / .chr_patch_hapl_scaff variants are skipped).
    prefix <- paste0(toupper(substr(sp$ensembl_dir, 1, 1)), substring(sp$ensembl_dir, 2))
    pat <- sprintf("%s\\.([A-Za-z0-9_.-]+?)\\.([0-9]+)\\.gtf\\.gz", gsub(".", "\\.", prefix, fixed = TRUE))
    txt <- tryCatch(listing(dir), error = function(e) "", warning = function(w) "")
    hits <- unique(regmatches(txt, gregexpr(pat, txt, ignore.case = TRUE))[[1]])
    if (!length(hits)) stop("Could not find a GTF in ", dir,
                            " -- check the organism name, division and release.")
    rels <- as.integer(sub(pat, "\\2", hits, ignore.case = TRUE))
    f <- hits[which.max(rels)]
    return(list(url = paste0(dir, f), file = f,
                label = sprintf("%s release %d, %s", div_label, max(rels), sp$label)))
  }

  if (identical(provider, "gencode")) {
    if (!species %in% c("human", "mouse")) stop("GENCODE covers human and mouse only; use Ensembl for ", sp$label, ".")
    base <- sprintf("https://ftp.ebi.ac.uk/pub/databases/gencode/%s/", sp$gencode_dir)
    if (latest) {
      dir <- paste0(base, "latest_release/")
      pat <- sprintf("gencode\\.v%s([0-9]+)\\.annotation\\.gtf\\.gz", sp$gencode_prefix)
      txt <- listing(dir)
      hits <- regmatches(txt, gregexpr(pat, txt))[[1]]
      if (!length(hits)) stop("Could not find the latest GENCODE GTF in ", dir)
      v <- max(as.integer(sub(pat, "\\1", hits)))
      f <- sprintf("gencode.v%s%d.annotation.gtf.gz", sp$gencode_prefix, v)
      return(list(url = paste0(dir, f), file = f,
                  label = sprintf("GENCODE %s%d, %s", sp$gencode_prefix, v, sp$label)))
    }
    v <- suppressWarnings(as.integer(sub("^[Mm]", "", release)))
    if (is.na(v)) stop("A GENCODE release is a number, e.g. 48 (human) or M37 / 37 (mouse).")
    f <- sprintf("gencode.v%s%d.annotation.gtf.gz", sp$gencode_prefix, v)
    return(list(url = sprintf("%srelease_%s%d/%s", base, sp$gencode_prefix, v, f), file = f,
                label = sprintf("GENCODE %s%d, %s", sp$gencode_prefix, v, sp$label)))
  }
  stop("Unknown provider: ", provider)
}


#' Folder where downloaded GTFs are kept, so each is fetched only once.
annotation_cache_dir <- function() {
  d <- tryCatch(tools::R_user_dir("DEGexplorer", "cache"), error = function(e) NULL)
  if (is.null(d)) d <- file.path(tempdir(), "DEGexplorer_cache")
  d <- file.path(d, "annotation")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Download a GTF (or reuse the cached copy) and return its local path.
#'
#' @param progress optional function(mb_done, mb_total) called as the file
#'   arrives (needs the curl package; mb_total is NA when the server does not
#'   say). Ensembl's server can be slow, so this keeps the progress bar honest.
download_annotation <- function(src, cache_dir = annotation_cache_dir(), progress = NULL) {
  dest <- file.path(cache_dir, src$file)
  if (file.exists(dest) && file.size(dest) > 1e6) return(dest)
  tmp <- paste0(dest, ".part")
  unlink(tmp)
  ok <- if (is.function(progress) && requireNamespace("curl", quietly = TRUE)) {
    last <- 0
    h <- curl::new_handle(noprogress = FALSE, progressfunction = function(down, up) {
      now <- down[2] / 1e6
      if (now - last >= 1) {            # report every MB
        last <<- now
        progress(now, if (down[1] > 0) down[1] / 1e6 else NA_real_)
      }
      TRUE
    })
    tryCatch({
      r <- curl::curl_fetch_disk(src$url, tmp, handle = h)
      if (r$status_code != 200) paste("the server answered HTTP", r$status_code) else TRUE
    }, error = function(e) conditionMessage(e))
  } else {
    old <- options(timeout = max(1800, getOption("timeout"))); on.exit(options(old), add = TRUE)
    tryCatch({ utils::download.file(src$url, tmp, mode = "wb", quiet = TRUE); TRUE },
             error = function(e) conditionMessage(e), warning = function(w) conditionMessage(w))
  }
  if (!isTRUE(ok) || !file.exists(tmp) || file.size(tmp) < 1e6) {
    unlink(tmp)
    stop("Download failed from ", src$url, if (is.character(ok)) paste0(": ", ok) else "")
  }
  file.rename(tmp, dest)
  dest
}

#' Gene lengths from a GTF: the total length of the union of each gene's exons.
#'
#' This is the length featureCounts reports, and what TPM / FPKM need when the
#' pipeline's per-sample effective lengths (e.g. rsem.merged.gene_lengths.tsv)
#' are not at hand.
#' @return named numeric vector (names = gene ids)
gtf_gene_lengths <- function(path) {
  cols <- readr::cols_only(X3 = readr::col_character(), X4 = readr::col_double(),
                           X5 = readr::col_double(), X9 = readr::col_character())
  cb <- readr::DataFrameCallback$new(function(x, pos) {
    x <- x[!is.na(x$X3) & x$X3 == "exon", , drop = FALSE]
    if (!nrow(x)) return(NULL)
    gff3 <- .is_gff3(x$X9)
    key <- if (gff3) "Parent" else "gene_id"
    data.frame(gene = .gtf_attr(x$X9, key, gff3), start = x$X4, end = x$X5,
               stringsAsFactors = FALSE)
  })
  ex <- suppressWarnings(readr::read_tsv_chunked(path, cb, chunk_size = 2e5, col_names = FALSE,
                                                 col_types = cols, comment = "#", progress = FALSE))
  if (is.null(ex) || !nrow(ex)) stop("No exon records found in ", basename(path), ".")
  ex <- ex[!is.na(ex$gene), , drop = FALSE]
  if (all(grepl("^transcript:|^rna-", ex$gene[seq_len(min(100, nrow(ex)))]))) {
    stop("This GFF3 links exons to transcripts, not genes; use the GTF version of the annotation.")
  }
  ex <- ex[order(ex$gene, ex$start), , drop = FALSE]
  # Running end of the exons seen so far in the same gene; each exon adds only
  # the part that extends past it, which merges overlaps from all transcripts.
  run_end <- stats::ave(ex$end, ex$gene, FUN = cummax)
  prev <- stats::ave(run_end, ex$gene, FUN = function(z) c(-Inf, utils::head(z, -1)))
  add <- pmax(0, ex$end - pmax(ex$start - 1, prev))
  len <- tapply(add, ex$gene, sum)
  stats::setNames(as.numeric(len), names(len))
}

#' Transcripts from a GTF: transcript -> gene, names, biotype and length.
#'
#' Needed for transcript-level analysis (gene symbols for each transcript, the
#' transcript-to-gene map for EBSeq's Ng vector and stage-wise testing).
#' Length is the sum of the transcript's exons.
#' @return data.frame(transcript_id, gene_id, gene_name, transcript_name,
#'   transcript_biotype, tx_length)
gtf_transcript_table <- function(path) {
  cols <- readr::cols_only(X3 = readr::col_character(), X4 = readr::col_double(),
                           X5 = readr::col_double(), X9 = readr::col_character())
  cb <- readr::DataFrameCallback$new(function(x, pos) {
    x <- x[!is.na(x$X3) & x$X3 == "exon", , drop = FALSE]
    if (!nrow(x)) return(NULL)
    gff3 <- .is_gff3(x$X9)
    if (gff3) return(data.frame(transcript_id = sub("^transcript:", "", .gtf_attr(x$X9, "Parent", TRUE)),
                                gene_id = NA_character_, gene_name = NA_character_,
                                transcript_name = NA_character_, transcript_biotype = NA_character_,
                                w = x$X5 - x$X4 + 1, stringsAsFactors = FALSE))
    bt <- .gtf_attr(x$X9, "transcript_biotype", FALSE)
    if (all(is.na(bt))) bt <- .gtf_attr(x$X9, "transcript_type", FALSE)
    data.frame(transcript_id = .gtf_attr(x$X9, "transcript_id", FALSE),
               gene_id = .gtf_attr(x$X9, "gene_id", FALSE),
               gene_name = .gtf_attr(x$X9, "gene_name", FALSE),
               transcript_name = .gtf_attr(x$X9, "transcript_name", FALSE),
               transcript_biotype = bt, w = x$X5 - x$X4 + 1, stringsAsFactors = FALSE)
  })
  ex <- suppressWarnings(readr::read_tsv_chunked(path, cb, chunk_size = 2e5, col_names = FALSE,
                                                 col_types = cols, comment = "#", progress = FALSE))
  if (is.null(ex) || !nrow(ex)) stop("No exon records found in ", basename(path), ".")
  ex <- ex[!is.na(ex$transcript_id), , drop = FALSE]
  len <- tapply(ex$w, ex$transcript_id, sum)
  first <- ex[!duplicated(ex$transcript_id), setdiff(colnames(ex), "w"), drop = FALSE]
  first$tx_length <- as.numeric(len[first$transcript_id])
  rownames(first) <- NULL
  first
}

#' Symbols, biotypes and gene lengths from a cached GTF, parsed once.
#'
#' Parsing a full human GTF takes a few minutes, so the result is saved next to
#' the downloaded file and reused.
load_annotation_bundle <- function(path, progress = NULL) {
  tick <- function(f, d) if (is.function(progress)) progress(f, d)
  rds <- sub("\\.gtf(\\.gz)?$", ".parsed.rds", path)
  b <- NULL
  if (!identical(rds, path) && file.exists(rds)) {
    b <- tryCatch(readRDS(rds), error = function(e) NULL)
    if (!is.null(b$ann) && !is.null(b$len) && !is.null(b$tx)) return(b)
  }
  if (is.null(b$ann) || is.null(b$len)) {
    tick(0.6, "Reading gene symbols and biotypes (first time only, a few minutes)")
    ann <- read_gtf_annotation(path)
    tick(0.75, "Computing gene lengths (first time only)")
    len <- gtf_gene_lengths(path)
    b <- list(ann = ann, len = len)
  }
  tick(0.85, "Reading transcripts (first time only)")
  b$tx <- tryCatch(gtf_transcript_table(path), error = function(e) NULL)
  if (!identical(rds, path)) try(saveRDS(b, rds), silent = TRUE)
  b
}

#' Re-key a gene-length vector to a matrix's ids (tolerating version suffixes).
#' @return named vector over the matrix ids that have a length
match_gene_lengths <- function(lengths, ids) {
  i <- annotation_match(ids, names(lengths))
  out <- stats::setNames(unname(lengths[i]), ids)
  out[!is.na(out)]
}

# ---- shared UI ------------------------------------------------------------------

#' "Download a GTF" controls; ids are "<prefix>_fetch_*".
annotation_fetch_ui <- function(p, summary = "No GTF? Download one from Ensembl or GENCODE (any organism)") {
  id <- function(x) paste0(p, "_fetch_", x)
  tags$details(
    tags$summary(tags$strong(summary)),
    div(class = "ctl-block",
      uiOutput(id("guess")),
      fluidRow(
        column(4, selectInput(id("species"), "Organism", organism_choices())),
        column(4, selectInput(id("provider"), "Source", ANNOTATION_PROVIDERS)),
        column(4, textInput(id("release"), "Release", "current"))),
      conditionalPanel(sprintf("input['%s'] == 'other'", id("species")),
        fluidRow(
          column(6, textInput(id("other_dir"), "Ensembl name of the organism", placeholder = "e.g. ovis_aries")),
          column(6, selectInput(id("other_div"), "Ensembl division",
                                c("Vertebrates (ensembl.org)" = "", "Plants" = "plants", "Metazoa" = "metazoa",
                                  "Fungi" = "fungi", "Protists" = "protists")))),
        helpText("The name is the organism's folder on Ensembl's FTP site: lower case with underscores.")),
      helpText(HTML(paste0(
        "Best is the release your counts were made with (a pipeline's run report usually records it in its ",
        "parameters, e.g. <code>gencode.v48</code> or <code>release-115</code>); <b>current</b> takes ",
        "the newest, which is fine for symbols and close enough for lengths. GENCODE covers human and mouse only. Human files are ",
        "about 125 MB (GENCODE) to 140 MB (Ensembl): several minutes on a slow connection, so the ",
        "progress bar shows the megabytes. Each file is downloaded once and kept for next time. ",
        "It provides gene symbols and biotypes, gene lengths (union of each gene's exons) for TPM / FPKM, and the transcript-to-gene map for transcript-level analysis."))),
      actionButton(id("go"), "Download and use", class = "btn-block")))
}

#' Server half: suggests species / source from `ids()`, downloads, parses, and
#' hands the result to `on_loaded(annotation, lengths, label, transcripts, organism)`.
annotation_fetch_server <- function(input, output, session, p, ids, on_loaded) {
  id <- function(x) paste0(p, "_fetch_", x)
  guess <- reactive({ i <- ids(); if (is.null(i) || !length(i)) NULL else guess_annotation_source(i) })

  observe({
    g <- guess(); req(g)
    if (!is.na(g$species)) updateSelectInput(session, id("species"), selected = g$species)
    updateSelectInput(session, id("provider"), selected = g$provider)
  })

  output[[id("guess")]] <- renderUI({
    g <- guess(); if (is.null(g)) return(NULL)
    div(style = "font-size:.88em;color:#42525f;margin-bottom:4px", g$note, " Choices set to match.")
  })

  observeEvent(input[[id("go")]], {
    res <- withProgress(message = "Getting the annotation", value = 0.05, {
      tryCatch({
        setProgress(0.1, detail = "Finding the file")
        src <- annotation_url(input[[id("provider")]], input[[id("species")]], input[[id("release")]],
                              other = list(ensembl_dir = input[[id("other_dir")]], division = input[[id("other_div")]]))
        setProgress(0.2, detail = paste("Downloading", src$file, "(first time only)"))
        path <- download_annotation(src, progress = function(done, total) {
          setProgress(0.2 + if (is.na(total)) 0 else 0.4 * done / total,
                      detail = sprintf("Downloading %s: %.0f%s MB", src$file, done,
                                       if (is.na(total)) "" else sprintf(" of %.0f", total)))
        })
        b <- load_annotation_bundle(path, progress = function(f, d) setProgress(f, detail = d))
        list(ann = b$ann, len = b$len, tx = b$tx, label = src$label)
      }, error = function(e) e)
    })
    if (inherits(res, "error")) {
      showNotification(paste("Could not get the annotation:", conditionMessage(res)), type = "error", duration = 20)
      return()
    }
    on_loaded(res$ann, res$len, res$label, res$tx, input[[id("species")]])
  })
}
