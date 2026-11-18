# =============================================================================
# make_example_enrichment.R -- (re)build data/example_enrichment/ from the
# current official Gene Ontology release (and, for local use only, current KEGG).
#
#   Rscript data-raw/make_example_enrichment.R                 # latest release
#   Rscript data-raw/make_example_enrichment.R --release 2026-06-19
#   Rscript data-raw/make_example_enrichment.R --evidence experimental
#   Rscript data-raw/make_example_enrichment.R --go-obo go-basic.obo \
#           --gaf HUMAN-uniprot.gaf.gz                         # files on disk
#   Rscript data-raw/make_example_enrichment.R --kegg          # also pathway.all.csv (local only)
#
# Run from the repository root. Base R only -- no Bioconductor needed.
#
# Sources of truth
#   GO     go-basic.obo + the human annotation file (GAF), both taken from the
#          SAME dated Gene Ontology Consortium release. The release date and
#          the ontology data-version are read from the files themselves and
#          written to data/example_enrichment/SOURCE.txt.
#   KEGG   only with --kegg: the KEGG REST API (pathway names, gene-pathway links,
#          gene symbols), with the KEGG release string recorded alongside. KEGG's
#          REST API is provided for academic use and KEGG data may not be
#          redistributed (https://www.kegg.jp/kegg/legal.html), so pathway.all.csv
#          is NOT shipped with DEGexplorer (it is listed in .gitignore). Build it
#          locally if you want the KEGG demo table and the KEGG README figure.
#
# Method
#   gene set   DEGs of data/example_deseq2_results.tsv.gz with adj. P <= 0.05
#              and |log2FC| >= 1
#   universe   every gene of that table with an adjusted p-value
#   GO sets    direct annotations propagated up the ontology (is_a + part_of by
#              default), so a gene annotated to a child term also counts for
#              every ancestor. NOT-qualified annotations are dropped, alternate
#              ids are mapped to their primary term, obsolete terms are removed.
#   evidence   --evidence all (default: every code, IEA and IBA included)
#                         no-iea        (drop electronic IEA)
#                         experimental  (EXP IDA IPI IMP IGI IEP + HTP codes)
#   test       one-sided hypergeometric, gene sets of 10-500 genes, BH-adjusted
#              -- the test clusterProfiler runs
#
# Options
#   --release YYYY-MM-DD   use that archived GO release instead of the current one
#   --evidence MODE        all | no-iea | experimental
#   --relations LIST       comma separated; default is_a,part_of. Add
#                          regulates,positively_regulates,negatively_regulates
#                          to reproduce org.Hs.eg.db / clusterProfiler's GOALL
#   --go-obo FILE --gaf FILE     use files already on disk (must be one release)
#   --kegg-dir DIR         use saved KEGG REST replies: pathway.txt, link.txt,
#                          genes.txt, info.txt
#   --kegg                 also build pathway.all.csv from KEGG (local use only)
#   --out DIR              output directory (default data/example_enrichment)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
has <- function(flag) flag %in% args

release   <- opt("--release")
evidence  <- opt("--evidence", "all")
relations <- strsplit(opt("--relations", "is_a,part_of"), ",")[[1]]
out_dir   <- opt("--out", "data/example_enrichment")
stopifnot(evidence %in% c("all", "no-iea", "experimental"))

source("R/utils_io.R")
options(timeout = 600)
cache <- file.path(tempdir(), "go_release"); dir.create(cache, showWarnings = FALSE)

fetch <- function(urls, dest) {
  for (u in urls) {
    ok <- tryCatch({ utils::download.file(u, dest, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && file.exists(dest) && file.size(dest) > 1000) { message("  downloaded ", u); return(dest) }
  }
  stop("Could not download any of:\n  ", paste(urls, collapse = "\n  "))
}

# ---- 1. the demo study ------------------------------------------------------

de  <- standardise_de_table(read_table_auto("data/example_deseq2_results.tsv.gz"))$data
n_deg_all <- sum(grepl("^ENSG", de$gene_id) & !is.na(de$padj) & de$padj <= 0.05 &
                   abs(de$log2FoldChange) >= 1)
# GO and KEGG files identify genes by symbol, so genes without one drop out.
de  <- de[grepl("^ENSG", de$gene_id) & !is.na(de$display_gene) & !grepl("^ENSG", de$display_gene), ]
de  <- de[!duplicated(de$display_gene), ]
uni <- de$display_gene[!is.na(de$padj)]
sig <- de$display_gene[!is.na(de$padj) & de$padj <= 0.05 & abs(de$log2FoldChange) >= 1]
message(sprintf("Demo study: %d of %d DEGs carry a gene symbol; %d background genes.",
                length(sig), n_deg_all, length(uni)))

ora <- function(t2g, term_names, min_gs = 10, max_gs = 500) {
  t2g <- unique(t2g[t2g$gene %in% uni, ])
  N <- length(unique(t2g$gene)); g <- intersect(sig, t2g$gene); n <- length(g)
  sets <- split(t2g$gene, t2g$term)
  sets <- sets[lengths(sets) >= min_gs & lengths(sets) <= max_gs]
  k <- vapply(sets, function(s) sum(s %in% g), 1L)
  sets <- sets[k > 0]; k <- k[k > 0]; M <- lengths(sets)
  p <- stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE)
  padj <- stats::p.adjust(p, "BH")
  out <- data.frame(
    ID = names(sets), Description = unname(term_names[names(sets)]),
    GeneRatio = paste0(k, "/", n), BgRatio = paste0(M, "/", N),
    pvalue = signif(p, 6), p.adjust = signif(padj, 6), qvalue = signif(padj, 6),
    geneID = vapply(sets, function(s) paste(intersect(sig, s), collapse = "/"), ""),
    Count = k, stringsAsFactors = FALSE)
  out[order(out$pvalue, out$ID), ]
}

# ---- 2. Gene Ontology: ontology + annotations from one release --------------

base <- if (is.null(release)) "https://current.geneontology.org" else
  paste0("https://release.geneontology.org/", release)
obo_file <- opt("--go-obo"); gaf_file <- opt("--gaf")
if (is.null(obo_file)) obo_file <- fetch(paste0(base, "/ontology/go-basic.obo"),
                                         file.path(cache, "go-basic.obo"))
# The human file was renamed in the 2026-06-19 release; older archived
# releases still carry the previous name.
if (is.null(gaf_file)) gaf_file <- fetch(
  c(paste0(base, "/annotations/gaf/HUMAN-uniprot.gaf.gz"),
    paste0(base, "/annotations/goa_human.gaf.gz")),
  file.path(cache, "human.gaf.gz"))

parse_obo <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  version <- sub("^data-version:\\s*", "", grep("^data-version:", x[1:50], value = TRUE)[1])
  # A line belongs to a [Term] stanza when the nearest stanza header above it
  # is "[Term]" (and not [Typedef] or [Instance]).
  hdr <- which(startsWith(x, "["))
  start <- hdr[x[hdr] == "[Term]"]
  near <- findInterval(seq_along(x), hdr)
  in_term <- near > 0 & x[hdr[pmax(near, 1)]] == "[Term]"
  stanza <- match(hdr[pmax(near, 1)], start)            # index of the term a line is in
  field <- function(tag) {
    i <- which(in_term & startsWith(x, paste0(tag, ": ")))
    data.frame(stanza = stanza[i], value = substring(x[i], nchar(tag) + 3), stringsAsFactors = FALSE)
  }
  id <- field("id"); id <- id[!duplicated(id$stanza), ]
  ids <- character(length(start)); ids[id$stanza] <- id$value
  one <- function(tag) { f <- field(tag); v <- rep(NA_character_, length(start))
                         v[f$stanza[!duplicated(f$stanza)]] <- f$value[!duplicated(f$stanza)]; v }
  obsolete <- !is.na(one("is_obsolete"))
  terms <- data.frame(id = ids, name = one("name"), namespace = one("namespace"),
                      obsolete = obsolete, stringsAsFactors = FALSE)
  isa <- field("is_a")
  edges <- data.frame(child = ids[isa$stanza], parent = sub("\\s.*$", "", isa$value),
                      rel = "is_a", stringsAsFactors = FALSE)
  r <- field("relationship")
  if (nrow(r)) {
    parts <- strsplit(r$value, "\\s+")
    edges <- rbind(edges, data.frame(child = ids[r$stanza],
                                     parent = vapply(parts, `[`, "", 2),
                                     rel = vapply(parts, `[`, "", 1), stringsAsFactors = FALSE))
  }
  alt <- field("alt_id")
  list(version = version, terms = terms[grepl("^GO:", terms$id), ],
       edges = edges[grepl("^GO:", edges$child) & grepl("^GO:", edges$parent), ],
       alt = stats::setNames(ids[alt$stanza], alt$value))
}

message("Reading the ontology ...")
go <- parse_obo(obo_file)
live <- go$terms[!go$terms$obsolete, ]
edges <- go$edges[go$edges$rel %in% relations & go$edges$child %in% live$id &
                    go$edges$parent %in% live$id, ]
parents <- split(edges$parent, edges$child)

# every term's ancestors (itself included), memoised
anc_env <- new.env(hash = TRUE, parent = emptyenv())
ancestors <- function(id) {
  hit <- anc_env[[id]]
  if (!is.null(hit)) return(hit)
  a <- id
  for (p in parents[[id]]) a <- c(a, ancestors(p))
  a <- unique(a); assign(id, a, envir = anc_env); a
}

message("Reading the annotations ...")
hdr <- readLines(gaf_file, n = 200, warn = FALSE)
hdr <- hdr[startsWith(hdr, "!")]
gaf_date <- sub("^![^:]*:\\s*", "", grep("^!\\s*(date-generated|generated-date|date)\\s*:", hdr,
                                         ignore.case = TRUE, value = TRUE)[1])
gaf_go_version <- sub("^![^:]*:\\s*", "", grep("^!\\s*go-version", hdr, ignore.case = TRUE, value = TRUE)[1])
gaf <- utils::read.delim(gaf_file, header = FALSE, comment.char = "!", quote = "", sep = "\t",
                         stringsAsFactors = FALSE, fill = TRUE,
                         colClasses = "character")[, c(3, 4, 5, 7, 9, 13)]
colnames(gaf) <- c("symbol", "qualifier", "go", "evidence", "aspect", "taxon")
n_raw <- nrow(gaf)
gaf <- gaf[grepl("taxon:9606", gaf$taxon, fixed = TRUE) & !grepl("NOT", gaf$qualifier, fixed = TRUE), ]
EXPERIMENTAL <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP")
if (evidence == "no-iea")       gaf <- gaf[gaf$evidence != "IEA", ]
if (evidence == "experimental") gaf <- gaf[gaf$evidence %in% EXPERIMENTAL, ]
is_alt <- gaf$go %in% names(go$alt)
gaf$go[is_alt] <- go$alt[gaf$go[is_alt]]
gaf <- unique(gaf[gaf$go %in% live$id & gaf$symbol %in% uni, c("symbol", "go")])
message(sprintf("  %s annotation lines; %s direct gene-term pairs kept for the demo background.",
                format(n_raw, big.mark = ","), format(nrow(gaf), big.mark = ",")))

message("Propagating annotations to ancestor terms (", paste(relations, collapse = " + "), ") ...")
anc <- lapply(unique(gaf$go), ancestors); names(anc) <- unique(gaf$go)
a <- anc[gaf$go]
t2g <- unique(data.frame(term = unlist(a, use.names = FALSE),
                         gene = rep(gaf$symbol, lengths(a)), stringsAsFactors = FALSE))
term_name <- stats::setNames(live$name, live$id)
term_ns   <- stats::setNames(live$namespace, live$id)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ns <- c(BP = "biological_process", CC = "cellular_component", MF = "molecular_function")
report <- character(0)
for (o in names(ns)) {
  res <- ora(t2g[term_ns[t2g$term] == ns[[o]], ], term_name)
  f <- sprintf("go.%s.all.csv", o)
  utils::write.csv(res, file.path(out_dir, f), row.names = FALSE)
  report <- c(report, sprintf("%-16s %5d terms, %4d with p < 0.05, %4d with adj. P < 0.05", f,
                              nrow(res), sum(res$pvalue < 0.05), sum(res$p.adjust < 0.05)))
}

# ---- 3. KEGG, current -------------------------------------------------------

kegg_info <- paste("not included: KEGG data are not redistributed with DEGexplorer (KEGG licence);",
                   "run with --kegg to build pathway.all.csv locally")
if (has("--kegg")) {
  kd <- opt("--kegg-dir")
  kget <- function(name, url) if (!is.null(kd)) file.path(kd, name) else
    fetch(url, file.path(cache, name))
  rd <- function(f) utils::read.delim(f, header = FALSE, quote = "", sep = "\t",
                                      stringsAsFactors = FALSE, colClasses = "character")
  pw   <- rd(kget("pathway.txt", "https://rest.kegg.jp/list/pathway/hsa"))
  link <- rd(kget("link.txt",    "https://rest.kegg.jp/link/pathway/hsa"))
  gene <- rd(kget("genes.txt",   "https://rest.kegg.jp/list/hsa"))
  # The release reply is only a few hundred bytes, below fetch()'s 1000-byte
  # sanity check, so read it directly and never let it stop the run.
  info <- if (!is.null(kd) && file.exists(file.path(kd, "info.txt")))
    readLines(file.path(kd, "info.txt"), warn = FALSE) else
    tryCatch(readLines("https://rest.kegg.jp/info/pathway", warn = FALSE),
             error = function(e) character(0))
  kegg_info <- trimws(sub("^pathway\\s*", "", grep("Release", info, value = TRUE)[1]))
  if (is.na(kegg_info)) kegg_info <- paste("KEGG REST API, release not reported; retrieved", format(Sys.Date()))
  pw_name <- stats::setNames(sub(" - Homo sapiens \\(human\\)$", "", pw[[2]]), sub("^path:", "", pw[[1]]))
  # "hsa:3434 <tab> ... <tab> IFIT1, G10P1, ...; description": first symbol is the official one
  desc <- gene[[ncol(gene)]]
  g_sym <- stats::setNames(trimws(sub("[,;].*$", "", desc)), gene[[1]])
  k2g <- data.frame(term = sub("^path:", "", link[[2]]), gene = unname(g_sym[link[[1]]]),
                    stringsAsFactors = FALSE)
  k2g <- k2g[!is.na(k2g$gene) & grepl("^hsa[0-9]{5}$", k2g$term), ]
  res <- ora(k2g, pw_name)
  utils::write.csv(res, file.path(out_dir, "pathway.all.csv"), row.names = FALSE)
  report <- c(report, sprintf("%-16s %5d terms, %4d with p < 0.05, %4d with adj. P < 0.05",
                              "pathway.all.csv", nrow(res), sum(res$pvalue < 0.05),
                              sum(res$p.adjust < 0.05)))
}

# ---- 4. provenance ----------------------------------------------------------

writeLines(c(
  "Provenance of the example enrichment tables in this folder",
  "===========================================================",
  sprintf("Built on:               %s", format(Sys.Date())),
  sprintf("Built by:               data-raw/make_example_enrichment.R %s", paste(args, collapse = " ")),
  "",
  "Gene Ontology",
  sprintf("  release requested:    %s", if (is.null(release)) "current" else release),
  sprintf("  ontology:             go-basic.obo, data-version %s", go$version),
  sprintf("  annotations:          %s", basename(gaf_file)),
  sprintf("  annotation file date: %s", gaf_date),
  sprintf("  annotation go-version:%s", if (is.na(gaf_go_version)) " (not stated in the file)" else paste0(" ", gaf_go_version)),
  sprintf("  evidence codes:       %s", evidence),
  sprintf("  propagated over:      %s", paste(relations, collapse = ", ")),
  "",
  "KEGG",
  sprintf("  %s", kegg_info),
  "",
  "Demo study",
  sprintf("  data/example_deseq2_results.tsv.gz: %d DEGs (adj. P <= 0.05, |log2FC| >= 1),", n_deg_all),
  sprintf("  of which %d carry a gene symbol; %d background genes. Genes are matched", length(sig), length(uni)),
  "  to annotations by symbol.",
  "", "Tables", paste0("  ", report)),
  file.path(out_dir, "SOURCE.txt"))
cat(readLines(file.path(out_dir, "SOURCE.txt")), sep = "\n")
