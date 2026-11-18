# =============================================================================
# organisms.R -- one registry of the organisms the app knows by name, and the
# helpers that let annotation download, enrichment and gene look-ups follow the
# organism of the loaded study.
#
# Differential expression itself never depends on the organism. What does:
#   - downloading a GTF        (Ensembl / Ensembl Genomes directory per species)
#   - GO / KEGG enrichment     (the organism's Bioconductor OrgDb + KEGG code)
#   - MSigDB collections       (msigdbr species; human genes mapped by orthology
#                               for other organisms)
#   - gene look-up links       (species name, KEGG code, human-only resources)
# Any other organism still works: pick "Other" and give its Ensembl name, or
# upload your own GTF; enrichment then needs an OrgDb package name.
# =============================================================================

ORGANISMS <- list(
  human = list(label = "Human (Homo sapiens)", species = "Homo sapiens",
               orgdb = "org.Hs.eg.db", kegg = "hsa", msigdbr = "Homo sapiens",
               ensembl_dir = "homo_sapiens", division = "", ensembl_file = "Homo_sapiens.GRCh38",
               gencode_dir = "Gencode_human", gencode_prefix = "",
               gene_prefix = "ENSG", tx_prefix = "ENST", mygene = "human"),
  mouse = list(label = "Mouse (Mus musculus)", species = "Mus musculus",
               orgdb = "org.Mm.eg.db", kegg = "mmu", msigdbr = "Mus musculus",
               ensembl_dir = "mus_musculus", division = "", ensembl_file = "Mus_musculus.GRCm39",
               gencode_dir = "Gencode_mouse", gencode_prefix = "M",
               gene_prefix = "ENSMUSG", tx_prefix = "ENSMUST", mygene = "mouse"),
  rat = list(label = "Rat (Rattus norvegicus)", species = "Rattus norvegicus",
             orgdb = "org.Rn.eg.db", kegg = "rno", msigdbr = "Rattus norvegicus",
             ensembl_dir = "rattus_norvegicus", division = "",
             gene_prefix = "ENSRNOG", tx_prefix = "ENSRNOT", mygene = "rat"),
  zebrafish = list(label = "Zebrafish (Danio rerio)", species = "Danio rerio",
                   orgdb = "org.Dr.eg.db", kegg = "dre", msigdbr = "Danio rerio",
                   ensembl_dir = "danio_rerio", division = "",
                   gene_prefix = "ENSDARG", tx_prefix = "ENSDART", mygene = "zebrafish"),
  chicken = list(label = "Chicken (Gallus gallus)", species = "Gallus gallus",
                 orgdb = "org.Gg.eg.db", kegg = "gga", msigdbr = "Gallus gallus",
                 ensembl_dir = "gallus_gallus", division = "",
                 gene_prefix = "ENSGALG", tx_prefix = "ENSGALT", mygene = "9031"),
  pig = list(label = "Pig (Sus scrofa)", species = "Sus scrofa",
             orgdb = "org.Ss.eg.db", kegg = "ssc", msigdbr = "Sus scrofa",
             ensembl_dir = "sus_scrofa", division = "",
             gene_prefix = "ENSSSCG", tx_prefix = "ENSSSCT", mygene = "pig"),
  cow = list(label = "Cow (Bos taurus)", species = "Bos taurus",
             orgdb = "org.Bt.eg.db", kegg = "bta", msigdbr = "Bos taurus",
             ensembl_dir = "bos_taurus", division = "",
             gene_prefix = "ENSBTAG", tx_prefix = "ENSBTAT", mygene = "9913"),
  dog = list(label = "Dog (Canis lupus familiaris)", species = "Canis lupus familiaris",
             orgdb = "org.Cf.eg.db", kegg = "cfa", msigdbr = "Canis lupus familiaris",
             ensembl_dir = "canis_lupus_familiaris", division = "",
             gene_prefix = "ENSCAFG", tx_prefix = "ENSCAFT", mygene = "9615"),
  rhesus = list(label = "Rhesus macaque (Macaca mulatta)", species = "Macaca mulatta",
                orgdb = "org.Mmu.eg.db", kegg = "mcc", msigdbr = "Macaca mulatta",
                ensembl_dir = "macaca_mulatta", division = "",
                gene_prefix = "ENSMMUG", tx_prefix = "ENSMMUT", mygene = "9544"),
  fly = list(label = "Fruit fly (Drosophila melanogaster)", species = "Drosophila melanogaster",
             orgdb = "org.Dm.eg.db", kegg = "dme", msigdbr = "Drosophila melanogaster",
             ensembl_dir = "drosophila_melanogaster", division = "",
             gene_prefix = "FBgn", tx_prefix = "FBtr", mygene = "fruitfly"),
  worm = list(label = "Nematode (Caenorhabditis elegans)", species = "Caenorhabditis elegans",
              orgdb = "org.Ce.eg.db", kegg = "cel", msigdbr = "Caenorhabditis elegans",
              ensembl_dir = "caenorhabditis_elegans", division = "",
              gene_prefix = "WBGene", tx_prefix = "", mygene = "nematode"),
  yeast = list(label = "Yeast (Saccharomyces cerevisiae)", species = "Saccharomyces cerevisiae",
               orgdb = "org.Sc.sgd.db", kegg = "sce", msigdbr = "Saccharomyces cerevisiae",
               ensembl_dir = "saccharomyces_cerevisiae", division = "",
               gene_prefix = "", tx_prefix = "", mygene = "yeast"),
  arabidopsis = list(label = "Arabidopsis (Arabidopsis thaliana)", species = "Arabidopsis thaliana",
                     orgdb = "org.At.tair.db", kegg = "ath", msigdbr = "",
                     ensembl_dir = "arabidopsis_thaliana", division = "plants",
                     gene_prefix = "AT", tx_prefix = "", mygene = "thale-cress")
)

#' The organism entry for a key; "other" builds one from the user's settings.
#'
#' @param other list(ensembl_dir, division, orgdb, kegg, species) for organisms
#'   not in ORGANISMS (any Ensembl or Ensembl Genomes species).
organism_info <- function(key, other = NULL) {
  if (!is.null(key) && key %in% names(ORGANISMS)) return(c(ORGANISMS[[key]], key = key))
  o <- other %||% list()
  dir <- tolower(gsub("[^A-Za-z_]+", "_", trimws(o$ensembl_dir %||% "")))
  sp <- trimws(o$species %||% "")
  if (!nzchar(sp) && nzchar(dir)) {
    sp <- sub("_", " ", dir); substr(sp, 1, 1) <- toupper(substr(sp, 1, 1))
  }
  list(key = "other", label = if (nzchar(sp)) sp else "Other organism", species = sp,
       orgdb = trimws(o$orgdb %||% ""), kegg = trimws(o$kegg %||% ""), msigdbr = sp,
       ensembl_dir = dir, division = trimws(o$division %||% ""),
       gene_prefix = "", tx_prefix = "", mygene = sp)
}

#' Choices for an organism selector: the known organisms, then "Other".
organism_choices <- function() {
  c(stats::setNames(names(ORGANISMS), vapply(ORGANISMS, `[[`, "", "label")),
    "Other (any Ensembl species)" = "other")
}

#' Guess the organism from gene or transcript ids (Ensembl-style prefixes).
#' @return an ORGANISMS key, or NA when the ids do not say.
guess_organism <- function(ids) {
  ids <- utils::head(as.character(ids[!is.na(ids) & nzchar(ids)]), 5000)
  if (!length(ids)) return(NA_character_)
  ids <- sub("\\.[0-9]+$", "", ids)
  best <- NA_character_; best_frac <- 0.5
  # longest prefixes first, so ENSMUSG is not mistaken for something shorter
  pre <- unlist(lapply(names(ORGANISMS), function(k) {
    p <- c(ORGANISMS[[k]]$gene_prefix, ORGANISMS[[k]]$tx_prefix)
    stats::setNames(p[nzchar(p)], rep(k, sum(nzchar(p))))
  }))
  pre <- pre[order(-nchar(pre))]
  taken <- rep(FALSE, length(ids))
  hits <- list()
  for (i in seq_along(pre)) {
    p <- pre[[i]]
    rx <- if (p == "AT") "^AT[1-5CM]G[0-9]{5}" else paste0("^", p, "[0-9]")
    m <- !taken & grepl(rx, ids)
    taken <- taken | m
    hits[[names(pre)[i]]] <- (hits[[names(pre)[i]]] %||% 0) + sum(m)
  }
  if (!length(hits)) return(NA_character_)
  fr <- unlist(hits) / length(ids)
  if (max(fr) > best_frac) best <- names(fr)[which.max(fr)]
  # yeast systematic names (YAL001C) and plain ORFs
  if (is.na(best) && mean(grepl("^Y[A-P][LR][0-9]{3}[WC]", ids)) > 0.5) best <- "yeast"
  best
}

#' Whether ids look like transcripts rather than genes (Ensembl-style).
looks_like_transcripts <- function(ids) {
  ids <- utils::head(as.character(ids[!is.na(ids)]), 5000)
  if (!length(ids)) return(FALSE)
  base <- sub("\\.[0-9]+$", "", ids)
  tx <- unlist(lapply(ORGANISMS, `[[`, "tx_prefix")); tx <- tx[nzchar(tx)]
  frac <- mean(vapply(base, function(x) any(startsWith(x, tx)), logical(1)))
  frac > 0.5 || mean(grepl("^AT[1-5CM]G[0-9]{5}\\.[0-9]+$", ids)) > 0.5
}

#' The Bioconductor install line for an organism's OrgDb, or NULL when installed.
organism_orgdb_missing <- function(org) {
  db <- org$orgdb %||% ""
  if (!nzchar(db)) return("No OrgDb annotation package is set for this organism.")
  if (requireNamespace(db, quietly = TRUE)) return(NULL)
  sprintf('BiocManager::install("%s")', db)
}

#' Pick the OrgDb key type that matches most of the ids (ENSEMBL, TAIR, SYMBOL ...).
#' @return list(keytype, ids) where ids are the inputs re-keyed if symbols were used
orgdb_keytype <- function(orgdb, ids, symbols = NULL) {
  db <- getExportedValue(orgdb, orgdb)
  kts <- AnnotationDbi::keytypes(db)
  base <- sub("\\.[0-9]+$", "", ids)
  try_kt <- function(kt, x) {
    if (!kt %in% kts) return(0)
    k <- tryCatch(AnnotationDbi::keys(db, kt), error = function(e) character(0))
    mean(x %in% k)
  }
  cand <- c("ENSEMBL", "TAIR", "FLYBASE", "WORMBASE", "ORF", "SGD", "ENTREZID", "SYMBOL")
  sc <- vapply(cand, try_kt, numeric(1), x = base)
  best <- names(sc)[which.max(sc)]
  if (max(sc) >= 0.3) return(list(keytype = best, map = stats::setNames(base, ids)))
  if (!is.null(symbols)) {
    s <- symbols[ids]
    sc_s <- try_kt("SYMBOL", s[!is.na(s)])
    if (sc_s >= 0.3) return(list(keytype = "SYMBOL", map = stats::setNames(unname(s), ids)))
  }
  stop(sprintf("Few of the gene ids match %s (best: %s, %.0f%%). Check the organism setting.",
               orgdb, best, 100 * max(sc)))
}
