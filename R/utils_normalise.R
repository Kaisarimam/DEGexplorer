# =============================================================================
# utils_normalise.R -- expression normalisation from raw counts
#
# All methods are implemented in base R so the app does not require DESeq2 or
# edgeR to be installed just to obtain their normalisations. Where a reference
# implementation exists the algorithm follows it exactly:
#
#   CPM      counts per million (library-size scaling only)
#   TPM      transcript per million   (length-corrected, then scaled)
#   FPKM/RPKM fragments per kb per million
#   DESeq2   median-of-ratios size factors (Anders & Huber 2010)
#   edgeR    TMM effective library sizes (Robinson & Oshlack 2010)
#   VST-like variance-stabilising surrogate: log2(normalised + 1)
#
# Every function takes a gene x sample numeric matrix and returns one.
# =============================================================================

NORM_METHODS <- c(
  "Raw counts"                    = "raw",
  "CPM (counts per million)"      = "cpm",
  "TPM (transcripts per million)" = "tpm",
  "FPKM / RPKM"                   = "fpkm",
  "DESeq2 median-of-ratios"       = "deseq2",
  "edgeR TMM"                     = "tmm",
  "log2 CPM"                      = "logcpm",
  "log2 TPM"                      = "logtpm",
  "log2 FPKM"                     = "logfpkm",
  "log2 DESeq2-normalised"        = "logdeseq2",
  "log2 TMM-normalised"           = "logtmm"
)

# Methods that need a gene-length matrix / vector
LENGTH_DEPENDENT <- c("tpm", "fpkm", "logtpm", "logfpkm")

#' Are the inputs sufficient for a given method?
norm_requirements_met <- function(method, has_lengths) {
  !(method %in% LENGTH_DEPENDENT) || isTRUE(has_lengths)
}

.clean_counts <- function(counts) {
  counts[!is.finite(counts)] <- 0
  counts[counts < 0] <- 0
  counts
}

#' Align a gene-length matrix (or vector) to the counts matrix.
#' RSEM emits per-sample effective lengths; a plain 1-column table or a named
#' vector is also accepted and recycled across samples.
align_lengths <- function(lengths, counts) {
  if (is.null(lengths)) return(NULL)
  if (is.vector(lengths) && !is.matrix(lengths)) {
    lengths <- matrix(lengths, ncol = 1,
                      dimnames = list(names(lengths), "length"))
  }
  common <- intersect(rownames(counts), rownames(lengths))
  if (length(common) == 0) {
    stop("The gene-length matrix shares no gene ids with the count matrix.")
  }
  L <- lengths[common, , drop = FALSE]
  if (ncol(L) == 1) {
    L <- matrix(rep(L[, 1], ncol(counts)), nrow = length(common),
                dimnames = list(common, colnames(counts)))
  } else {
    shared_cols <- intersect(colnames(counts), colnames(L))
    if (length(shared_cols) == ncol(counts)) {
      L <- L[, colnames(counts), drop = FALSE]
    } else {
      # column names don't line up: fall back to the per-gene mean length
      m <- rowMeans(L, na.rm = TRUE)
      L <- matrix(rep(m, ncol(counts)), nrow = length(common),
                  dimnames = list(common, colnames(counts)))
    }
  }
  L[!is.finite(L) | L <= 0] <- NA_real_
  L
}

# ---- individual methods -----------------------------------------------------

calc_cpm <- function(counts) {
  counts <- .clean_counts(counts)
  lib <- colSums(counts)
  lib[lib == 0] <- 1
  sweep(counts, 2, lib, "/") * 1e6
}

calc_tpm <- function(counts, lengths) {
  counts <- .clean_counts(counts)
  L <- align_lengths(lengths, counts)
  common <- rownames(L)
  cnt <- counts[common, , drop = FALSE]
  rate <- cnt / L                      # reads per base
  rate[!is.finite(rate)] <- 0
  denom <- colSums(rate, na.rm = TRUE)
  denom[denom == 0] <- 1
  sweep(rate, 2, denom, "/") * 1e6
}

calc_fpkm <- function(counts, lengths) {
  counts <- .clean_counts(counts)
  L <- align_lengths(lengths, counts)
  common <- rownames(L)
  cnt <- counts[common, , drop = FALSE]
  lib <- colSums(counts)               # library size from the full matrix
  lib[lib == 0] <- 1
  out <- sweep(cnt, 2, lib / 1e6, "/") / (L / 1000)
  out[!is.finite(out)] <- 0
  out
}

#' DESeq2 median-of-ratios size factors.
deseq2_size_factors <- function(counts) {
  counts <- .clean_counts(counts)
  logc <- log(counts)
  logc[!is.finite(logc)] <- NA_real_
  # geometric mean across samples, using only genes positive in every sample
  loggeo <- rowMeans(logc, na.rm = FALSE)
  keep <- is.finite(loggeo)
  if (!any(keep)) {
    # fall back to a poised-median approach when no gene is universally positive
    loggeo <- rowMeans(logc, na.rm = TRUE)
    keep <- is.finite(loggeo)
  }
  if (!any(keep)) return(setNames(rep(1, ncol(counts)), colnames(counts)))
  sf <- apply(logc[keep, , drop = FALSE], 2, function(x) {
    exp(median((x - loggeo[keep])[is.finite(x)], na.rm = TRUE))
  })
  sf[!is.finite(sf) | sf <= 0] <- 1
  sf
}

calc_deseq2 <- function(counts) {
  sf <- deseq2_size_factors(counts)
  sweep(.clean_counts(counts), 2, sf, "/")
}

#' edgeR TMM normalisation factors (Robinson & Oshlack 2010).
#' Returns factors that multiply to 1, as edgeR's calcNormFactors does.
tmm_norm_factors <- function(counts, ref_column = NULL,
                             logratio_trim = 0.3, sum_trim = 0.05,
                             a_cutoff = -1e10) {
  counts <- .clean_counts(counts)
  lib <- colSums(counts)
  lib[lib == 0] <- 1
  if (is.null(ref_column)) {
    # edgeR picks the sample whose upper-quartile is closest to the mean UQ
    f75 <- apply(sweep(counts, 2, lib, "/"), 2, function(x) quantile(x, 0.75, na.rm = TRUE))
    ref_column <- if (all(!is.finite(f75)) || all(f75 == 0)) 1L else
      which.min(abs(f75 - mean(f75, na.rm = TRUE)))
  }
  ref <- counts[, ref_column]
  ref_lib <- lib[ref_column]

  f <- vapply(seq_len(ncol(counts)), function(j) {
    obs <- counts[, j]; n_obs <- lib[j]
    keep <- obs > 0 & ref > 0
    if (sum(keep) < 2) return(1)
    o <- obs[keep]; r <- ref[keep]
    logR <- log2((o / n_obs) / (r / ref_lib))                    # M
    absE <- (log2(o / n_obs) + log2(r / ref_lib)) / 2            # A
    v <- (n_obs - o) / n_obs / o + (ref_lib - r) / ref_lib / r   # approx var

    fin <- is.finite(logR) & is.finite(absE) & absE > a_cutoff
    logR <- logR[fin]; absE <- absE[fin]; v <- v[fin]
    if (length(logR) < 2) return(1)
    if (max(abs(logR)) < 1e-6) return(1)

    n <- length(logR)
    loL <- floor(n * logratio_trim) + 1; hiL <- n + 1 - loL
    loS <- floor(n * sum_trim) + 1;      hiS <- n + 1 - loS
    rkR <- rank(logR); rkA <- rank(absE)
    keep2 <- (rkR >= loL & rkR <= hiL) & (rkA >= loS & rkA <= hiS)
    if (!any(keep2)) return(1)
    w <- 1 / v[keep2]
    w[!is.finite(w) | w <= 0] <- NA_real_
    est <- sum(w * logR[keep2], na.rm = TRUE) / sum(w, na.rm = TRUE)
    if (!is.finite(est)) est <- mean(logR[keep2], na.rm = TRUE)
    2^est
  }, numeric(1))

  f[!is.finite(f) | f <= 0] <- 1
  f <- f / exp(mean(log(f)))   # centre so the factors multiply to 1
  setNames(f, colnames(counts))
}

calc_tmm <- function(counts) {
  counts <- .clean_counts(counts)
  nf  <- tmm_norm_factors(counts)
  lib <- colSums(counts); lib[lib == 0] <- 1
  eff <- lib * nf
  sweep(counts, 2, eff, "/") * mean(eff)
}

# ---- dispatcher -------------------------------------------------------------

#' Produce a normalised expression matrix.
#'
#' @param counts  raw (or RSEM expected) counts, genes x samples
#' @param lengths optional gene-length matrix/vector for TPM & FPKM
#' @param method  one of `NORM_METHODS`
#' @param prior   pseudo-count added before log2 (default 1)
normalise_counts <- function(counts, lengths = NULL, method = "cpm", prior = 1) {
  method <- match.arg(method, unname(NORM_METHODS))
  if (method %in% LENGTH_DEPENDENT && is.null(lengths)) {
    stop("A gene-length matrix is required for ", toupper(sub("^log", "", method)),
         ". Upload the gene-lengths file (e.g. rsem.merged.gene_lengths.tsv) on the Data tab.")
  }
  base <- switch(method,
    raw       = .clean_counts(counts),
    cpm       = calc_cpm(counts),
    logcpm    = calc_cpm(counts),
    tpm       = calc_tpm(counts, lengths),
    logtpm    = calc_tpm(counts, lengths),
    fpkm      = calc_fpkm(counts, lengths),
    logfpkm   = calc_fpkm(counts, lengths),
    deseq2    = calc_deseq2(counts),
    logdeseq2 = calc_deseq2(counts),
    tmm       = calc_tmm(counts),
    logtmm    = calc_tmm(counts),
    stop("Unknown normalisation method: ", method)
  )
  if (grepl("^log", method)) log2(base + prior) else base
}

#' Human-readable axis label for a normalisation method.
norm_label <- function(method, log_prefix = TRUE) {
  lab <- switch(method,
    raw = "raw counts", cpm = "CPM", tpm = "TPM", fpkm = "FPKM",
    deseq2 = "DESeq2-normalised counts", tmm = "TMM-normalised counts",
    logcpm = "CPM", logtpm = "TPM", logfpkm = "FPKM",
    logdeseq2 = "DESeq2-normalised counts", logtmm = "TMM-normalised counts",
    method
  )
  if (log_prefix && grepl("^log", method)) paste0("log2 ", lab) else lab
}

is_log_method <- function(method) grepl("^log", method)

#' Filter lowly expressed genes (edgeR-style "at least `min_count` in at least
#' `min_samples` samples"), used before PCA / clustering.
filter_low_expression <- function(mat, min_count = 1, min_samples = 1) {
  if (min_samples < 1) return(mat)
  keep <- rowSums(mat >= min_count, na.rm = TRUE) >= min_samples
  mat[keep, , drop = FALSE]
}
