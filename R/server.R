# =============================================================================
# server.R -- reactive logic
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#' Collect the style inputs for a plot prefix into a style list.
collect_style <- function(input, p) {
  g <- function(suffix, default) {
    v <- input[[paste0(p, "_", suffix)]]
    if (is.null(v)) default else v
  }
  resolve_style(list(
    base_size        = g("base_size", 13),
    title_size       = g("title_size", 15),
    axis_title_size  = g("axis_title_size", 13),
    axis_text_size   = g("axis_text_size", 11),
    legend_text_size = g("legend_text_size", 11),
    legend_title_size= g("legend_text_size", 11) + 1,
    subtitle_size    = g("title_size", 15) - 3,
    label_size       = g("label_size", 3.2),
    label_bold       = g("label_bold", TRUE),
    label_colour     = g("label_colour", "#000000"),
    point_size       = g("point_size", 1.6),
    point_alpha      = g("point_alpha", 0.7),
    highlight_size   = g("highlight_size", 2.8),
    grid             = g("grid", TRUE)
  ))
}

#' Read a colour input, falling back to the plot's own default.
#' Without this a NULL input (colourpicker not installed, or a control that has
#' not been rendered yet) collapses scale_colour_manual() and every render *and
#' export* of that figure fails with "Insufficient values in manual scale".
collect_colours <- function(input, p, defaults) {
  pick <- function(suffix, default) {
    v <- input[[paste0(p, "_", suffix)]]
    if (is.null(v) || length(v) != 1 || is.na(v) || !nzchar(v)) default else v
  }
  list(up   = pick("up_colour",   defaults$up),
       down = pick("down_colour", defaults$down),
       ns   = pick("ns_colour",   defaults$ns))
}

collect_labels <- function(input, p) {
  list(title = input[[paste0(p, "_title")]],
       subtitle = input[[paste0(p, "_subtitle")]],
       xlab = input[[paste0(p, "_xlab")]],
       ylab = input[[paste0(p, "_ylab")]])
}

app_server <- function(input, output, session) {

  options(shiny.maxRequestSize = 2048 * 1024^2)  # 2 GB: full uncompressed GTFs

  rv <- reactiveValues(
    de = NULL, de_raw = NULL, de_notes = character(0), de_pipeline = NULL,
    counts = NULL, counts_ann = NULL, lengths = NULL, norm_supplied = NULL,
    gene_ann = NULL, universe_ids = NULL, de_info = NULL,
    meta = NULL, loaded = FALSE, source_label = "",
    enrich_single = NULL, enrich_up = NULL, enrich_down = NULL,
    enrich_split = FALSE, enrich_db = NULL, enrich_sets = list(),
    gene_report_html = NULL, gene_report_msg = NULL,
    gal = NULL, gal_objs = list(),
    gal_kind = "", gal_source = "", gal_notes = character(0),
    pv_png = NULL, pv_msg = NULL,
    messages = NULL
  )

  # =========================================================================
  # Data loading
  # =========================================================================

  # Each bundled example is one result-table layout for the same comparison
  # (stim vs mock), sharing the same counts / lengths / VST matrices and sample
  # sheet. data-raw/make_example_de_tables.R builds them from the example counts.
  EXAMPLE_FILES <- list(
    deseq2 = list(
      de  = "data/example_deseq2_results.tsv.gz",
      lab = "bundled example: stim vs mock (DESeq2, annotated)"),
    limma_annotated = list(
      de  = "data/example_limma_results.tsv.gz",
      lab = "bundled example: stim vs mock (limma, annotated)"),
    limma_results = list(
      de  = "data/example_limma_results_plain.tsv.gz",
      ann = "data/example_limma_results.tsv.gz",
      lab = "bundled example: stim vs mock (limma, plain results.tsv)"),
    limma_filtered = list(
      de  = "data/example_limma_results_filtered.tsv.gz",
      ann = "data/example_limma_results.tsv.gz",
      uni = "data/example_limma_results_plain.tsv.gz",
      lab = "bundled example: stim vs mock (limma, results_filtered.tsv)")
  )

  load_example <- function(which = "deseq2") {
    spec <- EXAMPLE_FILES[[which]] %||% EXAMPLE_FILES[["deseq2"]]

    de_raw <- read_table_auto(app_file(spec$de))
    info   <- describe_de_table(de_raw, spec$de)
    std    <- standardise_de_table(de_raw)

    ann <- if (!is.null(spec$ann)) read_gene_annotation(app_file(spec$ann)) else NULL
    fill <- apply_gene_annotation(std$data, ann)

    uni <- if (!is.null(spec$uni)) {
      u <- read_table_auto(app_file(spec$uni))
      um <- detect_de_columns(u)
      if (is.na(um$gene_id)) NULL else unique(as.character(u[[um$gene_id]]))
    } else NULL

    cnt  <- read_matrix_auto(app_file("data/example_rsem_gene_counts.tsv.gz"))
    len  <- read_matrix_auto(app_file("data/example_rsem_gene_lengths.tsv.gz"))
    nrm  <- read_matrix_auto(app_file("data/example_vst.tsv.gz"))
    meta <- read_sample_metadata(app_file("data/example_samplesheet.csv"))

    rv$de <- fill$data; rv$de_raw <- de_raw
    rv$de_notes <- c(std$notes, info$notes); rv$de_pipeline <- info$label
    rv$de_info <- info; rv$gene_ann <- ann; rv$universe_ids <- uni
    rv$counts <- cnt$matrix; rv$counts_ann <- cnt$annotation
    rv$lengths <- len$matrix; rv$norm_supplied <- nrm$matrix
    rv$meta <- meta; rv$loaded <- TRUE
    rv$source_label <- spec$lab
    rv$messages <- list(type = "ok", text = paste0(
      sprintf("Loaded %s \u2014 %s rows. ", spec$lab, format(nrow(fill$data), big.mark = ",")),
      if (fill$n_filled > 0) sprintf("%s gene symbol(s) filled in from the annotated table. ",
                                     format(fill$n_filled, big.mark = ",")) else "",
      if (!is.null(uni)) sprintf("Enrichment background set from the unfiltered results (%s genes). ",
                                 format(length(uni), big.mark = ",")) else "",
      "Also loaded: RSEM counts, gene lengths, VST matrix and sample sheet."))
    reset_enrichment()
  }

  reset_enrichment <- function() {
    rv$enrich_single <- NULL; rv$enrich_up <- NULL; rv$enrich_down <- NULL
    rv$enrich_split <- FALSE; rv$enrich_sets <- list()
    # Gallery results computed from the previous DE table no longer apply.
    # Imported tables are independent of it and stay.
    if (identical(rv$gal_kind, "run")) {
      rv$gal <- NULL; rv$gal_objs <- list(); rv$gal_kind <- ""; rv$gal_source <- ""
      rv$gal_notes <- character(0); rv$pv_png <- NULL; rv$pv_msg <- NULL
    }
  }

  # Analysis tab: run DE from raw counts and hand the result to the Explorer.
  analysis_state <- analysis_server(input, output, session, rv, reset_enrichment)   # kept for tests

  observeEvent(input$load_example, {
    res <- tryCatch({ load_example(input$example_choice %||% "deseq2"); NULL },
                    error = function(e) conditionMessage(e))
    if (!is.null(res)) rv$messages <- list(type = "err",
      text = paste("That bundled example could not be loaded:", res))
  }, ignoreInit = TRUE)

  observeEvent(input$reset_all, {
    for (nm in c("de","de_raw","counts","counts_ann","lengths","norm_supplied","meta",
                 "gene_ann","universe_ids","de_info"))
      rv[[nm]] <- NULL
    rv$loaded <- FALSE; rv$de_notes <- character(0); rv$source_label <- ""
    rv$messages <- list(type = "ok", text = "Cleared. Upload files or load the example.")
    reset_enrichment()
  }, ignoreInit = TRUE)

  # column-mapping UI appears once a DE file is uploaded
  uploaded_de_raw <- reactive({
    req(input$up_de)
    tryCatch(read_table_auto(input$up_de$datapath), error = function(e) e)
  })

  output$de_mapping_ui <- renderUI({
    if (is.null(input$up_de)) return(NULL)
    draw <- uploaded_de_raw()
    if (inherits(draw, "error")) {
      return(div(class = "note-err", strong("Could not read that file: "), conditionMessage(draw)))
    }
    det  <- detect_de_columns(draw)
    info <- describe_de_table(draw, input$up_de$name)
    cols <- c("(auto)", colnames(draw))
    tagList(
      div(class = "note-ok",
          sprintf("Detected format: %s \u2014 %s rows, %d columns.",
                  info$label, format(nrow(draw), big.mark = ","), ncol(draw))),
      if (length(info$notes))
        div(class = "note-warn", lapply(info$notes, function(n) div(style = "margin:2px 0", n))),
      tags$details(
        tags$summary(tags$strong("Column mapping (override if needed)")),
        div(class = "ctl-block",
          lapply(names(DE_SYNONYMS), function(f) {
            selectInput(paste0("map_", f), f, choices = cols,
                        selected = if (is.na(det[[f]])) "(auto)" else det[[f]])
          })))
    )
  })

  output$companion_note <- renderUI({
    bits <- list()
    if (!is.null(input$up_annot)) {
      a <- tryCatch(read_gene_annotation(input$up_annot$datapath, input$up_annot$name),
                    error = function(e) e)
      bits[[length(bits) + 1]] <- if (inherits(a, "error"))
        div(class = "note-err", conditionMessage(a))
      else div(class = "note-ok", sprintf(
        "%s parsed: %s genes, %s with symbols%s.",
        attr(a, "source") %||% "Annotation table",
        format(nrow(a), big.mark = ","),
        format(attr(a, "n_symbols") %||% sum(!is.na(a$gene_name)), big.mark = ","),
        {
          ex <- setdiff(colnames(a), c("gene_id", "gene_name"))
          if (length(ex)) paste0(" \u2014 also carrying ", paste(utils::head(ex, 5), collapse = ", ")) else ""
        }))
    }
    if (!is.null(input$up_universe)) {
      u <- tryCatch(read_table_auto(input$up_universe$datapath), error = function(e) e)
      bits[[length(bits) + 1]] <- if (inherits(u, "error"))
        div(class = "note-err", conditionMessage(u))
      else div(class = "note-ok",
               sprintf("Background: %s genes will be used as the enrichment universe.",
                       format(nrow(u), big.mark = ",")))
    }
    if (!length(bits)) return(NULL)
    tagList(bits)
  })

  # sample subsetting: crucial when a matrix holds several studies
  output$sample_subset_ui <- renderUI({
    cols <- NULL
    if (!is.null(input$up_counts)) {
      cols <- tryCatch(colnames(read_matrix_auto(input$up_counts$datapath)$matrix),
                       error = function(e) NULL)
    }
    if (is.null(cols)) return(NULL)
    tagList(
      div(class = "note-warn",
          sprintf("The count matrix has %d sample columns. Tick only the samples belonging to this comparison.", length(cols))),
      checkboxGroupInput("sample_subset", "Samples to keep", choices = cols, selected = cols)
    )
  })

  observeEvent(input$apply_upload, {
    res <- tryCatch({
      if (is.null(input$up_de)) stop("A differential-expression results file is required.")
      draw <- read_table_auto(input$up_de$datapath)
      info <- describe_de_table(draw, input$up_de$name)
      ov <- list()
      for (f in names(DE_SYNONYMS)) ov[[f]] <- input[[paste0("map_", f)]]
      std <- standardise_de_table(draw, ov)

      # optional companion files (gene annotation, unfiltered results)
      gene_ann <- if (!is.null(input$up_annot))
        read_gene_annotation(input$up_annot$datapath, input$up_annot$name) else d_fetch$ann
      fill <- apply_gene_annotation(std$data, gene_ann)
      std$data <- fill$data

      universe_ids <- NULL
      if (!is.null(input$up_universe)) {
        u  <- read_table_auto(input$up_universe$datapath)
        um <- detect_de_columns(u)
        if (is.na(um$gene_id)) stop("No gene id column found in the unfiltered results file.")
        universe_ids <- unique(as.character(u[[um$gene_id]]))
        universe_ids <- universe_ids[!is.na(universe_ids) & universe_ids != ""]
      }

      cnt <- ann <- len <- nrm <- NULL
      if (!is.null(input$up_counts)) {
        r <- read_matrix_auto(input$up_counts$datapath); cnt <- r$matrix; ann <- r$annotation
      }
      if (!is.null(input$up_lengths)) len <- read_matrix_auto(input$up_lengths$datapath)$matrix
      if (!is.null(input$up_norm))    nrm <- read_matrix_auto(input$up_norm$datapath)$matrix

      keep <- input$sample_subset
      if (!is.null(cnt) && !is.null(keep) && length(keep)) {
        cnt <- cnt[, intersect(colnames(cnt), keep), drop = FALSE]
        if (!is.null(len)) len <- len[, intersect(colnames(len), keep), drop = FALSE]
        if (!is.null(nrm)) nrm <- nrm[, intersect(colnames(nrm), keep), drop = FALSE]
      }
      # No lengths file: take gene lengths from a GTF, uploaded or downloaded.
      len_note <- NULL
      if (is.null(len) && !is.null(cnt)) {
        gl <- if (!is.null(input$up_annot) && is_gtf_file(input$up_annot$datapath, input$up_annot$name)) {
          len_note <- paste("gene lengths from", input$up_annot$name)
          gtf_gene_lengths(input$up_annot$datapath)
        } else if (!is.null(d_fetch$len)) {
          len_note <- paste("gene lengths from", d_fetch$label)
          d_fetch$len
        }
        if (!is.null(gl)) {
          len <- match_gene_lengths(gl, rownames(cnt))
          if (!length(len)) { len <- NULL; len_note <- NULL }
        }
      }

      meta <- if (!is.null(input$up_meta)) read_sample_metadata(input$up_meta$datapath) else {
        ids <- colnames(cnt %||% nrm)
        if (is.null(ids)) NULL else scaffold_metadata(ids)
      }
      if (!is.null(meta) && !is.null(cnt)) {
        meta <- meta[meta$sample %in% colnames(cnt), , drop = FALSE]
      }
      list(std = std, cnt = cnt, ann = ann, len = len, nrm = nrm, meta = meta,
           name = input$up_de$name, info = info, gene_ann = gene_ann,
           universe_ids = universe_ids, n_filled = fill$n_filled, len_note = len_note)
    }, error = function(e) e)

    if (inherits(res, "error")) {
      rv$messages <- list(type = "err", text = conditionMessage(res)); return(invisible())
    }
    rv$de <- res$std$data; rv$de_raw <- NULL
    rv$de_notes <- c(res$std$notes, res$info$notes); rv$de_pipeline <- res$info$label
    rv$de_info <- res$info; rv$gene_ann <- res$gene_ann
    rv$universe_ids <- res$universe_ids
    rv$counts <- res$cnt; rv$counts_ann <- res$ann
    rv$lengths <- res$len; rv$norm_supplied <- res$nrm
    rv$meta <- res$meta; rv$loaded <- TRUE
    rv$source_label <- sprintf("%s (%s)", res$name, res$info$label)
    reset_enrichment()

    msg <- sprintf("Loaded %s: %s rows detected as %s.", res$name,
                   format(nrow(res$std$data), big.mark = ","), res$info$label)
    extra <- c()
    if (!is.null(res$gene_ann) && res$n_filled > 0)
      extra <- c(extra, sprintf("%s gene symbol(s) from the %s",
                                format(res$n_filled, big.mark = ","),
                                attr(res$gene_ann, "source") %||% "annotation table"))
    if (!is.null(res$universe_ids))
      extra <- c(extra, sprintf("enrichment background of %s genes",
                                format(length(res$universe_ids), big.mark = ",")))
    if (!is.null(res$cnt)) extra <- c(extra, sprintf("counts %d x %d", nrow(res$cnt), ncol(res$cnt)))
    if (!is.null(res$len)) extra <- c(extra, res$len_note %||% "gene lengths")
    if (!is.null(res$nrm)) extra <- c(extra, "pre-normalised matrix")
    if (length(extra)) msg <- paste0(msg, " Also loaded: ", paste(extra, collapse = ", "), ".")
    rv$messages <- list(type = "ok", text = msg)
  }, ignoreInit = TRUE)

  output$load_messages <- renderUI({
    m <- rv$messages; if (is.null(m)) return(NULL)
    cls <- switch(m$type, ok = "note-ok", warn = "note-warn", "note-err")
    div(class = cls, m$text)
  })

  output$data_status <- renderUI({
    pill <- function(ok, label) span(class = paste("status-pill", if (ok) "pill-ok" else "pill-no"),
                                     paste0(label, if (ok) " \u2713" else " \u2014"))
    tagList(
      div(style = "margin-bottom:10px",
        pill(!is.null(rv$de), "DE results"),
        pill(!is.null(rv$counts), "raw counts"),
        pill(!is.null(rv$lengths), "gene lengths"),
        pill(!is.null(rv$norm_supplied), "normalised matrix"),
        pill(!is.null(rv$meta), "sample sheet")),
      if (nzchar(rv$source_label %||% "")) div(style = "color:#66737f;font-size:.9em;margin-bottom:8px",
                                               paste("Current source:", rv$source_label)),
      if (length(rv$de_notes)) div(class = "note-warn",
        strong("Notes on the loaded table: "), paste(rv$de_notes, collapse = " ")),
      if (is.null(rv$de)) div(class = "note-warn",
        "No data loaded yet. Upload a results file, or press ",
        strong("Load bundled example"), " to explore the app.")
    )
  })

  # =========================================================================
  # Derived data
  # =========================================================================

  de_classified <- reactive({
    req(rv$de)
    classify_degs(rv$de, input$padj_cutoff, input$lfc_cutoff, input$basemean_cutoff)
  })

  output$threshold_summary <- renderUI({
    if (is.null(rv$de)) return(NULL)
    d <- de_classified()
    div(style = "font-size:.88em;color:#42525f;margin-top:6px",
        sprintf("%s tested \u2014 ", format(sum(!is.na(d$padj)), big.mark = ",")),
        span(style = "color:#D6604D;font-weight:600", sprintf("%d up", sum(d$status == "Up"))),
        " / ",
        span(style = "color:#4393C3;font-weight:600", sprintf("%d down", sum(d$status == "Down"))))
  })

  # ---- Organism (app-wide) ---------------------------------------------------
  # Guessed from the ids whenever new data arrives; drives enrichment (OrgDb,
  # KEGG code, MSigDB species) and the gene look-up links.
  organism <- reactive(organism_info(input$organism %||% "human",
    list(species = input$org_other_species, orgdb = input$org_other_orgdb, kegg = input$org_other_kegg)))
  org_guess <- reactive({
    ids <- if (!is.null(rv$counts)) rownames(rv$counts) else rv$de$gene_id
    if (is.null(ids)) NA_character_ else guess_organism(ids)
  })
  observeEvent(org_guess(), {
    g <- org_guess()
    if (!is.na(g) && !identical(g, isolate(input$organism))) updateSelectInput(session, "organism", selected = g)
  })
  observeEvent(organism(), {
    if (nzchar(organism()$species)) updateTextInput(session, "g_species", value = organism()$species)
  }, ignoreInit = TRUE)
  output$organism_note <- renderUI({
    o <- organism(); g <- org_guess()
    msgs <- list(
      if (!is.na(g)) span(sprintf("Ids look like %s.", ORGANISMS[[g]]$label))
      else if (!is.null(rv$de)) span("The ids do not say which organism this is; set it here."),
      if (!is.null(m <- organism_orgdb_missing(o))) span(" GO / KEGG need ", tags$code(if (startsWith(m, "BiocManager")) m else "an OrgDb package"), "."))
    msgs <- Filter(Negate(is.null), msgs)
    if (!length(msgs)) return(NULL)
    div(style = "font-size:.82em;color:#42525f", msgs)
  })

  # available expression matrices, by key
  matrix_choices <- reactive({
    ch <- c()
    if (!is.null(rv$counts)) ch <- c(ch, "Computed from raw counts" = "counts")
    if (!is.null(rv$norm_supplied)) ch <- c(ch, "Supplied normalised matrix" = "supplied")
    ch
  })

  #' Return an expression matrix for a given source + normalisation choice.
  get_expression <- function(source, method) {
    if (identical(source, "supplied")) {
      m <- rv$norm_supplied
      if (is.null(m)) stop("No pre-normalised matrix has been loaded.")
      return(m)
    }
    if (is.null(rv$counts)) {
      stop("No raw count matrix loaded. Upload one on the Data tab (or use the supplied matrix).")
    }
    normalise_counts(rv$counts, rv$lengths, method)
  }

  observe({
    ch <- matrix_choices()
    sel_h <- if ("counts" %in% ch) "counts" else if (length(ch)) ch[1] else NULL
    updateSelectInput(session, "h_matrix", choices = ch, selected = sel_h)
    updateSelectInput(session, "p_matrix", choices = ch, selected = sel_h)
  })

  observe({
    req(rv$meta)
    grps <- sort(unique(rv$meta$condition))
    updateSelectInput(session, "s_group_x", choices = grps, selected = grps[1])
    updateSelectInput(session, "s_group_y", choices = grps,
                      selected = if (length(grps) > 1) grps[2] else grps[1])
    updateSelectInput(session, "h_reference", choices = grps, selected = grps[1])
    cols <- setdiff(colnames(rv$meta), "sample")
    updateSelectInput(session, "p_colour_by", choices = cols,
                      selected = if ("condition" %in% cols) "condition" else cols[1])
    updateSelectInput(session, "p_shape_by", choices = c("(none)" = "(none)", cols),
                      selected = "(none)")
    updateCheckboxGroupInput(session, "h_annotate_by", choices = cols, selected = cols)
  })

  observe({
    req(rv$de)
    ch <- stats::setNames(rv$de$gene_id, paste0(rv$de$display_gene, "  (", rv$de$gene_id, ")"))
    for (p in c("v", "m", "s", "h", "g")) {
      updateSelectizeInput(session, paste0(p, "_gene_dropdown"), choices = ch, server = TRUE)
    }
  })

  observe({
    gs <- rv$enrich_sets
    ch <- if (!length(gs)) c("Run enrichment first" = "") else stats::setNames(names(gs), names(gs))
    for (p in c("v", "m", "s", "h", "g")) updateSelectInput(session, paste0(p, "_pathway_term"), choices = ch)
  })

  # generic gene-selection resolver per prefix
  genes_for <- function(p) {
    d <- de_classified()
    mode <- input[[paste0(p, "_gene_mode")]] %||% "none"
    resolve_gene_selection(
      mode, d,
      top_n = input[[paste0(p, "_top_n")]] %||% 15,
      dropdown = input[[paste0(p, "_gene_dropdown")]] %||% character(0),
      pasted = input[[paste0(p, "_gene_paste")]] %||% "",
      pathway = rv$enrich_sets[[input[[paste0(p, "_pathway_term")]] %||% ""]] %||% character(0))
  }

  # feedback on unmatched pasted genes
  lapply(c("v", "m", "s", "h", "g"), function(p) {
    output[[paste0(p, "_paste_feedback")]] <- renderUI({
      txt <- input[[paste0(p, "_gene_paste")]]
      if (is.null(txt) || !nzchar(trimws(txt)) || is.null(rv$de)) return(NULL)
      r <- match_user_genes(txt, rv$de)
      tagList(
        div(style = "font-size:.85em;color:#1d6f42", sprintf("%d gene(s) matched.", length(r$ids))),
        if (length(r$missing))
          div(class = "note-warn", style = "font-size:.85em",
              sprintf("Not found: %s", paste(utils::head(r$missing, 12), collapse = ", ")),
              if (length(r$missing) > 12) sprintf(" (+%d more)", length(r$missing) - 12))
      )
    })
  })
  output$h_paste_feedback <- renderUI({
    txt <- input$h_gene_paste
    if (is.null(txt) || !nzchar(trimws(txt)) || is.null(rv$de)) return(NULL)
    r <- match_user_genes(txt, rv$de)
    tagList(div(style = "font-size:.85em;color:#1d6f42", sprintf("%d gene(s) matched.", length(r$ids))),
            if (length(r$missing)) div(class = "note-warn", style = "font-size:.85em",
                                       sprintf("Not found: %s", paste(utils::head(r$missing, 12), collapse = ", "))))
  })

  # =========================================================================
  # Data tab: sample sheet + computed matrices
  # =========================================================================

  # Annotation fetched from Ensembl / GENCODE: kept for the next Load, and
  # applied straight away to data already loaded.
  d_fetch <- reactiveValues(ann = NULL, len = NULL, label = NULL)
  annotation_fetch_server(input, output, session, "d",
    ids = reactive(if (!is.null(rv$counts)) rownames(rv$counts) else rv$de$gene_id),
    on_loaded = function(ann, len, label, tx = NULL, organism_key = NULL) {
      attr(ann, "source") <- label
      d_fetch$ann <- ann; d_fetch$len <- len; d_fetch$label <- label
      if (!is.null(organism_key)) updateSelectInput(session, "organism", selected = organism_key)
      done <- character(0)
      if (!is.null(rv$counts)) {
        l <- match_gene_lengths(len, rownames(rv$counts))
        if (length(l)) {
          rv$lengths <- l
          done <- c(done, sprintf("gene lengths for %s of %s genes",
                                  format(length(l), big.mark = ","), format(nrow(rv$counts), big.mark = ",")))
        }
      }
      if (!is.null(rv$de)) {
        f <- apply_gene_annotation(rv$de, ann)
        rv$de <- f$data
        if (is.null(rv$gene_ann)) rv$gene_ann <- ann
        if (f$n_filled > 0) done <- c(done, sprintf("%s gene symbols filled in", format(f$n_filled, big.mark = ",")))
      }
      showNotification(paste0(label, " ready. ",
        if (length(done)) paste0("Applied: ", paste(done, collapse = "; "), ".")
        else "It will be used when you press 'Load uploaded files'."), type = "message", duration = 12)
    })

  # Sample sheet: uploaded, or built here from the matrix column names.
  sheet_editor_server(
    input, output, session, "d",
    get_sheet = function() rv$meta,
    set_sheet = function(x) rv$meta <- x,
    get_types = function() {
      v <- setdiff(colnames(rv$meta %||% data.frame()), "sample")
      stats::setNames(as.list(rep("factor", length(v))), v)
    },
    protected = "condition")

  output$d_sheet_note <- renderUI({
    if (is.null(rv$meta)) return(div(class = "note-warn",
      "Load a count or normalised matrix and the sample sheet is built here from its column names."))
    g <- rv$meta$condition
    if (all(is.na(g) | g %in% c("", "ungrouped"))) return(div(class = "note-warn",
      "Every sample is still 'ungrouped'. Fill the condition column (e.g. 'Fill from names' or",
      "'Set for selected rows') so Scatter, PCA and Heatmap know the groups."))
    div(class = "note-ok", sprintf("%d samples in %d group(s): %s.", nrow(rv$meta),
        length(unique(g[nzchar(g)])), paste(utils::head(sort(unique(g[nzchar(g)])), 8), collapse = ", ")))
  })

  computed_matrix <- reactive({
    req(rv$counts)
    m <- input$dl_norm_method %||% "tpm"
    if (!norm_requirements_met(m, !is.null(rv$lengths))) return(NULL)
    normalise_counts(rv$counts, rv$lengths, m)
  })

  output$norm_preview_note <- renderUI({
    if (is.null(rv$counts)) return(div(class = "note-warn",
      "Upload a raw count matrix to compute TPM, FPKM, CPM, DESeq2 and TMM values."))
    m <- input$dl_norm_method %||% "tpm"
    if (!norm_requirements_met(m, !is.null(rv$lengths))) {
      return(div(class = "note-warn",
        sprintf("%s needs a gene-length matrix. Upload one (e.g. rsem.merged.gene_lengths.tsv) to enable it.",
                toupper(sub("^log", "", m)))))
    }
    NULL
  })

  output$norm_preview <- DT::renderDT({
    mat <- computed_matrix(); req(mat)
    head_mat <- utils::head(mat, 50)
    df <- data.frame(gene_id = rownames(head_mat), round(as.data.frame(head_mat), 3),
                     check.names = FALSE)
    if (!is.null(rv$counts_ann)) {
      df <- cbind(gene_name = rv$counts_ann$gene_name[match(df$gene_id, rv$counts_ann$gene_id)], df)
    }
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  dl_matrix <- function(fmt) {
    downloadHandler(
      filename = function() stamp_name(paste0(input$dl_norm_method %||% "matrix", "_matrix"), fmt),
      content = function(file) {
        mat <- computed_matrix()
        if (is.null(mat)) stop("This matrix cannot be computed with the files currently loaded.")
        write_matrix_file(mat, file, fmt, annotation = rv$counts_ann)
      })
  }
  output$dl_norm_csv <- dl_matrix("csv")
  output$dl_norm_tsv <- dl_matrix("tsv")

  # =========================================================================
  # Volcano
  # =========================================================================

  volcano_obj <- reactive({
    d <- de_classified(); lb <- collect_labels(input, "v")
    cl <- collect_colours(input, "v", DEG_COLOURS)
    make_volcano(d, input$padj_cutoff, input$lfc_cutoff, input$basemean_cutoff,
                 label_ids = genes_for("v"),
                 up_colour = cl$up, down_colour = cl$down, ns_colour = cl$ns,
                 style = collect_style(input, "v"),
                 show_cutoff_lines = isTRUE(input$v_cutoff_lines),
                 title = lb$title, subtitle = lb$subtitle, xlab = lb$xlab, ylab = lb$ylab,
                 max_overlaps = input$v_max_overlaps %||% 30)
  })
  output$volcano_plot <- renderPlot(volcano_obj())
  output$volcano_table <- DT::renderDT({
    DT::datatable(deg_display_table(significant_degs(de_classified())),
                  options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  output$v_dl_table <- downloadHandler(
    filename = function() stamp_name("significant_genes", "csv"),
    content = function(f) utils::write.csv(deg_display_table(significant_degs(de_classified())),
                                            f, row.names = FALSE))
  output$v_dl_all <- downloadHandler(
    filename = function() stamp_name("all_genes", "csv"),
    content = function(f) utils::write.csv(de_full_table(de_classified()), f, row.names = FALSE, na = ""))
  output$v_dl_xlsx <- downloadHandler(
    filename = function() stamp_name("de_results", "xlsx"),
    content = function(f) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) stop('Excel export needs install.packages("openxlsx")')
      write_de_workbook(de_classified(), f, input$padj_cutoff, input$lfc_cutoff, input$basemean_cutoff)
    })

  # =========================================================================
  # MA
  # =========================================================================

  ma_obj <- reactive({
    d <- de_classified(); lb <- collect_labels(input, "m")
    cl <- collect_colours(input, "m", DEG_COLOURS)
    make_ma(d, input$padj_cutoff, input$lfc_cutoff, input$basemean_cutoff,
            label_ids = genes_for("m"),
            up_colour = cl$up, down_colour = cl$down, ns_colour = cl$ns,
            style = collect_style(input, "m"),
            show_cutoff_lines = isTRUE(input$m_cutoff_lines),
            title = lb$title, subtitle = lb$subtitle, xlab = lb$xlab, ylab = lb$ylab,
            max_overlaps = input$m_max_overlaps %||% 30)
  })
  output$ma_plot <- renderPlot(ma_obj())
  output$ma_table <- DT::renderDT({
    DT::datatable(deg_display_table(significant_degs(de_classified())),
                  options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  # =========================================================================
  # Scatter
  # =========================================================================

  output$scatter_note <- renderUI({
    if (is.null(rv$counts) && is.null(rv$norm_supplied))
      return(div(class = "note-warn", "The scatter plot needs an expression matrix. Upload raw counts or a normalised matrix on the Data tab."))
    if (is.null(rv$meta))
      return(div(class = "note-warn", "A sample sheet is required so the app knows which samples form each group."))
    if (identical(input$s_value, "fpkm") || identical(input$s_value, "tpm")) {
      if (is.null(rv$lengths) && !is.null(rv$counts))
        return(div(class = "note-warn", "TPM and FPKM need a gene-length matrix. Upload one, or pick CPM / DESeq2 / TMM instead."))
    }
    NULL
  })

  # FPKM / TPM need gene lengths. When none are loaded (e.g. a run from the
  # Analysis tab), move the scatter off them so it draws instead of erroring.
  observe({
    if (!is.null(rv$counts) && is.null(rv$lengths) &&
        isolate(input$s_value %||% "fpkm") %in% c("fpkm", "tpm")) {
      updateSelectInput(session, "s_value", selected = "deseq2")
    }
  })

  scatter_obj <- reactive({
    req(rv$meta, input$s_group_x, input$s_group_y)
    d <- de_classified()
    expr <- if (identical(input$s_value, "asis")) {
      m <- rv$norm_supplied; if (is.null(m)) stop("No supplied matrix is loaded."); m
    } else get_expression("counts", input$s_value)

    sdf <- group_mean_expression(expr, rv$meta, input$s_group_x, input$s_group_y,
                                 log_transform = isTRUE(input$s_log))
    ot <- match_user_genes(input$s_offtarget %||% "", rv$de)$ids
    fc <- match_user_genes(input$s_focus %||% "", rv$de)$ids
    vl <- paste0(if (isTRUE(input$s_log)) "log2 " else "",
                 if (identical(input$s_value, "asis")) "expression" else norm_label(input$s_value, FALSE))
    lb <- collect_labels(input, "s")
    cl <- collect_colours(input, "s", list(up = "#E8A33D", down = "#2C5AA0", ns = "grey80"))
    make_expression_scatter(sdf, d, input$padj_cutoff, input$lfc_cutoff, input$basemean_cutoff,
                            label_ids = genes_for("s"), offtarget_ids = ot, focus_ids = fc,
                            up_colour = cl$up, down_colour = cl$down, ns_colour = cl$ns,
                            style = collect_style(input, "s"),
                            show_identity = isTRUE(input$s_show_identity),
                            identity_colour = input$s_identity_colour %||% "grey40",
                            show_counts = isTRUE(input$s_show_counts),
                            group_x_name = input$s_group_x, group_y_name = input$s_group_y,
                            value_label = vl,
                            title = lb$title, subtitle = lb$subtitle,
                            xlab = lb$xlab, ylab = lb$ylab,
                            max_overlaps = input$s_max_overlaps %||% 30)
  })
  output$scatter_plot <- renderPlot(scatter_obj())

  # =========================================================================
  # PCA + scree
  # =========================================================================

  output$pca_note <- renderUI({
    if (is.null(rv$counts) && is.null(rv$norm_supplied))
      return(div(class = "note-warn", "PCA needs an expression matrix. Upload one on the Data tab."))
    m <- input$p_norm %||% "logcpm"
    if (identical(input$p_matrix, "counts") && !norm_requirements_met(m, !is.null(rv$lengths)))
      return(div(class = "note-warn",
                 sprintf("%s needs a gene-length matrix.", toupper(sub("^log", "", m)))))
    NULL
  })

  pca_result <- reactive({
    req(rv$meta)
    expr <- get_expression(input$p_matrix %||% "counts", input$p_norm %||% "logcpm")
    keep <- intersect(colnames(expr), rv$meta$sample)
    validate(need(length(keep) >= 2, "At least 2 samples from the sample sheet must be present in the matrix."))
    run_pca(expr[, keep, drop = FALSE], ntop = input$p_ntop %||% 500,
            scale_genes = isTRUE(input$p_scale_genes))
  })

  pca_obj <- reactive({
    lb <- collect_labels(input, "p")
    make_pca_plot(pca_result(), rv$meta,
                  colour_by = input$p_colour_by, shape_by = input$p_shape_by,
                  pc_x = input$p_pc_x %||% 1, pc_y = input$p_pc_y %||% 2,
                  style = collect_style(input, "p"),
                  show_labels = isTRUE(input$p_show_labels),
                  show_ellipse = isTRUE(input$p_show_ellipse),
                  title = lb$title, subtitle = lb$subtitle, xlab = lb$xlab, ylab = lb$ylab,
                  max_overlaps = input$p_max_overlaps %||% 30)
  })
  scree_obj <- reactive({
    make_scree_plot(pca_result(), n_pc = input$p_scree_n %||% 10,
                    style = collect_style(input, "p"),
                    show_cumulative = isTRUE(input$p_scree_cumulative))
  })
  output$pca_plot <- renderPlot(pca_obj())
  output$scree_plot <- renderPlot(scree_obj())
  output$pca_table <- DT::renderDT({
    p <- pca_result()
    n <- min(length(p$percent), 20)
    DT::datatable(data.frame(
      PC = paste0("PC", seq_len(n)),
      `Variance explained (%)` = round(p$percent[seq_len(n)], 2),
      `Cumulative (%)` = round(cumsum(p$percent)[seq_len(n)], 2),
      check.names = FALSE), options = list(pageLength = 10, dom = "tip"), rownames = FALSE)
  })

  # =========================================================================
  # Heatmap
  # =========================================================================

  heatmap_gene_ids <- reactive({
    d <- de_classified()
    switch(input$h_gene_mode,
      top_n     = head(significant_degs(d)$gene_id, input$h_top_n %||% 30),
      all_sig   = significant_degs(d)$gene_id,
      all_genes = character(0),
      dropdown  = input$h_gene_dropdown %||% character(0),
      user_list = match_user_genes(input$h_gene_paste %||% "", rv$de)$ids,
      pathway   = rv$enrich_sets[[input$h_pathway_term %||% ""]] %||% character(0),
      character(0))
  })

  heatmap_obj <- reactive({
    req(rv$meta)
    expr <- get_expression(input$h_matrix %||% "counts", input$h_norm %||% "logcpm")
    all_genes <- identical(input$h_gene_mode, "all_genes")
    ids <- heatmap_gene_ids()
    if (!all_genes) validate(need(length(ids) >= 2, "Select at least 2 genes for the heatmap."))
    make_heatmap(ids, expr, rv$meta, rv$de, all_genes = all_genes,
                 scale_rows = isTRUE(input$h_scale_rows),
                 cluster_rows = isTRUE(input$h_cluster_rows),
                 cluster_cols = isTRUE(input$h_cluster_cols),
                 show_rownames = isTRUE(input$h_show_rownames),
                 show_colnames = isTRUE(input$h_show_colnames),
                 annotate_by = input$h_annotate_by %||% "condition",
                 reference = input$h_reference,
                 palette = input$h_palette %||% "RdBu",
                 style = collect_style(input, "h"),
                 max_genes = if (all_genes) (input$h_max_genes %||% 2000) else 5000,
                 title = axis_label(input$h_title, "Heatmap"),
                 gap_between_groups = isTRUE(input$h_gap))
  })

  output$heatmap_note <- renderUI({
    if (is.null(rv$counts) && is.null(rv$norm_supplied))
      return(div(class = "note-warn", "The heatmap needs an expression matrix. Upload one on the Data tab."))
    h <- tryCatch(heatmap_obj(), error = function(e) e)
    if (inherits(h, "error")) return(NULL)
    msgs <- list()
    msgs[[length(msgs) + 1]] <- sprintf("%s gene(s) plotted across %d sample(s); %s on the left.",
      format(h$n_genes, big.mark = ","), length(h$samples), input$h_reference %||% "reference group")
    if (isTRUE(h$truncated)) msgs[[length(msgs) + 1]] <-
      "Row cap reached: the most variable genes were kept."
    if (isTRUE(h$rownames_hidden)) msgs[[length(msgs) + 1]] <-
      "Gene names hidden automatically (more than 120 rows)."
    if (isTRUE(input$h_cluster_cols)) msgs[[length(msgs) + 1]] <-
      "Sample clustering is ON, so the dendrogram controls column order rather than the control/treated split."
    div(class = if (isTRUE(input$h_cluster_cols)) "note-warn" else "note-ok",
        paste(unlist(msgs), collapse = " "))
  })

  output$heatmap_plot <- renderPlot({
    h <- heatmap_obj(); grid::grid.newpage(); grid::grid.draw(h$plot$gtable)
  })

  # =========================================================================
  # Enrichment
  # =========================================================================

  output$enrich_pkg_note <- renderUI({
    miss <- enrichment_missing_packages(input$e_database %||% "BP", input$e_plot_type %||% "dot", organism())
    if (!length(miss)) {
      if ((input$e_database %||% "BP") %in% GO_RELEASE_DATABASES) {
        return(div(class = "note-ok", sprintf(paste(
          "Gene sets from the official Gene Ontology release %s (go-basic.obo + UniProt",
          "human GAF, propagated over is_a and part_of), matched to your table by gene symbol."),
          go_release_version())))
      }
      return(NULL)
    }
    div(class = "note-warn",
        if (all(startsWith(miss, "#"))) span(paste(sub("^# ", "", miss), collapse = " "))
        else tagList(strong("Optional packages needed for this selection are not installed. Run:"),
                     tags$pre(paste(miss, collapse = "\n"))))
  })

  observeEvent(input$e_run, {
    req(rv$de)
    d <- de_classified()
    g_up <- enrichment_ids(rv$de, d$gene_id[d$status == "Up"])
    g_dn <- enrichment_ids(rv$de, d$gene_id[d$status == "Down"])
    universe <- enrichment_universe(rv$de, rv$universe_ids, rownames(rv$counts))
    if (length(universe) <= nrow(rv$de) * 1.1 && isTRUE(rv$de_info$filtered)) {
      showNotification(paste(
        "This results table is pre-filtered, so enrichment is testing the DEGs against themselves.",
        "Load the unfiltered '*.results.tsv' as the enrichment background on the Data tab."),
        type = "warning", duration = 15)
    }
    do_split <- identical(input$e_gene_source, "all") && isTRUE(input$e_split)
    run1 <- function(g) tryCatch(
      run_enrichment(g, universe, input$e_database, input$e_pvalue_cutoff, input$e_qvalue_cutoff,
                     symbol_map = enrichment_symbol_map(rv$de), organism = organism()),
      error = function(e) e)

    if (do_split) {
      withProgress(message = "Running enrichment (up + down)...", value = 0.3, {
        ru <- run1(g_up); incProgress(0.4); rd <- run1(g_dn)
      })
      ok_u <- !inherits(ru, "error"); ok_d <- !inherits(rd, "error")
      if (!ok_u) showNotification(paste("Up-regulated:", conditionMessage(ru)), type = "error", duration = 12)
      if (!ok_d) showNotification(paste("Down-regulated:", conditionMessage(rd)), type = "error", duration = 12)
      if (!ok_u && !ok_d) return(invisible())
      rv$enrich_single <- NULL
      rv$enrich_up <- if (ok_u) ru else NULL
      rv$enrich_down <- if (ok_d) rd else NULL
      rv$enrich_split <- TRUE; rv$enrich_db <- input$e_database
      rv$enrich_sets <- enrichment_sets_to_rows(c(
        if (ok_u) enrichment_gene_sets(ru, input$e_database, "Up") else list(),
        if (ok_d) enrichment_gene_sets(rd, input$e_database, "Down") else list()), rv$de)
      showNotification(sprintf("Enrichment complete: %d up, %d down term(s).",
        if (ok_u) nrow(as.data.frame(ru)) else 0,
        if (ok_d) nrow(as.data.frame(rd)) else 0), type = "message")
    } else {
      g <- switch(input$e_gene_source, up = g_up, down = g_dn, all = c(g_up, g_dn))
      res <- withProgress(message = "Running enrichment...", value = 0.4, run1(g))
      if (inherits(res, "error")) {
        showNotification(paste("Enrichment failed:", conditionMessage(res)), type = "error", duration = 12)
        return(invisible())
      }
      n <- nrow(as.data.frame(res))
      showNotification(if (n == 0) "Enrichment ran but found no significant terms."
                       else sprintf("Enrichment complete: %d term(s).", n),
                       type = if (n == 0) "warning" else "message")
      rv$enrich_single <- res; rv$enrich_up <- NULL; rv$enrich_down <- NULL
      rv$enrich_split <- FALSE; rv$enrich_db <- input$e_database
      rv$enrich_sets <- enrichment_sets_to_rows(enrichment_gene_sets(res, input$e_database), rv$de)
    }
  }, ignoreInit = TRUE)

  # Plain-text results: exactly what the CSV / Excel downloads contain.
  enrich_table_df <- reactive({
    if (isTRUE(rv$enrich_split)) {
      rbind(enrich_result_table(rv$enrich_up,   rv$de, "Up"),
            enrich_result_table(rv$enrich_down, rv$de, "Down"))
    } else {
      req(rv$enrich_single)
      enrich_result_table(rv$enrich_single, rv$de)
    }
  })

  # Same table dressed for the screen: numbers rounded, and the two gene
  # columns collapsed to a click-to-expand summary so a term carrying 200 genes
  # no longer stretches its row down the page.
  enrich_display_df <- reactive({
    df <- enrich_table_df()
    if (is.null(df) || !nrow(df)) return(df)
    df <- format_enrich_numbers(df)
    for (cl in intersect(ENRICH_GENE_COLS, colnames(df))) {
      df[[cl]] <- collapsible_gene_cell(df[[cl]])
    }
    df
  })

  output$enrich_table <- DT::renderDT({
    df <- enrich_display_df()
    validate(need(!is.null(df) && nrow(df) > 0,
                  "No enrichment results yet \u2014 choose your options and press Run enrichment."))
    gene_cols <- which(colnames(df) %in% ENRICH_GENE_COLS) - 1L
    DT::datatable(
      df, escape = FALSE, rownames = FALSE,
      colnames = enrich_table_colnames(df),
      options = list(
        pageLength = 10, scrollX = TRUE, autoWidth = FALSE,
        columnDefs = list(
          list(targets = as.list(gene_cols), className = "gene-col", orderable = FALSE),
          list(targets = "_all", className = "dt-head-nowrap"))
      ))
  })

  fc_named_symbol <- reactive({
    req(rv$de); stats::setNames(rv$de$log2FoldChange, rv$de$display_gene)
  })

  enrich_plot_obj <- reactive({
    validate(need(!is.null(rv$enrich_single) || !is.null(rv$enrich_up) || !is.null(rv$enrich_down),
                  "No enrichment results yet \u2014 choose your options and press Run enrichment."))
    lb <- collect_labels(input, "e")
    ttl <- axis_label(lb$title, "Enrichment results")
    build_enrichment_figure(
      plot_type = input$e_plot_type, database = rv$enrich_db,
      res_single = rv$enrich_single, res_up = rv$enrich_up, res_down = rv$enrich_down,
      split = isTRUE(rv$enrich_split), fc_named = fc_named_symbol(),
      top_n = input$e_top_terms, style = collect_style(input, "e"),
      title = ttl, xlab = lb$xlab, ylab = lb$ylab)
  })

  output$enrich_plot_ui <- renderUI({
    h <- if (isTRUE(rv$enrich_split)) 880 else 620
    if (input$e_plot_type %in% c("cnet", "heat")) h <- max(h, 720)
    plotOutput("enrich_plot", height = paste0(h, "px"))
  })
  output$enrich_plot <- renderPlot(enrich_plot_obj())

  output$e_dl_csv <- downloadHandler(
    filename = function() stamp_name(paste0("enrichment_", rv$enrich_db %||% "results"), "csv"),
    content = function(f) utils::write.csv(enrich_table_df(), f, row.names = FALSE))
  output$e_dl_xlsx <- downloadHandler(
    filename = function() stamp_name(paste0("enrichment_", rv$enrich_db %||% "results"), "xlsx"),
    content = function(f) openxlsx::write.xlsx(enrich_table_df(), f))

  # =========================================================================
  # Enrichment gallery
  # =========================================================================

  set_gallery <- function(tables, kind, source, notes = character(0), objs = list()) {
    rv$gal <- tables; rv$gal_objs <- objs; rv$gal_kind <- kind
    rv$gal_source <- source; rv$gal_notes <- notes
    rv$pv_png <- NULL; rv$pv_msg <- NULL
  }

  do_gallery_import <- function(paths, names, hint, label) {
    imp <- tryCatch(import_enrichment_files(paths, names, hint, enrichment_gene_view(rv$de)), error = function(e) e)
    if (inherits(imp, "error")) {
      showNotification(paste("Import failed:", conditionMessage(imp)), type = "error", duration = 15)
      return(invisible(FALSE))
    }
    if (!length(imp$tables)) {
      showNotification(paste(c("Nothing could be imported.", imp$notes), collapse = " "),
                       type = "error", duration = 20)
      return(invisible(FALSE))
    }
    set_gallery(imp$tables, "import", paste(label, paste(names, collapse = ", ")), imp$notes)
    showNotification(sprintf("Imported %d database table(s).", length(imp$tables)), type = "message")
    invisible(TRUE)
  }

  observeEvent(input$x_import, {
    if (is.null(input$x_files)) {
      showNotification("Choose one or more enrichment tables first.", type = "warning")
      return(invisible())
    }
    do_gallery_import(input$x_files$datapath, input$x_files$name, input$x_hint %||% "auto",
                      "imported:")
  }, ignoreInit = TRUE)

  # The demo tables go through exactly the same import path as a user's files.
  demo_files <- function() example_enrichment_files(app_file("data", "example_enrichment"))

  observeEvent(input$x_import_demo, {
    f <- demo_files()
    if (!length(f)) {
      showNotification("The example files are missing from data/example_enrichment/.",
                       type = "error"); return(invisible())
    }
    do_gallery_import(f, basename(f), "auto", "example files (computed from the bundled demo data):")
  }, ignoreInit = TRUE)

  output$x_dl_demo <- downloadHandler(
    filename = function() "example_enrichment_tables.zip",
    contentType = "application/zip",
    content = function(f) {
      tmp <- tempfile("demo_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
      src <- c(demo_files(), app_file("data", "example_enrichment", "README.md"),
               app_file("data", "example_enrichment", "SOURCE.txt"))
      file.copy(src[file.exists(src)], tmp)
      r <- tryCatch({ zip_directory(tmp, f); NULL }, error = function(e) conditionMessage(e))
      if (!is.null(r)) showNotification(paste("Could not build the archive:", r), type = "error")
    })

  observeEvent(input$x_run, {
    req(rv$de)
    miss <- enrichment_missing_packages("BP", organism = organism())
    if (length(miss)) {
      showNotification(paste("Required package(s) not installed. Run:", paste(miss, collapse = "; ")),
                       type = "error", duration = 20)
      return(invisible())
    }
    d <- de_classified()
    genes <- enrichment_ids(rv$de, switch(input$x_gene_source %||% "all",
                    up = d$gene_id[d$status == "Up"], down = d$gene_id[d$status == "Down"],
                    d$gene_id[d$status %in% c("Up", "Down")]))
    universe <- enrichment_universe(rv$de, rv$universe_ids, rownames(rv$counts))
    results <- list(); notes <- character(0)
    withProgress(message = "Running enrichment", value = 0, {
      for (db in GALLERY_DATABASES) {
        incProgress(0.24, detail = GALLERY_DB_LABELS[[db]])
        # Cutoffs of 1 keep every tested term, which is what the ".all" tables
        # and the adjustable significance filter both need.
        r <- tryCatch(run_enrichment(genes, universe, db, 1, 1, symbol_map = enrichment_symbol_map(rv$de),
                                     organism = organism()), error = function(e) e)
        if (inherits(r, "error")) notes <- c(notes, sprintf("%s: %s", db, conditionMessage(r)))
        else if (is.null(r) || !nrow(as.data.frame(r))) notes <- c(notes, sprintf("%s: no terms returned.", db))
        else results[[db]] <- r
      }
    })
    tables <- tryCatch(enrich_results_to_gallery(results, enrichment_gene_view(rv$de)), error = function(e) e)
    if (inherits(tables, "error") || !length(tables)) {
      showNotification(paste(c("Enrichment produced no results.",
                               if (inherits(tables, "error")) conditionMessage(tables), notes),
                             collapse = " "), type = "error", duration = 20)
      return(invisible())
    }
    objs <- lapply(results[intersect(c("BP", "CC", "MF"), names(results))],
                   function(r) make_readable(r, "BP"))
    set_gallery(tables, "run",
                sprintf("%s DEGs (%d genes) from %s",
                        switch(input$x_gene_source %||% "all", up = "up-regulated",
                               down = "down-regulated", "all"),
                        length(genes), rv$source_label),
                notes, objs = objs)
    showNotification(sprintf("Enrichment complete for %s.", paste(names(tables), collapse = ", ")),
                     type = "message")
  }, ignoreInit = TRUE)

  # Fold changes used to colour genes come from the loaded DE table, matched
  # by gene symbol.
  gallery_fc <- reactive({
    if (is.null(rv$de)) return(NULL)
    gv <- enrichment_gene_view(rv$de)      # gene symbols, also for transcript-level tables
    fc <- stats::setNames(gv$log2FoldChange, gv$display_gene)
    fc[!is.na(fc) & !duplicated(names(fc))]
  })

  gallery_opts <- reactive({
    list(top_n = input$x_top_n %||% 10, rank_by = input$x_rank_by %||% "p.adjust",
         cutoff = input$x_cutoff %||% 0.05,
         palette = { v <- input$x_palette %||% "auto"; if (identical(v, "auto")) NULL else v },
         bar_colour = { v <- input$x_bar_colour; if (is.null(v) || !nzchar(v)) "#ADD8E6" else v },
         dot_x = input$x_dot_x %||% "score",
         show_values = isTRUE(input$x_show_values),
         digits = as.integer(input$x_digits %||% 2), angle = input$x_angle %||% 70,
         y_breaks = { t <- input$x_y_breaks %||% ""
                      if (nzchar(trimws(t))) unlist(strsplit(t, "[,; ]+")) else NULL },
         kegg_value = input$x_kegg_value %||% "count",
         sankey_ratio = input$x_sankey_ratio %||% 7, dot_ratio = input$x_dot_ratio %||% 5,
         max_genes = input$x_max_genes %||% 100, min_similarity = input$x_min_sim %||% 0.2,
         label_chars = input$x_label_chars %||% 60)
  })

  gallery_style <- reactive({
    st <- collect_style(input, "x")
    st$base_family <- input$x_font %||% "sans"
    st
  })

  # Database / pathway pickers follow whatever is loaded.
  observe({
    tabs <- rv$gal
    fig <- input$x_fig %||% "score_bar"
    avail <- intersect(GALLERY_DATABASES, names(tabs))
    if (identical(gallery_figure_scope(fig), "go1")) avail <- setdiff(avail, "KEGG")
    ch <- if (!length(avail)) c("Load results first" = "") else
      stats::setNames(avail, paste0(if_else_kegg(avail), GALLERY_DB_LABELS[avail]))
    keep <- isolate(input$x_db)
    updateSelectInput(session, "x_db", choices = ch,
                      selected = if (!is.null(keep) && keep %in% avail) keep else
                        if (fig == "sankey" && "KEGG" %in% avail) "KEGG" else avail[1])
  })
  if_else_kegg <- function(db) ifelse(db == "KEGG", "", "GO: ")

  observe({
    k <- rv$gal[["KEGG"]]
    sig <- if (is.null(k)) NULL else gallery_sig(k, input$x_rank_by %||% "p.adjust", input$x_cutoff %||% 0.05)
    if (!is.null(sig)) sig <- sig[.is_kegg_id(sig$ID), , drop = FALSE]
    ch <- if (is.null(sig) || !nrow(sig)) c("No significant KEGG pathways loaded" = "") else
      stats::setNames(sig$ID, sprintf("%s (%s)", sig$Description, sig$ID))
    updateSelectInput(session, "x_pathway", choices = ch)
  })

  # Export size follows the figure type (the user can still override it).
  observeEvent(input$x_fig, {
    sz <- gallery_figure_size(input$x_fig)
    updateNumericInput(session, "x_w", value = sz[1])
    updateNumericInput(session, "x_h", value = sz[2])
  })

  observeEvent(input$x_pv_render, {
    id <- input$x_pathway %||% ""
    if (!nzchar(id)) {
      showNotification("Choose a pathway first.", type = "warning"); return(invisible())
    }
    fc <- gallery_fc()
    lim <- suppressWarnings(as.numeric(input$x_pv_limit)); if (is.na(lim) || lim <= 0) lim <- NULL
    r <- withProgress(message = paste("Rendering", id, "with pathview..."), value = 0.5,
      tryCatch(render_pathview_map(id, fc, "SYMBOL", limit = lim),
               error = function(e) e))
    if (inherits(r, "error")) {
      rv$pv_png <- NULL; rv$pv_msg <- conditionMessage(r)
    } else {
      rv$pv_png <- r; rv$pv_msg <- NULL
    }
  }, ignoreInit = TRUE)

  gallery_obj <- reactive({
    fig <- input$x_fig %||% "score_bar"
    if (fig == "pathview") {
      if (is.null(rv$pv_png)) stop("Render a pathway map first.")
      return(png_draw_function(rv$pv_png))
    }
    if (is.null(rv$gal) || !length(rv$gal)) stop("No enrichment results are loaded yet.")
    build_gallery_figure(fig, rv$gal, input$x_db, gallery_opts(), gallery_style(),
                         collect_labels(input, "x"), gallery_fc(), rv$gal_objs)
  })

  output$gallery_status <- renderUI({
    tabs <- rv$gal
    if (is.null(tabs) || !length(tabs))
      return(div(class = "note-warn",
                 "No enrichment results loaded yet. In the sidebar, import your own enrichment ",
                 "tables or run GO + KEGG on the current DEGs. Example files in the import ",
                 "layout are under ", strong("Example files"), " in the import option."))
    rb <- input$x_rank_by %||% "p.adjust"; co <- input$x_cutoff %||% 0.05
    pills <- lapply(names(tabs), function(db) span(class = "status-pill pill-ok", sprintf(
      "%s: %s terms, %s significant", db, format(nrow(tabs[[db]]), big.mark = ","),
      format(nrow(gallery_sig(tabs[[db]], rb, co)), big.mark = ","))))
    tagList(
      div(style = "margin-bottom:6px", pills),
      div(style = "color:#66737f;font-size:.9em;margin-bottom:6px", paste("Source:", rv$gal_source)),
      if (identical(rv$gal_kind, "import")) div(style = "color:#66737f;font-size:.9em;margin-bottom:6px",
        "Genes are coloured by the log2 fold changes of the table on the Data tab, matched by symbol."),
      if (length(rv$gal_notes)) div(class = "note-ok", style = "font-size:.88em",
                                    lapply(rv$gal_notes, function(n) div(n))))
  })

  output$gallery_plot_ui <- renderUI({
    fig <- input$x_fig %||% "score_bar"
    if (fig == "pathview") {
      return(tagList(
        if (!HAS_PATHVIEW) div(class = "note-warn",
          strong("pathview is not installed. Run:"),
          tags$pre('BiocManager::install("pathview")')),
        if (!is.null(rv$pv_msg)) div(class = "note-err", rv$pv_msg),
        if (is.null(rv$pv_png)) div(class = "note-warn",
          "Choose a pathway in the sidebar and press ", strong("Render pathway map"),
          ". Genes are coloured by log2 fold change.")
        else div(style = "overflow:auto;border:1px solid #dfe5ea;border-radius:4px;padding:6px",
                 imageOutput("gallery_pathview", height = "auto"))))
    }
    h <- switch(fig, go3_h = 720, go3_v = 700, sankey = 760, cnet = 780, emap = 720,
                goplot = 860, kegg_class = 640, 600)
    plotOutput("gallery_plot", height = paste0(h, "px"))
  })

  output$gallery_plot <- renderPlot({
    obj <- tryCatch(gallery_obj(), error = function(e) e)
    validate(need(!inherits(obj, "error"),
                  if (inherits(obj, "error")) conditionMessage(obj) else ""))
    print(obj)
  })

  output$gallery_pathview <- renderImage({
    req(rv$pv_png)
    list(src = rv$pv_png, contentType = "image/png", alt = "KEGG pathway map",
         style = "max-width:100%;height:auto")
  }, deleteFile = FALSE)

  # The table under the figure: the database(s) the figure draws from.
  gallery_table_dbs <- reactive({
    avail <- intersect(GALLERY_DATABASES, names(rv$gal))
    switch(gallery_figure_scope(input$x_fig %||% "score_bar"),
           go3 = intersect(c("BP", "CC", "MF"), avail),
           kegg = intersect("KEGG", avail),
           intersect(input$x_db %||% "", avail))
  })

  gallery_tables_for_export <- function(scope) {
    tabs <- rv$gal
    if (is.null(tabs) || !length(tabs)) stop("No enrichment results are loaded.")
    out <- lapply(tabs, function(d) gallery_export_table(
      if (scope == "sig") gallery_sig(d, input$x_rank_by %||% "p.adjust", input$x_cutoff %||% 0.05) else d))
    names(out) <- ifelse(names(tabs) == "KEGG", "pathway", paste0("go.", names(tabs)))
    out
  }

  output$gallery_table_title <- renderText({
    dbs <- gallery_table_dbs()
    if (!length(dbs)) "Enriched terms" else
      paste("Enriched terms \u2014", paste(GALLERY_DB_LABELS[dbs], collapse = ", "))
  })

  output$gallery_table <- DT::renderDT({
    dbs <- gallery_table_dbs()
    validate(need(length(dbs) > 0, "No table to show for this figure yet."))
    parts <- lapply(dbs, function(db) {
      d <- rv$gal[[db]]
      if (identical(input$x_table_scope %||% "sig", "sig"))
        d <- gallery_sig(d, input$x_rank_by %||% "p.adjust", input$x_cutoff %||% 0.05)
      d <- gallery_export_table(d)
      if (length(dbs) > 1 && nrow(d)) d <- cbind(Ontology = db, d)
      d
    })
    df <- do.call(rbind, parts)
    validate(need(!is.null(df) && nrow(df) > 0, "No terms pass the current cutoff."))
    df <- format_enrich_numbers(df)
    df$geneID <- collapsible_gene_cell(df$geneID)
    gene_col <- which(colnames(df) == "geneID") - 1L
    DT::datatable(df, escape = FALSE, rownames = FALSE,
      colnames = sub("^geneID$", "Genes", enrich_table_colnames(df)),
      options = list(pageLength = 10, scrollX = TRUE, autoWidth = FALSE,
        columnDefs = list(
          list(targets = as.list(gene_col), className = "gene-col", orderable = FALSE),
          list(targets = "_all", className = "dt-head-nowrap"))))
  })

  # One workbook, one sheet per database (CSV of everything without openxlsx).
  gallery_table_handler <- function(scope) {
    has_xlsx <- requireNamespace("openxlsx", quietly = TRUE)
    downloadHandler(
      filename = function() stamp_name(paste0("enrichment_", scope), if (has_xlsx) "xlsx" else "csv"),
      content = function(f) {
        tabs <- tryCatch(gallery_tables_for_export(scope), error = function(e) e)
        if (inherits(tabs, "error")) {
          showNotification(conditionMessage(tabs), type = "error")
          tabs <- list(note = data.frame(message = conditionMessage(tabs)))
        }
        if (has_xlsx) openxlsx::write.xlsx(tabs, f)
        else utils::write.csv(do.call(rbind, Map(function(d, n) if (nrow(d)) cbind(Database = n, d),
                                                 tabs, names(tabs))), f, row.names = FALSE)
      })
  }
  output$x_dl_all <- gallery_table_handler("all")
  output$x_dl_sig <- gallery_table_handler("sig")

  output$x_dl_bundle <- downloadHandler(
    filename = function() stamp_name("enrichment_results", "zip"),
    contentType = "application/zip",
    content = function(f) {
      out <- tempfile("bundle_"); dir.create(out)
      on.exit(unlink(out, recursive = TRUE), add = TRUE)
      res <- withProgress(message = "Building the results bundle", value = 0, {
        tryCatch({
          write_gallery_bundle(
            rv$gal %||% list(), out, gallery_opts(), gallery_style(),
            fc = gallery_fc(), objects = rv$gal_objs,
            dpi = input$x_dpi %||% 300,
            pathview = isTRUE(input$x_bundle_pathview),
            pathview_fc = gallery_fc(), pathview_id_type = "SYMBOL",
            progress = function(msg) incProgress(0.12, detail = msg))
          zip_directory(out, f); NULL
        }, error = function(e) conditionMessage(e))
      })
      if (!is.null(res)) {
        showNotification(paste("The bundle could not be built:", res), type = "error", duration = 15)
        writeLines(paste("The bundle could not be built:", res), file.path(out, "ERROR.txt"))
        tryCatch(zip_directory(out, f), error = function(e) writeLines(res, f))
      }
    })

  # =========================================================================
  # Gene intelligence
  # =========================================================================

  gene_info_ids <- reactive(genes_for("g"))

  output$gene_link_table <- DT::renderDT({
    ids <- gene_info_ids()
    validate(need(length(ids) > 0, "Choose a gene set in the sidebar to see resources."))
    DT::datatable(
      gene_link_table(ids, rv$de, input$g_species %||% organism()$species, input$g_context %||% ""),
      escape = FALSE, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  observeEvent(input$g_build, {
    ids <- gene_info_ids()
    if (!length(ids)) {
      showNotification("Select at least one gene first.", type = "warning"); return(invisible())
    }
    sym <- rv$de$display_gene[match(ids, rv$de$gene_id)]
    summaries <- NULL; pubs <- NULL
    if (isTRUE(input$g_live)) {
      withProgress(message = "Fetching gene information...", value = 0.2, {
        lim <- min(length(ids), input$g_max_live %||% 25)
        summaries <- fetch_gene_summaries(utils::head(sym, lim), species = organism()$mygene %||% organism()$species)
        incProgress(0.4, detail = "publications")
        pubs <- list()
        for (s in utils::head(sym, lim)) {
          pubs[[s]] <- fetch_publications(s, input$g_context %||% "", input$g_n_pubs %||% 5)
        }
      })
      if (is_fetch_error(summaries)) {
        rv$gene_report_msg <- list(type = "warn", text = paste0(
          "Live lookups failed (", fetch_error_message(summaries),
          "). The report was built from database links only."))
      } else {
        rv$gene_report_msg <- list(type = "ok", text = "Live gene descriptions and publications retrieved.")
      }
    } else {
      rv$gene_report_msg <- list(type = "ok",
        text = "Report built with database links. Tick 'Fetch live...' to add descriptions and publications.")
    }
    d <- de_classified()
    rv$gene_report_html <- build_gene_report(
      ids, rv$de, summaries = summaries, publications = pubs,
      species = input$g_species %||% organism()$species,
      context_terms = input$g_context %||% "",
      meta_info = c(
        Source = rv$source_label,
        Thresholds = sprintf("adj. P <= %s, |log2FC| >= %s", input$padj_cutoff, input$lfc_cutoff),
        `Significant genes` = sprintf("%d up / %d down",
                                      sum(d$status == "Up"), sum(d$status == "Down"))),
      title = "Differentially expressed gene report")
  }, ignoreInit = TRUE)

  output$gene_report_note <- renderUI({
    m <- rv$gene_report_msg; if (is.null(m)) return(NULL)
    div(class = if (m$type == "ok") "note-ok" else "note-warn", m$text)
  })

  output$gene_report_preview <- renderUI({
    if (is.null(rv$gene_report_html))
      return(div(class = "note-warn", "Press ", strong("Build report"), " to generate the dossier."))
    tags$iframe(srcdoc = rv$gene_report_html, style = "width:100%;height:640px;border:1px solid #dfe5ea;border-radius:4px;")
  })

  output$g_dl_html <- downloadHandler(
    filename = function() stamp_name("gene_report", "html"),
    content = function(f) {
      html <- rv$gene_report_html
      if (is.null(html)) {
        ids <- gene_info_ids()
        validate(need(length(ids) > 0, "Select genes first."))
        html <- build_gene_report(ids, rv$de, species = input$g_species %||% organism()$species,
                                  context_terms = input$g_context %||% "")
      }
      writeLines(html, f)
    })

  output$g_dl_csv <- downloadHandler(
    filename = function() stamp_name("gene_resources", "csv"),
    content = function(f) {
      ids <- gene_info_ids()
      validate(need(length(ids) > 0, "Select genes first."))
      tb <- gene_link_table(ids, rv$de, input$g_species %||% organism()$species, input$g_context %||% "")
      tb$Resources <- gsub("<[^>]+>", " ", tb$Resources)
      tb$Literature <- gsub("<[^>]+>", " ", tb$Literature)
      utils::write.csv(tb, f, row.names = FALSE)
    })

  # =========================================================================
  # Figure downloads
  # =========================================================================

  #' Build one TIFF + one PDF download handler for a figure.
  #'
  #' Both formats go through exactly the same code path, and any failure is
  #' reported to the user (notification plus an explanatory page inside the
  #' downloaded file) rather than leaving the button looking inert.
  make_figure_handler <- function(obj_fun, stem, fmt, w_in, h_in, dpi_in) {
    downloadHandler(
      filename    = function() stamp_name(stem, fmt),
      contentType = switch(fmt, pdf = "application/pdf", png = "image/png", "image/tiff"),
      content = function(f) {
        obj <- tryCatch(obj_fun(), error = function(e) e)
        if (inherits(obj, "error")) {
          showNotification(paste0("Could not build the figure: ", conditionMessage(obj)),
                           type = "error", duration = 15)
          export_plot_safe(NULL, f, fmt, w_in(), h_in(), dpi_in())
          return(invisible())
        }
        err <- export_plot_safe(obj, f, fmt, w_in(), h_in(), dpi_in())
        if (!is.null(err)) {
          showNotification(paste0(toupper(fmt), " export failed: ", err),
                           type = "error", duration = 15)
        }
      })
  }

  register_export <- function(p, obj_fun, stem,
                              h_scale = 1) {
    w   <- function() input[[paste0(p, "_w")]]
    h   <- function() {
      v <- suppressWarnings(as.numeric(input[[paste0(p, "_h")]]))
      if (length(v) != 1 || is.na(v)) v <- 6
      v * h_scale
    }
    dpi <- function() input[[paste0(p, "_dpi")]]
    output[[paste0(p, "_dl_tiff")]] <- make_figure_handler(obj_fun, stem, "tiff", w, h, dpi)
    output[[paste0(p, "_dl_pdf")]]  <- make_figure_handler(obj_fun, stem, "pdf",  w, h, dpi)
    output[[paste0(p, "_dl_png")]]  <- make_figure_handler(obj_fun, stem, "png",  w, h, dpi)
  }
  register_export("v", volcano_obj, "volcano_plot")
  register_export("m", ma_obj, "ma_plot")
  register_export("s", scatter_obj, "expression_scatter")
  register_export("h", function() heatmap_obj()$plot, "heatmap")
  register_export("e", enrich_plot_obj, "enrichment")
  register_export("x", gallery_obj, "enrichment_gallery")

  # Gene explorer: every gene in the matrix (defined here so it can use
  # register_export, get_expression and matrix_choices).
  gene_explorer_server(input, output, session, rv, de_classified, get_expression,
                       matrix_choices, register_export)

  # PCA exports both panels: the scree plot is stacked underneath so one file
  # carries the full QC figure, hence the taller page.
  register_export("p", function() combine_stacked(pca_obj(), scree_obj()),
                  "pca_scree", h_scale = 1.6)

  # Load the bundled example at start-up so the app is never empty. Failure is
  # non-fatal: the Data tab then simply invites the user to upload files.
  tryCatch(load_example(), error = function(e) {
    rv$messages <- list(type = "warn",
      text = paste("The bundled example could not be loaded:", conditionMessage(e),
                   "-- upload your own files on this tab."))
  })
}
