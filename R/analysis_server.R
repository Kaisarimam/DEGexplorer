# =============================================================================
# analysis_server.R -- server logic for the Analysis tab
#
# Called from app_server(). A finished run is handed to the Explorer by filling
# the same reactive values an uploaded results file fills, so every other tab
# works on it unchanged.
# =============================================================================

analysis_server <- function(input, output, session, rv, reset_enrichment) {

  av <- reactiveValues(
    counts = NULL, counts_ann = NULL, annot = NULL, source = "",
    len_file = NULL, len_gtf = NULL, len_label = "",
    sheet = NULL, types = list(), contrasts = empty_contrasts(),
    result = NULL, explored = NULL,
    # transcript level
    tx_counts = NULL,            # transcript -> gene columns that came with the counts
    tx_file = NULL,              # an uploaded tx2gene table
    tx_gtf = NULL, gtf_path = NULL,  # transcripts of a GTF (parsed on demand)
    overdisp = NULL, quant_note = NULL, ngvec = NULL
  )
  level <- reactive(input$a_level %||% "gene")
  tx_level <- reactive(level() %in% c("transcript", "dtu"))

  # DTU has its own methods; the engine list follows the analysis type.
  observeEvent(level(), {
    dtu <- identical(level(), "dtu")
    ch <- if (dtu) DTU_ENGINES else DE_ENGINES
    cur <- isolate(input$a_engine)
    updateSelectInput(session, "a_engine", choices = ch,
                      selected = if (!is.null(cur) && cur %in% ch) cur else if (dtu) "edger_dtu" else "deseq2")
  }, ignoreInit = TRUE)

  variables <- reactive(names(av$types))
  factor_vars <- reactive({
    v <- variables(); v[vapply(v, function(x) identical(av$types[[x]], "factor"), logical(1))]
  })
  included <- reactive({
    req(av$sheet)
    intersect(av$sheet$sample, input$a_include %||% character(0))
  })

  # ---- 1. counts and annotation ---------------------------------------------

  set_counts <- function(m, label, ann = NULL) {
    if (is.null(m) || !ncol(m)) stop("No numeric sample columns found in that file.")
    if (any(m < 0, na.rm = TRUE)) stop("The matrix has negative values; raw counts are needed, not normalised data.")
    av$counts <- m; av$counts_ann <- ann; av$source <- label
    av$len_file <- NULL   # a lengths file belongs to its own matrix; a GTF still applies
    av$overdisp <- NULL; av$quant_note <- NULL; av$ngvec <- NULL; av$tx_counts <- NULL
    av$sheet <- design_scaffold(colnames(m)); av$types <- list()
    av$contrasts <- empty_contrasts(); av$result <- NULL; av$explored <- NULL
  }

  load_counts <- function(path, label, gene_level = FALSE) {
    r <- read_matrix_auto(path)
    is_tx <- !gene_level && (tx_level() || looks_like_transcripts(rownames(r$matrix)))
    if (is_tx) {
      t <- read_transcript_matrix(path)
      set_counts(t$matrix, label)
      av$tx_counts <- t$tx2gene
      if (!tx_level()) {
        updateRadioButtons(session, "a_level", selected = "transcript")
        showNotification("The ids look like transcripts, so the feature level was set to Transcripts.",
                         type = "message", duration = 8)
      }
    } else set_counts(r$matrix, label, r$annotation)
  }

  # Salmon / kallisto / RSEM folders: counts, effective lengths, overdispersion.
  observeEvent(input$a_quant_read, {
    root <- NULL
    if (!is.null(input$a_quant_zip)) {
      root <- file.path(tempdir(), paste0("quant_", as.integer(Sys.time())))
      dir.create(root); utils::unzip(input$a_quant_zip$datapath, exdir = root)
    } else if (nzchar(trimws(input$a_quant_path %||% ""))) {
      root <- path.expand(trimws(input$a_quant_path))
    }
    if (is.null(root)) { showNotification("Give a folder path or a .zip of the quantification folders.", type = "warning"); return() }
    q <- withProgress(message = "Reading quantifications", value = 0.3,
                      tryCatch(read_quant_dirs(root), error = function(e) e))
    if (inherits(q, "error")) {
      showNotification(paste("Could not read the quantifications:", conditionMessage(q)), type = "error", duration = 20)
      return()
    }
    set_counts(q$counts, sprintf("%s quantifications (%d samples)", tools::toTitleCase(q$tool), q$n_samples))
    av$len_file <- q$lengths; av$overdisp <- q$overdispersion; av$quant_note <- q$note
    av$tx_counts <- q$tx2gene
    if (!tx_level()) updateRadioButtons(session, "a_level", selected = "transcript")
    showNotification(paste0(sprintf("%s transcripts \u00d7 %d samples read, with effective lengths",
                                    format(nrow(q$counts), big.mark = ","), q$n_samples),
                            if (!is.null(q$overdispersion)) " and quantification overdispersion." else ".",
                            if (!is.null(q$note)) paste0(" ", q$note) else ""),
                     type = if (is.null(q$note)) "message" else "warning", duration = 15)
  })

  observeEvent(input$a_tx2gene, {
    t <- tryCatch(withProgress(message = "Reading the transcript-to-gene map", value = 0.5,
                               read_tx2gene(input$a_tx2gene$datapath)), error = function(e) e)
    if (inherits(t, "error")) {
      showNotification(paste("Could not read the transcript-to-gene map:", conditionMessage(t)), type = "error", duration = 12)
      return()
    }
    av$tx_file <- t
    n <- if (is.null(av$counts)) NA else sum(rownames(av$counts) %in% t$transcript_id |
                                              sub("\\.[0-9]+$", "", rownames(av$counts)) %in% sub("\\.[0-9]+$", "", t$transcript_id))
    showNotification(sprintf("Transcript-to-gene map: %s transcripts%s.", format(nrow(t), big.mark = ","),
                             if (is.na(n)) "" else sprintf(", %s of the matrix's", format(n, big.mark = ","))),
                     type = "message")
  })

  observeEvent(input$a_ngvec, {
    req(av$counts)
    v <- tryCatch(read_ngvector(input$a_ngvec$datapath, rownames(av$counts)), error = function(e) e)
    if (inherits(v, "error")) {
      showNotification(paste("Could not read the Ng vector:", conditionMessage(v)), type = "error", duration = 15)
      return()
    }
    av$ngvec <- v
    showNotification(sprintf("Ng vector: %s.", paste(sprintf("Ng=%s: %s", names(table(v)), table(v)), collapse = "; ")),
                     type = "message")
  })

  #' Everything a transcript-level run needs, merged from the available sources.
  transcript_inputs <- function() {
    if (is.null(av$counts)) return(NULL)
    if (is.null(av$tx_file) && is.null(av$tx_counts) && is.null(av$tx_gtf) && !is.null(av$gtf_path)) {
      av$tx_gtf <- withProgress(message = "Reading transcripts from the GTF", value = 0.5,
                                tryCatch(gtf_transcript_table(av$gtf_path), error = function(e) NULL))
    }
    list(tx2gene = merge_tx2gene(rownames(av$counts), av$tx_file, av$tx_counts, av$tx_gtf),
         overdispersion = av$overdisp, ngvec = av$ngvec)
  }

  observeEvent(input$a_counts, {
    r <- tryCatch({ load_counts(input$a_counts$datapath, input$a_counts$name); NULL },
                  error = function(e) conditionMessage(e))
    if (!is.null(r)) showNotification(paste("Could not read the count matrix:", r), type = "error", duration = 12)
  })

  observeEvent(input$a_annot, {
    a <- tryCatch(read_gene_annotation(input$a_annot$datapath, input$a_annot$name),
                  error = function(e) e)
    if (inherits(a, "error")) {
      showNotification(paste("Could not read the annotation:", conditionMessage(a)), type = "error", duration = 12)
      return()
    }
    av$annot <- a
    if (is_gtf_file(input$a_annot$datapath, input$a_annot$name)) {
      av$gtf_path <- input$a_annot$datapath; av$tx_gtf <- NULL
      l <- withProgress(message = "Computing gene lengths from the GTF", value = 0.5,
                        tryCatch(gtf_gene_lengths(input$a_annot$datapath), error = function(e) NULL))
      if (!is.null(l)) { av$len_gtf <- l; av$len_label <- paste("exon union from", input$a_annot$name) }
    }
  })

  observeEvent(input$a_lengths, {
    m <- tryCatch(read_matrix_auto(input$a_lengths$datapath)$matrix, error = function(e) e)
    if (inherits(m, "error") || is.null(m)) {
      showNotification(paste("Could not read the gene lengths:",
                             if (inherits(m, "error")) conditionMessage(m) else "no numeric columns"),
                       type = "error", duration = 12)
      return()
    }
    av$len_file <- m
  })

  annotation_fetch_server(input, output, session, "a",
    ids = reactive(rownames(av$counts)),
    on_loaded = function(ann, len, label, tx = NULL, organism_key = NULL) {
      attr(ann, "source") <- label
      av$annot <- ann; av$len_gtf <- len; av$len_label <- paste("exon union from", label)
      if (!is.null(tx)) av$tx_gtf <- tx
      if (!is.null(organism_key)) updateSelectInput(session, "organism", selected = organism_key)
      showNotification(sprintf("%s: symbols for %s genes and lengths for %s.", label,
                               format(sum(!is.na(ann$gene_name)), big.mark = ","),
                               format(length(len), big.mark = ",")), type = "message", duration = 10)
    })

  #' Gene lengths for the Explorer: an uploaded lengths file wins over a GTF.
  lengths_for_explorer <- function() {
    if (!is.null(av$len_file)) return(av$len_file)
    if (!is.null(av$len_gtf) && !is.null(av$counts)) {
      l <- match_gene_lengths(av$len_gtf, rownames(av$counts))
      if (length(l)) return(l)
    }
    NULL
  }

  observeEvent(input$a_use_example, {
    r <- tryCatch({
      load_counts(app_file("data", "example_rsem_gene_counts.tsv.gz"),
                  "bundled example (6 samples, RSEM counts)", gene_level = TRUE)
      updateRadioButtons(session, "a_level", selected = "gene")
      av$len_file <- read_matrix_auto(app_file("data", "example_rsem_gene_lengths.tsv.gz"))$matrix
      meta <- read_sample_metadata(app_file("data", "example_samplesheet.csv"))
      sh <- av$sheet
      for (v in setdiff(colnames(meta), "sample")) {
        sh[[v]] <- as.character(meta[[v]][match(sh$sample, meta$sample)])
        sh[[v]][is.na(sh[[v]])] <- ""
        av$types[[v]] <- "factor"
      }
      av$sheet <- sh
      if ("condition" %in% colnames(sh)) {
        lv <- sort(unique(sh$condition[nzchar(sh$condition)]))
        ref <- if ("control" %in% lv) "control" else lv[1]
        tg <- setdiff(lv, ref)
        if (length(tg)) av$contrasts <- data.frame(
          id = contrast_id("condition", tg, ref), variable = "condition",
          reference = ref, target = tg, stringsAsFactors = FALSE)
      }
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(r)) showNotification(paste("Could not load the example:", r), type = "error", duration = 12)
  })

  gene_annotation <- reactive({
    if (!is.null(av$annot)) return(av$annot)
    a <- av$counts_ann
    if (is.null(a) || !all(c("gene_id", "gene_name") %in% colnames(a))) return(NULL)
    a
  })

  output$a_counts_note <- renderUI({
    if (is.null(av$counts)) return(div(class = "note-warn",
      "Upload a raw count matrix, or use the bundled example."))
    if (tx_level()) {
      t2g <- merge_tx2gene(rownames(av$counts), av$tx_file, av$tx_counts, av$tx_gtf)
      n_map <- sum(!is.na(t2g$gene_id))
      return(div(class = if (n_map == 0) "note-warn" else "note-ok", sprintf(
        "%s: %s transcripts \u00d7 %d samples; transcript-to-gene map for %s%s; %s; %s.",
        av$source, format(nrow(av$counts), big.mark = ","), ncol(av$counts),
        format(n_map, big.mark = ","),
        if (n_map == 0) (if (!is.null(av$gtf_path)) " (read from the GTF at run time)" else " (load a GTF or a tx2gene table)") else "",
        if (is.null(av$len_file)) "no effective lengths" else "effective lengths loaded",
        if (!is.null(av$overdisp)) sprintf("quantification overdispersion loaded (median %.2f)", stats::median(av$overdisp))
        else av$quant_note %||% "no quantification resamples (read the Salmon / kallisto folders to add them)")))
    }
    ann <- gene_annotation()
    len <- if (!is.null(av$len_file)) "gene lengths from the uploaded file" else
           if (!is.null(av$len_gtf)) paste("gene lengths:", av$len_label) else
           "no gene lengths (TPM / FPKM unavailable; add a GTF or a lengths file)"
    div(class = "note-ok", sprintf("%s: %s genes \u00d7 %d samples%s; %s.",
      av$source, format(nrow(av$counts), big.mark = ","), ncol(av$counts),
      if (is.null(ann)) "; no gene symbols yet (add or download a GTF)" else
        sprintf("; gene symbols from %s", if (!is.null(av$annot)) attr(av$annot, "source") %||% "the annotation file"
                                           else "the matrix's gene_name column"),
      len))
  })

  output$a_tx_sop <- renderUI({
    e <- input$a_engine %||% "deseq2"
    txt <- switch(e,
      deseq2 = paste("DESeq2 on estimated transcript counts with each sample's effective transcript lengths",
                     "as offsets (tximport txOut): tick Gene-length offsets below. Fold-change shrinkage helps",
                     "with the many low-count transcripts."),
      limma  = paste("limma-voom on counts divided by their quantification overdispersion (Baldoni et al. 2024);",
                     "the filter works on the divided counts."),
      edger  = paste("edgeR quasi-likelihood on counts divided by their quantification overdispersion",
                     "(Baldoni et al. 2024, edgeR 4); the filter works on the divided counts."),
      ebseq  = paste("EBSeq isoform mode: transcripts are grouped by their Ng (isoforms per gene, or",
                     "rsem-generate-ngvector), each group with its own variance."))
    helpText(style = "font-size:.82em", txt)
  })
  output$a_dtu_sop <- renderUI({
    e <- input$a_engine %||% "edger_dtu"
    txt <- switch(e,
      edger_dtu = paste("edgeR (Baldoni et al. 2025, NAR): counts divided by their quantification overdispersion,",
                        "filterByExpr, glmQLFit and diffSplice; gene-level Simes p-values. Controlled the FDR in every",
                        "benchmark scenario of that paper and was the most powerful. Recommended."),
      limma_dtu = paste("limma (Baldoni et al. 2025): divided counts, voomLmFit and diffSplice. As good as edgeR, and",
                        "suits complex designs."),
      dexseq = paste("DEXSeq as in Love, Soneson & Patro 2018 (F1000Research): scaledTPM counts, DRIMSeq's filter,",
                     "~sample + exon + condition:exon, per-gene q-values. Each comparison is fitted on its two groups.",
                     "Can be liberal with three samples per group."),
      drimseq = paste("DRIMSeq (Dirichlet-multinomial) as in Love et al. 2018, with its post-hoc proportion-SD filter.",
                      "Each comparison is fitted on its two groups. Liberal in recent benchmarks; slow."),
      saturn = paste("satuRn (Gilis et al. 2022): quasi-binomial models, empirical p-values. Scales to large",
                     "datasets; conservative with fewer than about ten samples per group."))
    helpText(style = "font-size:.82em", txt,
             "All need a transcript-to-gene map. Results are per transcript, with a gene-level FDR (gene_padj).")
  })
  output$a_dtu_filter_note <- renderUI({
    if (!identical(input$a_dtu_filter %||% "auto", "auto")) return(NULL)
    e <- input$a_engine %||% "edger_dtu"
    helpText(style = "font-size:.82em", switch(e,
      dexseq = , drimseq = "Recommended here: DRIMSeq dmFilter (values below).",
      saturn = "Recommended here: filterByExpr with satuRn's settings, single-transcript genes removed.",
      "Recommended here: filterByExpr on the design (defaults)."))
  })
  output$a_recommended_help <- renderUI({
    helpText(style = "font-size:.82em", if (identical(level(), "dtu"))
      "edgeR diffSplice on divided counts, filterByExpr, Simes gene-level test, stageR." else
      "DESeq2; at least 10 counts in 1 sample; ashr-shrunken fold changes; gene-length offsets when sample-specific lengths are loaded.")
  })
  output$a_stager_note <- renderUI({
    if (requireNamespace("stageR", quietly = TRUE)) return(NULL)
    div(class = "note-warn", style = "font-size:.82em", "Needs stageR:", tags$code('BiocManager::install("stageR")'))
  })

  output$a_engine_note <- renderUI({
    miss <- de_engine_missing(input$a_engine %||% "deseq2")
    if (is.null(miss)) return(NULL)
    div(class = "note-warn", strong("Not installed. Run:"), tags$pre(miss))
  })

  # ---- 2. sample sheet ---------------------------------------------------------

  output$a_include_ui <- renderUI({
    req(av$counts)   # only a new matrix resets the selection
    s <- isolate(av$sheet$sample)
    checkboxGroupInput("a_include", sprintf("Samples to include (%d in the matrix)", length(s)),
                       choices = s, selected = s, inline = TRUE)
  })

  sheet_editor_server(
    input, output, session, "a",
    get_sheet = function() av$sheet, set_sheet = function(x) av$sheet <- x,
    get_types = function() av$types, set_types = function(x) av$types <- x,
    on_remove = function(v) av$contrasts <- av$contrasts[av$contrasts$variable != v, , drop = FALSE])

  observe({
    f <- factor_vars()
    sel <- isolate(input$a_main)
    updateSelectInput(session, "a_main", choices = f,
                      selected = if (!is.null(sel) && sel %in% f) sel else f[1])
  })

  pending_covars <- reactiveVal(NULL)   # covariates named by a loaded contrasts file
  observe({
    others <- setdiff(variables(), input$a_main %||% "")
    pc <- pending_covars()
    sel <- intersect(if (!is.null(pc)) pc else isolate(input$a_covars) %||% character(0), others)
    updateCheckboxGroupInput(session, "a_covars", choices = others, selected = sel)
  })

  # ---- design files: load and save (CSV) -------------------------------

  # A loaded contrasts file's covariates hold until the browser reports the new ticks.
  observeEvent(input$a_covars, pending_covars(NULL), ignoreInit = TRUE, ignoreNULL = FALSE)

  observeEvent(input$a_load_sheet, {
    if (is.null(av$sheet)) { showNotification("Load a count matrix first.", type = "warning"); return() }
    r <- tryCatch(read_design_samplesheet(input$a_load_sheet$datapath, av$sheet$sample),
                  error = function(e) e)
    if (inherits(r, "error")) {
      showNotification(paste("Could not read the sample sheet:", conditionMessage(r)), type = "error", duration = 12)
      return()
    }
    sh <- av$sheet; ty <- av$types
    for (v in setdiff(colnames(r$sheet), "sample")) { sh[[v]] <- r$sheet[[v]]; ty[[v]] <- r$types[[v]] }
    av$sheet <- sh; av$types <- ty
    showNotification(sprintf("Sample sheet loaded: %d of %d samples matched; variables %s.%s",
      r$matched, nrow(sh), paste(setdiff(colnames(r$sheet), "sample"), collapse = ", "),
      if (length(r$unmatched)) paste0(" Not in the file: ", paste(utils::head(r$unmatched, 5), collapse = ", "),
                                     if (length(r$unmatched) > 5) ", ..." else "", ".") else ""),
      type = if (length(r$unmatched)) "warning" else "message", duration = 12)
  })

  observeEvent(input$a_load_contrasts, {
    r <- tryCatch(read_design_contrasts(input$a_load_contrasts$datapath), error = function(e) e)
    if (inherits(r, "error")) {
      showNotification(paste("Could not read the comparisons:", conditionMessage(r)), type = "error", duration = 12)
      return()
    }
    vars <- unique(r$contrasts$variable)
    miss <- setdiff(c(vars, r$blocking), variables())
    if (length(miss)) {
      showNotification(sprintf("The comparisons use %s, which the sample sheet does not have. Load or add %s first.",
                               paste(miss, collapse = ", "), if (length(miss) > 1) "these variables" else "this variable"),
                       type = "error", duration = 15)
      return()
    }
    if (length(vars) > 1) showNotification(sprintf(
      "The file compares several variables (%s); the app tests one main variable at a time, so only '%s' comparisons were loaded.",
      paste(vars, collapse = ", "), vars[1]), type = "warning", duration = 15)
    av$contrasts <- r$contrasts[r$contrasts$variable == vars[1], , drop = FALSE]
    updateSelectInput(session, "a_main", selected = vars[1])
    pending_covars(setdiff(r$blocking, vars[1]))
    showNotification(sprintf("%d comparison(s) loaded%s.", nrow(av$contrasts),
                             if (length(r$blocking)) paste0("; covariates ticked: ", paste(r$blocking, collapse = ", ")) else ""),
                     type = "message")
  })

  output$a_dl_sheet <- downloadHandler(
    filename = function() "samplesheet.csv",
    content = function(file) {
      req(av$sheet)
      sh <- av$sheet[av$sheet$sample %in% included(), , drop = FALSE]
      utils::write.csv(sh, file, row.names = FALSE)
    })

  output$a_dl_contrasts <- downloadHandler(
    filename = function() "contrasts.csv",
    content = function(file) {
      utils::write.csv(de_contrasts_table(av$contrasts, setdiff(input$a_covars %||% character(0), input$a_main)),
                       file, row.names = FALSE)
    })

  # ---- 3. model and comparisons -----------------------------------------------

  main_levels <- reactive({
    m <- input$a_main
    req(av$sheet, m, m %in% colnames(av$sheet))
    x <- trimws(av$sheet[[m]][av$sheet$sample %in% included()])
    sort(unique(x[!is.na(x) & nzchar(x)]))
  })

  observe({
    lv <- main_levels()
    if (!length(lv)) {
      updateSelectInput(session, "a_ref", choices = character(0))
      updateSelectInput(session, "a_target", choices = character(0))
      return()
    }
    r0 <- isolate(input$a_ref); t0 <- isolate(input$a_target)
    ref <- if (!is.null(r0) && r0 %in% lv) r0 else lv[1]
    tg  <- if (!is.null(t0) && t0 %in% lv && t0 != ref) t0 else c(setdiff(lv, ref), ref)[1]
    updateSelectInput(session, "a_ref", choices = lv, selected = ref)
    updateSelectInput(session, "a_target", choices = lv, selected = tg)
  })

  output$a_formula <- renderUI({
    req(input$a_main)
    f <- paste(deparse(design_formula(input$a_main, setdiff(input$a_covars %||% character(0), input$a_main))),
               collapse = "")
    div(style = "margin-top:6px", "Model: ", tags$code(f))
  })

  add_contrasts <- function(targets, ref) {
    m <- input$a_main
    new <- data.frame(id = contrast_id(m, targets, ref), variable = m, reference = ref,
                      target = targets, stringsAsFactors = FALSE)
    new <- new[new$target != new$reference, , drop = FALSE]
    keep <- !(new$id %in% av$contrasts$id)
    if (!any(keep)) { showNotification("That comparison is already in the list.", type = "warning"); return() }
    av$contrasts <- rbind(av$contrasts, new[keep, , drop = FALSE])
  }

  observeEvent(input$a_add_cmp, {
    req(input$a_main, input$a_ref, input$a_target)
    if (identical(input$a_ref, input$a_target)) {
      showNotification("Choose two different levels.", type = "warning"); return()
    }
    nm <- trimws(input$a_cmp_name %||% "")
    if (nzchar(nm)) {
      id <- make.names(nm)
      if (id %in% av$contrasts$id) { showNotification("That name is taken.", type = "warning"); return() }
      av$contrasts <- rbind(av$contrasts, data.frame(id = id, variable = input$a_main,
        reference = input$a_ref, target = input$a_target, stringsAsFactors = FALSE))
      updateTextInput(session, "a_cmp_name", value = "")
    } else add_contrasts(input$a_target, input$a_ref)
  })

  observeEvent(input$a_all_vs_ref, {
    req(input$a_main, input$a_ref)
    add_contrasts(setdiff(main_levels(), input$a_ref), input$a_ref)
  })

  observeEvent(input$a_clear_cmp, av$contrasts <- empty_contrasts())

  observeEvent(input$a_remove_cmp, {
    rows <- input$a_cmp_table_rows_selected
    if (!length(rows)) { showNotification("Select comparisons in the table first.", type = "warning"); return() }
    av$contrasts <- av$contrasts[-rows, , drop = FALSE]
  })

  output$a_cmp_table <- DT::renderDT({
    df <- av$contrasts
    colnames(df) <- c("Name", "Variable", "Reference", "Target")
    DT::datatable(df, rownames = FALSE, selection = "multiple",
                  options = list(dom = "t", pageLength = 100, ordering = FALSE,
                                 language = list(emptyTable = "No comparisons yet.")))
  })

  checks <- reactive({
    req(av$sheet)
    ck <- validate_design(av$sheet, included(), input$a_main,
                          input$a_covars %||% character(0), av$types, av$contrasts)
    ck$errors <- c(ck$errors, engine_design_errors(input$a_engine %||% "deseq2",
                                                   input$a_covars %||% character(0), av$types,
                                                   analysis_options()))
    if (identical(level(), "dtu")) {
      e <- input$a_engine %||% "edger_dtu"
      if (!is_dtu_engine(e)) ck$errors <- c(ck$errors, "Choose a transcript-usage method.")
      if (!is.null(av$counts) && is.null(av$gtf_path) &&
          all(is.na(merge_tx2gene(rownames(av$counts), av$tx_file, av$tx_counts, av$tx_gtf)$gene_id)))
        ck$errors <- c(ck$errors, "Transcript usage needs a transcript-to-gene map: load a GTF or a tx2gene table.")
      if (is_dtu_engine(e)) {
        len <- lengths_for_explorer()
        rc <- dtu_resolve_counts(e, input$a_dtu_counts %||% "auto", !is.null(dim(len)), !is.null(av$overdisp))
        if (!is.null(rc$note)) ck$warnings <- c(ck$warnings, rc$note)
      }
    }
    if (identical(level(), "transcript")) {
      e <- input$a_engine %||% "deseq2"
      if (e == "ebseq" && identical(input$a_ebseq_ng, "file") && is.null(av$ngvec))
        ck$errors <- c(ck$errors, "EBSeq: choose the Ng vector file, or use 'Isoforms per gene'.")
      if (e %in% c("edger", "limma") && isTRUE(input$a_divide_od) && is.null(av$overdisp))
        ck$warnings <- c(ck$warnings, paste(
          "No quantification resamples are loaded, so transcript counts will be used without dividing out",
          "read-to-transcript ambiguity. Read the Salmon / kallisto folders (quantified with Gibbs or bootstrap samples) to add them."))
      if (!is.null(av$counts) && is.null(av$tx_file) && is.null(av$tx_counts) && is.null(av$tx_gtf) && is.null(av$gtf_path))
        ck$warnings <- c(ck$warnings, "No transcript-to-gene map: load a GTF or a tx2gene table for gene symbols, EBSeq isoform groups and stage-wise testing.")
    }
    ck
  })

  analysis_options <- reactive({
    o <- de_default_options()
    o$feature_level <- level()
    o$divide_overdispersion <- isTRUE(input$a_divide_od %||% TRUE)
    o$ebseq_ng <- input$a_ebseq_ng %||% o$ebseq_ng
    o$ebseq_trunthre <- input$a_ebseq_trunthre %||% o$ebseq_trunthre
    o$stager <- isTRUE(input$a_stager %||% TRUE)
    o$stager_alpha <- input$a_stager_alpha %||% o$stager_alpha
    o$dtu_counts <- input$a_dtu_counts %||% o$dtu_counts
    o$dtu_filter <- input$a_dtu_filter %||% o$dtu_filter
    o$dtu_min_gene_expr <- input$a_dtu_min_gene %||% o$dtu_min_gene_expr
    o$dtu_min_feature_expr <- input$a_dtu_min_tx %||% o$dtu_min_feature_expr
    o$dtu_min_feature_prop <- input$a_dtu_min_prop %||% o$dtu_min_feature_prop
    o$dtu_gene_test <- input$a_dtu_gene_test %||% o$dtu_gene_test
    o$dtu_robust <- !isFALSE(input$a_dtu_robust)
    o$drimseq_sd_filter <- input$a_drimseq_sd %||% o$drimseq_sd_filter
    o$saturn_pvalues <- input$a_saturn_p %||% o$saturn_pvalues
    o$min_abundance <- input$a_min_abundance %||% o$min_abundance
    o$use_cpm <- isTRUE(input$a_use_cpm)
    o$filter_mode <- input$a_filter_mode %||% o$filter_mode
    o$min_samples <- input$a_min_samples %||% o$min_samples
    o$min_proportion <- input$a_min_proportion %||% o$min_proportion
    o$fit_type <- input$a_fit_type %||% o$fit_type
    o$alpha <- input$a_alpha %||% o$alpha
    o$lfc_threshold <- input$a_lfc_threshold %||% o$lfc_threshold
    o$independent_filtering <- isTRUE(input$a_indep)
    o$cooks_cutoff <- isTRUE(input$a_cooks)
    o$transform <- input$a_transform %||% o$transform
    o$blind <- isTRUE(input$a_blind)
    o$norm_method <- input$a_norm_method %||% o$norm_method
    o$limma_voom <- input$a_limma_voom %||% o$limma_voom
    o$limma_trend <- isTRUE(input$a_limma_trend)
    o$limma_robust <- isTRUE(input$a_limma_robust)
    o$limma_block <- input$a_limma_block %||% "(none)"
    o$limma_lfc <- input$a_limma_lfc %||% 0
    o$edger_test <- input$a_edger_test %||% o$edger_test
    o$edger_robust <- !isFALSE(input$a_edger_robust)
    o$edger_lfc <- input$a_edger_lfc %||% 0
    o$ebseq_maxround <- input$a_ebseq_maxround %||% o$ebseq_maxround
    o$deseq2_test <- input$a_deseq2_test %||% o$deseq2_test
    o$shrink <- input$a_shrink %||% o$shrink
    o$length_offset <- isTRUE(input$a_length_offset)
    o
  })

  # One button sets every control to the recommended DESeq2 settings.
  observeEvent(input$a_recommended, {
    if (identical(level(), "dtu")) {
      updateSelectInput(session, "a_engine", selected = "edger_dtu")
      updateSelectInput(session, "a_dtu_counts", selected = "auto")
      updateSelectInput(session, "a_dtu_filter", selected = "auto")
      updateRadioButtons(session, "a_dtu_gene_test", selected = "simes")
      updateCheckboxInput(session, "a_dtu_robust", value = TRUE)
      updateCheckboxInput(session, "a_stager", value = TRUE)
      updateNumericInput(session, "a_stager_alpha", value = 0.05)
      showNotification(paste0("Recommended transcript-usage settings applied (edgeR diffSplice, Baldoni et al. 2025)",
                              if (is.null(av$overdisp)) "; read the Salmon / kallisto folders with resamples to divide the counts" else "",
                              "."), type = "message", duration = 10)
      return()
    }
    o <- recommended_options()
    updateSelectInput(session, "a_engine", selected = "deseq2")
    updateNumericInput(session, "a_min_abundance", value = o$min_abundance)
    updateCheckboxInput(session, "a_use_cpm", value = o$use_cpm)
    updateRadioButtons(session, "a_filter_mode", selected = o$filter_mode)
    updateNumericInput(session, "a_min_samples", value = o$min_samples)
    updateSelectInput(session, "a_fit_type", selected = o$fit_type)
    updateNumericInput(session, "a_alpha", value = o$alpha)
    updateNumericInput(session, "a_lfc_threshold", value = o$lfc_threshold)
    updateCheckboxInput(session, "a_indep", value = o$independent_filtering)
    updateCheckboxInput(session, "a_cooks", value = o$cooks_cutoff)
    updateRadioButtons(session, "a_deseq2_test", selected = o$deseq2_test)
    shrink <- if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal"
    updateSelectInput(session, "a_shrink", selected = shrink)
    len <- lengths_for_explorer()
    updateCheckboxInput(session, "a_length_offset", value = !is.null(dim(len)))
    updateRadioButtons(session, "a_transform", selected = o$transform)
    updateCheckboxInput(session, "a_blind", value = o$blind)
    showNotification(paste0("Recommended DESeq2 settings applied",
      if (shrink != "ashr") " (normal shrinkage; install.packages(\"ashr\") for ashr)" else "",
      if (is.null(dim(len))) "; gene-length offsets are off because no sample-specific lengths are loaded" else "",
      "."), type = "message", duration = 10)
  })

  output$a_shrink_note <- renderUI({
    sh <- input$a_shrink %||% "none"
    if (!sh %in% c("apeglm", "ashr") || requireNamespace(sh, quietly = TRUE)) return(NULL)
    div(class = "note-warn", strong("Not installed. Run:"),
        tags$pre(if (sh == "apeglm") 'BiocManager::install("apeglm")' else 'install.packages("ashr")'))
  })

  output$a_offset_note <- renderUI({
    len <- lengths_for_explorer()
    msg <- if (is.null(len)) "No gene lengths loaded: load a lengths file such as rsem.merged.gene_lengths.tsv." else
      if (is.null(dim(len))) "The loaded lengths are one per gene (from a GTF); offsets need sample-specific lengths (RSEM / Salmon) and would have no effect." else
      sprintf("Sample-specific lengths loaded (%s genes). Corrects for differences in transcript usage between samples.",
              format(nrow(len), big.mark = ","))
    if (!isTRUE(input$a_length_offset) && !is.null(dim(len))) return(helpText(style = "font-size:.82em", msg))
    if (!isTRUE(input$a_length_offset)) return(NULL)
    div(class = if (is.null(dim(len))) "note-warn" else "note-ok", style = "font-size:.82em", msg)
  })

  # limma's random block is chosen among the ticked factor covariates.
  observe({
    cv <- input$a_covars %||% character(0)
    f <- cv[vapply(cv, function(x) identical(av$types[[x]], "factor"), logical(1))]
    sel <- isolate(input$a_limma_block)
    updateSelectInput(session, "a_limma_block", choices = c("(none)", f),
                      selected = if (!is.null(sel) && sel %in% f) sel else "(none)")
  })

  output$a_checks <- renderUI({
    if (is.null(av$sheet)) return(NULL)
    ck <- checks()
    parts <- list()
    if (length(ck$errors)) parts[[1]] <- div(class = "note-err", strong("Fix before running: "),
      tags$ul(lapply(ck$errors, tags$li)))
    if (length(ck$warnings)) parts[[length(parts) + 1]] <- div(class = "note-warn",
      tags$ul(lapply(ck$warnings, tags$li)))
    if (!length(ck$errors) && identical(level(), "dtu")) {
      e <- input$a_engine %||% "edger_dtu"
      parts[[length(parts) + 1]] <- div(class = "note-ok", sprintf(
        "Design OK: %d samples, %d comparison(s). %s filters the transcripts when it runs (%s).",
        length(included()), nrow(av$contrasts), engine_label(e),
        if (identical(dtu_resolve_filter(e, input$a_dtu_filter %||% "auto"), "drimseq")) "DRIMSeq dmFilter" else "edgeR filterByExpr"))
    } else if (!length(ck$errors)) {
      o <- analysis_options()
      cd <- design_coldata(av$sheet, included(), input$a_main, av$types)
      flt <- tryCatch(filter_counts(round(av$counts[, cd$sample, drop = FALSE]), o$min_abundance,
                                    o$filter_mode, o$min_samples, o$min_proportion,
                                    group = cd[[input$a_main]], use_cpm = o$use_cpm),
                      error = function(e) NULL)
      parts[[length(parts) + 1]] <- div(class = "note-ok", sprintf(
        "Design OK: %d samples, %d comparison(s).%s", length(included()), nrow(av$contrasts),
        if (is.null(flt)) "" else sprintf(" Filtering keeps %s of %s genes (%s).",
          format(sum(flt$keep), big.mark = ","), format(length(flt$keep), big.mark = ","), flt$description)))
    }
    tagList(parts)
  })

  # ---- 4. run -------------------------------------------------------------------

  observeEvent(input$a_run, {
    if (is.null(av$counts)) { showNotification("Load a count matrix first.", type = "warning"); return() }
    ck <- checks()
    if (length(ck$errors)) {
      showNotification("The design has problems; see step 3.", type = "error", duration = 10); return()
    }
    engine <- input$a_engine %||% "deseq2"
    res <- withProgress(message = paste("Running", engine_label(engine)), value = 0, {
      tryCatch(run_de_analysis(av$counts, av$sheet, included(), input$a_main,
                               input$a_covars %||% character(0), av$types, av$contrasts,
                               engine, analysis_options(),
                               progress = function(f, d) setProgress(value = f, detail = d),
                               lengths = if (isTRUE(input$a_length_offset) || identical(level(), "dtu")) lengths_for_explorer() else NULL,
                               transcripts = if (tx_level()) transcript_inputs() else NULL),
               error = function(e) e)
    })
    if (inherits(res, "error")) {
      showNotification(paste("The analysis failed:", conditionMessage(res)), type = "error", duration = 20)
      return()
    }
    av$result <- res
    ids <- names(res$results)
    updateSelectInput(session, "a_explore_cmp", choices = ids, selected = ids[1])
    explore_comparison(ids[1])
    showNotification(sprintf("%s finished: %d comparison(s). '%s' is now loaded in the Explorer tabs.",
                             res$engine, length(ids), ids[1]), type = "message", duration = 10)
  })

  #' Load one comparison into the Explorer.
  explore_comparison <- function(id) {
    res <- av$result; req(res, id %in% names(res$results))
    df <- res$results[[id]]
    std <- standardise_de_table(df)
    tx <- isTRUE(res$level %in% c("transcript", "dtu"))
    ann <- if (tx) NULL else gene_annotation()
    fill <- apply_gene_annotation(std$data, ann)
    if (tx) {
      # transcripts are labelled by transcript name (IFIT1-201), else "symbol (transcript id)"
      d <- fill$data
      lab <- ifelse(!is.na(d$transcript_name) & nzchar(d$transcript_name), d$transcript_name,
             ifelse(!is.na(d$gene_name) & nzchar(d$gene_name), paste0(d$gene_name, " (", d$gene_id, ")"), d$gene_id))
      d$display_gene <- lab
      fill$data <- d
      ann <- data.frame(gene_id = d$gene_id, gene_name = lab, stringsAsFactors = FALSE)
    }
    samples <- res$coldata$sample

    meta <- av$sheet[av$sheet$sample %in% samples, , drop = FALSE]
    meta[] <- lapply(meta, as.character)
    if (res$main != "condition") meta$condition <- meta[[res$main]]
    rownames(meta) <- NULL

    rv$de <- fill$data; rv$de_raw <- df
    rv$de_notes <- c(std$notes, res$notes)
    rv$de_pipeline <- sprintf("%s (run in DEGexplorer)", res$engine_label %||% res$engine)
    rv$de_info <- list(pipeline = res$engine, flavour = "results", nfcore = FALSE,
                       annotated = !is.null(ann), filtered = FALSE,
                       label = rv$de_pipeline, notes = character(0))
    rv$gene_ann <- ann
    rv$universe_ids <- res$kept
    rv$counts <- av$counts[, samples, drop = FALSE]; rv$counts_ann <- av$counts_ann
    rv$lengths <- lengths_for_explorer(); rv$norm_supplied <- res$transformed
    rv$meta <- meta; rv$loaded <- TRUE
    rv$source_label <- sprintf("Analysis tab: %s, %s (%s)", res$engine, id, res$formula)
    rv$messages <- list(type = "ok", text = sprintf(
      "Loaded comparison %s from the Analysis tab: %s %s tested; %s matrix loaded for PCA and heatmaps.",
      id, format(nrow(df), big.mark = ","), if (tx) "transcripts" else "genes", res$transform_label))
    reset_enrichment()
    av$explored <- id
  }

  observeEvent(input$a_explore, {
    req(input$a_explore_cmp)
    explore_comparison(input$a_explore_cmp)
    showNotification(sprintf("'%s' is now loaded in the Explorer tabs.", input$a_explore_cmp),
                     type = "message")
  })

  output$a_result_note <- renderUI({
    res <- av$result
    if (is.null(res)) return(div(class = "note-warn", "No analysis run yet."))
    div(class = "note-ok",
      sprintf("%s (%s), run %s. Model %s. Filter: %s; %s of %s %s tested. %s matrix for PCA and heatmaps.",
              res$engine_label %||% res$engine,
              paste(names(res$versions), res$versions, collapse = ", "), res$run_at, res$formula, res$filter,
              format(length(res$kept), big.mark = ","), format(res$n_input, big.mark = ","),
              if (isTRUE(res$level %in% c("transcript", "dtu"))) "transcripts" else "genes",
              res$transform_label),
      if (!is.null(av$explored)) div(strong("In the Explorer now: "), av$explored),
      if (length(res$notes)) tags$ul(style = "margin:6px 0 0 0", lapply(res$notes, tags$li)))
  })

  output$a_result_table <- DT::renderDT({
    res <- av$result; req(res)
    p <- input$padj_cutoff %||% 0.05; l <- input$lfc_cutoff %||% 1
    tab <- do.call(rbind, lapply(names(res$results), function(id) {
      n <- de_result_counts(res$results[[id]], p, l)
      cr <- res$contrasts[res$contrasts$id == id, ]
      out <- data.frame(Comparison = id, Target = cr$target, Reference = cr$reference,
                        Tested = n[["tested"]], Up = n[["up"]], Down = n[["down"]],
                        stringsAsFactors = FALSE)
      if (identical(res$level, "dtu")) out$`Genes with DTU` <- dtu_gene_count(res$results[[id]], p)
      out
    }))
    DT::datatable(tab, rownames = FALSE, selection = "none",
                  caption = if (identical(res$level, "dtu"))
                    sprintf("Transcripts whose share of their gene rises (Up) or falls (Down) at adjusted p \u2264 %g and |log2 FC| \u2265 %g; genes at gene-level FDR \u2264 %g.", p, l, p)
                  else sprintf("Counts at adjusted p \u2264 %g and |log2 FC| \u2265 %g (sidebar thresholds).", p, l),
                  options = list(dom = "t", pageLength = 100, ordering = FALSE))
  })

  output$a_dl_result <- downloadHandler(
    filename = function() de_result_filename(input$a_explore_cmp %||% "results", result_file_tag(av$result)),
    content = function(file) {
      res <- av$result; req(res, input$a_explore_cmp)
      df <- de_result_for_export(res, input$a_explore_cmp, gene_annotation())
      utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    })

  output$a_dl_sig <- downloadHandler(
    filename = function() sub("\\.results\\.tsv$", ".significant.tsv",
                              de_result_filename(input$a_explore_cmp %||% "results", result_file_tag(av$result))),
    content = function(file) {
      res <- av$result; req(res, input$a_explore_cmp)
      df <- res$results[[input$a_explore_cmp]]
      p <- input$padj_cutoff %||% 0.05; l <- input$lfc_cutoff %||% 1
      df <- df[!is.na(df$padj) & df$padj <= p & !is.na(df$log2FoldChange) & abs(df$log2FoldChange) >= l, , drop = FALSE]
      df$direction <- ifelse(df$log2FoldChange > 0, "Up", "Down")
      df <- de_result_for_export(res, input$a_explore_cmp, gene_annotation(), df = df)
      utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    })

  # ---- run record, script and bundle ------------------------------------------------

  record_lines <- reactive({
    res <- av$result; req(res)
    de_run_record(res, input$padj_cutoff %||% 0.05, input$lfc_cutoff %||% 1, av$source)
  })
  output$a_record <- renderUI({
    if (is.null(av$result)) return(helpText("Run an analysis first."))
    div(style = "font-size:.88em;max-height:480px;overflow:auto;border:1px solid #e3e7ea;padding:8px 12px;border-radius:4px",
        shiny::markdown(paste(record_lines(), collapse = "\n")))
  })
  output$a_dl_record <- downloadHandler(
    filename = function() stamp_name("DEGexplorer_run_record", "md"),
    content = function(file) writeLines(record_lines(), file))
  output$a_dl_script <- downloadHandler(
    filename = function() "reproduce.R",
    content = function(file) { res <- av$result; req(res); writeLines(de_reproduce_script(res), file) })
  output$a_dl_bundle <- downloadHandler(
    filename = function() stamp_name(paste0("DEGexplorer_", result_file_tag(av$result)), "zip"),
    content = function(file) {
      res <- av$result; req(res)
      withProgress(message = "Building the export bundle", value = 0.3, {
        r <- tryCatch(write_de_bundle(res, file, input$padj_cutoff %||% 0.05, input$lfc_cutoff %||% 1,
                                      annotation = gene_annotation(), source = av$source),
                      error = function(e) e)
      })
      if (inherits(r, "error")) showNotification(paste("Could not build the bundle:", conditionMessage(r)),
                                                 type = "error", duration = 15)
    })

  # ---- 5. quality control ---------------------------------------------------------

  QC_HINTS <- c(
    library = "Very unequal library sizes are handled by normalisation, but a sample with far fewer reads or detected genes than the rest deserves a look.",
    distribution = "After normalisation the boxes should line up; a sample that stays shifted may be a technical outlier.",
    pca = "Samples should group by the condition you compare. If they group by donor or batch instead, add that variable as a covariate.",
    distance = "Replicates of a group should sit close together (small distances, joined early in the tree).",
    dispersion = "Gene-wise estimates (grey) should scatter around the red trend and fall with expression; blue points are the values used for testing.",
    pvalues = "Expect a flat background with a peak at 0. A peak at 1 or a hump in the middle suggests a mis-specified model or outliers.")
  QC_HINTS_ENGINE <- list(
    limma = c(dispersion = "voom: the square-root standard deviation should fall smoothly with expression along the red trend; the weights used in testing come from this trend."),
    edgeR = c(dispersion = "Biological coefficient of variation per gene (grey) around the trend (red) and the common value (blue line); a BCV of 0.1-0.4 is typical for replicate cell lines or mice, higher for human donors."),
    EBSeq = c(dispersion = "EBSeq fits no dispersion trend to show.",
              pvalues = "PPDE per gene: most genes sit near 0, differentially expressed ones near 1. Genes right of the dashed line (0.95) are called at FDR 5%."),
    `edgeR diffSplice` = c(dispersion = "Quasi-likelihood dispersions per transcript (quarter-root): raw (grey), squeezed (blue) and the trend (red), as edgeR's plotQLDisp draws them."),
    `limma diffSplice` = c(dispersion = "voomLmFit: the square-root standard deviation should fall smoothly with expression along the red trend."),
    DEXSeq = c(dispersion = "DEXSeq dispersions per transcript (first comparison): gene-wise estimates (grey) around the fitted trend (red); blue are the values used."),
    DRIMSeq = c(dispersion = "DRIMSeq precision per gene (first comparison) against mean gene expression; the red line is the common precision. Higher precision = less variable transcript shares."),
    satuRn = c(dispersion = "satuRn shrinks its quasi-binomial dispersions internally and has no plot to show here."))

  observe({
    res <- av$result; req(res)
    vars <- setdiff(colnames(qc_sample_table(res)), "sample")
    updateSelectInput(session, "a_qc_colour", choices = vars, selected = res$main)
    updateSelectInput(session, "a_qc_shape", choices = c("(none)", vars),
                      selected = if (length(res$covariates)) res$covariates[1] else "(none)")
    ids <- names(res$results)
    updateSelectizeInput(session, "a_qc_cmps", choices = ids, selected = ids)
  })

  output$a_qc_note <- renderUI({
    if (is.null(av$result)) return(div(class = "note-warn", "Run an analysis to see its quality-control plots."))
    w <- input$a_qc_which %||% "library"
    h <- QC_HINTS_ENGINE[[av$result$engine]][w]
    div(style = "font-size:.9em;color:#42525f;margin-bottom:6px",
        if (!is.null(h) && !is.na(h)) h else QC_HINTS[[w]])
  })

  qc_obj <- reactive({
    res <- av$result; req(res, input$a_qc_which)
    shp <- input$a_qc_shape
    build_qc_figure(input$a_qc_which, res,
                    colour_by = input$a_qc_colour %||% res$main,
                    shape_by = if (is.null(shp) || identical(shp, "(none)")) NULL else shp,
                    ntop = input$a_qc_ntop %||% 500,
                    annotate_by = unique(c(input$a_qc_colour %||% res$main, res$covariates)),
                    ids = input$a_qc_cmps)
  })

  output$a_qc_plot <- renderPlot({
    o <- qc_obj()
    if (inherits(o, "pheatmap")) { grid::grid.newpage(); grid::grid.draw(o$gtable) } else print(o)
  })

  for (fmt in c("png", "pdf", "tiff")) local({
    f <- fmt
    output[[paste0("a_qc_dl_", f)]] <- downloadHandler(
      filename = function() stamp_name(paste0("qc_", input$a_qc_which %||% "plot"), f),
      content = function(file) {
        o <- tryCatch(qc_obj(), error = function(e) e)
        if (inherits(o, "error")) {
          export_plot_safe(NULL, file, f, input$a_qc_w, input$a_qc_h, input$a_qc_dpi)
          showNotification(paste("Could not build the figure:", conditionMessage(o)), type = "error", duration = 12)
          return(invisible())
        }
        err <- export_plot_safe(o, file, f, input$a_qc_w, input$a_qc_h, input$a_qc_dpi)
        if (!is.null(err)) showNotification(paste("Export failed:", err), type = "error", duration = 12)
      })
  })

  output$a_qc_table <- DT::renderDT({
    res <- av$result; req(res)
    DT::datatable(qc_sample_summary(res), rownames = FALSE, selection = "none",
                  options = list(dom = "t", pageLength = 200, scrollX = TRUE))
  })

  # ---- transcript usage plot (DTU runs) ------------------------------------------

  output$a_usage_ui <- renderUI({
    res <- av$result
    if (is.null(res) || !identical(res$level, "dtu")) return(NULL)
    tagList(
      h5("Transcript usage of one gene"),
      fluidRow(
        column(5, selectizeInput("a_usage_gene", "Gene (most significant first)", choices = NULL)),
        column(7, helpText(style = "margin-top:26px",
          "Each point is one sample's share of the gene's counts carried by that transcript; bars are group means. * marks transcripts significant at the sidebar's adjusted p."))),
      plotOutput("a_usage_plot", height = "460px"),
      fluidRow(
        column(2, numericInput("a_usage_w", "Width (in)", 8, min = 2, max = 30, step = 0.5)),
        column(2, numericInput("a_usage_h", "Height (in)", 5, min = 2, max = 30, step = 0.5)),
        column(2, numericInput("a_usage_dpi", "DPI", 300, min = 72, max = 1200, step = 1)),
        column(2, br(), downloadButton("a_usage_dl_png", "PNG", class = "btn-block")),
        column(2, br(), downloadButton("a_usage_dl_pdf", "PDF", class = "btn-block")),
        column(2, br(), downloadButton("a_usage_dl_tiff", "TIFF", class = "btn-block"))))
  })
  observe({
    res <- av$result; id <- input$a_explore_cmp
    req(res, identical(res$level, "dtu"), id %in% names(res$results))
    ch <- dtu_gene_choices(res, id)
    updateSelectizeInput(session, "a_usage_gene", choices = ch, selected = if (length(ch)) ch[[1]] else NULL,
                         server = TRUE)
  })
  usage_obj <- reactive({
    res <- av$result; req(res, identical(res$level, "dtu"), input$a_usage_gene, input$a_explore_cmp)
    dtu_usage_plot(res, input$a_explore_cmp, input$a_usage_gene, padj_cutoff = input$padj_cutoff %||% 0.05)
  })
  output$a_usage_plot <- renderPlot(print(usage_obj()))
  for (fmt in c("png", "pdf", "tiff")) local({
    f <- fmt
    output[[paste0("a_usage_dl_", f)]] <- downloadHandler(
      filename = function() stamp_name(paste0("usage_", input$a_usage_gene %||% "gene"), f),
      content = function(file) {
        o <- tryCatch(usage_obj(), error = function(e) e)
        if (inherits(o, "error")) { export_plot_safe(NULL, file, f, input$a_usage_w, input$a_usage_h, input$a_usage_dpi); return(invisible()) }
        err <- export_plot_safe(o, file, f, input$a_usage_w, input$a_usage_h, input$a_usage_dpi)
        if (!is.null(err)) showNotification(paste("Export failed:", err), type = "error", duration = 12)
      })
  })

  # Returned for tests; the app ignores it.
  invisible(list(av = av, checks = checks, pending_covars = pending_covars))
}
