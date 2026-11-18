# =============================================================================
# gene_explorer.R -- every gene in the count matrix: expression plots, a full
# table with DE statistics, and whole-gene-list downloads
# =============================================================================

#' Annotation columns worth carrying into a full results table, when present.
DE_ANNOTATION_COLS <- c("gene_biotype", "gene_type", "biotype", "chromosome", "start", "end",
                        "strand", "description", "entrezid", "entrez_id")

#' The complete DE table at full precision: statistics, status against the
#' current thresholds, then any annotation columns the table carries.
de_full_table <- function(d) {
  num <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  out <- data.frame(gene_id = d$gene_id,
                    gene_name = if (!is.null(d$display_gene)) d$display_gene else d$gene_id,
                    stringsAsFactors = FALSE)
  for (cl in intersect(num, colnames(d))) out[[cl]] <- d[[cl]]
  out$status <- as.character(d$status %||% NA_character_)
  for (cl in intersect(DE_ANNOTATION_COLS, colnames(d))) out[[cl]] <- d[[cl]]
  out <- out[order(out$padj, out$pvalue, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' One workbook: all genes, significant, up and down, plus the thresholds used.
write_de_workbook <- function(d, path, padj_cutoff, lfc_cutoff, basemean_cutoff = 0) {
  full <- de_full_table(d)
  sheets <- list(
    `All genes`   = full,
    `Significant` = full[full$status %in% c("Up", "Down"), , drop = FALSE],
    `Up`          = full[full$status %in% "Up", , drop = FALSE],
    `Down`        = full[full$status %in% "Down", , drop = FALSE],
    `Thresholds`  = data.frame(
      setting = c("adjusted p-value <=", "|log2 fold change| >=", "minimum mean expression",
                  "genes", "up", "down", "written"),
      value = c(padj_cutoff, lfc_cutoff, basemean_cutoff, nrow(full),
                sum(full$status %in% "Up"), sum(full$status %in% "Down"),
                format(Sys.time(), "%Y-%m-%d %H:%M")),
      stringsAsFactors = FALSE))
  openxlsx::write.xlsx(sheets, path)
  invisible(path)
}

# ---- gene catalogue and expression table ----------------------------------------

#' Every gene the app knows about, with the best available symbol.
#'
#' Count-matrix genes come first; genes only in the DE table are added so the
#' table never hides a tested gene. Symbols come from the DE table, then the
#' gene annotation, then the count matrix's own gene_name column.
gene_catalog <- function(count_ids = NULL, de = NULL, gene_ann = NULL, counts_ann = NULL) {
  ids <- unique(c(as.character(count_ids), if (!is.null(de)) as.character(de$gene_id)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  sym <- rep(NA_character_, length(ids))
  fill <- function(sym, map_ids, map_sym) {
    if (is.null(map_ids)) return(sym)
    i <- annotation_match(ids, map_ids)
    v <- as.character(map_sym)[i]
    todo <- is.na(sym) | !nzchar(sym) | sym == ids
    sym[todo & !is.na(v) & nzchar(v)] <- v[todo & !is.na(v) & nzchar(v)]
    sym
  }
  if (!is.null(de)) sym <- fill(sym, de$gene_id, de$display_gene)
  if (!is.null(gene_ann)) sym <- fill(sym, gene_ann$gene_id, gene_ann$gene_name)
  if (!is.null(counts_ann) && all(c("gene_id", "gene_name") %in% colnames(counts_ann)))
    sym <- fill(sym, counts_ann$gene_id, counts_ann$gene_name)
  sym[is.na(sym) | !nzchar(sym)] <- ids[is.na(sym) | !nzchar(sym)]
  data.frame(gene_id = ids, symbol = sym, in_matrix = ids %in% count_ids,
             stringsAsFactors = FALSE)
}

#' Mean of each group for every gene.
group_means <- function(expr, meta, group_col = "condition") {
  smp <- intersect(colnames(expr), meta$sample)
  g <- as.character(meta[[group_col]][match(smp, meta$sample)])
  g[is.na(g) | !nzchar(g)] <- "(no group)"
  lv <- unique(g[order(g)])
  out <- vapply(lv, function(l) rowMeans(expr[, smp[g == l], drop = FALSE], na.rm = TRUE),
                numeric(nrow(expr)))
  if (!is.matrix(out)) out <- matrix(out, ncol = length(lv), dimnames = list(rownames(expr), lv))
  colnames(out) <- paste0("mean: ", lv)
  out
}

#' All genes in one table: ids, symbol, DE statistics, group means and every
#' sample's value.
#'
#' @param expr  gene x sample matrix (any normalisation), or NULL
#' @param d     DE table already classified (de_classified()), or NULL
gene_expression_table <- function(catalog, expr = NULL, meta = NULL, d = NULL,
                                  group_col = "condition", digits = 3) {
  out <- data.frame(`Gene ID` = catalog$gene_id, Symbol = catalog$symbol,
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(d)) {
    i <- match(catalog$gene_id, d$gene_id)
    out$`log2 FC` <- round(d$log2FoldChange[i], digits)
    out$`P value` <- signif(d$pvalue[i], 3)
    out$`adj. P` <- signif(d$padj[i], 3)
    st <- as.character(d$status[i])
    st[is.na(i)] <- "Not tested"
    out$Status <- st
  }
  if (!is.null(expr)) {
    j <- match(catalog$gene_id, rownames(expr))
    if (!is.null(meta) && group_col %in% colnames(meta)) {
      gm <- group_means(expr, meta, group_col)
      out <- cbind(out, round(gm[j, , drop = FALSE], digits))
    }
    smp <- if (!is.null(meta)) c(intersect(meta$sample, colnames(expr)), setdiff(colnames(expr), meta$sample))
           else colnames(expr)
    out <- cbind(out, round(expr[j, smp, drop = FALSE], digits))
  }
  rownames(out) <- NULL
  out
}

# ---- plot -------------------------------------------------------------------------

#' Expression of one or more genes across samples.
#'
#' @param layout "samples": one bar per sample, coloured by group;
#'               "groups": group mean +/- SE with the samples as points.
gene_bar_plot <- function(expr, meta, genes, catalog, d = NULL, group_col = "condition",
                          layout = "samples", log_scale = FALSE, show_points = TRUE,
                          value_label = "Expression", style = default_style(),
                          title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL) {
  st <- resolve_style(style)
  genes <- intersect(genes, rownames(expr))
  if (!length(genes)) stop("None of the selected genes is in this expression matrix.")
  smp <- intersect(meta$sample, colnames(expr))
  if (!length(smp)) smp <- colnames(expr)
  grp <- as.character(meta[[group_col]][match(smp, meta$sample)])
  grp[is.na(grp) | !nzchar(grp)] <- "(no group)"
  ord <- smp[order(grp, smp)]

  lab <- function(g) {
    s <- catalog$symbol[match(g, catalog$gene_id)]
    s <- ifelse(is.na(s), g, s)
    if (is.null(d)) return(ifelse(s == g, g, sprintf("%s (%s)", s, g)))
    k <- match(g, d$gene_id)
    stat <- ifelse(is.na(k), "not tested",
                   sprintf("log2FC %.2f \u00b7 adj. P %s \u00b7 %s", d$log2FoldChange[k],
                           ifelse(is.na(d$padj[k]), "NA", formatC(d$padj[k], format = "g", digits = 2)),
                           as.character(d$status[k])))
    sprintf("%s\n%s", ifelse(s == g, g, sprintf("%s (%s)", s, g)), stat)
  }
  df <- do.call(rbind, lapply(genes, function(g) data.frame(
    gene = lab(g), sample = smp, group = grp, value = as.numeric(expr[g, smp]),
    stringsAsFactors = FALSE)))
  df$gene <- factor(df$gene, levels = lab(genes))
  df$sample <- factor(df$sample, levels = ord)
  if (log_scale) df$value <- log2(pmax(df$value, 0) + 1)
  ylab_default <- if (log_scale) sprintf("log2(%s + 1)", value_label) else value_label
  pal <- stats::setNames(group_palette(length(unique(grp))), sort(unique(grp)))

  if (identical(layout, "groups")) {
    sm <- stats::aggregate(value ~ gene + group, df, function(v)
      c(mean = mean(v), se = if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else 0))
    sm <- data.frame(gene = sm$gene, group = sm$group, mean = sm$value[, "mean"], se = sm$value[, "se"])
    p <- ggplot2::ggplot(sm, ggplot2::aes(x = group, y = mean, fill = group)) +
      ggplot2::geom_col(width = 0.65, alpha = 0.85) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.2)
    if (show_points) p <- p + ggplot2::geom_point(
      data = df, ggplot2::aes(x = group, y = value), inherit.aes = FALSE,
      position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1),
      size = st$point_size + 0.6, alpha = 0.8)
    default_sub <- "bars: group mean ± SE; points: samples"
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = value, fill = group)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    default_sub <- NULL
  }
  p <- p + ggplot2::facet_wrap(~gene, scales = "free_y", ncol = if (length(genes) > 3) 2 else 1) +
    ggplot2::scale_fill_manual(values = pal, name = group_col) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    deg_theme(st)
  if (!identical(layout, "groups")) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  apply_labels(p, title, subtitle %||% default_sub, xlab, ylab,
               default_title = if (length(genes) == 1) "Gene expression" else "Gene expression",
               default_x = NULL, default_y = ylab_default)
}

# ---- UI ------------------------------------------------------------------------------

gene_explorer_sidebar_ui <- function() {
  tagList(
    h4("Gene explorer", class = "section"),
    selectInput("ge_source", "Expression matrix", choices = NULL),
    conditionalPanel("input.ge_source == 'counts'",
      selectInput("ge_method", "Normalisation", choices = NORM_METHODS, selected = "deseq2")),
    selectInput("ge_group_by", "Group samples by", choices = NULL),
    radioButtons("ge_layout", "Plot",
                 choices = c("A bar per sample" = "samples", "Group mean ± SE" = "groups")),
    checkboxInput("ge_log", "log2 scale (value + 1)", FALSE),
    conditionalPanel("input.ge_layout == 'groups'",
      checkboxInput("ge_points", "Show the samples as points", TRUE)),
    label_controls("ge", "Gene expression"),
    style_controls("ge"),
    export_controls("ge", w = 8, h = 5, png = TRUE))
}

gene_explorer_tab_ui <- function() {
  tagList(
    uiOutput("ge_note"),
    fluidRow(
      column(8, selectizeInput("ge_gene", "Genes to plot (search by symbol or Ensembl id; up to 12)",
                               choices = NULL, multiple = TRUE,
                               options = list(maxItems = 12, placeholder = "type a gene, e.g. IFIT1",
                                              maxOptions = 200))),
      column(4, br(), actionButton("ge_to_geneinfo", "Open in Gene intelligence", class = "btn-block"))),
    uiOutput("ge_links"),
    uiOutput("ge_plot_ui"),
    hr(),
    h4("All genes"),
    fluidRow(
      column(6, radioButtons("ge_scope", NULL, inline = TRUE,
                             choices = c("All" = "all", "Tested" = "tested", "Significant" = "sig",
                                         "Up" = "up", "Down" = "down"))),
      column(3, downloadButton("ge_dl_csv", "This table (CSV)", class = "btn-block")),
      column(3, downloadButton("ge_dl_xlsx", "This table (Excel)", class = "btn-block"))),
    helpText(paste("Every gene in the count matrix, with its DE statistics (Status uses the sidebar thresholds;",
                   "'Not tested' = not in the DE results, e.g. removed by filtering), group means and each",
                   "sample's value in the matrix chosen in the sidebar. Click a row to plot that gene;",
                   "the search box filters by symbol or id.")),
    DT::DTOutput("ge_table"))
}

# ---- server ---------------------------------------------------------------------------

#' @param get_expression,matrix_choices,de_classified,register_export helpers
#'   from app_server, passed in so the module uses the app's own matrices.
gene_explorer_server <- function(input, output, session, rv, de_classified, get_expression,
                                 matrix_choices, register_export) {

  observe({
    ch <- matrix_choices()
    sel <- isolate(input$ge_source)
    updateSelectInput(session, "ge_source", choices = ch,
                      selected = if (!is.null(sel) && sel %in% ch) sel else
                        if ("counts" %in% ch) "counts" else if (length(ch)) ch[1] else character(0))
  })

  # Gene lengths decide whether TPM / FPKM can be offered.
  observe({
    m <- if (is.null(rv$lengths)) NORM_METHODS[!NORM_METHODS %in% LENGTH_DEPENDENT] else NORM_METHODS
    sel <- isolate(input$ge_method) %||% "deseq2"
    updateSelectInput(session, "ge_method", choices = m, selected = if (sel %in% m) sel else "deseq2")
  })

  observe({
    req(rv$meta)
    cols <- setdiff(colnames(rv$meta), "sample")
    sel <- isolate(input$ge_group_by)
    updateSelectInput(session, "ge_group_by", choices = cols,
                      selected = if (!is.null(sel) && sel %in% cols) sel else
                        if ("condition" %in% cols) "condition" else cols[1])
  })

  catalog <- reactive({
    ids <- if (!is.null(rv$counts)) rownames(rv$counts) else
           if (!is.null(rv$norm_supplied)) rownames(rv$norm_supplied) else NULL
    req(!is.null(ids) || !is.null(rv$de))
    gene_catalog(ids, rv$de, rv$gene_ann, rv$counts_ann)
  })
  gene_choices <- reactive({
    cg <- catalog()
    stats::setNames(cg$gene_id, ifelse(cg$symbol == cg$gene_id, cg$gene_id,
                                       paste0(cg$symbol, "  (", cg$gene_id, ")")))
  })
  de_now <- reactive(if (is.null(rv$de)) NULL else de_classified())

  # New data: offer every gene and start on the most significant one.
  observe({
    ch <- gene_choices()
    d <- isolate(de_now())
    first <- if (!is.null(d) && any(!is.na(d$padj))) d$gene_id[which.min(d$padj)] else unname(ch[1])
    updateSelectizeInput(session, "ge_gene", choices = ch, selected = first, server = TRUE)
  })

  expr <- reactive({
    src <- input$ge_source
    validate(need(!is.null(src) && nzchar(src),
                  "Load a count matrix (Data or Analysis tab) to see expression values."))
    tryCatch(get_expression(src, input$ge_method %||% "deseq2"),
             error = function(e) validate(need(FALSE, conditionMessage(e))))
  })
  value_label <- reactive({
    if (identical(input$ge_source, "supplied")) return("Supplied normalised values")
    nm <- names(NORM_METHODS)[match(input$ge_method %||% "deseq2", NORM_METHODS)]
    if (is.na(nm)) "Expression" else nm
  })
  meta_now <- reactive({
    m <- rv$meta
    if (is.null(m)) m <- scaffold_metadata(colnames(expr()))
    m
  })
  group_col <- reactive({
    g <- input$ge_group_by
    if (is.null(g) || !g %in% colnames(meta_now())) "condition" else g
  })

  output$ge_note <- renderUI({
    if (is.null(rv$counts) && is.null(rv$norm_supplied) && is.null(rv$de)) {
      return(div(class = "note-warn", "Load data on the Data tab, or run an analysis on the Analysis tab."))
    }
    cg <- catalog()
    div(class = "note-ok", sprintf("%s genes in the matrix%s.", format(sum(cg$in_matrix), big.mark = ","),
        if (is.null(rv$de)) "; no DE results loaded, so the table shows expression only" else
          sprintf("; %s of them in the DE results", format(sum(cg$gene_id %in% rv$de$gene_id), big.mark = ","))))
  })

  ge_obj <- reactive({
    g <- input$ge_gene
    validate(need(length(g) > 0, "Choose one or more genes above, or click a row in the table."))
    lb <- collect_labels(input, "ge")
    gene_bar_plot(expr(), meta_now(), g, catalog(), de_now(), group_col = group_col(),
                  layout = input$ge_layout %||% "samples", log_scale = isTRUE(input$ge_log),
                  show_points = !isFALSE(input$ge_points), value_label = value_label(),
                  style = collect_style(input, "ge"),
                  title = if (nzchar(lb$title %||% "")) lb$title else NULL,
                  subtitle = if (nzchar(lb$subtitle %||% "")) lb$subtitle else NULL,
                  xlab = lb$xlab, ylab = lb$ylab)
  })

  output$ge_plot_ui <- renderUI({
    n <- max(1, length(input$ge_gene))
    h <- if (n == 1) 420 else if (n <= 3) 300 * n else 280 * ceiling(n / 2) + 80
    plotOutput("ge_plot", height = paste0(h, "px"))
  })
  output$ge_plot <- renderPlot(ge_obj())
  register_export("ge", ge_obj, "gene_expression")

  output$ge_links <- renderUI({
    g <- utils::head(input$ge_gene, 4); req(length(g))
    cg <- catalog(); sp <- input$g_species %||% "Homo sapiens"  # set from the app-wide Organism
    div(style = "font-size:.88em;margin:-4px 0 8px 0", lapply(g, function(id) {
      s <- cg$symbol[match(id, cg$gene_id)]
      L <- gene_resource_links(s, id, sp)
      div(strong(s), ": ", HTML(paste(sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>',
                                             unlist(L), names(L)), collapse = " &middot; ")))
    }))
  })

  observeEvent(input$ge_to_geneinfo, {
    g <- input$ge_gene
    if (!length(g)) { showNotification("Choose a gene first.", type = "warning"); return() }
    if (is.null(rv$de)) { showNotification("Gene intelligence needs DE results loaded.", type = "warning"); return() }
    ok <- intersect(g, rv$de$gene_id)
    if (!length(ok)) {
      showNotification("These genes are not in the DE results (filtered out), so Gene intelligence cannot report on them; the links above still work.",
                       type = "warning", duration = 10)
      return()
    }
    updateSelectInput(session, "g_gene_mode", selected = "dropdown")
    ch <- stats::setNames(rv$de$gene_id, paste0(rv$de$display_gene, "  (", rv$de$gene_id, ")"))
    updateSelectizeInput(session, "g_gene_dropdown", choices = ch, selected = ok, server = TRUE)
    updateTabsetPanel(session, "tabs", selected = "geneinfo")
    if (length(ok) < length(g)) showNotification(sprintf(
      "%d gene(s) are not in the DE results and were left out.", length(g) - length(ok)), type = "message")
  })

  full_table <- reactive({
    e <- tryCatch(expr(), error = function(e) NULL)
    gene_expression_table(catalog(), e, if (is.null(e)) NULL else meta_now(), de_now(), group_col())
  })
  scoped_table <- reactive({
    t <- full_table(); sc <- input$ge_scope %||% "all"
    if (sc == "all" || is.null(t$Status)) return(t)
    keep <- switch(sc, tested = t$Status != "Not tested", sig = t$Status %in% c("Up", "Down"),
                   up = t$Status == "Up", down = t$Status == "Down")
    t[keep, , drop = FALSE]
  })

  output$ge_table <- DT::renderDT({
    t <- scoped_table()
    if (!is.null(t$`adj. P`)) t <- t[order(t$`adj. P`, na.last = TRUE), , drop = FALSE]
    DT::datatable(t, rownames = FALSE, selection = "single",
                  options = list(pageLength = 15, scrollX = TRUE, deferRender = TRUE,
                                 search = list(regex = FALSE, caseInsensitive = TRUE)))
  }, server = TRUE)

  observeEvent(input$ge_table_rows_selected, {
    t <- scoped_table()
    if (!is.null(t$`adj. P`)) t <- t[order(t$`adj. P`, na.last = TRUE), , drop = FALSE]
    id <- t$`Gene ID`[input$ge_table_rows_selected]
    req(length(id) == 1)
    updateSelectizeInput(session, "ge_gene", choices = gene_choices(), selected = id, server = TRUE)
  })

  scope_name <- function() sprintf("genes_%s", input$ge_scope %||% "all")
  output$ge_dl_csv <- downloadHandler(
    filename = function() stamp_name(scope_name(), "csv"),
    content = function(f) utils::write.csv(scoped_table(), f, row.names = FALSE, na = ""))
  output$ge_dl_xlsx <- downloadHandler(
    filename = function() stamp_name(scope_name(), "xlsx"),
    content = function(f) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) stop('Excel export needs install.packages("openxlsx")')
      openxlsx::write.xlsx(list(Genes = scoped_table()), f)
    })

  invisible(list(catalog = catalog, full_table = full_table, ge_obj = ge_obj))
}
