# =============================================================================
# ui.R -- user interface
# =============================================================================

# ---- small builders ---------------------------------------------------------

#' Style controls (fonts, points, labels) -- one copy per plot tab so each
#' panel keeps its own look.
style_controls <- function(p, label_size = 3.2, point_size = 1.6) {
  tagList(
    tags$details(
      tags$summary(tags$strong("Fonts, points & labels")),
      div(class = "ctl-block",
        sliderInput(paste0(p, "_base_size"), "Base font size (pt)",
                    min = 6, max = 28, value = 13, step = 0.5),
        sliderInput(paste0(p, "_title_size"), "Title size",
                    min = 6, max = 36, value = 15, step = 0.5),
        sliderInput(paste0(p, "_axis_title_size"), "Axis title size",
                    min = 6, max = 30, value = 13, step = 0.5),
        sliderInput(paste0(p, "_axis_text_size"), "Axis tick-label size",
                    min = 4, max = 26, value = 11, step = 0.5),
        sliderInput(paste0(p, "_legend_text_size"), "Legend text size",
                    min = 4, max = 26, value = 11, step = 0.5),
        tags$hr(),
        sliderInput(paste0(p, "_label_size"), "Gene tag size",
                    min = 1, max = 12, value = label_size, step = 0.1),
        checkboxInput(paste0(p, "_label_bold"), "Bold gene tags", TRUE),
        colourpicker::colourInput(paste0(p, "_label_colour"), "Gene tag colour", "#000000"),
        sliderInput(paste0(p, "_max_overlaps"), "Max overlapping tags",
                    min = 5, max = 200, value = 30, step = 5),
        tags$hr(),
        sliderInput(paste0(p, "_point_size"), "Point size",
                    min = 0.2, max = 8, value = point_size, step = 0.1),
        sliderInput(paste0(p, "_point_alpha"), "Point opacity",
                    min = 0.05, max = 1, value = 0.7, step = 0.05),
        sliderInput(paste0(p, "_highlight_size"), "Highlighted point size",
                    min = 0.5, max = 12, value = 2.8, step = 0.1),
        checkboxInput(paste0(p, "_grid"), "Show grid lines", TRUE)
      )
    )
  )
}

#' Title / axis-label overrides.
label_controls <- function(p, default_title = "") {
  tags$details(
    tags$summary(tags$strong("Titles & axis labels")),
    div(class = "ctl-block",
      textInput(paste0(p, "_title"), "Plot title", value = default_title),
      textInput(paste0(p, "_subtitle"), "Subtitle", value = ""),
      textInput(paste0(p, "_xlab"), "X-axis label", value = "",
                placeholder = "leave blank for the default"),
      textInput(paste0(p, "_ylab"), "Y-axis label", value = "",
                placeholder = "leave blank for the default")
    )
  )
}

#' Gene-selection controls shared by the DE plot tabs.
gene_select_controls <- function(p, default_mode = "top_n", default_n = 15) {
  tagList(
    selectInput(paste0(p, "_gene_mode"), "Genes to label",
                choices = GENE_SELECT_MODES, selected = default_mode),
    conditionalPanel(
      sprintf("['top_n','top_fc','top_up','top_down'].indexOf(input.%s_gene_mode) > -1", p),
      numericInput(paste0(p, "_top_n"), "How many (N)", value = default_n, min = 1, max = 500)),
    conditionalPanel(
      sprintf("input.%s_gene_mode == 'dropdown'", p),
      selectizeInput(paste0(p, "_gene_dropdown"), "Search genes", choices = NULL,
                     multiple = TRUE,
                     options = list(placeholder = "type a gene symbol or Ensembl ID...",
                                    maxOptions = 500))),
    conditionalPanel(
      sprintf("input.%s_gene_mode == 'user_list'", p),
      textAreaInput(paste0(p, "_gene_paste"), "Paste gene list", rows = 4,
                    placeholder = "IFIT1, MX1, OAS1\nor Ensembl IDs, one per line"),
      uiOutput(paste0(p, "_paste_feedback"))),
    conditionalPanel(
      sprintf("input.%s_gene_mode == 'pathway'", p),
      selectInput(paste0(p, "_pathway_term"), "Enriched term",
                  choices = c("Run enrichment first" = "")))
  )
}

#' Export block (size / DPI / TIFF / PDF).
export_controls <- function(p, w = 7, h = 6, note = NULL, png = FALSE) {
  tagList(
    tags$hr(),
    tags$strong("Export (publication quality)"),
    if (!is.null(note)) helpText(note),
    fluidRow(
      column(6, numericInput(paste0(p, "_w"), "Width (in)", w, min = 2, max = 30, step = 0.5)),
      column(6, numericInput(paste0(p, "_h"), "Height (in)", h, min = 2, max = 30, step = 0.5))),
    numericInput(paste0(p, "_dpi"), "Resolution (DPI)", 300, min = 72, max = 1200, step = 1),
    if (png) fluidRow(
      column(4, downloadButton(paste0(p, "_dl_tiff"), "TIFF", class = "btn-block")),
      column(4, downloadButton(paste0(p, "_dl_pdf"), "PDF", class = "btn-block")),
      column(4, downloadButton(paste0(p, "_dl_png"), "PNG", class = "btn-block")))
    else fluidRow(
      column(6, downloadButton(paste0(p, "_dl_tiff"), "TIFF", class = "btn-block")),
      column(6, downloadButton(paste0(p, "_dl_pdf"), "PDF", class = "btn-block")))
  )
}

# ---- the UI -----------------------------------------------------------------

app_ui <- function() {
  fluidPage(
    theme = shinythemes::shinytheme("flatly"),
    shinyjs::useShinyjs(),
    tags$head(tags$style(HTML("
      .btn-block{width:100%;margin-bottom:4px}
      .well{background:#f7f9fa}
      h4.section{margin-top:6px;color:#2c3e50}
      details{margin:8px 0;border:1px solid #dfe5ea;border-radius:4px;padding:6px 9px;background:#fff}
      details summary{cursor:pointer;color:#2c3e50;outline:none}
      .ctl-block{padding-top:8px}
      .shiny-notification{font-size:.92em}
      .note-ok{color:#1d6f42;background:#eaf6ef;border-left:4px solid #4caf82;padding:8px 12px;margin:8px 0;border-radius:3px}
      .note-warn{color:#8a6d3b;background:#fcf8e3;border-left:4px solid #dcc27a;padding:8px 12px;margin:8px 0;border-radius:3px}
      .note-err{color:#a4372a;background:#fdecea;border-left:4px solid #e0796b;padding:8px 12px;margin:8px 0;border-radius:3px}
      .status-pill{display:inline-block;padding:2px 9px;border-radius:11px;font-size:.8em;margin-right:5px}
      .pill-ok{background:#d6f0e0;color:#1d6f42}.pill-no{background:#f2d7d3;color:#a4372a}
      table.dataTable td{word-break:break-word;vertical-align:top}
      table.dataTable td.gene-col{max-width:280px;min-width:170px}
      table.dataTable th.dt-head-nowrap{white-space:nowrap}
      details.gene-cell{margin:0;border:none;padding:0;background:transparent}
      details.gene-cell > summary{list-style:none;font-size:.9em;line-height:1.35}
      details.gene-cell > summary::-webkit-details-marker{display:none}
      details.gene-cell > summary::before{content:\"\\25B8\";color:#7f8c8d;margin-right:5px}
      details.gene-cell[open] > summary::before{content:\"\\25BE\"}
      details.gene-cell .gene-more{color:#2980b9;font-size:.9em;white-space:nowrap}
      details.gene-cell .gene-list{margin-top:5px;padding:5px 7px;max-height:170px;overflow:auto;
        font-size:.86em;line-height:1.45;background:#f7f9fa;border:1px solid #e3e8ec;border-radius:3px}
      .gene-cell-plain{font-size:.9em}
    "))),
    titlePanel(div(strong("DEGexplorer"),
                   span(" \u2014 differential expression, pathways & gene intelligence",
                        style = "font-weight:300;color:#556;font-size:.8em;"))),

    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Significance thresholds", class = "section"),
        helpText("Applied everywhere: point colouring, gene sets, heatmaps and enrichment input."),
        numericInput("padj_cutoff", "Adjusted p-value \u2264", 0.05, min = 0, max = 1, step = 0.01),
        numericInput("lfc_cutoff", "|log2 fold change| \u2265", 1, min = 0, max = 20, step = 0.1),
        numericInput("basemean_cutoff", "Minimum mean expression", 0, min = 0, step = 10),
        uiOutput("threshold_summary"),
        tags$details(open = NA,
          tags$summary(tags$strong("Organism")),
          selectInput("organism", NULL, organism_choices()),
          conditionalPanel("input.organism == 'other'",
            textInput("org_other_species", "Species name", placeholder = "e.g. Ovis aries"),
            textInput("org_other_orgdb", "OrgDb package (for GO / KEGG)", placeholder = "e.g. org.Xl.eg.db"),
            textInput("org_other_kegg", "KEGG organism code", placeholder = "e.g. oas")),
          uiOutput("organism_note")),
        tags$hr(),

        # ---- per-tab control blocks ----
        conditionalPanel("input.tabs == 'volcano'",
          h4("Volcano plot", class = "section"),
          gene_select_controls("v"),
          tags$details(tags$summary(tags$strong("Colours & cutoff lines")),
            div(class = "ctl-block",
              fluidRow(
                column(4, colourpicker::colourInput("v_up_colour", "Up", "#D6604D")),
                column(4, colourpicker::colourInput("v_down_colour", "Down", "#4393C3")),
                column(4, colourpicker::colourInput("v_ns_colour", "NS", "#BFBFBF"))),
              checkboxInput("v_cutoff_lines", "Show cutoff lines", TRUE))),
          label_controls("v", "Volcano plot"),
          style_controls("v"),
          export_controls("v")),

        conditionalPanel("input.tabs == 'ma'",
          h4("MA plot", class = "section"),
          gene_select_controls("m"),
          tags$details(tags$summary(tags$strong("Colours & cutoff lines")),
            div(class = "ctl-block",
              fluidRow(
                column(4, colourpicker::colourInput("m_up_colour", "Up", "#D6604D")),
                column(4, colourpicker::colourInput("m_down_colour", "Down", "#4393C3")),
                column(4, colourpicker::colourInput("m_ns_colour", "NS", "#BFBFBF"))),
              checkboxInput("m_cutoff_lines", "Show cutoff lines", TRUE))),
          label_controls("m", "MA plot"),
          style_controls("m"),
          export_controls("m")),

        conditionalPanel("input.tabs == 'scatter'",
          h4("Expression scatter", class = "section"),
          helpText("Mean expression per group, one axis each, with the identity line."),
          selectInput("s_value", "Expression measure",
                      choices = c("FPKM / RPKM" = "fpkm", "TPM" = "tpm",
                                  "CPM" = "cpm", "DESeq2-normalised" = "deseq2",
                                  "TMM-normalised" = "tmm", "Raw counts" = "raw",
                                  "Supplied matrix as-is" = "asis"),
                      selected = "fpkm"),
          checkboxInput("s_log", "log2 transform (recommended)", TRUE),
          fluidRow(
            column(6, selectInput("s_group_x", "X axis group", choices = NULL)),
            column(6, selectInput("s_group_y", "Y axis group", choices = NULL))),
          gene_select_controls("s", default_mode = "top_n", default_n = 8),
          tags$details(tags$summary(tags$strong("Off-target & focus genes")),
            div(class = "ctl-block",
              textAreaInput("s_offtarget", "Off-target genes (open diamonds)", rows = 2,
                            placeholder = "PPIB, PPIC, PPIE, PPIF, PPIH"),
              textAreaInput("s_focus", "Focus / on-target genes (filled dot)", rows = 2,
                            placeholder = "PPIA"))),
          tags$details(tags$summary(tags$strong("Colours & guides")),
            div(class = "ctl-block",
              fluidRow(
                column(4, colourpicker::colourInput("s_up_colour", "Up", "#E8A33D")),
                column(4, colourpicker::colourInput("s_down_colour", "Down", "#2C5AA0")),
                column(4, colourpicker::colourInput("s_ns_colour", "NS", "#CCCCCC"))),
              colourpicker::colourInput("s_identity_colour", "Identity line", "#2E7D32"),
              checkboxInput("s_show_identity", "Show identity line", TRUE),
              checkboxInput("s_show_counts", "Show up/down gene counts", TRUE))),
          label_controls("s"),
          style_controls("s", label_size = 3.4, point_size = 1.2),
          export_controls("s", w = 7.5, h = 6.5)),

        conditionalPanel("input.tabs == 'pca'",
          h4("PCA & scree", class = "section"),
          selectInput("p_matrix", "Expression source", choices = NULL),
          selectInput("p_norm", "Normalisation", choices = NORM_METHODS, selected = "logcpm"),
          numericInput("p_ntop", "Most variable genes to use", 500, min = 10, max = 50000, step = 100),
          checkboxInput("p_scale_genes", "Scale genes to unit variance", FALSE),
          fluidRow(
            column(6, numericInput("p_pc_x", "PC on X", 1, min = 1, max = 20)),
            column(6, numericInput("p_pc_y", "PC on Y", 2, min = 1, max = 20))),
          selectInput("p_colour_by", "Colour by", choices = NULL),
          selectInput("p_shape_by", "Shape by", choices = c("(none)" = "(none)")),
          checkboxInput("p_show_labels", "Label sample names", TRUE),
          checkboxInput("p_show_ellipse", "Draw group ellipses", FALSE),
          numericInput("p_scree_n", "PCs in scree plot", 10, min = 2, max = 50),
          checkboxInput("p_scree_cumulative", "Show cumulative line", TRUE),
          label_controls("p", "Principal component analysis"),
          style_controls("p", point_size = 1.8),
          export_controls("p", w = 7.5, h = 6)),

        conditionalPanel("input.tabs == 'heatmap'",
          h4("Heatmap", class = "section"),
          radioButtons("h_gene_mode", "Genes to show",
                       choices = c("Top N DEGs" = "top_n",
                                   "All significant DEGs" = "all_sig",
                                   "ALL genes in the matrix" = "all_genes",
                                   "Select genes (searchable)" = "dropdown",
                                   "Paste my own list" = "user_list",
                                   "Genes from an enriched pathway" = "pathway")),
          conditionalPanel("input.h_gene_mode == 'top_n'",
            numericInput("h_top_n", "How many (N)", 30, min = 2, max = 2000)),
          conditionalPanel("input.h_gene_mode == 'all_genes'",
            div(class = "note-warn",
                "All-genes mode uses the full expression matrix, so genes that were never ",
                "significant are included. Rows are capped below for responsiveness."),
            numericInput("h_max_genes", "Row cap (most variable kept)", 2000, min = 50, max = 20000, step = 50)),
          conditionalPanel("input.h_gene_mode == 'dropdown'",
            selectizeInput("h_gene_dropdown", "Search genes", choices = NULL, multiple = TRUE,
                           options = list(placeholder = "type a gene...", maxOptions = 500))),
          conditionalPanel("input.h_gene_mode == 'user_list'",
            textAreaInput("h_gene_paste", "Paste gene list", rows = 4),
            uiOutput("h_paste_feedback")),
          conditionalPanel("input.h_gene_mode == 'pathway'",
            selectInput("h_pathway_term", "Enriched term", choices = c("Run enrichment first" = ""))),
          tags$hr(),
          selectInput("h_matrix", "Expression source", choices = NULL),
          selectInput("h_norm", "Normalisation", choices = NORM_METHODS, selected = "logcpm"),
          selectInput("h_reference", "Control / reference group (drawn on the LEFT)", choices = NULL),
          checkboxInput("h_scale_rows", "Row Z-score", TRUE),
          fluidRow(
            column(6, checkboxInput("h_cluster_rows", "Cluster genes", TRUE)),
            column(6, checkboxInput("h_cluster_cols", "Cluster samples", FALSE))),
          helpText("Leave sample clustering off to keep control on the left and treated on the right."),
          fluidRow(
            column(6, checkboxInput("h_show_rownames", "Gene names", TRUE)),
            column(6, checkboxInput("h_show_colnames", "Sample names", TRUE))),
          checkboxInput("h_gap", "Gap between groups", TRUE),
          checkboxGroupInput("h_annotate_by", "Annotate samples by", choices = NULL),
          selectInput("h_palette", "Colour palette",
                      choices = c("Red-Blue" = "RdBu", "Red-Yellow-Blue" = "RdYlBu",
                                  "Viridis" = "Viridis", "Magma" = "Magma",
                                  "Navy-White-Firebrick" = "navy",
                                  "Green-White-Purple" = "GreenRed")),
          label_controls("h", "Heatmap"),
          style_controls("h"),
          export_controls("h", w = 7.5, h = 8)),

        conditionalPanel("input.tabs == 'enrich'",
          h4("Enrichment", class = "section"),
          radioButtons("e_gene_source", "Gene set",
                       choices = c("Up-regulated DEGs" = "up",
                                   "Down-regulated DEGs" = "down",
                                   "All DEGs (up + down)" = "all")),
          conditionalPanel("input.e_gene_source == 'all'",
            checkboxInput("e_split", "Separate Up / Down panels", TRUE)),
          selectInput("e_database", "Database", choices = ENRICH_DATABASES),
          selectInput("e_plot_type", "Plot type",
                      choices = c("Dot plot" = "dot", "Bar plot" = "bar",
                                  "Gene-concept network" = "cnet",
                                  "Enrichment map" = "emap",
                                  "Gene-term heatmap" = "heat")),
          numericInput("e_pvalue_cutoff", "p-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
          numericInput("e_qvalue_cutoff", "q-value cutoff", 0.2, min = 0, max = 1, step = 0.01),
          numericInput("e_top_terms", "Terms per panel", 15, min = 3, max = 60),
          actionButton("e_run", "Run enrichment", class = "btn-primary btn-block"),
          label_controls("e", "Enrichment results"),
          style_controls("e"),
          export_controls("e", w = 8, h = 7),
          downloadButton("e_dl_csv", "Results table (CSV)", class = "btn-block"),
          downloadButton("e_dl_xlsx", "Results table (Excel)", class = "btn-block")),

        conditionalPanel("input.tabs == 'gallery'",
          h4("Enrichment gallery", class = "section"),
          helpText("GO (BP, CC, MF) and KEGG results together, drawn as a set of ready-made figures."),
          radioButtons("x_source", "Results to draw",
                       choices = c("Import my enrichment tables" = "import",
                                   "Run GO + KEGG on the current DEGs" = "run"),
                       selected = "import"),
          conditionalPanel("input.x_source == 'import'",
            fileInput("x_files", "Enrichment tables (.xlsx, .csv, .tsv)", multiple = TRUE,
                      accept = c(".xlsx", ".csv", ".tsv", ".txt", ".gz")),
            helpText(HTML(paste0(
              "clusterProfiler layout (<code>ID, Description, GeneRatio, BgRatio, pvalue, ",
              "p.adjust, qvalue, geneID, Count</code>); DAVID and g:Profiler tables are ",
              "recognised too. Select several files at once, e.g. <code>go.BP.all.xlsx</code>, ",
              "<code>go.CC.all.xlsx</code>, <code>go.MF.all.xlsx</code>, ",
              "<code>pathway.all.xlsx</code>."))),
            selectInput("x_hint", "If a file does not say which database it is",
                        choices = c("Work it out from the file name" = "auto",
                                    "GO: Biological process" = "BP",
                                    "GO: Cellular component" = "CC",
                                    "GO: Molecular function" = "MF",
                                    "KEGG pathway" = "KEGG")),
            actionButton("x_import", "Import tables", class = "btn-primary btn-block"),
            tags$details(tags$summary(tags$strong("Example files")),
              div(class = "ctl-block",
                helpText(HTML(paste0(
                  "Four tables in the import layout, computed from this app's own demo ",
                  "data (the bundled example's DESeq2 results): <code>go.BP.all.csv</code>, ",
                  "<code>go.CC.all.csv</code>, <code>go.MF.all.csv</code>, ",
                  "<code>pathway.all.csv</code>. Download them as a template, or import ",
                  "them directly to try the figures."))),
                downloadButton("x_dl_demo", "Download example files (.zip)", class = "btn-block"),
                actionButton("x_import_demo", "Import the example files", class = "btn-block")))),
          conditionalPanel("input.x_source == 'run'",
            radioButtons("x_gene_source", "Gene set",
                         choices = c("All DEGs (up + down)" = "all",
                                     "Up-regulated DEGs" = "up",
                                     "Down-regulated DEGs" = "down")),
            helpText("Uses the significance thresholds above and the Organism setting. Needs clusterProfiler and the organism's OrgDb package; KEGG needs an internet connection."),
            actionButton("x_run", "Run GO (BP, CC, MF) + KEGG", class = "btn-primary btn-block")),
          tags$hr(),
          selectInput("x_fig", "Figure", choices = GALLERY_FIGURES, selected = "score_bar"),
          conditionalPanel(
            "['score_bar','score_dot','sankey','cnet','emap','goplot'].indexOf(input.x_fig) > -1",
            selectInput("x_db", "Database", choices = c("Load results first" = ""))),
          conditionalPanel("input.x_fig == 'pathview'",
            selectInput("x_pathway", "Pathway", choices = c("Load KEGG results first" = "")),
            numericInput("x_pv_limit", "Colour scale limit (|log2FC|; 0 = automatic)", 0,
                         min = 0, max = 10, step = 0.5),
            actionButton("x_pv_render", "Render pathway map", class = "btn-primary btn-block"),
            helpText("Rendered by pathview, which downloads the map from KEGG (internet needed).")),
          conditionalPanel("input.x_fig != 'pathview'",
            fluidRow(
              column(6, selectInput("x_rank_by", "Significance by",
                                    choices = c("adjusted p" = "p.adjust", "p-value" = "pvalue"))),
              column(6, numericInput("x_cutoff", "Cutoff \u2264", 0.05, min = 0, max = 1, step = 0.01))),
            numericInput("x_top_n", "Terms per database", 10, min = 1, max = 60),
            tags$details(tags$summary(tags$strong("Figure options")),
              div(class = "ctl-block",
                conditionalPanel("['go3_h','go3_v','sankey','kegg_class','cnet'].indexOf(input.x_fig) > -1",
                  selectInput("x_palette", "Colour set",
                              choices = c("Match the figure type" = "auto", ENRICH_PALETTE_CHOICES),
                              selected = "auto")),
                conditionalPanel("input.x_fig == 'score_bar'",
                  colourpicker::colourInput("x_bar_colour", "Bar colour", "#ADD8E6")),
                conditionalPanel("input.x_fig == 'score_dot'",
                  radioButtons("x_dot_x", "X axis", inline = TRUE,
                               choices = c("Enrichment score" = "score", "Gene ratio" = "ratio"))),
                conditionalPanel("['score_bar','go3_h','go3_v','kegg_class'].indexOf(input.x_fig) > -1",
                  checkboxInput("x_show_values", "Print the value on each bar", FALSE),
                  radioButtons("x_digits", "Decimal places", choices = c(0, 1, 2), selected = 2,
                               inline = TRUE)),
                conditionalPanel("input.x_fig == 'go3_v'",
                  sliderInput("x_angle", "Label angle (degrees)", min = 40, max = 90, value = 70, step = 1),
                  textInput("x_y_breaks", "Y ticks (comma separated; blank = automatic)", "",
                            placeholder = "0,1,2,3,4")),
                conditionalPanel("input.x_fig == 'kegg_class'",
                  radioButtons("x_kegg_value", "Bar length", inline = TRUE,
                               choices = c("Gene count" = "count", "Enrichment score" = "score"))),
                conditionalPanel("input.x_fig == 'sankey'",
                  fluidRow(
                    column(6, numericInput("x_sankey_ratio", "Sankey width", 7, min = 1, max = 20, step = 0.5)),
                    column(6, numericInput("x_dot_ratio", "Dot-panel width", 5, min = 1, max = 20, step = 0.5)))),
                conditionalPanel("['sankey','cnet'].indexOf(input.x_fig) > -1",
                  numericInput("x_max_genes", "Maximum genes drawn", 100, min = 5, max = 400, step = 5)),
                conditionalPanel("input.x_fig == 'emap'",
                  sliderInput("x_min_sim", "Join terms sharing at least (Jaccard)", min = 0.05,
                              max = 0.9, value = 0.2, step = 0.05)),
                sliderInput("x_label_chars", "Longest term label (characters)", min = 15, max = 120,
                            value = 60, step = 5),
                radioButtons("x_font", "Font family", inline = TRUE,
                             choices = c("Arial / Helvetica" = "sans", "Times" = "serif"))))),
          label_controls("x", ""),
          style_controls("x", label_size = 3, point_size = 1.6),
          export_controls("x", w = 8, h = 6, png = TRUE,
                          note = "Width and height follow the figure type until you change them."),
          tags$hr(),
          tags$strong("Tables"),
          fluidRow(
            column(6, downloadButton("x_dl_all", "All terms", class = "btn-block")),
            column(6, downloadButton("x_dl_sig", "Significant", class = "btn-block"))),
          tags$hr(),
          tags$strong("Complete results bundle"),
          helpText("Every table and figure for every loaded database in one .zip: go/ and pathway/ folders, each figure as PNG, PDF and TIFF."),
          checkboxInput("x_bundle_pathview", "Include KEGG pathway maps (pathview; slow, needs internet)", FALSE),
          downloadButton("x_dl_bundle", "Download bundle (.zip)", class = "btn-block")),

        conditionalPanel("input.tabs == 'geneinfo'",
          h4("Gene intelligence", class = "section"),
          gene_select_controls("g", default_mode = "top_n", default_n = 20),
          textInput("g_species", "Species", "Homo sapiens"),
          textInput("g_context", "Literature context terms",
                    placeholder = "T cell, interferon"),
          helpText("Comma-separated. Narrows PubMed / Europe PMC searches to your study."),
          tags$hr(),
          checkboxInput("g_live", "Fetch live gene descriptions & publications", FALSE),
          div(class = "note-warn",
              "Live lookups need an internet connection (MyGene.info and Europe PMC). ",
              "Without them the report still contains every database link."),
          conditionalPanel("input.g_live",
            numericInput("g_n_pubs", "Publications per gene", 5, min = 1, max = 25),
            numericInput("g_max_live", "Max genes to fetch live", 25, min = 1, max = 200)),
          actionButton("g_build", "Build report", class = "btn-primary btn-block"),
          tags$hr(),
          downloadButton("g_dl_html", "Download report (HTML)", class = "btn-block"),
          downloadButton("g_dl_csv", "Gene table (CSV)", class = "btn-block")),

        conditionalPanel("input.tabs == 'analysis'", analysis_sidebar_ui()),

        conditionalPanel("input.tabs == 'genes'", gene_explorer_sidebar_ui()),

        conditionalPanel("input.tabs == 'data'",
          h4("Data", class = "section"),
          helpText("Load your own study, or explore the bundled example."),
          helpText("Starting from raw counts? The Analysis tab runs the differential expression here."))
      ),

      mainPanel(
        width = 9,
        tabsetPanel(
          id = "tabs",

          tabPanel("Data", value = "data", br(),
            uiOutput("data_status"),
            fluidRow(
              column(6,
                h4("1. Differential expression results"),
                helpText("DESeq2, limma / limma-voom, edgeR or any table with fold-change and p-value columns. Column names are detected automatically."),
                div(class = "note-ok", style = "font-size:.86em",
                    strong("nf-core/differentialabundance: "),
                    "any of the per-contrast outputs works here \u2014 ",
                    tags$code("*.limma.results.tsv"), ", ",
                    tags$code("*.limma.results_filtered.tsv"), ", ",
                    tags$code("*.limma.annotated.tsv"),
                    " and the matching ", tags$code("*.deseq2.*"), " files."),
                fileInput("up_de", NULL, accept = c(".tsv", ".txt", ".csv", ".gz")),
                uiOutput("de_mapping_ui"),
                tags$details(
                  tags$summary(tags$strong("Companion files (optional)")),
                  div(class = "ctl-block",
                    helpText(HTML(paste0(
                      "Supply these when the results file above is a bare ",
                      "<code>results.tsv</code> or a <code>results_filtered.tsv</code>."))),
                    fileInput("up_annot",
                              "Gene annotation \u2014 a GTF/GFF3, or a table of ids and symbols",
                              accept = c(".gtf", ".gff", ".gff3", ".tsv", ".txt", ".csv", ".gz")),
                    helpText(HTML(paste0(
                      "Drop in the same <code>.gtf</code> / <code>.gtf.gz</code> your pipeline ",
                      "used and gene symbols, biotypes and loci are read straight out of it \u2014 ",
                      "no id/symbol table needed. GENCODE version suffixes ",
                      "(<code>ENSG\u2026.5</code>) are matched automatically. GFF3 works too."))),
                    fileInput("up_universe", "Unfiltered results \u2014 enrichment background (e.g. *.results.tsv)",
                              accept = c(".tsv", ".txt", ".csv", ".gz")),
                    uiOutput("companion_note")))),
              column(6,
                h4("2. Expression matrices"),
                helpText("Raw counts (e.g. rsem.merged.gene_counts.tsv) let the app compute TPM, FPKM, CPM, DESeq2 and TMM values. Gene lengths are needed for TPM and FPKM."),
                fileInput("up_counts", "Raw / expected counts", accept = c(".tsv", ".txt", ".csv", ".gz")),
                fileInput("up_lengths", "Gene lengths (for TPM / FPKM)", accept = c(".tsv", ".txt", ".csv", ".gz")),
                annotation_fetch_ui("d", "No lengths file? Get gene lengths (and symbols) from Ensembl or GENCODE"),
                helpText("A GTF uploaded as the gene annotation (companion files, left) also supplies gene lengths when no lengths file is given."),
                fileInput("up_norm", "Pre-normalised matrix (optional, e.g. VST)", accept = c(".tsv", ".txt", ".csv", ".gz")))),
            fluidRow(
              column(6,
                h4("3. Sample sheet (optional)"),
                helpText("First column = sample id (matching the matrix columns), plus a condition / group column. No file? Leave this empty and build the sheet in the app, below, after loading."),
                fileInput("up_meta", NULL, accept = c(".csv", ".tsv", ".txt")),
                uiOutput("sample_subset_ui")),
              column(6,
                h4("4. Load"),
                actionButton("apply_upload", "Load uploaded files", class = "btn-primary btn-block"),
                tags$hr(),
                selectInput("example_choice", "Bundled example (stim vs mock, 3 donors)",
                            choices = EXAMPLE_DATASETS, selected = "deseq2"),
                actionButton("load_example", "Load bundled example", class = "btn-block"),
                helpText("The limma entries are three common layouts of the same comparison (annotated, plain results, DEGs only), so you can see exactly how each behaves."),
                actionButton("reset_all", "Clear everything", class = "btn-block"),
                br(), uiOutput("load_messages"))),
            hr(),
            h4("Sample sheet"),
            helpText(HTML(paste0("Uploaded in step 3, or built here from the matrix column names. ",
              "The <code>condition</code> column defines the groups used by Scatter, PCA and Heatmap; ",
              "add any other variables (donor, batch\u2026) to colour or annotate by them."))),
            uiOutput("d_sheet_note"),
            sheet_editor_ui("d", show_type = FALSE),
            hr(),
            h4("Computed expression matrices"),
            helpText("Derived from the raw counts you supplied. Download any of them for use elsewhere."),
            fluidRow(
              column(4, selectInput("dl_norm_method", "Matrix", choices = NORM_METHODS, selected = "tpm")),
              column(4, br(), downloadButton("dl_norm_csv", "Download CSV", class = "btn-block")),
              column(4, br(), downloadButton("dl_norm_tsv", "Download TSV", class = "btn-block"))),
            uiOutput("norm_preview_note"),
            DT::DTOutput("norm_preview")),

          tabPanel("Analysis", value = "analysis", br(), analysis_tab_ui()),

          tabPanel("Volcano", value = "volcano", br(),
            plotOutput("volcano_plot", height = "640px"),
            hr(), h4("Significant genes"), DT::DTOutput("volcano_table"),
            fluidRow(
              column(3, downloadButton("v_dl_table", "Significant genes (CSV)", class = "btn-block")),
              column(3, downloadButton("v_dl_all", "All genes (CSV)", class = "btn-block")),
              column(4, downloadButton("v_dl_xlsx", "All + significant + up + down (Excel)", class = "btn-block"))),
            helpText("Downloads use the sidebar thresholds for the Status column; the full lists keep every tested gene at full precision.")),

          tabPanel("MA", value = "ma", br(),
            plotOutput("ma_plot", height = "640px"),
            hr(), h4("Significant genes"), DT::DTOutput("ma_table")),

          tabPanel("Scatter", value = "scatter", br(),
            uiOutput("scatter_note"),
            plotOutput("scatter_plot", height = "680px")),

          tabPanel("PCA", value = "pca", br(),
            uiOutput("pca_note"),
            plotOutput("pca_plot", height = "520px"),
            hr(), h4("Scree plot"),
            plotOutput("scree_plot", height = "400px"),
            hr(), h4("Variance explained"), DT::DTOutput("pca_table")),

          tabPanel("Heatmap", value = "heatmap", br(),
            uiOutput("heatmap_note"),
            plotOutput("heatmap_plot", height = "760px")),

          tabPanel("Gene explorer", value = "genes", br(), gene_explorer_tab_ui()),

          tabPanel("Enrichment", value = "enrich", br(),
            uiOutput("enrich_pkg_note"),
            h4("Enriched terms"), DT::DTOutput("enrich_table"),
            hr(), uiOutput("enrich_plot_ui")),

          tabPanel("Enrichment gallery", value = "gallery", br(),
            uiOutput("gallery_status"),
            uiOutput("gallery_plot_ui"),
            hr(),
            fluidRow(
              column(8, h4(textOutput("gallery_table_title", inline = TRUE))),
              column(4, radioButtons("x_table_scope", NULL, inline = TRUE,
                                     choices = c("Significant" = "sig", "All terms" = "all")))),
            DT::DTOutput("gallery_table")),

          tabPanel("Gene intelligence", value = "geneinfo", br(),
            uiOutput("gene_report_note"),
            h4("Genes & resources"),
            helpText("Every link opens in a new tab. Build the report for a combined, downloadable dossier."),
            DT::DTOutput("gene_link_table"),
            hr(), h4("Report preview"), uiOutput("gene_report_preview"))
        )
      )
    )
  )
}
