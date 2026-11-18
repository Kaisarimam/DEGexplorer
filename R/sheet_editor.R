# =============================================================================
# sheet_editor.R -- in-app sample sheet editing, shared by the Data and
# Analysis tabs
#
# A sample sheet is a data frame with a `sample` column plus one character
# column per variable. Each tab keeps its sheet wherever it likes and hands
# this module accessor functions, so the same controls work on both.
# Input ids are "<prefix>_<name>", e.g. a_add_var on the Analysis tab.
# =============================================================================

#' Controls and table for editing a sample sheet.
#'
#' @param p          id prefix
#' @param show_type  offer factor / numeric for new variables
sheet_editor_ui <- function(p, show_type = TRUE) {
  id <- function(x) paste0(p, "_", x)
  tagList(
    fluidRow(
      column(4,
        tags$strong("Add a variable"),
        textInput(id("var_name"), NULL, placeholder = "e.g. condition, donor, batch"),
        if (show_type) radioButtons(id("var_type"), NULL, inline = TRUE,
                                    choices = c("Factor (groups)" = "factor", "Numeric covariate" = "numeric")),
        actionButton(id("add_var"), "Add variable", class = "btn-block")),
      column(4,
        tags$strong("Set a value for the selected rows"),
        selectInput(id("assign_var"), NULL, choices = NULL),
        textInput(id("assign_value"), NULL, placeholder = "value, e.g. treated"),
        actionButton(id("assign"), "Set for selected rows", class = "btn-block")),
      column(4,
        tags$strong("Fill a variable from the sample names"),
        selectInput(id("guess_var"), NULL, choices = NULL),
        fluidRow(
          column(6, textInput(id("guess_sep"), "Split on", "_")),
          column(6, numericInput(id("guess_pos"), "Take piece", 1, step = 1))),
        helpText("1 = first piece, -1 = last piece."),
        uiOutput(id("guess_preview")),
        actionButton(id("guess_apply"), "Fill from names", class = "btn-block"))),
    helpText("Double-click a cell to edit it. Click rows to select them for ‘Set for selected rows’."),
    DT::DTOutput(id("sheet")),
    fluidRow(
      column(4, selectInput(id("var_remove"), "Remove a variable", choices = NULL)),
      column(4, br(), actionButton(id("remove_var"), "Remove", class = "btn-block")))
  )
}

#' Server half of the sample sheet editor.
#'
#' @param get_sheet,set_sheet  read / replace the sheet (a data frame)
#' @param get_types,set_types  read / replace the named list of variable types;
#'                             set_types may be NULL when every variable is a factor
#' @param protected            variables that cannot be removed
#' @param on_remove            optional function(variable) called after a removal
sheet_editor_server <- function(input, output, session, p, get_sheet, set_sheet,
                                get_types, set_types = NULL, protected = character(0),
                                on_remove = NULL) {
  id <- function(x) paste0(p, "_", x)
  vars <- reactive({ s <- get_sheet(); if (is.null(s)) character(0) else setdiff(colnames(s), "sample") })

  # Re-draw the table only when its shape changes (new data, a variable added
  # or removed). Value edits are pushed with replaceData so the page, the
  # scroll position and the row selection survive.
  shape <- reactiveVal(NULL)
  observe({
    s <- get_sheet()
    shape(if (is.null(s)) NULL else paste(c(colnames(s), s$sample), collapse = "\r"))
  })

  output[[id("sheet")]] <- DT::renderDT({
    req(shape())
    DT::datatable(isolate(get_sheet()), rownames = FALSE, selection = "multiple",
                  editable = list(target = "cell", disable = list(columns = 0)),
                  options = list(pageLength = 200, dom = "t", ordering = FALSE, scrollX = TRUE))
  })
  proxy <- DT::dataTableProxy(id("sheet"))
  push <- function() DT::replaceData(proxy, get_sheet(), resetPaging = FALSE, rownames = FALSE,
                                     clearSelection = "none")

  observeEvent(input[[id("sheet_cell_edit")]], {
    sh <- DT::editData(get_sheet(), input[[id("sheet_cell_edit")]], rownames = FALSE)
    for (v in vars()) sh[[v]] <- as.character(sh[[v]])
    set_sheet(sh)
  })

  observe({
    v <- vars(); rem <- setdiff(v, protected)
    for (x in c("assign_var", "guess_var")) {
      sel <- isolate(input[[id(x)]])
      updateSelectInput(session, id(x), choices = v, selected = if (!is.null(sel) && sel %in% v) sel else v[1])
    }
    sel <- isolate(input[[id("var_remove")]])
    updateSelectInput(session, id("var_remove"), choices = rem,
                      selected = if (!is.null(sel) && sel %in% rem) sel else rem[1])
  })

  observeEvent(input[[id("add_var")]], {
    sh <- get_sheet()
    if (is.null(sh)) { showNotification("Load data first.", type = "warning"); return() }
    nm <- tryCatch(clean_variable_name(input[[id("var_name")]], colnames(sh)), error = function(e) e)
    if (inherits(nm, "error")) { showNotification(conditionMessage(nm), type = "warning"); return() }
    sh[[nm]] <- rep("", nrow(sh)); set_sheet(sh)
    if (is.function(set_types)) {
      ty <- get_types(); ty[[nm]] <- input[[id("var_type")]] %||% "factor"; set_types(ty)
    }
    updateTextInput(session, id("var_name"), value = "")
  })

  observeEvent(input[[id("remove_var")]], {
    v <- input[[id("var_remove")]]
    req(v, v %in% vars(), !v %in% protected)
    sh <- get_sheet(); sh[[v]] <- NULL; set_sheet(sh)
    if (is.function(set_types)) { ty <- get_types(); ty[[v]] <- NULL; set_types(ty) }
    if (is.function(on_remove)) on_remove(v)
  })

  observeEvent(input[[id("assign")]], {
    v <- input[[id("assign_var")]]; rows <- input[[id("sheet_rows_selected")]]
    if (is.null(v) || !v %in% vars()) { showNotification("Add a variable first.", type = "warning"); return() }
    if (!length(rows)) { showNotification("Click one or more rows in the table first.", type = "warning"); return() }
    sh <- get_sheet(); sh[[v]][rows] <- trimws(input[[id("assign_value")]] %||% ""); set_sheet(sh)
    push()
  })

  guess <- reactive({
    sh <- get_sheet(); req(sh, input[[id("guess_sep")]])
    tryCatch(guess_from_names(sh$sample, input[[id("guess_sep")]], input[[id("guess_pos")]] %||% 1),
             error = function(e) NULL)
  })

  output[[id("guess_preview")]] <- renderUI({
    g <- guess(); if (is.null(g)) return(NULL)
    lv <- unique(g[!is.na(g)])
    div(style = "font-size:.85em;color:#42525f;margin-bottom:6px",
        if (!length(lv)) "No piece at that position." else
          sprintf("Would give %d value(s): %s%s", length(lv),
                  paste(utils::head(lv, 6), collapse = ", "), if (length(lv) > 6) ", …" else ""),
        if (any(is.na(g))) sprintf(" (%d sample(s) left blank)", sum(is.na(g))))
  })

  observeEvent(input[[id("guess_apply")]], {
    v <- input[[id("guess_var")]]; g <- guess()
    if (is.null(v) || !v %in% vars()) { showNotification("Add a variable first.", type = "warning"); return() }
    req(g)
    sh <- get_sheet(); sh[[v]] <- ifelse(is.na(g), "", g); set_sheet(sh)
    push()
  })

  invisible(list(vars = vars))
}
