#' The three-step selector - the ui function
#' @description Function should only be used for the developers
#' @param id id
#' @param right_margin margin on the right side, in px. For example, "20" translates to "20px".
#' @importFrom shinyWidgets pickerInput updatePickerInput
#' @export
#' @examples
#' if (interactive()) {
#'   library(shiny)
#'   library(Biobase)
#'   
#'   file <- system.file("extdata/demo.RDS", package = "omicsViewer")
#'   dat <- readRDS(file)
#'   fData <- fData(dat)
#'   triset <- stringr::str_split_fixed(colnames(fData), '\\|', n= 3)
#'   
#'   ui <- fluidPage(
#'     triselector_ui("tres"),
#'     triselector_ui("tres2")
#'   )
#'   server <- function(input, output, session) {
#'     v1 <- triselector_module("tres", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("mean.diff")
#'     )
#'     v2 <- triselector_module("tres2", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("log.fdr"))
#'     observe({
#'       print("/////////////////////////")
#'       print(v1())
#'     })
#'   }
#'   
#'   shinyApp(ui, server)
#' }
#' @return a tagList of UI components

triselector_ui <- function(id, right_margin = "20") {
  rmar <- sprintf("padding-left:2px; padding-right:%spx; padding-top:2px; padding-bottom:2px", right_margin)
  ns <- NS(id)
  tagList(
    fluidRow(
      column(2, offset = 0, align = "right",
             style='padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px',
             uiOutput(ns("groupLabel"))
      ),
      column(3, offset = 0, style='padding:2px;',
        selectInput(inputId = ns("analysis"), label = NULL, choices = NULL, selectize = TRUE, width = "100%") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-analysis-selector"))),
      column(4, offset = 0, style='padding:2px;',
        selectInput(inputId = ns("subset"), label = NULL, choices = NULL, selectize = TRUE, width = "100%") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-subset-selector"))),
      column(3, offset = 0, style=rmar, #'padding:2px;',
        selectInput(inputId = ns("variable"), label = NULL, choices = NULL, selectize = TRUE, width = "100%") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-variable-selector")))
      # column(3, offset = 0, style='padding:2px;'
        # pickerInput(inputId = ns("analysis"), label = NULL, choices = NULL, options = list(`live-search` = TRUE))),
      # column(4, offset = 0, style='padding:2px;',
        # pickerInput(inputId = ns("subset"), label = NULL, choices = NULL, options = list(`live-search` = TRUE))),
      # column(3, offset = 0, style="padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px", #'padding:2px;',
        # pickerInput(inputId = ns("variable"), label = NULL, choices = NULL, options = list(`live-search` = TRUE)))
    )
  )
}

#' The three-step selector - the module function
#' @description The selector is used to select columns of phenotype and feature data.
#' The cascade is self-repairing: when an upstream component changes (by the
#' user or by a store write) and the downstream values are no longer valid,
#' the first valid ones are pushed automatically. The module reports only
#' settled, coherent triples (all components acknowledged, naming a real
#' column); while a server push is in flight or the triple is incomplete,
#' the last committed triple is held, so consumers never observe mid-cascade
#' echoes.
#' Function should only be used for the developers.
#' @param id module id
#' @param reactive_x an nx3 matrix
#' @param reactive_selector1 default value for selector 1
#' @param reactive_selector2 default value for selector 2
#' @param reactive_selector3 default value for selector 3
#' @param reactive_axis_request version bump that forces re-derivation of the
#'   cascade (e.g. a widget-store epoch); unchanged values are never re-sent
#' @param label of the triselector
#' @export
#' @examples
#' if (interactive()) {
#'   library(shiny)
#'   library(Biobase)
#'
#'   file <- system.file("extdata/demo.RDS", package = "omicsViewer")
#'   dat <- readRDS(file)
#'   fData <- fData(dat)
#'   triset <- stringr::str_split_fixed(colnames(fData), '\\|', n= 3)
#'
#'   ui <- fluidPage(
#'     triselector_ui("tres"),
#'     triselector_ui("tres2")
#'   )
#'   server <- function(input, output, session) {
#'     v1 <- triselector_module("tres", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("mean.diff")
#'     )
#'     v2 <- triselector_module("tres2", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("log.fdr"))
#'     observe({
#'       print("/////////////////////////")
#'       print(v1())
#'     })
#'   }
#'
#'   shinyApp(ui, server)
#' }
#' @return an reactive object containing the selected values

triselector_module <- function(id,
                               reactive_x,
                               reactive_selector1 = reactive(NULL),
                               reactive_selector2 = reactive(NULL),
                               reactive_selector3 = reactive(NULL),
                               reactive_axis_request = reactive(NULL),
                               label = "Group Label:") {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # Validate reactive_x input - silently handle NULL during startup
  validated_x <- reactive({
    # Allow NULL during initial app startup
    if (is.null(reactive_x())) {
      return(NULL)
    }

    # Perform validation
    validation <- validate_triselector_input(
      reactive_x(),
      name = "triselector input",
      allow_null = TRUE
    )

    # Only show error notifications after app has initialized
    # (prevents error spam during startup)
    if (!validation$valid && length(names(input)) > 0) {
      showNotification(
        validation$message,
        type = "error",
        duration = 5
      )
      return(NULL)
    }

    # Return validated data or NULL
    if (validation$valid) {
      reactive_x()
    } else {
      NULL
    }
  })

  output$groupLabel <- renderUI({
    h5(HTML(sprintf("<b>%s</b>", label)))
  })

  # ------------------------------------------------------------------
  # One cascade implementation (WP2, HANDOVER section 8/F3).
  #
  # Per-component pending: a server->client push that the browser has not
  # acknowledged yet. While a push is in flight the widget may still report
  # the previous value; that report is a stale echo, not a user choice, and
  # must not re-derive the cascade (the historical revert bug). A report
  # equal to the pushed value retires the pending entry; an empty/NULL/
  # "--select--" report means the widget was re-bound or reset, which
  # invalidates any in-flight push (it must be re-sent).
  #
  # Inputs are read UNCONDITIONALLY: `reactive_selector1() %||%
  # input$analysis` never registered the input dependency while a store
  # value existed, which froze the cascade to the store's analysis and made
  # a manual analysis change impossible (RC6).
  #
  # The cascade repairs downstream itself: the store validates but never
  # repairs dependent keys, so when an upstream change invalidates the
  # downstream values the triselector pushes the first valid ones.
  #
  # Only coherent, settled triples are emitted: while any pending is set or
  # the current inputs name no real column, the last committed triple is
  # held, so echoes never reach the store or the plot.
  #
  # All observers are kept referenced (observer-GC rule).
  # ------------------------------------------------------------------
  .tri_keep <- list()
  .tri_keep_obs <- function(o) {
    .tri_keep[[length(.tri_keep) + 1L]] <<- o
    invisible(o)
  }
  pend <- list(analysis = reactiveVal(NULL),
               subset = reactiveVal(NULL),
               variable = reactiveVal(NULL))
  # last (choices, selected) sent per widget - skips no-op re-sends, which
  # is what keeps store epoch bumps from fanning updates out to widgets
  # that already show the requested state
  .last_sent <- list(analysis = NULL, subset = NULL, variable = NULL)
  # last selector value seen per component; a CHANGING non-NULL selector
  # value marks a fresh canonical (store-side) intent that supersedes any
  # in-flight push for the same component (e.g. a second store write before
  # the browser acknowledged the first)
  .last_sel <- list(analysis = NULL, subset = NULL, variable = NULL)
  committed <- reactiveVal(NULL)

  .live <- function(comp)
    switch(comp, analysis = input$analysis, subset = input$subset,
           variable = input$variable)

  for (comp in c("analysis", "subset", "variable")) local({
    cc <- comp
    .tri_keep_obs(observeEvent(.live(cc), {
      # Any report retires an in-flight push: an equal value is the
      # acknowledgement; a diverging value is a genuine user action - a
      # real browser fires no input event for an update that has not been
      # applied to the DOM, so divergence can only mean the user picked
      # (and the user always wins over the push). An empty/NULL/
      # "--select--" report is a re-bound or reset widget.
      pend[[cc]](NULL)
    }, ignoreInit = TRUE))
  })

  .first_valid <- function(cands, choices) {
    for (v in cands)
      if (!is.null(v) && length(v) == 1L && v %in% choices)
        return(v)
    NULL
  }

  .update_widget <- function(comp, choices, selected) {
    last <- .last_sent[[comp]]
    livev <- isolate(.live(comp))
    in_flight <- !is.null(isolate(pend[[comp]]()))
    # no-op skip: identical choices AND the widget already shows the
    # requested selection (acknowledged), or this exact push is still in
    # flight. A diverging live value (re-bound/cleared widget) re-sends.
    if (identical(choices, last$choices) &&
        (identical(selected, livev) ||
         (in_flight && identical(selected, last$selected))))
      return(invisible(FALSE))
    if (comp == "variable")
      updateSelectizeInput(session, inputId = "variable", choices = choices,
                           selected = selected,
                           server = length(choices) > 1000L)
    else
      updateSelectInput(session, inputId = comp, choices = choices,
                        selected = selected)
    .last_sent[[comp]] <<- list(choices = choices, selected = selected)
    # arming pending only for a value the widget does not already show: an
    # input whose value does not change fires no event, so a pending entry
    # set there could never be acknowledged. The "--select--" placeholder
    # is "no selection", not a value to protect, and must not arm pending
    # either - it would block the settle gate without an ack to clear it
    if (identical(selected, livev))
      pend[[comp]](NULL)
    else if (!is.null(selected) && !identical(selected, "--select--"))
      pend[[comp]](selected)
    invisible(TRUE)
  }

  .tri_keep_obs(observe({
    reactive_axis_request()  # version bumps force re-derivation after restores
    req(vx <- validated_x())
    if (length(names(input)) == 0L)
      return(NULL)

    # read everything unconditionally - every source must be a dependency
    live <- list(analysis = input$analysis, subset = input$subset,
                 variable = input$variable)
    sel <- list(analysis = reactive_selector1(),
                subset = reactive_selector2(),
                variable = reactive_selector3())
    pd <- list(analysis = pend$analysis(), subset = pend$subset(),
               variable = pend$variable())
    fresh <- list(
      analysis = !identical(sel$analysis, .last_sel$analysis),
      subset = !identical(sel$subset, .last_sel$subset),
      variable = !identical(sel$variable, .last_sel$variable))
    on.exit({
      .last_sel <<- list(analysis = sel$analysis, subset = sel$subset,
                         variable = sel$variable)
    })
    # once a triple has settled, an invalidated downstream component is
    # REPAIRED to the first valid choice; before that (cold start, unset
    # cascades such as the attr4 aesthetic selectors) no default selection
    # is forced, matching the historical behaviour
    ever_settled <- !is.null(isolate(committed()))

    # ---- analysis ----
    cc1 <- unique(vx[, 1])
    w1 <- .first_valid(c(
      if (fresh$analysis) sel$analysis,
      pd$analysis,
      live$analysis,
      sel$analysis), cc1) %||% cc1[1]

    # ---- subset ----
    cc2 <- unique(vx[vx[, 1] == w1, 2])
    w2 <- .first_valid(c(
      if (fresh$subset) sel$subset,
      pd$subset,
      live$subset,
      sel$subset), cc2)
    if (is.null(w2) && (ever_settled || !is.null(sel$subset)))
      w2 <- cc2[1]

    # ---- variable ----
    cc3 <- if (is.null(w2)) character() else
      vx[, 3][vx[, 1] == w1 & vx[, 2] == w2]
    w3 <- .first_valid(c(
      if (fresh$variable) sel$variable,
      pd$variable,
      live$variable,
      sel$variable), cc3)
    if (is.null(w3) && (ever_settled || !is.null(sel$variable)))
      w3 <- cc3[1]

    .update_widget("analysis", cc1, w1)
    .update_widget("subset", cc2, w2)
    .update_widget("variable", c("--select--", cc3), w3 %||% "--select--")
  }))


  # emit only coherent, settled triples; hold the last committed otherwise
  .tri_keep_obs(observe({
    if (!is.null(pend$analysis()) || !is.null(pend$subset()) ||
        !is.null(pend$variable()))
      return(NULL)
    vx <- tryCatch(validated_x(),
                   shiny.silent.error = function(e) NULL,
                   error = function(e) NULL)
    if (is.null(vx) || !nrow(vx))
      return(NULL)
    tr <- list(analysis = input$analysis, subset = input$subset,
               variable = input$variable)
    ok3 <- function(v) !is.null(v) && nzchar(v) && !identical(v, "--select--")
    if (!ok3(tr$analysis) || !ok3(tr$subset) || !ok3(tr$variable))
      return(NULL)
    if (!paste(tr$analysis, tr$subset, tr$variable, sep = "|") %in%
        paste(vx[, 1], vx[, 2], vx[, 3], sep = "|"))
      return(NULL)
    committed(tr)  # reactiveVal: identical triples do not re-invalidate
  }))


  reactive({
    committed()
  })




  }) # end moduleServer
}
