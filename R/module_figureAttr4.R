#' @description Utility - extended figure control shiny ui
#' @param id id
#' @param circle circle icon for dropdown manu
#' @importFrom shinyWidgets dropdown textInputIcon updateTextInputIcon updateSwitchInput
#' 
attr4selector_ui <- function(id, circle = TRUE, right = FALSE) {
  ns <- NS(id)
  dropdown(
    margin = "25px",
    circle = circle, right = right,
    status = "default", icon = icon("cog"), width = "788px",
    tooltip = tooltipOptions(title = "Click to modify figure!"),
    br(),
    triselector_ui(ns("selectColorUI")),
    triselector_ui(ns("selectShapeUI")),
    triselector_ui(ns("selectSizeUI")),
    triselector_ui(ns("selectTooltipUI")),
    triselector_ui(ns("selectSearchCol")),
    conditionalPanel(
      "1 == 2",
      checkboxInput(ns("showSearchBox"), label = "show", value = FALSE)
      ),
    conditionalPanel(
      "input.showSearchBox == true",
      ns = ns,    
      div(
        style='padding-left:100px; padding-right:0px; padding-top:0px; padding-bottom:0px',
        selectInput(ns("searchon"), label = NULL, choices = NULL, multiple = TRUE, width = "100%")
        )
      ),
    fluidRow(      
      column(2),
      column(
        4, offset = 0, style='padding:2px;', 
        textInputIcon(inputId = ns("xcut"), label = "Select points by x/y cutoffs", value = "log10(2)", placeholder = "e.g. -1 or -log10(2)", icon = list("x-cut"))),
      column(
        4, offset = 0, style='padding-left:5px; padding-right:2px; padding-top:27px; padding-bottom:2px;',         
        textInputIcon(inputId = ns("ycut"), label = NULL, value = "-log10(0.05)", placeholder = "e.g 2 or -log10(0.05)", icon = list("y-cut"))),
      column(
        2, offset = 0, style='padding-left:5px; padding-right:20px; padding-top:4px; padding-bottom:2px;', 
        selectInput(inputId = ns("scorner"), label = "Area", choices = "None", selectize = TRUE))
    )
  )
}

#' @description Utility - extended figure control shiny module
#' @param id module id
#' @param reactive_meta reactive meta info, usually phenotype data or feature data of ExpressoinSet
#' @param reactive_expr expression matrix
#' @param reactive_triset reactive value, a matrix of nx3 use for triselector
#' @param pre_volcano logical; whether select areas using volcano cutoff
#' @param reactive_status the status of scatter plot, e.g. color variable, shape variable, etc.
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for the owning module. When given, the
#'   panel's user-editable surface (the five attribute cascades plus the
#'   x/y cutoffs and the selected area) is registered under a nested
#'   \code{<prefix>.attr4.*} namespace and driven through the store
#'   protocol; when NULL the legacy status-restore path is kept.
#' @param corner_apply_gate Reactive logical, default \code{reactive(TRUE)}.
#'   The volcano corner auto-selection resolves only while the gate is TRUE;
#'   the owning scatter passes its "displayed axes have converged to the
#'   canonical store axes" reactive so the corner rectangles are never
#'   applied to an outgoing figure mid-switch.
#' @examples
#' #' # library(shiny)
#' # library(shinyjs)
#' # library(Biobase)
#' # library(shinyWidgets)
#' # source("Git/R/module_triselector.R")
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # pd <- pData(dat)
#' # fd <- fData(dat)
#' # expr <- exprs(dat)
#' # ts <- trisetter(meta = pd, expr = expr, combine = "pheno")
#' # # ts <- stringr::str_split_fixed(colnames(pd), pattern = "\\|", n = 3)
#' # ui <- fluidPage(
#' #   attr4selector_ui("a4test")
#' # )
#' # server <- function(input, output, session) {
#   k <- attr4selector_module("a4test", reactive_meta=reactive(pd), #' reactive_expr=reactive(expr), reactive_triset = reactive(ts))
#' #   # observe(
#' #   #   print(k())
#' #   # )
#' # }#'
#' # shinyApp(ui, server)
#'
attr4selector_module <- function(
  id, reactive_meta=reactive(NULL), reactive_expr=reactive(NULL),
  reactive_triset = reactive(NULL), pre_volcano = reactive(FALSE),
  reactive_status = reactive(NULL), store = NULL,
  corner_apply_gate = reactive(TRUE)
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  params <- reactiveValues(highlight = NULL, highlightName = NULL, color = NULL, shape = NULL, size = NULL, tooltips = NULL, cutoff = NULL)

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # The panel is shared by several modules; each store-backed owner passes
  # its child view and the keys register under <prefix>.attr4.*. The five
  # triselector cascades (color/shape/size/tooltip/search) are driven by
  # store_watch selectors; cutoffs and the selected area sync/push through
  # the standard observers. The hidden search-value select (searchon) is
  # deliberately NOT registered: its highlight wiring is dormant in the
  # current UI (searchValue is never fed from input$searchon), so it is
  # not a functional user-editable widget.
  # ------------------------------------------------------------------
  store4 <- NULL
  if (!is.null(store)) {
    .a4_root_store <- if (is.null(store$parent)) store else store$parent
    .a4_base <- if (is.null(store$parent)) "attr4" else paste0(store$prefix, ".attr4")
    store4 <- widget_store_child(.a4_root_store, .a4_base)
  }
  .a4_key <- function(k) paste0(store4$prefix, ".", k)

  # Defensive choice sources for the store validators: they read module
  # reactives WITHOUT req() (a shiny.validation condition inside a
  # choices_provider aborts the whole store_apply transaction) and mirror
  # exactly what the triselector cascade offers in the UI.
  .a4_ts <- function() {
    ts <- tryCatch(reactive_triset(),
                   shiny.silent.error = function(e) NULL,
                   error = function(e) NULL)
    if (is.null(ts) || !is.matrix(ts) || nrow(ts) == 0L) NULL else ts
  }
  .a4_lv1 <- function() {
    ts <- .a4_ts()
    if (is.null(ts)) character(0) else unique(ts[, 1])
  }
  .a4_lv2 <- function(a) {
    ts <- .a4_ts()
    if (is.null(ts) || is.null(a) || !nzchar(a)) character(0)
    else unique(ts[ts[, 1] %in% a, 2])
  }
  .a4_lv3 <- function(a, b) {
    ts <- .a4_ts()
    if (is.null(ts)) character(0)
    else {
      i <- rep(TRUE, nrow(ts))
      if (!is.null(a) && nzchar(a)) i <- i & ts[, 1] %in% a
      if (!is.null(b) && nzchar(b)) i <- i & ts[, 2] %in% b
      unique(ts[i, 3])
    }
  }
  .a4_corner_choices <- function(xcut, ycut) {
    x <- if (is.null(xcut)) NULL else
      suppressWarnings(tryCatch(text2num(xcut), error = function(e) NULL))
    y <- if (is.null(ycut)) NULL else
      suppressWarnings(tryCatch(text2num(ycut), error = function(e) NULL))
    if (is.numeric(x) && is.null(y)) c("None", "left", "right")
    else if (is.null(x) && is.numeric(y)) c("None", "top", "bottom")
    else if (is.numeric(x) && is.numeric(y))
      c("None", "volcano", "left", "right", "top", "bottom",
        "topleft", "topright", "bottomleft", "bottomright")
    else "None"
  }

  if (!is.null(store4)) {
    .a4_groups <- c(color = "Color", shape = "Shape", size = "Size",
                    tooltip = "Tooltips", search = "Search")
    .a4_binds <- unlist(lapply(names(.a4_groups), function(g) {
      ka <- paste0(g, "_analysis"); ks <- paste0(g, "_subset"); kv <- paste0(g, "_variable")
      lab <- .a4_groups[[g]]
      list(
        widget_binding(ka, "select",
          label = paste(lab, "attribute category"),
          help = paste0("First level of the ", tolower(lab),
                        " attribute cascade (annotation category)"),
          choices_provider = function(v) .a4_lv1()),
        widget_binding(ks, "select",
          label = paste(lab, "attribute subcategory"),
          help = paste0("Second level of the ", tolower(lab),
                        " attribute cascade"),
          depends_on = ka,
          choices_provider = function(v) .a4_lv2(v[[.a4_key(ka)]])),
        widget_binding(kv, "select_cascaded",
          label = paste(lab, "attribute variable"),
          help = paste0("Annotation variable mapped to the ", tolower(lab),
                        " aesthetic"),
          depends_on = c(ka, ks),
          choices_provider = function(v) .a4_lv3(v[[.a4_key(ka)]], v[[.a4_key(ks)]]))
      )
    }), recursive = FALSE)
    .a4_binds <- c(.a4_binds, list(
      widget_binding("xcut", "string", label = "X cutoff",
        help = paste("X-axis cutoff for area selection; a number or an R",
                     "expression such as -log10(0.05)")),
      widget_binding("ycut", "string", label = "Y cutoff",
        help = paste("Y-axis cutoff for area selection; a number or an R",
                     "expression such as log10(2)")),
      widget_binding("scorner", "select", label = "Selected area",
        help = "Plot region selected by the cutoffs (corner or volcano mode)",
        choices_provider = function(v)
          .a4_corner_choices(v[[.a4_key("xcut")]], v[[.a4_key("ycut")]]))))
    do.call(store_register, c(list(store4), .a4_binds))
  }

  selectColor_s1 <- reactiveVal()
  selectColor_s2 <- reactiveVal()
  selectColor_s3 <- reactiveVal()
  
  selectShape_s1 <- reactiveVal()
  selectShape_s2 <- reactiveVal()
  selectShape_s3 <- reactiveVal()
  
  selectSize_s1 <- reactiveVal()
  selectSize_s2 <- reactiveVal()
  selectSize_s3 <- reactiveVal()
  
  selectTooltip_s1 <- reactiveVal()
  selectTooltip_s2 <- reactiveVal()
  selectTooltip_s3 <- reactiveVal()
  
  searchOnCol_s1 <- reactiveVal()
  searchOnCol_s2 <- reactiveVal()
  searchOnCol_s3 <- reactiveVal()

  # store-backed owners drive the cascades from the canonical store;
  # store-less owners keep the restore-only reactive selectors
  .a4_make_selector <- function(ui_id, key, label, s1, s2, s3) {
    if (!is.null(store4))
      triselector_module(ui_id, reactive_x = reactive_triset, label = label,
        reactive_selector1 = store_watch(store4, paste0(key, "_analysis")),
        reactive_selector2 = store_watch(store4, paste0(key, "_subset")),
        reactive_selector3 = store_watch(store4, paste0(key, "_variable")),
        reactive_axis_request = store_epoch(store4))
    else
      triselector_module(ui_id, reactive_x = reactive_triset, label = label,
        reactive_selector1 = s1, reactive_selector2 = s2, reactive_selector3 = s3)
  }

  selectColor <- .a4_make_selector("selectColorUI", "color", "Color", selectColor_s1, selectColor_s2, selectColor_s3)
  selectShape <- .a4_make_selector("selectShapeUI", "shape", "Shape", selectShape_s1, selectShape_s2, selectShape_s3)
  selectSize <- .a4_make_selector("selectSizeUI", "size", "Size", selectSize_s1, selectSize_s2, selectSize_s3)
  selectTooltip <- .a4_make_selector("selectTooltipUI", "tooltip", "Tooltips", selectTooltip_s1, selectTooltip_s2, selectTooltip_s3)
  searchOnCol <- .a4_make_selector("selectSearchCol", "search", "Search", searchOnCol_s1, searchOnCol_s2, searchOnCol_s3)

  # Store glue (observer-GC rule: keep every observer referenced)
  if (!is.null(store4)) {
    .a4_store_observers <- list()
    .a4_keep <- function(obs) {
      .a4_store_observers[[length(.a4_store_observers) + 1L]] <<- obs
      invisible(obs)
    }
    .a4_read_tris <- function(sel)
      tryCatch(sel(), shiny.silent.error = function(e) NULL,
               error = function(e) NULL)
    .a4_component_set <- function(sel)
      !is.null(sel) &&
        nzchar(sel$analysis %||% "") && !identical(sel$analysis, "--select--") &&
        nzchar(sel$subset %||% "") && !identical(sel$subset, "--select--") &&
        nzchar(sel$variable %||% "") && !identical(sel$variable, "--select--")

    # UI -> store: user edits and widget acknowledgements
    .a4_sync_group <- function(sel, key) {
      .a4_keep(observe({
        cv <- .a4_read_tris(sel)
        if (.a4_component_set(cv)) {
          store_sync_from_ui(store4, paste0(key, "_analysis"), cv$analysis)
          store_sync_from_ui(store4, paste0(key, "_subset"), cv$subset)
          store_sync_from_ui(store4, paste0(key, "_variable"), cv$variable)
        }
      }))
    }
    .a4_sync_group(selectColor, "color")
    .a4_sync_group(selectShape, "shape")
    .a4_sync_group(selectSize, "size")
    .a4_sync_group(selectTooltip, "tooltip")
    .a4_sync_group(searchOnCol, "search")
    .a4_keep(observeEvent(input$xcut, {
      store_sync_from_ui(store4, "xcut", input$xcut)
    }, ignoreInit = TRUE))
    .a4_keep(observeEvent(input$ycut, {
      store_sync_from_ui(store4, "ycut", input$ycut)
    }, ignoreInit = TRUE))
    .a4_keep(observeEvent(input$scorner, {
      store_sync_from_ui(store4, "scorner", input$scorner)
    }, ignoreInit = TRUE))

    # Seed unset cutoff keys with the widget defaults once inputs exist
    # (restore-first-wins; the cascades are genuinely unset by default).
    # mark_pending = FALSE: these values are copied FROM the live inputs,
    # so the widgets already display them - arming pending/re-assert here
    # can never be acknowledged and would later clobber unrelated direct
    # widget updates (the load-time volcano corner bug).
    .a4_seeded <- FALSE
    .a4_keep(observe({
      if (.a4_seeded) return(NULL)
      if (is.null(input$xcut) || is.null(input$ycut) || is.null(input$scorner))
        return(NULL)
      .a4_seeded <<- TRUE
      vals <- list(xcut = input$xcut, ycut = input$ycut, scorner = input$scorner)
      held <- store_read(store4, names(vals))
      patch <- vals[vapply(names(vals), function(k)
        is.null(held[[paste0(store4$prefix, ".", k)]]), logical(1))]
      if (length(patch))
        tryCatch(store_apply(store4, patch, origin = "system", strict = FALSE,
                             mark_pending = FALSE),
                 error = function(e) NULL)
    }))

    # Store -> UI push for external writes only; the cascades re-derive from
    # the store_watch selectors themselves (meta_scatter pattern)
    .a4_epoch <- store_epoch(store4)
    .a4_keep(observe({
      .a4_epoch()
      vals <- store_read(store4, c("xcut", "ycut", "scorner"))
      pending <- .a4_root_store$pending
      if (!is.null(vals[[.a4_key("scorner")]]) &&
          !is.null(pending[[.a4_key("scorner")]]))
        updateSelectInput(session, "scorner", selected = vals[[.a4_key("scorner")]])
      for (k in c("xcut", "ycut"))
        if (!is.null(vals[[.a4_key(k)]]) && !is.null(pending[[.a4_key(k)]]))
          updateTextInputIcon(session, k, value = vals[[.a4_key(k)]])
    }))
  }


  vv <- reactive( varSelector(searchOnCol(), expr = reactive_expr(), meta = reactive_meta()) )

  pre_search <- reactiveVal()
  observe({    
    updateSelectInput(session, "searchon", choices = vv(), selected = pre_search())
  })
  observe(
    updateCheckboxInput(session, "showSearchBox", value = !is.null(vv()))
  )

  # Debounced text inputs to avoid unnecessary reactive updates while typing
  val_xcut <- reactive({ text2num(input$xcut) }) %>% debounce(500)
  val_ycut <- reactive({ text2num(input$ycut) }) %>% debounce(500)

  observeEvent(list(val_xcut(), val_ycut()), {    
    if (is.numeric(val_xcut()) && is.null(val_ycut())) {
      ac <- c("None", "left", "right")
    } else if (is.null(val_xcut()) && is.numeric(val_ycut())) {
      ac <- c("None", "top", "bottom")
    } else if (is.numeric(val_xcut()) && is.numeric(val_ycut())) {
      ac <- c("None", "volcano", "left", "right", "top", "bottom", "topleft", "topright", "bottomleft", "bottomright")      
    } else 
      ac <- "None"
    ps <- ac[1]
    if (!is.null(input$scorner) && input$scorner %in% ac) 
      ps <- input$scorner
    # never clobber a corner we just pushed and the browser has not
    # acknowledged yet (input$scorner still reports the pre-push value)
    pushed <- .a4_scorner_pushed()
    if (!is.null(pushed) && pushed %in% ac)
      ps <- pushed
    updateSelectInput(session, inputId = "scorner", choices = ac, selected = ps)
  })  

  # Volcano corner auto-selection, structured as INTENT + reactively
  # resolved corner:
  #
  # - pre_volcano() (the canonical, atomic volcano detection owned by the
  #   scatter module) only fires on genuine volcano-ness TRANSITIONS, so
  #   switching between two volcano views never touches the corner at all.
  # - The intended corner resolves through corner_effective() ONLY when
  #   corner_apply_gate() is TRUE - the scatter passes its "displayed axes
  #   have converged to the canonical axes" reactive, so the volcano
  #   rectangles are never applied to the outgoing figure mid-switch.
  # - The RENDER path (rectval in the owning module) consumes
  #   params$cutoff_reactive, so the resolved corner lands in the SAME
  #   reactive recompute as the new axes - one paint, rects included. An
  #   observer-based apply landed one flush later and painted the new axes
  #   with the stale corner first (a visible extra flash).
  # - Side effects (scorner widget, store, params$cutoff status mirror) ride
  #   a content-deduped observer: the historical per-input seed timestamp
  #   forced an invalidation on every event even when nothing changed, each
  #   one a redundant full plot redraw. That observer is the ONLY writer of
  #   params$cutoff: corner_effective deliberately does not read params (a
  #   reactiveValues write always re-invalidates its readers), and the
  #   generic input observer no longer mirrors input$scorner into it - two
  #   disagreeing writers plus the corner -> params feedback edge
  #   ping-ponged forever, and while the loop spun shiny never drained the
  #   outgoing message queue, so the browser could never ack the scorner
  #   update and the two writers never agreed (the load-time infinite loop).
  pendingCorner <- reactiveVal(NULL)
  # the corner value last pushed to the scorner widget and not yet
  # acknowledged by the browser; while it is set, a widget report that
  # still shows the pre-push value is stale, not a user choice (same
  # pattern as pendingAnalysis in the triselector module)
  .a4_scorner_pushed <- reactiveVal(NULL)
  widget_corner <- reactive({
    sc <- input$scorner
    pushed <- .a4_scorner_pushed()
    if (!is.null(pushed) && !identical(sc, pushed))
      return(pushed)
    sc %||% "None"
  })
  corner_effective <- reactive({
    intent <- pendingCorner()
    if (!is.null(intent) &&
        isTRUE(tryCatch(corner_apply_gate(), error = function(e) FALSE)))
      return(intent)
    widget_corner()
  })
  cutoff_effective <- reactive(
    list(x = val_xcut(), y = val_ycut(), corner = corner_effective())
  )
  .a4_cutoff_key <- reactiveVal(NULL)
  observe({
    l <- cutoff_effective()
    key <- paste(l$corner %||% "", l$x %||% "", l$y %||% "", sep = "\r")
    if (identical(key, .a4_cutoff_key()))
      return(NULL)
    .a4_cutoff_key(key)
    params$cutoff <- l
    if (!is.null(store4))
      tryCatch(
        store_apply(store4, list(scorner = l$corner), origin = "system",
                    strict = FALSE, mark_pending = FALSE),
        error = function(e) NULL
      )
    if (!identical(l$corner, isolate(input$scorner))) {
      .a4_scorner_pushed(l$corner)
      # choices ride along so the selected value always exists in the
      # selectize options even if this lands before the debounced choices
      # rebuild below (a bare setValue for a value not yet among the options
      # is silently dropped by selectize)
      updateSelectInput(session, inputId = "scorner",
        choices = .a4_corner_choices(l$x, l$y), selected = l$corner)
    }
  })
  observeEvent(pre_volcano(), {
    pendingCorner(if (isTRUE(pre_volcano())) "volcano" else "None")
  })
  # Any scorner widget report retires the auto-selection intent: either it
  # acknowledges our own push (the value then survives via widget_corner)
  # or it is a genuine user/restore choice, which always wins over the
  # automatic volcano corner.
  observeEvent(input$scorner, {
    .a4_scorner_pushed(NULL)
    pendingCorner(NULL)
  }, ignoreInit = TRUE)
    
  searchValue <- reactiveVal()
  observe({
    foo <- function() searchValue(input$searchon)
    debounce(foo, 1000) 
    })  

  observe({    
    if (is.null(vv())) {
      updateSelectInput(session, "searchon", choices = NULL, selected = NULL)
      searchValue(NULL)
      pre_search(NULL)
    }
    })
  observe({    
    req(vv())
    if ( is.null(searchValue()) )
      params$highlight <- NULL else
        params$highlight <- which(vv() %in% searchValue())
    isolate( params$highlightName <- searchOnCol()$variable )
  })

  observe(
    params$color <- varSelector(selectColor(), reactive_expr(), reactive_meta(), alternative = selectShape()$variable)
  )
  observe(
    params$shape <- varSelector(selectShape(), reactive_expr(), reactive_meta(), alternative = selectColor()$variable)
  )
  observe(
    params$size <- varSelector(selectSize(), reactive_expr(), reactive_meta())
  )
  observe(
    params$tooltips <- varSelector(selectTooltip(), reactive_expr(), reactive_meta())
  )

  acorner <- reactiveVal()    
  i_xcut <- reactiveVal()    
  i_ycut <- reactiveVal()    
  # Widget state mirror for the status snapshot ONLY (reactiveVal writes
  # are deduped). The cutoff/corner truth flows through cutoff_effective /
  # the side-effect observer above; writing params here as well raced that
  # observer with stale widget reports.
  observe({
    req( !is.null(input$scorner) && nchar(input$scorner) != 0 )

    acorner( input$scorner )
    i_xcut( input$xcut )
    i_ycut( input$ycut )
  })

  observe({
    params$status <- list(
      selectColor = selectColor(),
      selectShape = selectShape(),
      selectSize = selectSize(),
      selectTooltip = selectTooltip(),
      searchOnCol = searchOnCol(),
      searchValue = searchValue(),
      xcut = i_xcut(),
      ycut = i_ycut(),
      acorner = acorner()
      # ,
      # nClickCorner = input$actSelect
    )
  })
  ############### restore status ##############

  # observe({
  #   if (grepl("feature_space", ns("xx")))
  #     print(reactive_status())
  #   })

  # Store-backed owners restore through the canonical store (single
  # transactional path, per-key resilient); store-less owners keep the
  # legacy reactive-selector writes below.
  if (!is.null(store4)) {
    .a4_status_names <- c(color = "selectColor", shape = "selectShape",
                          size = "selectSize", tooltip = "selectTooltip",
                          search = "searchOnCol")
    .a4_keep(observe({
      if (is.null(s <- reactive_status()))
        return(NULL)
      patch <- list()
      for (g in names(.a4_status_names)) {
        tr <- s[[.a4_status_names[[g]]]]
        if (is.list(tr) && length(tr) == 3L &&
            !identical(tr$variable, "--select--") && nzchar(tr$variable %||% "")) {
          patch <- c(patch, stats::setNames(
            list(tr$analysis, tr$subset, tr$variable),
            paste0(g, c("_analysis", "_subset", "_variable"))))
        }
      }
      if (!is.null(s$xcut)) patch$xcut <- s$xcut
      if (!is.null(s$ycut)) patch$ycut <- s$ycut
      if (!is.null(s$acorner)) patch$scorner <- s$acorner
      if (length(patch))
        tryCatch(store_apply(store4, patch, origin = "restore", strict = FALSE),
                 error = function(e) NULL)
    }))
  } else {
  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)

    selectColor_s1( s$selectColor[[1]] )
    selectColor_s2( s$selectColor[[2]] )
    selectColor_s3( s$selectColor[[3]] )
    })
  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)
    selectShape_s1( s$selectShape[[1]] )
    selectShape_s2( s$selectShape[[2]] )
    selectShape_s3( s$selectShape[[3]] )
  })

  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)
    selectSize_s1( s$selectSize[[1]] )
    selectSize_s2( s$selectSize[[2]] )
    selectSize_s3( s$selectSize[[3]] )
  })

  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)
    selectTooltip_s1( s$selectTooltip[[1]] )
    selectTooltip_s2( s$selectTooltip[[2]] )
    selectTooltip_s3( s$selectTooltip[[3]] )
  })

  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)
    searchOnCol_s1( s$searchOnCol[[1]] )
    searchOnCol_s2( s$searchOnCol[[2]] )
    searchOnCol_s3( s$searchOnCol[[3]] )
  })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return(NULL)
    updateTextInputIcon(session, "xcut", value = s$xcut)
    updateTextInputIcon(session, "ycut", value = s$ycut) 
    updateSelectInput(session, "scorner", selected = s$acorner)    
  })
  }
  observe({
    if (is.null(s <- reactive_status()))
      return(NULL)
    pre_search(s$searchValue)
  })

  # Reactive cutoff for the owning module's render path: consuming this
  # (instead of the params$cutoff mirror written by observers) makes the
  # resolved corner land in the same reactive recompute as new axes.
  params$cutoff_reactive <- cutoff_effective

  params

  }) # end moduleServer
}

