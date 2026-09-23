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
  reactive_status = reactive(NULL), store = NULL
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
    updateSelectInput(session, inputId = "scorner", choices = ac, selected = ps)
  })  

  # Volcano corner auto-selection. When the owning scatter detects volcano
  # axes, the "volcano" area (both top corners) is selected automatically.
  # Store-backed panels route the scorner change through the canonical store
  # so it becomes the single source of truth: a raw updateSelectInput here
  # used to be clobbered by an unrelated in-flight store re-assert on the
  # same widget (load-time bug: corner ended at "None" despite this branch
  # running). The direct updateSelectInput is kept as an immediate visual
  # sync - it carries the same value, and the input ack closes the loop.
  observeEvent(pre_volcano(), {
    corner <- if (isTRUE(pre_volcano())) "volcano" else "None"
    if (!is.null(store4))
      tryCatch(
        store_apply(store4, list(scorner = corner), origin = "system",
                    strict = FALSE, mark_pending = FALSE),
        error = function(e) NULL
      )
    if (isTRUE(pre_volcano())) {
      l <- list(x = val_xcut(), y = val_ycut(), corner = "volcano")
      attr(l, "seed") <- Sys.time()
      params$cutoff <- l
    } else {
      params$cutoff <- list(x = val_xcut(), y = val_ycut(), corner = "None")
    }
    updateSelectInput(session, inputId = "scorner", selected = corner)
  })
    
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
  # observeEvent(input$actSelect, {
  # observeEvent(input$scorner, {    
  observe({
    req( !is.null(input$scorner) && nchar(input$scorner) != 0 )

    acorner( input$scorner )
    i_xcut( input$xcut )
    i_ycut( input$ycut )
    l <- list(x = val_xcut(), y = val_ycut(), corner =  input$scorner)
    attr(l, "seed") <- Sys.time()
    params$cutoff <- l
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

  params

  }) # end moduleServer
}

