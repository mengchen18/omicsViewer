#' Create a stable signature for a pair of scatter axes
#'
#' Selection emphasis is tied to the exact axes on which that selection was
#' made or restored. This prevents an opacity vector from being carried into a
#' different quick/custom figure.
#'
#' @keywords internal
.scatter_axis_signature <- function(x, y) {
  axis <- function(z) {
    if (is.null(z)) return(rep("", 3L))
    value <- unlist(z[c("analysis", "subset", "variable")], use.names = FALSE)
    if (length(value) != 3L) return(rep("", 3L))
    as.character(value)
  }
  paste(c(axis(x), axis(y)), collapse = "\r")
}

#' Meta Scatter Plot UI Function
#'
#' @description
#' Creates the user interface for the metadata scatter plot visualization module.
#' Provides interactive scatter plots with advanced selection tools including
#' corner selection for volcano plots.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{meta_scatter_module}}.
#'
#' @return
#' A \code{tagList} containing:
#' \itemize{
#'   \item Figure attribute selector (color, shape, size controls)
#'   \item Clear selection button
#'   \item A compact Quick view / Custom visualization mode switch
#'   \item Quick-view badges for common X/Y axis combinations
#'   \item X-axis and Y-axis variable selectors for custom visualizations
#'   \item Interactive plotly scatter plot with lasso/box selection
#' }
#'
#' @family visualization modules
#' @seealso
#' \code{\link{meta_scatter_module}} for the corresponding server logic.
#' \code{\link{plotly_scatter_module}} for the scatter plot implementation.
#'
#' @keywords internal
#' @importFrom shinyWidgets actionBttn radioGroupButtons updateRadioGroupButtons
#'
meta_scatter_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h3("Plot Controls and Variable Selection", class = "sr-only", `aria-label` = "Controls for customizing scatter plot appearance including color, shape, size mapping and selecting X and Y axis variables"),
    fluidRow(
      column(
        1,
        attr4selector_ui(ns("a4selector")),
        actionBttn(ns("clear"), "Clear figure selection", style = "minimal", color = "primary", size = "xs") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-clear-selection-button"))
      ), # style = "margin-top: 20px;",
      column(
        11,
        radioGroupButtons(
          inputId = ns("axisMode"),
          label = NULL,
          choices = c("Quick view" = "quick", "Custom visualization" = "custom"),
          selected = "quick",
          status = "primary",
          size = "sm",
          justified = TRUE,
          individual = FALSE,
          width = "100%"
        ),
        hr(style = "margin: 8px 0;"),
        tabsetPanel(
          id = ns("axisModeTabs"),
          type = "hidden",
          tabPanelBody(
            value = "quick",
            quick_badges_ui(ns("quickViews"))
          ),
          tabPanelBody(
            value = "custom",
            triselector_ui(ns("tris_main_scatter1")),
            triselector_ui(ns("tris_main_scatter2"))
          )
        )
      )
    ),
    tags$h3("Interactive Scatter Plot Visualization", class = "sr-only", `aria-label` = "Scatter plot with lasso and box selection tools, regression line option, and corner selection for volcano plots"),
    plotly_scatter_ui(ns("main_scatterOutput"), height = META_SCATTER_PLOT_HEIGHT),
    # Hidden text summary for AI browsers and screen readers
    div(
      class = "sr-only", `aria-live` = "polite", `aria-atomic` = "true",
      uiOutput(ns("plotSummary"))
    )
  )
}

#' @title Utility - scatter plot for meta shiny module
#' @param id Character. Namespace ID for the Shiny module.
#' @param reactive_meta reactive meta data, phenotype data or feature data
#' @param reactive_expr reactive expression data
#' @param combine how to combine the expression and meta data, pheno or feature?
#' @param source source id for plotly object
#' @param reactive_x reactive value for pre-selected x-aixs
#' @param reactive_y reactive value for pre-selected y-aixs
#' @param reactive_status the status of scatter plot, e.g. x-, y-axis, color variable, shape variable, etc.
#' @param store Child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this scatter space. Axis and
#'   mode state lives in the store; the hand-rolled xax/yax/axisRequest
#'   sync machinery was replaced by the store protocol (plan section 6,
#'   phase S2).
#'
meta_scatter_module <- function(
  id, reactive_meta = reactive(NULL), reactive_expr = reactive(NULL),
  combine = c("pheno", "feature"), source = "plotlyscattersource",
  reactive_x = reactive(NULL), reactive_y = reactive(NULL),
  reactive_status = reactive(NULL),
  store = NULL
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    notNullAndPositiveLength <- function(x) !is.null(x) && length(x) > 0

    # Helper: Get feature/sample names based on combine mode
    get_names <- function() {
      if (combine == "pheno") {
        colnames(reactive_expr())
      } else {
        rownames(reactive_expr())
      }
    }

    triset <- reactive({
      ts <- trisetter(expr = reactive_expr(), meta = reactive_meta(), combine = combine[1])
      ts[ts[, 1] != "Surv", ]
    })

    quick_views <- reactive({
      prepare_quick_views(reactive_meta(), triset())
    })

    activeQuickView <- reactive({
      active_quick_view(quick_views(), v1(), v2())
    })

    quickBadge <- quick_badges_module(
      "quickViews",
      views = quick_views,
      activeId = activeQuickView
    )

    # ------------------------------------------------------------------
    # Canonical widget-store bindings (control plane, plan section 6).
    # The store is the single source of truth for the axis and mode state;
    # the triselector cascade pushes store values to the UI, and the sync
    # observer below mirrors user edits (and acknowledgements) back.
    # ------------------------------------------------------------------
    stopifnot(!is.null(store))
    kx1 <- paste0(store$prefix, ".x_analysis")
    kx2 <- paste0(store$prefix, ".x_subset")
    kx3 <- paste0(store$prefix, ".x_variable")
    ky1 <- paste0(store$prefix, ".y_analysis")
    ky2 <- paste0(store$prefix, ".y_subset")
    ky3 <- paste0(store$prefix, ".y_variable")
    kmode <- paste0(store$prefix, ".axis_mode")
    store_register(
      store,
      widget_binding("x_analysis", "select", label = "X-axis analysis",
        help = "Analysis category for the X axis",
        choices_provider = function(v) unique(triset()[, 1])),
      widget_binding("x_subset", "select", label = "X-axis subset",
        help = "Subset within the X-axis analysis", depends_on = "x_analysis",
        choices_provider = function(v)
          unique(triset()[triset()[, 1] %in% v[[kx1]], 2])),
      widget_binding("x_variable", "select_cascaded", label = "X-axis variable",
        help = "Annotation variable plotted on the X axis",
        depends_on = c("x_analysis", "x_subset"),
        choices_provider = function(v)
          triset()[triset()[, 1] %in% v[[kx1]] & triset()[, 2] %in% v[[kx2]], 3]),
      widget_binding("y_analysis", "select", label = "Y-axis analysis",
        help = "Analysis category for the Y axis",
        choices_provider = function(v) unique(triset()[, 1])),
      widget_binding("y_subset", "select", label = "Y-axis subset",
        help = "Subset within the Y-axis analysis", depends_on = "y_analysis",
        choices_provider = function(v)
          unique(triset()[triset()[, 1] %in% v[[ky1]], 2])),
      widget_binding("y_variable", "select_cascaded", label = "Y-axis variable",
        help = "Annotation variable plotted on the Y axis",
        depends_on = c("y_analysis", "y_subset"),
        choices_provider = function(v)
          triset()[triset()[, 1] %in% v[[ky1]] & triset()[, 2] %in% v[[ky2]], 3]),
      widget_binding("axis_mode", "enum", label = "Axis mode",
        help = "Quick-view badges or custom triselectors",
        values = c("quick", "custom"))
    )

    # Reactive axis watchers: the canonical, atomic view of the current
    # axes (used by the volcano detection below). Created once; each call
    # inside a reactive registers the per-key dependency.
    store_watchers <- list(
      x_analysis = store_watch(store, "x_analysis"),
      x_subset = store_watch(store, "x_subset"),
      x_variable = store_watch(store, "x_variable"),
      y_analysis = store_watch(store, "y_analysis"),
      y_subset = store_watch(store, "y_subset"),
      y_variable = store_watch(store, "y_variable")
    )

    # Seed the store with the dataset's default axes whenever the defaults
    # change (initial load, dataset reload). A snapshot/agent restore that
    # lands first wins: seeding skips keys the store already holds.
    .scatter_axis_triple <- function(axis_string) {
      if (is.null(axis_string)) return(NULL)
      l <- strsplit(axis_string, "\\|")[[1]]
      if (length(l) != 3L || any(!nzchar(l))) return(NULL)
      l
    }
    # Store-glue observers are kept referenced: observers whose dependencies
    # are only weakly held by the reactive graph are garbage collected
    # between flushes, silently killing later sync/seed/restore/push work
    # (found while hardening the S3 generic tier; explains sporadic
    # restore/stress flakiness).
    .scatter_store_observers <- list()
    .scatter_keep <- function(obs) {
      .scatter_store_observers[[length(.scatter_store_observers) + 1L]] <<- obs
      invisible(obs)
    }
    last_seeded <- character()
    .scatter_keep(observe({
      dx <- reactive_x()
      dy <- reactive_y()
      req(nrow(ts <- triset()) > 0)
      stamp <- c(dx %||% "", dy %||% "")
      if (identical(stamp, last_seeded))
        return(NULL)
      last_seeded <<- stamp
      patch <- list()
      tx <- .scatter_axis_triple(dx)
      ty <- .scatter_axis_triple(dy)
      vals <- store_read(store, c("x_analysis", "y_analysis"))
      if (!is.null(tx) && is.null(vals[[kx1]]))
        patch <- c(patch, stats::setNames(as.list(tx), c("x_analysis", "x_subset", "x_variable")))
      if (!is.null(ty) && is.null(vals[[ky1]]))
        patch <- c(patch, stats::setNames(as.list(ty), c("y_analysis", "y_subset", "y_variable")))
      if (!length(patch))
        return(NULL)
      tryCatch(store_apply(store, patch, origin = "system", strict = FALSE),
               error = function(e) NULL)
    }))

    v1 <- triselector_module("tris_main_scatter1",
      reactive_x = triset, label = "X-axis",
      reactive_selector1 = store_watch(store, "x_analysis"),
      reactive_selector2 = store_watch(store, "x_subset"),
      reactive_selector3 = store_watch(store, "x_variable"),
      reactive_axis_request = store_epoch(store)
    )
    v2 <- triselector_module("tris_main_scatter2",
      reactive_x = triset, label = "Y-axis",
      reactive_selector1 = store_watch(store, "y_analysis"),
      reactive_selector2 = store_watch(store, "y_subset"),
      reactive_selector3 = store_watch(store, "y_variable"),
      reactive_axis_request = store_epoch(store)
    )

    # UI -> store synchronisation: user edits and widget acknowledgements.
    # store_sync_from_ui is acknowledgement-aware: confirming an in-flight
    # external write clears its pending entry; a diverging value is a user
    # override (user always wins). Reads are req-guarded, hence tryCatch.
    .scatter_read_tris <- function(sel) {
      tryCatch(sel(), shiny.silent.error = function(e) NULL,
               error = function(e) NULL)
    }
    .scatter_component_set <- function(sel) {
      # triselectors report the "--select--" placeholder while unset; that is
      # "no value", not a user choice, and must not enter the store
      !is.null(sel) &&
        nzchar(sel$analysis %||% "") && !identical(sel$analysis, "--select--") &&
        nzchar(sel$subset %||% "") && !identical(sel$subset, "--select--") &&
        nzchar(sel$variable %||% "") && !identical(sel$variable, "--select--")
    }
    .scatter_triple_coherent <- function(sel) {
      # A coherent triple names a real column in the current triset. While a
      # cascaded select catches up with a store push, the triselector briefly
      # reports MIXED components (new analysis + old subset/variable); those
      # are cascade echoes, not user choices, and must NOT be mirrored into
      # the store - the eager user-override there cleared the pending entries
      # of an in-flight quick-view apply and made the canonical axes (and the
      # volcano corner driven from them) oscillate after every switch.
      if (!.scatter_component_set(sel))
        return(FALSE)
      triple <- paste(sel$analysis, sel$subset, sel$variable, sep = "|")
      ts <- tryCatch(triset(), shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.null(ts) || !nrow(ts))
        return(FALSE)
      triple %in% paste(ts[, 1], ts[, 2], ts[, 3], sep = "|")
    }
    .scatter_keep(observe({
      xv <- .scatter_read_tris(v1)
      yv <- .scatter_read_tris(v2)
      if (.scatter_triple_coherent(xv)) {
        store_sync_from_ui(store, "x_analysis", xv$analysis)
        store_sync_from_ui(store, "x_subset", xv$subset)
        store_sync_from_ui(store, "x_variable", xv$variable)
      }
      if (.scatter_triple_coherent(yv)) {
        store_sync_from_ui(store, "y_analysis", yv$analysis)
        store_sync_from_ui(store, "y_subset", yv$subset)
        store_sync_from_ui(store, "y_variable", yv$variable)
      }
    }))
    .scatter_keep(observeEvent(input$axisMode, {
      updateTabsetPanel(session, "axisModeTabs", selected = input$axisMode)
      if (input$axisMode %in% c("quick", "custom"))
        store_sync_from_ui(store, "axis_mode", input$axisMode)
    }, ignoreInit = TRUE))

    # Store -> UI push for the axis mode radio: external writes only (a
    # pending entry marks them); user clicks sync back above. The epoch
    # reactive is created once and the observer is kept referenced: an
    # observer whose dependencies are only weakly held by the reactive
    # graph is garbage collected between flushes, silently killing later
    # pushes (same retention rule as the heatmap store glue).
    axisModeRoot <- if (is.null(store$parent)) store else store$parent
    axis_mode_epoch <- store_epoch(store)
    .scatter_store_observers$axis_mode_push <- observe({
      axis_mode_epoch()
      mode <- store_read(store, "axis_mode")[[1]]
      if (is.null(mode)) return(NULL)
      if (is.null(axisModeRoot$pending[[kmode]])) return(NULL)
      updateRadioGroupButtons(
        session, "axisMode",
        choices = c("Quick view" = "quick", "Custom visualization" = "custom"),
        selected = mode
      )
      updateTabsetPanel(session, "axisModeTabs", selected = mode)
    })

    # Plotly owns the visible box/lasso immediately after a user selection. We
    # therefore do not make selection emphasis a reactive dependency of the
    # plot. These values allow a later redraw to emphasize only the selection
    # belonging to the current axes; changing axes clears that emphasis.
    selectionDisplayAxes <- reactiveVal(NULL)
    pendingSelectionDisplayAxes <- reactiveVal(NULL)
    selectionDisplayTrigger <- reactiveVal(0L)
    observe({
      current_axes <- .scatter_axis_signature(v1(), v2())
      pending_axes <- pendingSelectionDisplayAxes()
      if (!is.null(pending_axes)) {
        if (identical(current_axes, pending_axes)) {
          selectionDisplayAxes(pending_axes)
          pendingSelectionDisplayAxes(NULL)
        }
        return(NULL)
      }

      displayed_axes <- isolate(selectionDisplayAxes())
      if (!is.null(displayed_axes) && !identical(current_axes, displayed_axes))
        selectionDisplayAxes(NULL)
    })

    .scatter_keep(observeEvent(quickBadge()$trigger, {
      qv <- quickBadge()$view
      req(nrow(qv) == 1)
      req(xx <- .quick_view_axis(qv$x))
      req(yy <- .quick_view_axis(qv$y))

      store_apply(store, c(
        stats::setNames(as.list(xx), c("x_analysis", "x_subset", "x_variable")),
        stats::setNames(as.list(yy), c("y_analysis", "y_subset", "y_variable")),
        list(axis_mode = "quick")
      ), origin = "system")
    }, ignoreInit = TRUE))

    # Keep the compact mode switch and the header-less tab panel in sync. The
    # tab panel switches content immediately without adding another tab bar.
    observeEvent(input$axisMode, {
      updateTabsetPanel(session, "axisModeTabs", selected = input$axisMode)
    }, ignoreInit = TRUE)

    # Detect volcano plot: x=mean.diff, y=log.fdr/log.pvalue (both from
    # ttest). Two invariants keep the corner auto-selection from flashing:
    # (1) the axes are read from the CANONICAL store watchers - the store
    # transaction lands atomically while the v1()/v2() triselector reports
    # lag and transiently revert during the cascade, so a v1/v2-based
    # detection oscillated FALSE->TRUE->FALSE->TRUE and re-fired the corner
    # chain on every oscillation; (2) the detection only turns TRUE once the
    # DISPLAYED axes have converged to the canonical axes - otherwise the
    # volcano corner (and its rectangles) was applied to the outgoing
    # figure before the axis switch landed (a visible pre-flash on the old
    # plot, e.g. volcano rectangles drawn over a correlation plot).
    pre_vol <- reactive({
      vals <- lapply(store_watchers, function(w) w())
      if (any(vapply(vals, is.null, logical(1))))
        return(FALSE)
      if (!identical(vals$x_analysis, "ttest") ||
          !identical(vals$y_analysis, "ttest") ||
          !identical(vals$x_variable, "mean.diff") ||
          !vals$y_variable %in% c("log.fdr", "log.pvalue"))
        return(FALSE)
      xv <- .scatter_read_tris(v1)
      yv <- .scatter_read_tris(v2)
      .scatter_component_set(xv) && .scatter_component_set(yv) &&
        identical(xv$analysis, vals$x_analysis) &&
        identical(xv$subset, vals$x_subset) &&
        identical(xv$variable, vals$x_variable) &&
        identical(yv$analysis, vals$y_analysis) &&
        identical(yv$subset, vals$y_subset) &&
        identical(yv$variable, vals$y_variable)
    })

    attr4select_status <- reactiveVal()
    attr4select <- attr4selector_module(
      "a4selector",
      reactive_meta = reactive_meta, reactive_expr = reactive_expr,
      reactive_triset = triset, pre_volcano = pre_vol, reactive_status = attr4select_status,
      store = store
    )

    xycoord <- reactive({
      req(v1()$variable)
      req(v2()$variable)
      req(!v1()$variable %in% c("--select--", ""))
      req(!v2()$variable %in% c("--select--", ""))
      x <- varSelector(v1(), reactive_expr(), reactive_meta())
      y <- varSelector(v2(), reactive_expr(), reactive_meta())
      req(x)
      req(y)
      req(is.numeric(x) || is.numeric(y))
      req(length(x) == length(y))
      list(x = x, y = y)
    })

    # Track clear button clicks
    clear_counter <- reactiveVal(0)
    observeEvent(input$clear, {
      clear_counter(clear_counter() + 1)
    })

    # Rectangle for corner selection (volcano plot)
    rectval <- reactive({
      # Recalculate when clear button clicked
      clear_counter()

      # Return NULL if no cutoff selected or "None" corner
      if (is.null(attr4select$cutoff) || attr4select$cutoff$corner == "None") {
        return(NULL)
      }

      # Force dependency on axis changes
      v1()
      v2()

      # Calculate rectangle based on coordinates and cutoff
      # Note: xycoord() already depends on v1() and v2(), so we don't need
      # to touch them explicitly - reactive graph handles transitive dependencies
      coords <- xycoord()
      if (is.null(coords)) {
        return(NULL)
      }

      line_rect(l = attr4select$cutoff, coords)$rect
    })

    scatter_vars <- reactive({
      req(l <- xycoord())
      l$source <- source
      l$xlab <- attr(l$x, "label")
      l$ylab <- attr(l$y, "label")
      l$color <- attr4select$color
      l$shape <- attr4select$shape
      l$size <- attr4select$size
      l$tooltips <- attr4select$tooltips
      l$highlight <- attr4select$highlight
      l$highlightName <- attr4select$highlightName
      l$rect <- rectval()
      # A restoration on unchanged axes needs one deliberate redraw. Ordinary
      # selection emphasis is deliberately isolated so a lasso/box event does
      # not immediately erase Plotly's browser-owned selection shape. When the
      # axes change, the triselector outputs already invalidate this reactive.
      selectionDisplayTrigger()

      # Do not carry an opacity vector from one figure into another. Emphasize
      # semantic IDs only when the current axes match either a restored axis
      # pair or the axis pair on which the user made the selection.
      current_axes <- .scatter_axis_signature(v1(), v2())
      restored_axes <- isolate(pendingSelectionDisplayAxes())
      displayed_axes <- isolate(selectionDisplayAxes())
      if (!is.null(restored_axes) && identical(restored_axes, current_axes)) {
        displayed_axes <- restored_axes
      }
      selected_ids <- if (
        !is.null(displayed_axes) && identical(displayed_axes, current_axes)
      ) {
        isolate(selVal()$selected)
      } else {
        character(0)
      }
      l$inSelection <- if (notNullAndPositiveLength(selected_ids)) {
        which(get_names() %in% selected_ids)
      } else {
        NA
      }
      l
    })

    showRegLine <- reactiveVal(FALSE)
    htestV1 <- reactiveVal()
    htestV2 <- reactiveVal()
    v_scatter <- plotly_scatter_module(
      "main_scatterOutput",
      reactive_param_plotly_scatter = scatter_vars,
      reactive_regLine = showRegLine, htest_var1 = htestV1, htest_var2 = htestV2
    )
    observe({
      showRegLine(v_scatter()$regline)
    })

    selVal <- reactiveVal(
      list(
        clicked = character(0),
        selected = character(0)
      )
    )
    sbc <- reactiveVal(FALSE)

    observeEvent(list(input$clear, reactive_expr()), {
      selVal(list(
        clicked = character(0),
        selected = character(0)
      ))
      selectionDisplayAxes(NULL)
      pendingSelectionDisplayAxes(NULL)
      sbc(FALSE)
    })

    # Workaround: Track previous selection to prevent redundant updates
    # Plotly events can fire even when selection hasn't actually changed,
    # causing unnecessary reactive chain invalidations. We store the previous
    # selection and only update selVal when it truly changes.
    clientSideSelection <- reactiveVal(character(0))
    observeEvent(v_scatter(), {
      l <- get_names()
      u_c <- l[v_scatter()$clicked]
      u_s <- l[v_scatter()$selected]

      # Only update if selection actually changed
      req(!identical(tmp <- c(u_c, u_s), clientSideSelection()))

      clientSideSelection(tmp)
      selVal(list(
        clicked = u_c,
        selected = u_s
      ))
      selectionDisplayAxes(
        if (notNullAndPositiveLength(tmp))
          .scatter_axis_signature(v1(), v2())
        else
          NULL
      )
      sbc(FALSE)
    })

    returnCornerSelection <- reactiveVal(TRUE)
    observeEvent(rectval(), {
      if (!returnCornerSelection()) {
        return(NULL)
      }

      rec <- rectval()
      if (is.null(rec)) {
        selVal(list(
          clicked = character(0),
          selected = character(0)
        ))
        selectionDisplayAxes(NULL)
        return(NULL)
      }
      req(cc <- xycoord())
      l <- get_names()

      i <- lapply(rec, function(r1) {
        which(cc$x > r1["x0"] & cc$x < r1["x1"] & cc$y > r1["y0"] & cc$y < r1["y1"])
      })
      i <- sort(unique(unlist(i)))
      selVal(list(
        clicked = character(0),
        selected = l[i]
      ))
      selectionDisplayAxes(.scatter_axis_signature(v1(), v2()))
      sbc(TRUE)
    })

    ############## status save ###############
    # Derive snapshot state when it is read. Besides avoiding the historical
    # circular reactiveVal update, this ensures a Save click evaluates the
    # current axis/selector controls rather than an observer's last write.
    scatter_status <- reactive({
      safe_state_value <- function(value) {
        tryCatch(value, shiny.silent.error = function(e) NULL,
                 error = function(e) NULL)
      }
      current <- isolate(selVal())
      vals <- store_read(store, c("x_analysis", "x_subset", "x_variable",
                                   "y_analysis", "y_subset", "y_variable"))
      list(
        axisMode = input$axisMode,
        xax = list(v1 = vals[[kx1]], v2 = vals[[kx2]], v3 = vals[[kx3]]),
        yax = list(v1 = vals[[ky1]], v2 = vals[[ky2]], v3 = vals[[ky3]]),
        showRegLine = showRegLine(),
        attr4 = safe_state_value(attr4select$status),
        selection_clicked = current$clicked,
        selection_selected = current$selected,
        selectByCorner = sbc()
      )
    })

    ############## status restore ###############
    # Consolidate all status restoration into single observer. Axis and mode
    # state goes through the canonical store (transactional, diff-only),
    # which replaces the former direct radio updates and xax/yax writes:
    # keys absent from the status are never touched, so restores have no
    # side effects beyond what the snapshot recorded.
    .scatter_keep(observeEvent(reactive_status(), {
      s <- reactive_status()
      if (is.null(s)) {
        return()
      }

      restored_axes <- .scatter_axis_signature(s$xax, s$yax)
      # v1()/v2() are req(input$variable)-guarded; while a triselector cascade
      # is still in flight they abort. current_axes only feeds a redraw-trigger
      # comparison, so NULL on failure is harmless.
      current_axes <- tryCatch(
        .scatter_axis_signature(isolate(v1()), isolate(v2())),
        error = function(e) NULL
      )
      selectionDisplayAxes(NULL)
      pendingSelectionDisplayAxes(restored_axes)

      patch <- list()
      if (identical(length(s$xax), 3L))
        patch <- c(patch, stats::setNames(as.list(s$xax),
                 c("x_analysis", "x_subset", "x_variable")))
      if (identical(length(s$yax), 3L))
        patch <- c(patch, stats::setNames(as.list(s$yax),
                 c("y_analysis", "y_subset", "y_variable")))
      if (!is.null(s$axisMode) && s$axisMode %in% c("quick", "custom"))
        patch$axis_mode <- s$axisMode
      tryCatch(store_apply(store, patch, origin = "restore", strict = FALSE),
               error = function(e) NULL)
      if (identical(restored_axes, current_axes))
        selectionDisplayTrigger(isolate(selectionDisplayTrigger()) + 1L)

      # Restore attribute selector status
      attr4select_status(NULL)
      attr4select_status(s$attr4)

      # Restore regression line setting
      showRegLine(s$showRegLine)

      # Computed hypothesis-test output is intentionally not restored; it is
      # recalculated from the restored widget and selection state.
      returnCornerSelection(s$selectByCorner)
      selVal(list(
        clicked = s$selection_clicked,
        selected = s$selection_selected
      ))
    }))
    #############################################

    # Generate hidden text summary for AI browsers and screen readers
    output$plotSummary <- renderUI({
      req(scatter_vars())
      vars <- scatter_vars()

      # Calculate basic statistics
      x_vals <- vars$x
      y_vals <- vars$y
      n_points <- length(x_vals)

      # Build summary text
      summary_parts <- c(
        sprintf("Scatter plot visualization with %d data points.", n_points),
        sprintf("X-axis: %s.", vars$xlab %||% "Variable"),
        sprintf("Y-axis: %s.", vars$ylab %||% "Variable")
      )

      # Add numeric range info for numeric axes
      if (is.numeric(x_vals)) {
        summary_parts <- c(
          summary_parts,
          sprintf("X-axis range: %.3g to %.3g.", min(x_vals, na.rm = TRUE), max(x_vals, na.rm = TRUE))
        )
      }
      if (is.numeric(y_vals)) {
        summary_parts <- c(
          summary_parts,
          sprintf("Y-axis range: %.3g to %.3g.", min(y_vals, na.rm = TRUE), max(y_vals, na.rm = TRUE))
        )
      }

      # Add correlation for numeric-numeric plots
      if (is.numeric(x_vals) && is.numeric(y_vals)) {
        cor_result <- tryCatch(
          {
            cor.test(x_vals, y_vals, use = "complete.obs")
          },
          error = function(e) NULL
        )

        if (!is.null(cor_result)) {
          summary_parts <- c(
            summary_parts,
            sprintf(
              "Pearson correlation: r = %.3f, p-value = %.3g.",
              cor_result$estimate, cor_result$p.value
            )
          )
        }
      }

      # Add selection info
      sel <- selVal()
      if (length(sel$selected) > 0) {
        summary_parts <- c(
          summary_parts,
          sprintf("%d points currently selected.", length(sel$selected))
        )
      }

      tags$p(paste(summary_parts, collapse = " "))
    })

    reactive({
      current <- selVal()
      sta <- scatter_status()
      # Keep an explicit list member as a robust fallback; some reactive
      # consumers historically lost attributes when forwarding module values.
      current$state <- sta
      attr(current, "status") <- sta
      # Derived quick views are not snapshot state, but the optional AI
      # assistant uses this metadata to describe and validate one-click views.
      attr(current, "quickViews") <- quick_views()
      current
    })
  }) # end moduleServer
}
