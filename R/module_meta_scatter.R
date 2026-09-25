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
#' @param selection Selection-bus port (\code{\link{selection_port}})
#'   bound to this scatter's space ("feature" or "sample"). The scatter
#'   is the figure-space writer of the unified selection: lasso/box/click
#'   events, the volcano corner auto-selection, clear and snapshot
#'   restore all report through the port (origins "figure", "corner",
#'   "clear", "restore"), and every consumer (tables, result space,
#'   snapshot) reads the bus instead of adopting module returns.
#'
meta_scatter_module <- function(
  id, reactive_meta = reactive(NULL), reactive_expr = reactive(NULL),
  combine = c("pheno", "feature"), source = "plotlyscattersource",
  reactive_x = reactive(NULL), reactive_y = reactive(NULL),
  reactive_status = reactive(NULL),
  store = NULL,
  selection = NULL
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    stopifnot(!is.null(store), !is.null(selection))
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
      tx <- .scatter_axis_triple(dx)
      ty <- .scatter_axis_triple(dy)
      vals <- list()
      if (!is.null(tx))
        vals <- c(vals, stats::setNames(as.list(tx),
                   c("x_analysis", "x_subset", "x_variable")))
      if (!is.null(ty))
        vals <- c(vals, stats::setNames(as.list(ty),
                   c("y_analysis", "y_subset", "y_variable")))
      if (length(vals))
        store_seed(store, vals)
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

    # UI -> store synchronisation: the binding mirrors the SETTLED triples
    # (the triselector holds its last committed triple while a cascade is
    # in flight), so mid-cascade echoes never revert an in-flight store
    # write, and acknowledgements of store pushes clear their pending
    # entries through the ack branch of store_sync_from_ui.
    .scatter_read_tris <- function(sel) {
      tryCatch(sel(), shiny.silent.error = function(e) NULL,
               error = function(e) NULL)
    }
    store_bind_triselector(store,
      keys = c(analysis = "x_analysis", subset = "x_subset", variable = "x_variable"),
      sel = v1, keep = .scatter_keep)
    store_bind_triselector(store,
      keys = c(analysis = "y_analysis", subset = "y_subset", variable = "y_variable"),
      sel = v2, keep = .scatter_keep)
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

    # ------------------------------------------------------------------
    # Selection display contract (unified).
    #
    # Emphasis (the inSelection opacity vector) is resolved IN-RENDER,
    # from the current display authority:
    #
    #   - "corner": the volcano corner owns the display. The emphasized
    #     ids are recomputed from the LIVE rectangles and coordinates, so
    #     the emphasis always lands in the SAME paint as an axis or cutoff
    #     change (a selection written by an observer lands one queue
    #     position after the render consumers and would otherwise paint a
    #     figure whose emphasis contradicts the propagated selection - the
    #     log.fdr -> log.pvalue switch left the corner genes selected in
    #     the tables but emphasized nothing).
    #   - "browser"/"restore": the anchor path - Plotly owns the visible
    #     box/lasso immediately after a user selection, so a selection is
    #     never a reactive dependency of the plot (a replot would erase
    #     the browser-owned selection shape); emphasis comes from the
    #     settled selVal, gated on the axes it was made/restored on.
    #   - "none": no emphasis (cleared).
    #
    # cornerAuthority is read ISOLATED: an authority flip alone must not
    # repaint (it would erase a browser-owned lasso). Writers that need a
    # repaint Plotly does not own (clear, restore on unchanged axes, the
    # first corner arm) bump selectionDisplayTrigger, whose seed rides
    # the params through an attribute the render barrier's identical()
    # compares but do.call drops.
    #
    # cornerEngaged (propagation gate) is a SEPARATE concept: whether
    # corner changes may CLAIM the selection at all. Disengaged by clear
    # and by a manual-selection restore; re-engaged by any genuine
    # cutoff/corner edit. The historic returnCornerSelection/sbc pair
    # collapsed into these two flags.
    selectionDisplayAxes <- reactiveVal(NULL)
    pendingSelectionDisplayAxes <- reactiveVal(NULL)
    selectionDisplayTrigger <- reactiveVal(0L)
    cornerEngaged <- reactiveVal(TRUE)
    cornerAuthority <- reactiveVal(FALSE)
    .scatter_ids_in_rects <- function(coords, rects) {
      i <- lapply(rects, function(r1) {
        which(coords$x > r1["x0"] & coords$x < r1["x1"] &
              coords$y > r1["y0"] & coords$y < r1["y1"])
      })
      sort(unique(unlist(i)))
    }
    .scatter_keep(observe({
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
    }))

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
    .scatter_keep(observeEvent(input$axisMode, {
      updateTabsetPanel(session, "axisModeTabs", selected = input$axisMode)
    }, ignoreInit = TRUE))

    # Detect volcano plot: x=mean.diff, y=log.fdr/log.pvalue (both from
    # ttest), read from the DISPLAYED (settled) triselector triples. The
    # store watchers lag one observer hop behind a user edit (the store
    # sync runs after the committed flip), which armed the volcano-exit
    # intent after the render had already painted the new axes with the
    # outgoing corner (two paints on every volcano exit edit). The
    # committed triples flip ATOMICALLY (the triselector holds its last
    # settled triple mid-cascade), so this predicate still transitions
    # exactly once per view change and never dips between two volcano
    # views - the property the historical store-watcher version was
    # written to guarantee. Timing of the corner APPLICATION is handled
    # separately, by the convergence gate passed to attr4selector.
    displayed_volcano <- reactive({
      xv <- .scatter_read_tris(v1)
      yv <- .scatter_read_tris(v2)
      if (is.null(xv) || is.null(yv))
        return(FALSE)
      identical(xv$analysis, "ttest") &&
        identical(yv$analysis, "ttest") &&
        identical(xv$variable, "mean.diff") &&
        yv$variable %in% c("log.fdr", "log.pvalue")
    })

    # The displayed axes have caught up with the canonical store axes. The
    # corner auto-selection waits for this so the volcano rectangles are
    # never applied to the outgoing figure mid-switch. A settled triselector
    # triple is complete by construction (never a placeholder / partial
    # cascade), so a non-NULL value is a component-complete triple.
    .scatter_axes_converged <- reactive({
      vals <- lapply(store_watchers, function(w) w())
      if (any(vapply(vals, is.null, logical(1))))
        return(FALSE)
      xv <- .scatter_read_tris(v1)
      yv <- .scatter_read_tris(v2)
      !is.null(xv) && !is.null(yv) &&
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
      reactive_triset = triset, pre_volcano = displayed_volcano, reactive_status = attr4select_status,
      store = store,
      corner_apply_gate = .scatter_axes_converged,
      corner_valid_on_axes = displayed_volcano
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
      cutoff <- attr4select$cutoff_reactive
      l <- if (is.function(cutoff)) cutoff() else attr4select$cutoff
      if (is.null(l) || l$corner == "None") {
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

      line_rect(l = l, coords)$rect
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
      rr0 <- rectval()
      l$rect <- rr0
      # A restoration on unchanged axes needs one deliberate redraw. The
      # seed rides as an attribute: the render barrier's identical()
      # compares attributes (so a seed bump forces exactly one commit and
      # repaint), while c() in plotly_scatter's do.call drops attributes,
      # so it never reaches the plotly call itself.
      attr(l, "redrawSeed") <- selectionDisplayTrigger()

      # ---- emphasis resolution (see the display contract above) ----
      # Two display sources, resolved by VALUE (no observer ordering):
      #
      # - the anchor path (browser/restore authority): ids from the
      #   settled selVal, gated on the axes the selection was made or
      #   restored on. Read isolated: a browser event must not repaint
      #   (Plotly owns the visible box/lasso; a replot would erase it).
      # - the corner path (server-owned rects): ids recomputed from the
      #   LIVE rectangles and coordinates while the corner is engaged, so
      #   the emphasis lands in the SAME paint as an axis or cutoff
      #   change (a selection written by an observer lands one queue
      #   position after the render consumers and would otherwise paint
      #   a figure whose emphasis contradicts the propagated selection -
      #   the log.fdr -> log.pvalue switch left the corner genes selected
      #   in the tables but emphasized nothing).
      #
      # The anchor path wins only when it holds a DIFFERENT, non-empty
      # selection on these exact axes (a browser lasso or a manual
      # restore overriding the corner). When the corner owned the display
      # and its rects vanished (disarmed corner, left the volcano), the
      # anchor is skipped: nothing else has claimed the display since,
      # and the settle observer is about to clear the residual state.
      current_axes <- .scatter_axis_signature(v1(), v2())
      restored_axes <- isolate(pendingSelectionDisplayAxes())
      displayed_axes <- isolate(selectionDisplayAxes())
      if (!is.null(restored_axes) && identical(restored_axes, current_axes)) {
        displayed_axes <- restored_axes
      }
      anchor_ids <- if (!is.null(displayed_axes) && identical(displayed_axes, current_axes)) {
        as.character(isolate(selVal()$selected))
      } else {
        character(0)
      }
      if (isTRUE(isolate(cornerAuthority())) && is.null(rr0))
        anchor_ids <- character(0)
      corner_ids <- if (isTRUE(isolate(cornerEngaged())) && notNullAndPositiveLength(rr0)) {
        get_names()[.scatter_ids_in_rects(l, rr0)]
      } else {
        character(0)
      }
      selected_ids <- if (notNullAndPositiveLength(anchor_ids) &&
                          !identical(anchor_ids, corner_ids))
        anchor_ids else corner_ids
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
      reactive_regLine = showRegLine, htest_var1 = htestV1, htest_var2 = htestV2,
      # render barrier gate: commits only when BOTH displayed axes have
      # caught up with the canonical store axes - blocks the stale first
      # paint of the outgoing figure and the mixed-axis frames the serial
      # acknowledgement model produces mid-cascade (RC1)
      reactive_ready = .scatter_axes_converged
    )
    .scatter_keep(observe({
      showRegLine(v_scatter()$regline)
    }))

    selVal <- reactiveVal(
      list(
        clicked = character(0),
        selected = character(0)
      )
    )

    # Clear (user button or dataset change): the selection is empty
    # everywhere - propagation through the bus (mirror = TRUE un-filters
    # the tables), emphasis gone (authority "none"), and the corner is
    # DISENGAGED so a later rectval echo cannot resurrect the selection
    # (observeEvent fires on invalidation, not on value change: the
    # clear_counter dependency re-armed the corner observer and it
    # re-selected the corner genes immediately after the clear). A genuine
    # cutoff/corner edit re-engages (observer below). The seed bump forces
    # the one repaint that drops the emphasis; the rects stay rendered -
    # they visualize the cutoff configuration, which clear does not touch.
    # The FIRST computation of reactive_expr() also fires this observer
    # (observeEvent does not ignore it); it is not a dataset change, so a
    # signature guard skips it - an init-time clear report would disarm
    # the corner before the load-time volcano intent arms it.
    .clear_expr_sig <- reactiveVal(NULL)
    .scatter_keep(observeEvent(list(input$clear, reactive_expr()), {
      if (is.null(input$clear) || input$clear == 0L) {
        e <- tryCatch(reactive_expr(), error = function(e) NULL)
        sig <- paste(dim(e), collapse = "x")
        prev <- .clear_expr_sig()
        .clear_expr_sig(sig)
        if (is.null(prev) || identical(sig, prev))
          return(NULL)
      }
      selVal(list(
        clicked = character(0),
        selected = character(0)
      ))
      selectionDisplayAxes(NULL)
      pendingSelectionDisplayAxes(NULL)
      cornerEngaged(FALSE)
      cornerAuthority(FALSE)
      selectionDisplayTrigger(isolate(selectionDisplayTrigger()) + 1L)
      selection$report(
        origin = "clear",
        report = list(src = "clear", n = isolate(selectionDisplayTrigger())),
        ids = character(0), anchor = NULL, mirror = TRUE)
    }))

    # Workaround: Track previous selection to prevent redundant updates
    # Plotly events can fire even when selection hasn't actually changed,
    # causing unnecessary reactive chain invalidations. We store the previous
    # selection and only update selVal when it truly changes.
    clientSideSelection <- reactiveVal(character(0))
    .scatter_keep(observeEvent(v_scatter(), {
      l <- get_names()
      u_c <- l[v_scatter()$clicked]
      u_s <- l[v_scatter()$selected]

      # Only update if selection actually changed
      req(!identical(tmp <- c(u_c, u_s), clientSideSelection()))
      clientSideSelection(tmp)
      axes <- if (notNullAndPositiveLength(tmp))
        .scatter_axis_signature(v1(), v2()) else NULL
      selVal(list(
        clicked = u_c,
        selected = u_s
      ))
      selectionDisplayAxes(axes)
      # a browser event overrides the corner display authority; the flip
      # is isolated in scatter_vars, so Plotly's own lasso/box shape is
      # never erased by a repaint
      cornerAuthority(FALSE)
      # unified propagation: selected wins over clicked (the effective
      # selection), anchored to the axes it was made on and mirrored into
      # the table row highlight
      eff <- if (notNullAndPositiveLength(u_s)) u_s else u_c
      selection$report(
        origin = "figure",
        report = list(clicked = u_c, selected = u_s),
        ids = eff,
        clicked = u_c, anchor = axes,
        mirror = if (notNullAndPositiveLength(eff)) eff else TRUE)
    }))

    # Re-engagement: any GENUINE cutoff/corner configuration change (user
    # edit of xcut/ycut/scorner) re-arms the corner auto-selection. The
    # historic returnCornerSelection flag was only ever cleared (by a
    # manual-selection restore) and never re-armed, leaving the corner
    # dead for the rest of the session. Value-guarded: cutoff_effective
    # recomputes on every axis write (corner_valid_on_axes reads the
    # displayed axes) without the cutoff VALUE changing.
    .cutoff_key <- reactiveVal(NULL)
    .scatter_keep(observe({
      cutoff <- attr4select$cutoff_reactive
      l <- if (is.function(cutoff)) cutoff() else attr4select$cutoff
      if (is.null(l))
        return(NULL)
      key <- paste(l$corner %||% "", l$x %||% "", l$y %||% "", sep = "\r")
      if (identical(key, .cutoff_key()))
        return(NULL)
      .cutoff_key(key)
      cornerEngaged(TRUE)
    }))

    # The corner selection consumes only SETTLED state, and the observer
    # is kept referenced (observer-GC rule). Mid-cascade evaluations of
    # rectval carry stale components - the triselector outputs hold the
    # outgoing triple while the store already holds the incoming one, and
    # the scorner widget lags the corner resolution by a round trip - so
    # an ungated fire re-selected the OUTGOING view's corner regions (or
    # a mixed axes/corner combination) nondeterministically, depending on
    # intra-flush ordering; and an unreferenced observer is garbage
    # collected between flushes, silently swallowing the clearing fire
    # and leaving the stale selection alive (both observed live: quick
    # view switches resurrected the previous volcano's corner genes).
    #
    # observeEvent fires on INVALIDATION, not on value change: rectval
    # also invalidates on clear_counter bumps and cutoff echoes, so the
    # value is deduped explicitly - a same-value re-fire must not claim
    # the selection again (that is what resurrected the corner genes
    # right after a clear).
    #
    # The observer triggers on BOTH the rects and the AXES-CONVERGED
    # reactive: in a real browser the committed triselector triple
    # settles one flush BEFORE the canonical store write lands (the
    # UI->store sync runs behind the client acknowledgement round
    # trips), so a rectval-only event fired exactly once mid-cascade,
    # bailed on the convergence gate, and never re-ran after the store
    # caught up - an axis switch (y log.fdr -> log.pvalue) re-derived the
    # corner rects in-render but never re-REPORTED the selection, and the
    # result space kept the stale ids. The convergence reactive provides
    # the second trigger; the value dedupe below absorbs the extra runs.
    .rectval_last <- reactiveVal(NULL)
    .scatter_keep(observe({
      conv <- tryCatch(.scatter_axes_converged(), error = function(e) FALSE)
      rec <- rectval()
      if (!isTRUE(cornerEngaged())) {
        return(NULL)
      }
      if (!isTRUE(conv)) {
        return(NULL)
      }

      if (identical(rec, .rectval_last())) {
        return(NULL)
      }
      .rectval_last(rec)

      if (is.null(rec)) {
        selVal(list(
          clicked = character(0),
          selected = character(0)
        ))
        selectionDisplayAxes(NULL)
        cornerAuthority(FALSE)
        selection$report(
          origin = "corner",
          report = list(corner = "None"),
          ids = character(0), mirror = TRUE)
        return(NULL)
      }
      req(cc <- xycoord())
      l <- get_names()

      i <- .scatter_ids_in_rects(cc, rec)
      axes <- .scatter_axis_signature(v1(), v2())
      selVal(list(
        clicked = character(0),
        selected = l[i]
      ))
      selectionDisplayAxes(axes)
      # the corner claims the display. No repaint seed is needed: every
      # state change that alters the emphasis (axis switch, cutoff edit)
      # already recomputes the params, and the emphasis is resolved
      # IN-RENDER from the live rects - the value-identical recompute a
      # flip would cause is absorbed by the render barrier's identical().
      cornerAuthority(TRUE)
      selection$report(
        origin = "corner",
        report = list(rects = rec),
        ids = l[i], anchor = axes,
        mirror = if (notNullAndPositiveLength(l[i])) l[i] else TRUE)
    }))

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
        # the authority flag is "the corner made the current selection"
        # (the historic sbc): restored as the corner engagement below
        selectByCorner = isTRUE(isolate(cornerAuthority()))
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
      #
      # Corner engagement: a corner-made snapshot re-engages the corner
      # (the rectval observer re-derives and re-claims it once the restored
      # axes settle); a manual-selection snapshot disengages it (the corner
      # must not override the restored selection) until the user edits a
      # cutoff (re-engagement observer above).
      cornerEngaged(isTRUE(s$selectByCorner))
      cornerAuthority(FALSE)
      selVal(list(
        clicked = s$selection_clicked %||% character(0),
        selected = s$selection_selected %||% character(0)
      ))
      # unified propagation: the restored selection is the bus record for
      # this space (ids only; the table-row mirror is restored separately
      # by the owning data-space module from its own status fields)
      selection$report(
        origin = "restore",
        report = list(clicked = s$selection_clicked %||% character(0),
                      selected = s$selection_selected %||% character(0)),
        ids = s$selection_selected %||% character(0),
        clicked = s$selection_clicked %||% character(0),
        anchor = restored_axes)
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
