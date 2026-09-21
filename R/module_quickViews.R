#' Quick-view shortcut helpers and badge module
#'
#' These internal utilities define one-click scatter-plot shortcuts. Custom
#' shortcuts are read from the \code{"quickViews"} or \code{"shortcut"}
#' attribute of feature/sample metadata. If no shortcuts are available, common
#' PCA and volcano views are detected from the standard omicsViewer column
#' naming conventions.
#'
#' @keywords internal
#' @name quickViews-internal
NULL

.empty_quick_views <- function() {
  data.frame(
    id = character(0),
    label = character(0),
    x = character(0),
    y = character(0),
    description = character(0),
    source = character(0),
    stringsAsFactors = FALSE
  )
}

.quick_view_value <- function(x, default = "") {
  if (is.null(x) || length(x) == 0)
    return(default)
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x))
    default else x
}

.quick_view_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x, perl = TRUE)
  x <- sub("^_+|_+$", "", x)
  if (!nzchar(x))
    x <- "view"
  x
}

.quick_view_axis <- function(axis) {
  if (!is.character(axis) || length(axis) != 1 || is.na(axis))
    return(NULL)
  v <- strsplit(axis, "|", fixed = TRUE)[[1]]
  if (length(v) != 3 || any(!nzchar(v)))
    return(NULL)
  list(v1 = v[1], v2 = v[2], v3 = v[3])
}

.quick_view_field <- function(entry, name, position) {
  if (!is.null(names(entry)) && name %in% names(entry))
    return(entry[[name]])
  if (length(entry) >= position)
    entry[[position]]
  else
    NULL
}

.axis_in_triset <- function(axis, triset) {
  v <- .quick_view_axis(axis)
  if (is.null(v) || is.null(triset) || !nrow(triset))
    return(FALSE)
  any(triset[, 1] == v$v1 & triset[, 2] == v$v2 & triset[, 3] == v$v3)
}

#' Parse and validate user-defined quick-view shortcuts
#'
#' Both the recommended list-of-records format and a compact named-vector
#' format are supported:
#' \code{list(pca = c(x, y, name = "PCA"))}.
#'
#' @param shortcuts A list or data.frame of shortcut definitions.
#' @param triset Matrix returned by \code{\link{trisetter}}. Used to ensure
#'   both axes can be represented by a triselector.
#' @param source Character vector of length one used to mark the origin.
#' @return A data.frame with id, label, x, y, description, and source columns.
#' @keywords internal
parse_quick_views <- function(shortcuts, triset = NULL, source = "custom") {
  if (is.null(shortcuts) || length(shortcuts) == 0)
    return(.empty_quick_views())

  if (is.data.frame(shortcuts)) {
    entries <- lapply(seq_len(nrow(shortcuts)), function(i) as.list(shortcuts[i, ]))
    entryNames <- rep("", length(entries))
  } else {
    entries <- shortcuts
    entryNames <- names(shortcuts)
    if (is.null(entryNames))
      entryNames <- rep("", length(entries))
  }

  parsed <- lapply(seq_along(entries), function(i) {
    e <- entries[[i]]
    if (is.atomic(e)) {
      x <- .quick_view_value(.quick_view_field(e, "x", 1), NA_character_)
      y <- .quick_view_value(.quick_view_field(e, "y", 2), NA_character_)
      label <- .quick_view_value(.quick_view_field(e, "name", 3))
      id <- .quick_view_value(.quick_view_field(e, "id", 4), if (!is.na(entryNames[i]) && nzchar(entryNames[i])) entryNames[i] else label)
      description <- .quick_view_value(.quick_view_field(e, "description", 5))
    } else {
      x <- .quick_view_value(.quick_view_field(e, "x", 1), NA_character_)
      y <- .quick_view_value(.quick_view_field(e, "y", 2), NA_character_)
      label <- .quick_view_value(.quick_view_field(e, "label", 3))
      if (!nzchar(label))
        label <- .quick_view_value(e$name)
      id <- .quick_view_value(.quick_view_field(e, "id", 4), if (!is.na(entryNames[i]) && nzchar(entryNames[i])) entryNames[i] else label)
      description <- .quick_view_value(.quick_view_field(e, "description", 5))
      if (!nzchar(description))
        description <- .quick_view_value(e$desc)
    }

    id <- .quick_view_id(.quick_view_value(id, "view"))
    label <- .quick_view_value(label, id)
    description <- .quick_view_value(description, sprintf("Use %s for X and %s for Y", x, y))
    data.frame(
      id = id, label = label, x = x, y = y,
      description = description, source = source,
      stringsAsFactors = FALSE, row.names = NULL
    )
  })

  out <- do.call(rbind, parsed)
  if (!nrow(out))
    return(.empty_quick_views())

  out <- out[!is.na(out$x) & !is.na(out$y), , drop = FALSE]
  validAxes <- vapply(
    seq_len(nrow(out)),
    function(i) {
      .axis_in_triset(out$x[i], triset) &&
        .axis_in_triset(out$y[i], triset)
    },
    logical(1)
  )
  out <- out[validAxes, , drop = FALSE]
  if (!nrow(out))
    return(.empty_quick_views())

  out$id <- make.unique(out$id, sep = "_")
  out <- out[!duplicated(out[c("x", "y")]), , drop = FALSE]
  out
}

#' Detect standard quick views from metadata column names
#'
#' @param meta Feature or sample metadata.
#' @param triset Matrix returned by \code{\link{trisetter}}.
#' @return A data.frame of quick-view definitions.
#' @keywords internal
detect_quick_views <- function(meta, triset = NULL) {
  if (is.null(meta))
    return(.empty_quick_views())
  cn <- colnames(meta)
  if (is.null(cn))
    return(.empty_quick_views())

  entries <- list()

  pcaPrefixes <- c("PCA\\|All", "PCA\\|removeMissing")
  for (pcaPrefix in pcaPrefixes) {
    px <- grep(paste0("^", pcaPrefix, "\\|PC1\\("), cn, value = TRUE)
    py <- grep(paste0("^", pcaPrefix, "\\|PC2\\("), cn, value = TRUE)
    if (length(px) && length(py)) {
      entries <- c(entries, list(list(
        id = "pca",
        label = "PCA",
        x = px[1],
        y = py[1],
        description = "First two principal components"
      )))
      break
    }
  }

  # Correlation analyses are stored in feature metadata as
  # Cor|phenotype variable|R and Cor|phenotype variable|logP.
  corR <- grep("^Cor\\|.+\\|R$", cn, value = TRUE)
  for (xx in corR) {
    variable <- sub("^Cor\\|(.*)\\|R$", "\\1", xx)
    yy <- paste0("Cor|", variable, "|logP")
    if (!yy %in% cn)
      next
    entries <- c(entries, list(list(
      id = paste0("cor_", .quick_view_id(variable)),
      label = paste("Cor:", variable),
      x = xx,
      y = yy,
      description = sprintf("Correlation coefficient versus -log10 p-value for %s", variable)
    )))
  }

  fx <- grep("^ttest\\|.+\\|mean\\.diff$", cn, value = TRUE)
  for (xx in fx) {
    contrast <- sub("^ttest\\|(.*)\\|mean\\.diff$", "\\1", xx)
    yy <- paste0("ttest|", contrast, "|log.fdr")
    yLabel <- "log.fdr"
    if (!yy %in% cn) {
      yy <- paste0("ttest|", contrast, "|log.pvalue")
      yLabel <- "log.pvalue"
    }
    if (!yy %in% cn)
      next
    entries <- c(entries, list(list(
      id = paste0("volcano_", .quick_view_id(contrast)),
      label = paste("Volcano", gsub("_vs_", " vs ", contrast, fixed = TRUE)),
      x = xx,
      y = yy,
      description = sprintf("mean.diff versus %s for %s", yLabel, contrast)
    )))
  }

  if (!length(entries))
    return(.empty_quick_views())
  parse_quick_views(entries, triset = triset, source = "auto")
}

#' Read custom quick-view metadata
#' @param meta Feature or sample metadata.
#' @return User-defined shortcut metadata, or NULL.
#' @keywords internal
get_quick_views <- function(meta) {
  if (is.null(meta))
    return(NULL)
  value <- attr(meta, "quickViews")
  if (is.null(value))
    value <- attr(meta, "shortcut")
  value
}

#' Combine custom and automatically detected quick views
#'
#' @param meta Feature or sample metadata.
#' @param triset Matrix returned by \code{\link{trisetter}}.
#' @return A data.frame of available quick-view definitions.
#' @keywords internal
prepare_quick_views <- function(meta, triset = NULL) {
  custom <- parse_quick_views(get_quick_views(meta), triset = triset, source = "custom")
  auto <- detect_quick_views(meta, triset = triset)
  if (nrow(custom) && nrow(auto)) {
    used <- paste(custom$x, custom$y, sep = "\r")
    auto <- auto[!paste(auto$x, auto$y, sep = "\r") %in% used, , drop = FALSE]
  }
  out <- rbind(custom, auto)
  if (!nrow(out))
    .empty_quick_views() else out
}

#' Determine the active quick-view badge
#'
#' @param views Quick-view data.frame returned by \code{prepare_quick_views}.
#' @param xAxis X-axis value returned by a triselector.
#' @param yAxis Y-axis value returned by a triselector.
#' @return The id of the matching view, or character(0).
#' @keywords internal
active_quick_view <- function(views, xAxis, yAxis) {
  if (is.null(views) || !nrow(views) || is.null(xAxis) || is.null(yAxis))
    return(character(0))
  x <- paste(xAxis$analysis, xAxis$subset, xAxis$variable, sep = "|")
  y <- paste(yAxis$analysis, yAxis$subset, yAxis$variable, sep = "|")
  i <- which(views$x == x & views$y == y)
  if (!length(i))
    character(0) else views$id[i[1]]
}

#' Quick-view badge UI
#'
#' @param id Module id.
#' @return A dynamic UI containing quick visualization buttons.
#' @keywords internal
quick_badges_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("badgeContainer")) %>%
    tagAppendAttributes(`data-testid` = paste0(id, "-quick-view-badges"))
}

#' Quick-view badge server logic
#'
#' @param id Module id.
#' @param views Reactive data.frame of quick views.
#' @param activeId Reactive id of the currently active view.
#' @return A reactive list containing the selected view and a monotonically
#'   increasing trigger counter.
#' @keywords internal
quick_badges_module <- function(id, views, activeId) {
  moduleServer(id, function(input, output, session) {
    selectedId <- reactiveVal(NULL)
    selectedTrigger <- reactiveVal(0)
    clickCounts <- reactiveVal(list())

    observeEvent(views(), {
      vv <- isolate(views())
      selectedId(NULL)
      clickCounts(list())
      if (is.null(vv) || !nrow(vv))
        return(NULL)

      invisible(lapply(seq_len(nrow(vv)), function(i) {
        buttonId <- paste0("badge_", vv$id[i])
        local({
          buttonIdLocal <- buttonId
          viewIdLocal <- vv$id[i]
          observeEvent(input[[buttonIdLocal]], {
            currentClick <- input[[buttonIdLocal]]
            counts <- isolate(clickCounts())
            previousClick <- if (is.null(counts[[buttonIdLocal]])) 0 else counts[[buttonIdLocal]]
            counts[[buttonIdLocal]] <- currentClick
            clickCounts(counts)

            # Re-rendering the active button resets its actionButton counter.
            # Only increasing counters represent genuine user clicks.
            if (currentClick > previousClick) {
              selectedId(viewIdLocal)
              selectedTrigger(selectedTrigger() + 1)
            }
          }, ignoreInit = TRUE)
        })
      }))
    }, ignoreInit = FALSE)

    output$badgeContainer <- renderUI({
      vv <- views()
      if (is.null(vv) || !nrow(vv))
        return(tags$p("No quick views detected; use Custom visualization.", style = "color:#777; font-size:12px; margin:0;"))
      # Never read activeId() while creating buttons: before the triselectors
      # initialize, its req() can suspend this renderUI permanently. Keep the
      # button DOM stable and update active styling in place after each flush.
      active <- character(0)
      ns <- session$ns

      buttons <- lapply(seq_len(nrow(vv)), function(i) {
        isActive <- identical(vv$id[i], active)
        actionButton(
          ns(paste0("badge_", vv$id[i])),
          vv$label[i],
          class = paste("btn", "btn-xs", if (isActive) "btn-primary" else "btn-default"),
          title = vv$description[i],
          `aria-pressed` = if (isActive) "true" else "false"
        ) %>%
          tagAppendAttributes(
            `data-testid` = paste0(session$ns(""), "badge-", vv$id[i]),
            `data-quick-view-id` = vv$id[i]
          )
      })

      badgeUi <- tags$div(
        class = "quick-view-badges",
        style = "display:flex; flex-wrap:wrap; align-items:center; gap:4px; margin-bottom:2px;",
        role = "group",
        `aria-label` = "Quick visualization shortcuts",
        buttons
      )
      badgeUi
    })

    observe({
      active <- activeId()
      session$onFlushed(function() {
        containerId <- session$ns("badgeContainer")
        activeLiteral <- if (length(active)) paste0("\x27", active, "\x27") else "null"
        shinyjs::runjs(paste0(
          "(function() {",
          "  var container = document.getElementById(\x27", containerId, "\x27);",
          "  if (!container) return;",
          "  var active = ", activeLiteral, ";",
          "  container.querySelectorAll(\x27[data-quick-view-id]\x27).forEach(function(button) {",
          "    var isActive = button.getAttribute(\x27data-quick-view-id\x27) === active;",
          "    button.classList.toggle(\x27btn-primary\x27, isActive);",
          "    button.classList.toggle(\x27btn-default\x27, !isActive);",
          "    button.setAttribute(\x27aria-pressed\x27, isActive ? \x27true\x27 : \x27false\x27);",
          "  });",
          "})();"
        ))
      }, once = TRUE)
    })

    reactive({
      trigger <- selectedTrigger()
      if (!trigger)
        return(NULL)
      vv <- isolate(views())
      view <- vv[vv$id == selectedId(), , drop = FALSE]
      if (!nrow(view))
        return(NULL)
      list(view = view, trigger = trigger)
    })
  })
}

#' Local null-coalescing helper for quick-view code
#'
#' @param lhs Left-hand value.
#' @param rhs Right-hand value returned when \code{lhs} is NULL.
#' @return \code{lhs}, or \code{rhs} if \code{lhs} is NULL.
#' @name quickViewsNullCoalesce
#' @keywords internal
'%||%' <- function(lhs, rhs) if (is.null(lhs)) rhs else lhs
