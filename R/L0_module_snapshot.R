#' Snapshot save/restore module (L0 split, todo 4.4c)
#'
#' Extracted verbatim from \code{app_module}'s server: the .ESS listing,
#' save observer (with the opt-in assistant-history payload), delete flow,
#' restore modal + confirm dialog, and the phased restore controller
#' (todo 2.8/4.1). Dependencies are explicit formals with the names the
#' inline code closed over, so the moved code is unchanged.
#'
#' @param input,output,session The CALLER's moduleServer scope: the
#'   snapshot modal UI lives in app_module's UI function, so this split
#'   keeps the app namespace verbatim instead of re-namespacing (a
#'   moduleServer wrapper would prefix every input id).
#' @param .dir Reactive data directory.
#' @param dataset_id Reactive dataset id (exact-id snapshot listing).
#' @param dataset_ready Reactive gate (vEset-equivalent).
#' @param eset Reactive ExpressionSet.
#' @param data_space Reactive data-space module return (selection ports,
#'   axes-convergence handles, and the status attribute).
#' @param result_status Reactive result-space status (\code{v2()}).
#' @param selection_features,selection_samples L0 mirror reactiveVals
#'   (\code{ri}/\code{rh}).
#' @param esv_status The L0 panel-status reactiveVal (written by the
#'   restore's phase-1 delivery).
#' @param data_space The data-space module return (selection ports +
#'   axes-convergence handles, todo 4.1).
#' @param store Canonical widget store.
#' @param assistant Optional WP11 assistant API.
#' @keywords internal
#' @rdname snapshotModule
snapshot_module <- function(input, output, session, .dir, dataset_id,
                            dataset_ready, eset, result_status,
                            selection_features, selection_samples,
                            esv_status, data_space, store,
                            assistant = NULL) {
  {
    ns <- session$ns
    vEset <- function() dataset_ready()
    reactive_eset <- function() eset()
    # the data-space module RETURN (selection ports + convergence handles
    # + the status attribute the save observer reads via attr(v1(), "status"))
    v1 <- function() data_space()
    v2 <- function() result_status()
    ri <- selection_features
    rh <- selection_samples
    app_store <- store
    assistant_api <- assistant
    current_dataset_id <- function() dataset_id()

  savedSS <- reactiveVal(
    data.frame(name = character(), link = character(), schema = integer(),
               created_at = character(), package_version = character(),
               stringsAsFactors = FALSE)
  )
  snapshot_refresh <- reactiveVal(0L)

  observe({
    req(.dir())
    snapshot_refresh()
    dsid_raw <- current_dataset_id()
    dsid <- sanitize_snapshot_name(dsid_raw, fallback = "ESVObj.RDS")
    prefix <- paste0("ESVSnapshot_", dsid, "_")
    ff <- list.files(.dir(), pattern = "\\.ESS$", ignore.case = TRUE)
    ff <- ff[startsWith(ff, prefix)]

    # Read compact metadata only; large panel payloads stay on disk. The
    # dataset id stored INSIDE the file decides membership (todo 2.7): a
    # filename prefix alone leaks across datasets whose sanitized ids
    # share a prefix ("demo.RDS" also listed "demo.RDS_v2.RDS"
    # snapshots). Files without a readable id (legacy) keep the prefix
    # match as a fallback.
    meta <- lapply(ff, function(f) tryCatch({
      x <- readRDS(file.path(.dir(), f))
      list(
        id = if (is.null(x$dataset$id)) NA_character_ else as.character(x$dataset$id),
        schema = if (is.null(x$schema_version)) NA_integer_ else as.integer(x$schema_version),
        created_at = x$created_at %.or_default% NA_character_,
        package_version = x$package_version %.or_default% NA_character_
      )
    }, error = function(e) list(id = NA_character_, schema = NA_integer_,
                                created_at = NA_character_,
                                package_version = NA_character_)))
    keep <- vapply(meta, function(m)
      is.na(m$id) || identical(m$id, dsid_raw), logical(1))
    ff <- ff[keep]
    meta <- meta[keep]

    if (length(ff) == 0) {
      savedSS(data.frame(name = character(), link = character(), schema = integer(),
                         created_at = character(), package_version = character(),
                         stringsAsFactors = FALSE))
      return(NULL)
    }

    savedSS(data.frame(
      name = sub("\\.ESS$", "", sub(prefix, "", ff)),
      link = ff,
      schema = vapply(meta, function(x) x$schema, integer(1)),
      created_at = vapply(meta, function(x) x$created_at, character(1)),
      package_version = vapply(meta, function(x) x$package_version, character(1)),
      stringsAsFactors = FALSE
    ))
  })

  shinyInput <- function(FUN, len, id, ...) {
    inputs <- c()
    for (i in len) {
      inputs <- c(inputs, as.character(FUN(paste0(id, i), ...)))
    }
    inputs
  }

  output$tab_saveSS <- renderDT({
    req(nrow(dt <- savedSS()) > 0)
    dt$delete <- shinyInput(
      actionButton, dt$name, "deletess_", label = "Delete",
      onclick = sprintf('Shiny.setInputValue("%s", this.id, {priority: "event"})', ns("deletess_button"))
    )
    dt$info <- ifelse(
      is.na(dt$schema),
      "legacy",
      paste0(
        "v", dt$schema,
        ifelse(is.na(dt$created_at), "", paste0(" | ", dt$created_at)),
        ifelse(is.na(dt$package_version), "", paste0(" | ", dt$package_version))
      )
    )
    DT::datatable(
      dt[, c("name", "info", "delete"), drop = FALSE],
      rownames = FALSE, colnames = c("Name", "Metadata", ""),
      selection = list(mode = "single", target = "cell",
                       selectable = -cbind(seq_len(nrow(dt)), 3)),
      escape = FALSE,
      options = list(
        dom = "t", autoWidth = FALSE, style = "compact-hover", scrollY = "450px",
        paging = FALSE,
        columns = list(list(width = "40%"), list(width = "42%"), list(width = "18%"))
      )
    )
  })

  selectedSS <- reactiveVal()
  observe({
    ss <- input$tab_saveSS_cells_selected
    if (length(ss) == 0 || ss[2] > 1)
      return(NULL)
    selectedSS(ss[1])
  })

  observeEvent(list(v1(), v2()), {
    selectedSS(NULL)
  })

  deleteSS <- reactiveVal()
  observeEvent(input$deletess_button, {
    selectedRow <- sub("deletess_", "", input$deletess_button, fixed = TRUE)
    deleteSS(selectedRow)
  })

  observeEvent(deleteSS(), {
    req(nrow(df <- savedSS()) > 0)
    req(i <- match(deleteSS(), df$name))
    showModal(modalDialog(
      title = tagList(icon("trash"), "Delete snapshot"),
      sprintf("Delete snapshot %s? This cannot be undone.", df$name[i]),
      footer = tagList(
        actionButton(ns("snapshot_delete_cancel"), "Cancel"),
        actionButton(ns("snapshot_delete_confirm"), label = tagList(icon("trash"), "Delete"), class = "btn-danger")
      ),
      easyClose = TRUE
    ) %>% tagAppendAttributes(class = "omicsviewer-modal"))
  })

  observeEvent(input$snapshot_delete_cancel, {
    removeModal()
    deleteSS(NULL)
  })

  observeEvent(input$snapshot_delete_confirm, {
    # Snapshot I/O failures (read-only data dir, full disk, stale listing)
    # must surface as notifications, not unhandled errors that close the
    # session (todo 1.5).
    tryCatch(
      {
        req(nrow(df <- savedSS()) > 0)
        req(i <- match(deleteSS(), df$name))
        unlink(file.path(.dir(), df$link[i]))
        removeModal()
        deleteSS(NULL)
        snapshot_refresh(snapshot_refresh() + 1L)
      },
      shiny.silent.error = function(e) NULL,
      error = function(e) {
        showNotification(
          sprintf("Could not delete snapshot: %s", conditionMessage(e)),
          type = "error", duration = 10
        )
      }
    )
  })

  # shared experimental-feature disclaimer for the snapshot popups
  .snapshot_experimental_note <- tags$div(
    style = paste("padding: 8px 12px; margin-bottom: 10px;",
                  "border-left: 4px solid #f0ad4e; background: #fcf8e3;",
                  "color: #8a6d3b; font-size: 90%;"),
    tags$strong("The snapshot function is experimental."),
    "Saved states may not restore exactly across package versions or",
    "revised datasets. Plotly-internal events (lasso/box shapes) are not",
    "restored - only the semantic selection of features/samples is.",
    "Do not rely on snapshots as the only record of an analysis."
  )

  observeEvent(input$snapshot, {
    # refresh the listing on modal open (todo 2.7): files changed outside
    # this session (another user, another tab) must be visible immediately
    snapshot_refresh(snapshot_refresh() + 1L)
    showModal(
      modalDialog(
        title = tagList(icon("camera-retro"), "Snapshots"),
        .snapshot_experimental_note,
        fluidRow(
          column(9, textInput(ns("snapshot_name"), label = "Save new snapshot", placeholder = "snapshot name", width = "100%")),
          column(3, style = "padding-top:31px", actionButton(ns("snapshot_save"), label = tagList(icon("save"), "Save")))
        ),
        # WP11: conversation persistence is opt-in per snapshot - transcripts
        # may contain sensitive dataset content, so nothing is saved silently.
        checkboxInput(
          ns("snapshot_include_chat"),
          "Include the AI assistant conversation (and its figure registry)",
          value = FALSE
        ) %>%
          tagAppendAttributes(
            title = "The chat transcript, tool results, and figure specifications are stored inside this snapshot. Credentials are never included."
          ),
        hr(),
        strong("Load saved snapshots:"),
        DTOutput(ns("tab_saveSS")),
        footer = NULL,
        easyClose = TRUE
      ) %>% tagAppendAttributes(class = "omicsviewer-modal")
    )
  })

  observeEvent(input$snapshot_save, {
    req(vEset())
    req(reactive_eset())
    name <- sanitize_snapshot_name(input$snapshot_name, fallback = paste0("snapshot-", format(Sys.time(), "%Y%m%d-%H%M%S")))

    df <- savedSS()
    if (name %in% df$name) {
      showNotification(sprintf("Snapshot name %s is already in use.", name), type = "error")
      return(NULL)
    }

    flink <- file.path(.dir(), snapshot_file_name(name, dataset_id = current_dataset_id(), fallback = name))
    if (file.exists(flink)) {
      showNotification("A snapshot file with this name already exists.", type = "error")
      return(NULL)
    }

    # Snapshot I/O failures (read-only data dir, full disk, name longer than
    # the file system allows) must surface as notifications, not unhandled
    # errors that close the session (todo 1.5).
    saved_ok <- tryCatch(
      {
        data_status <- tryCatch(
          attr(v1(), "status"),
          shiny.silent.error = function(e) list(),
          error = function(e) stop(e)
        )
        result_status <- tryCatch(
          v2(),
          shiny.silent.error = function(e) list(),
          error = function(e) stop(e)
        )
        obj <- build_app_state(
          dataset = reactive_eset(),
          dataset_id = current_dataset_id(),
          data_status = data_status,
          result_status = result_status,
          selected_features = ri(),
          selected_samples = rh(),
          # Canonical widget-store state rides along (schema 2): keeps every
          # registered widget's desired value in one authoritative snapshot,
          # with an explicit unset list so restores REPLACE instead of merge
          widget_store = store_snapshot(app_store),
          label = name
        )
        # WP11: the conversation rides along only when the user opted in and
        # there is a conversation to save (assistant unconfigured -> NULL).
        if (isTRUE(input$snapshot_include_chat)) {
          chat_payload <- tryCatch(
            assistant_api$snapshot_payload(),
            error = function(e) NULL
          )
          if (!is.null(chat_payload))
            obj$assistant <- chat_payload
        }
        write_app_state(obj, flink)
        TRUE
      },
      error = function(e) {
        showNotification(
          sprintf("Could not save snapshot %s: %s", name, conditionMessage(e)),
          type = "error", duration = 10
        )
        FALSE
      }
    )
    if (!isTRUE(saved_ok))
      return(NULL)
    snapshot_refresh(snapshot_refresh() + 1L)
    removeModal()
    showNotification(sprintf("Snapshot %s saved.", name), type = "message", duration = 3)
  })

  # Restore flow (todo 2.6/2.7): selecting a row opens a CONFIRM dialog
  # (a stray click must not silently reconfigure the session); confirming
  # runs the restore, surfaces a receipt notification (restored/rejected/
  # unknown keys, collected validation warnings) and resets the selection
  # state so re-clicking the same row works.
  restoreSS <- reactiveVal(NULL)

  observeEvent(selectedSS(), {
    req(vEset())
    req(nrow(df <- savedSS()) > 0)
    if (length(i <- selectedSS()) == 0)
      return(NULL)

    removeModal()
    ss <- tryCatch(readRDS(file.path(.dir(), df$link[i])), error = function(e) {
      showNotification("Could not read the selected snapshot.", type = "error")
      NULL
    })
    req(ss)

    ss <- tryCatch(
      validate_app_state(ss, dataset = reactive_eset(), dataset_id = current_dataset_id()),
      error = function(e) {
        showNotification(paste("Invalid snapshot:", conditionMessage(e)), type = "error")
        NULL
      }
    )
    req(ss)

    restoreSS(list(state = ss, name = df$name[i],
                   warnings = attr(ss, "warnings") %||% character()))
    showModal(modalDialog(
      title = tagList(icon("history"), "Restore snapshot"),
      .snapshot_experimental_note,
      sprintf("Restore snapshot \"%s\"? The current view state will be replaced.", df$name[i]),
      if (length(restoreSS()$warnings))
        tags$ul(tags$li(paste(restoreSS()$warnings, collapse = " "))) else NULL,
      footer = tagList(
        actionButton(ns("snapshot_restore_cancel"), "Cancel"),
        actionButton(ns("snapshot_restore_confirm"), label = tagList(icon("rotate-right"), "Restore"), class = "btn-primary")
      ),
      easyClose = TRUE
    ) %>% tagAppendAttributes(class = "omicsviewer-modal"))
  })

  observeEvent(input$snapshot_restore_cancel, {
    removeModal()
    restoreSS(NULL)
    selectedSS(NULL)
  })

  .snapshot_strip_store_owned <- function(state) {
    # The canonical widget store is the single write plane (todo 2.1/4.1):
    # when the snapshot carries a widget_store section, the panel copies
    # must not re-apply (they are legacy transport and may disagree with
    # the store after an agent write - the saved table status reports the
    # status-assembly-time column set, which can lag the store by a view)
    if (is.null(state$widget_store) || !is.list(state$widget_store$values))
      return(state)
    for (k in c("eset_fdata_fig", "eset_pdata_fig")) {
      fig <- state$panels$data_space[[k]]
      if (is.list(fig)) {
        fig$xax <- NULL
        fig$yax <- NULL
        fig$axisMode <- NULL
        state$panels$data_space[[k]] <- fig
      }
    }
    # Table columns/multi-selection are registered store keys
    # (dataspace.tab_*.columns / .multi_selection): the status copies
    # (showColumns/multiSelection) race the canonical store_restore and
    # the ack of the stale legacy push re-enters the store as a user
    # edit (observed live: a restored multi-column selection overwrote
    # the agent-set single column). DT one-shot state (page, filters,
    # ordering, rows) is NOT store-owned and keeps riding the status.
    for (k in c("eset_pdata_tab", "eset_fdata_tab", "eset_exprs_tab")) {
      tab <- state$panels$data_space[[k]]
      if (is.list(tab)) {
        tab$showColumns <- NULL
        tab$multiSelection <- NULL
        state$panels$data_space[[k]] <- tab
      }
    }
    state
  }

  # =====================================================================
  # Single restore controller (todo 2.8 deferred into 4.1): the restore
  # runs in phases so the selection lands exactly once, AFTER the axes it
  # belongs to. Phase 1 (confirm click): validate -> canonical
  # store_restore(replace) -> legacy panel-status delivery -> arm phase 2.
  # Phase 2 (observer): wait until BOTH scatter axes have converged onto
  # the restored store values, then apply the selection-bus records once
  # and show the receipt. The former code applied the selection from the
  # L1 status observer in the same flush as the axes cascade.
  # =====================================================================
  .restore_pending <- reactiveVal(NULL)

  .restore_valid_origin <- function(o)
    if (is.character(o) && length(o) == 1L && !is.na(o) &&
        o %in% c("figure", "corner", "clear", "table", "heatmap",
                 "cor_heatmap", "dyn_heatmap", "gslist", "restore", "system"))
      o else "restore"

  .restore_clean_ids <- function(x) {
    if (is.null(x) || length(x) == 0L) return(character(0))
    x <- as.character(x)
    x <- trimws(x[!is.na(x)])
    x[nzchar(x)]
  }

  .restore_selection_records <- function(ss) {
    # v2 snapshots carry the full bus record per space; v1/legacy fall
    # back to the plain id fields and the saved table mirrors.
    rec_f <- if (is.list(ss$selection$records)) ss$selection$records$feature else NULL
    rec_s <- if (is.list(ss$selection$records)) ss$selection$records$sample else NULL
    list(
      feature = if (is.list(rec_f) && !is.null(rec_f$ids)) list(
        ids = .restore_clean_ids(rec_f$ids),
        clicked = as.character(rec_f$clicked %||% character(0)),
        origin = .restore_valid_origin(rec_f$origin),
        anchor = if (is.null(rec_f$anchor)) NULL else as.character(rec_f$anchor),
        mirror = rec_f$mirror
      ) else list(
        ids = .restore_clean_ids(
          ss$selection$features %||%
            ss$panels$data_space$eset_selected_features),
        clicked = character(0),
        origin = "restore",
        anchor = NULL,
        mirror = ss$panels$data_space$eset_fdata_tabrows
      ),
      sample = if (is.list(rec_s) && !is.null(rec_s$ids)) list(
        ids = .restore_clean_ids(rec_s$ids),
        clicked = as.character(rec_s$clicked %||% character(0)),
        origin = .restore_valid_origin(rec_s$origin),
        anchor = if (is.null(rec_s$anchor)) NULL else as.character(rec_s$anchor),
        mirror = rec_s$mirror
      ) else list(
        ids = .restore_clean_ids(
          ss$selection$samples %||%
            ss$panels$data_space$eset_selected_samples),
        clicked = character(0),
        origin = "restore",
        anchor = NULL,
        mirror = ss$panels$data_space$eset_pdata_tabrows
      )
    )
  }

  .restore_receipt_notification <- function(name, receipt, n_warn) {
    n_applied <- if (is.null(receipt)) 0L else length(receipt$applied)
    n_reset <- if (is.null(receipt)) 0L else length(receipt$reset %||% character())
    n_rejected <- if (is.null(receipt)) 0L else length(receipt$rejected %||% list())
    n_unknown <- if (is.null(receipt)) 0L else length(receipt$unknown_ids %||% character())
    n_adj <- n_applied + n_reset + n_rejected + n_unknown + n_warn
    detail <- character()
    if (n_applied) detail <- c(detail, sprintf("%d widget%s restored", n_applied, if (n_applied == 1L) "" else "s"))
    if (n_reset) detail <- c(detail, sprintf("%d reset to default", n_reset))
    if (n_rejected) detail <- c(detail, sprintf("%d skipped (no longer valid)", n_rejected))
    if (n_unknown) detail <- c(detail, sprintf("%d unknown", n_unknown))
    if (n_warn) detail <- c(detail, sprintf("%d warning%s", n_warn, if (n_warn == 1L) "" else "s"))
    showNotification(
      sprintf("Snapshot %s restored%s%s", name,
              if (n_adj) sprintf(" (%d adjustment%s)", n_adj, if (n_adj == 1L) "" else "s") else "",
              if (length(detail)) paste0(": ", paste(detail, collapse = ", ")) else ""),
      type = if (n_warn || n_rejected) "warning" else "message",
      duration = 10
    )
  }

  observeEvent(input$snapshot_restore_confirm, {
    req(vEset())
    pend <- restoreSS()
    if (is.null(pend))
      return(NULL)
    # Graceful exit (todo 1.5 discipline extended to the restore path):
    # ANY failure in the restore pipeline (widget-store transaction,
    # assistant history revival, receipt assembly) surfaces as an error
    # notification with the session kept alive and the modal flow usable -
    # an escaping error would close the session via unhandledError.
    tryCatch({
      removeModal()
      restoreSS(NULL)
      selectedSS(NULL)

      ss <- .snapshot_strip_store_owned(pend$state)

      # Phase 1a: legacy panel-status delivery FIRST (DT state, heatmap
      # zoom, scatter local display state) - the transactional
      # NULL-then-state boundary child modules distinguish. Ordering
      # matters: the tables' one-shot status consumers must eat the legacy
      # transport BEFORE the canonical store values land, or the two
      # writers race the same selectize and the ack of the stale legacy
      # push re-enters the store as a user edit (observed live: restored
      # multi-column table selections reverted to the seeded single
      # column; store-owned axes/mode are stripped from the status below
      # and never race).
      esv_status(NULL)
      esv_status(ss)

      # Phase 1b: canonical widget-store restore (schema 2, replace
      # semantics). Per-key resilient: values invalidated by a revised
      # dataset are rejected and reported, not vetoed. Older snapshots
      # without widget_store skip this.
      receipt <- NULL
      if (!is.null(ss$widget_store) && is.list(ss$widget_store$values)) {
        receipt <- tryCatch(
          store_restore(app_store, ss$widget_store, replace = TRUE),
          error = function(e) {
            # degrade gracefully: the widget-store transaction failing must
            # not abort the panel-status restoration above
            showNotification(
              paste("Some widget settings could not be restored:", conditionMessage(e)),
              type = "warning", duration = 10)
            NULL
          }
        )
      }

      # WP11: revive the conversation + figure registry when the snapshot
      # carries one (opt-in at save time). Restored turns are inert context;
      # figure specs re-validate against the current dataset when reused.
      if (!is.null(ss$assistant)) {
        tryCatch(
          assistant_api$restore_history(ss$assistant),
          error = function(e)
            warning("Assistant snapshot restore failed: ", conditionMessage(e))
        )
      }

      # Phase 2 marker: selection applies once the restored axes converge.
      .restore_pending(list(
        records = .restore_selection_records(ss),
        receipt = receipt,
        name = pend$name,
        warnings = pend$warnings
      ))
    },
    shiny.silent.error = function(e) NULL,
    error = function(e) {
      removeModal()
      restoreSS(NULL)
      selectedSS(NULL)
      .restore_pending(NULL)
      showNotification(
        sprintf("Could not restore snapshot %s: %s", pend$name, conditionMessage(e)),
        type = "error", duration = 15
      )
    })
  })

  observe({
    pend <- .restore_pending()
    if (is.null(pend))
      return(NULL)
    handles <- tryCatch(isolate(v1())$axes_converged,
                        shiny.silent.error = function(e) NULL,
                        error = function(e) NULL)
    # Reactively wait for BOTH scatter axes to catch up with the restored
    # store values (calling the closures takes the dependencies; a
    # corner-origin selection then re-engages against the correct axes).
    # Handles unavailable (panel not yet warm) -> treat as converged so a
    # restore can never deadlock.
    converged <- if (is.null(handles)) TRUE else
      isTRUE(tryCatch(handles$feature(), error = function(e) FALSE)) &&
      isTRUE(tryCatch(handles$sample(), error = function(e) FALSE))
    if (!converged)
      return(NULL)
    .restore_pending(NULL)

    # Phase 2: the single selection apply, through the bus ports.
    tryCatch({
      ports <- tryCatch(isolate(v1())$selection,
                        shiny.silent.error = function(e) NULL,
                        error = function(e) NULL)
      if (is.null(ports))
        stop("Data-space panel unavailable for selection restore.")
      ports$feature$apply(
        ids = pend$records$feature$ids,
        clicked = pend$records$feature$clicked,
        origin = pend$records$feature$origin,
        anchor = pend$records$feature$anchor,
        mirror = pend$records$feature$mirror)
      ports$sample$apply(
        ids = pend$records$sample$ids,
        clicked = pend$records$sample$clicked,
        origin = pend$records$sample$origin,
        anchor = pend$records$sample$anchor,
        mirror = pend$records$sample$mirror)
      ri(pend$records$feature$ids)
      rh(pend$records$sample$ids)

      # Receipt (todo 2.6/2.7): every adjustment in one place
      .restore_receipt_notification(pend$name, pend$receipt,
                                    length(pend$warnings))
      if (length(pend$warnings))
        showNotification(paste(pend$warnings, collapse = "\n"),
                         type = "warning", duration = 15)
    }, error = function(e) {
      showNotification(
        paste("Selection could not be restored:", conditionMessage(e)),
        type = "warning", duration = 10
      )
    })
  })
  }
}
