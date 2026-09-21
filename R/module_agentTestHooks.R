#' Internal environment-gated agent bridge test hooks
#'
#' These hooks exist exclusively for headless UI testing (shinytest2 /
#' Playwright). They expose the exact \code{apply_agent_state} and
#' \code{apply_agent_scatter_view} callbacks used by the ellmer assistant
#' tools, so UI-effect tests drive the identical server-side pathway a model
#' tool call drives. The hooks are rendered and wired only when the
#' environment variable \code{OMICSVIEWER_TEST_HOOKS} is set to a true-like
#' value; ordinary deployments are unaffected. The hook UI is hidden from
#' users and never provides capabilities beyond the assistant tools
#' themselves.
#'
#' @return UI tags and a module server, plus the hook-enabled predicate.
#'
#' @keywords internal
#' @name agentTestHooksModule
NULL

#' Report whether agent test hooks are enabled
#'
#' @return TRUE when \code{OMICSVIEWER_TEST_HOOKS} is true, yes, on, or 1.
#' @keywords internal
#' @rdname agentTestHooksModule
agent_test_hooks_enabled <- function() {
  tolower(trimws(Sys.getenv("OMICSVIEWER_TEST_HOOKS")[1])) %in%
    c("true", "t", "yes", "y", "on", "1")
}

#' @param id Module ID.
#' @rdname agentTestHooksModule
#' @keywords internal
agent_test_hooks_ui <- function(id) {
  ns <- NS(id)
  tags$div(
    id = ns("container"),
    style = "display:none;",
    `data-testid` = "agent-test-hooks",
    selectInput(
      ns("op"), "Operation",
      choices = c(
        "state" = "state",
        "scatter" = "scatter",
        "overview" = "overview"
      ),
      selected = "state"
    ),
    textInput(
      ns("payload"), "JSON payload",
      value = "{}",
      width = "100%"
    ),
    actionButton(ns("run"), "Run hook"),
    verbatimTextOutput(ns("result"))
  )
}

#' @param id Module ID.
#' @param apply_state Callback applying a validated assistant state update
#'   (as used by the \code{set_omics_viewer_state} tool).
#' @param apply_scatter_view Callback applying a validated scatter view
#'   (as used by the \code{set_scatter_view} tool).
#' @param state Reactive compact assistant state (for the overview hook).
#' @rdname agentTestHooksModule
#' @keywords internal
agent_test_hooks_module <- function(id, apply_state, apply_scatter_view,
                                    state) {
  moduleServer(id, function(input, output, session) {
    last_result <- reactiveVal(NULL)

    run_hook <- function() {
      op <- input$op
      payload <- trimws(input$payload)
      if (!nzchar(payload))
        payload <- "{}"
      parsed <- tryCatch(
        jsonlite::fromJSON(payload, simplifyVector = FALSE),
        error = function(e) {
          list(hook_error = paste("Invalid JSON payload:", conditionMessage(e)))
        }
      )
      if (!is.null(parsed$hook_error))
        return(parsed)

      tryCatch(
        if (identical(op, "overview")) {
          isolate(state())
        } else if (identical(op, "scatter")) {
          do.call(
            apply_scatter_view,
            c(
              list(space = parsed$space),
              if (!is.null(parsed$quick_view_id)) list(quick_view_id = parsed$quick_view_id) else NULL,
              if (!is.null(parsed$x_axis)) list(x_axis = parsed$x_axis) else NULL,
              if (!is.null(parsed$y_axis)) list(y_axis = parsed$y_axis) else NULL
            )
          )
        } else {
          update <- parsed[c(
            "data_space_tab", "analysis_space_tab", "features", "samples"
          )]
          update <- update[!vapply(update, is.null, logical(1))]
          if (!length(update))
            stop("Payload contains no recognized state fields.")
          apply_state(update)
        },
        error = function(e) {
          list(hook_error = conditionMessage(e))
        }
      )
    }

    observeEvent(input$run, {
      result <- run_hook()
      if (!is.list(result))
        result <- list(value = result)
      last_result(c(list(hook_run = input$run), result))
    })

    output$result <- renderPrint({
      result <- last_result()
      if (is.null(result))
        cat("no result yet")
      else
        cat(jsonlite::toJSON(result, auto_unbox = TRUE, null = "null",
                             na = "null", pretty = 2))
    })

    invisible(NULL)
  })
}
