
survival_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Kaplan-Meier Survival Analysis"),
      tags$p("Kaplan-Meier survival analysis is a statistical method for analyzing time-to-event data, commonly used in clinical studies to estimate survival probability over time. This analysis accounts for censored observations (patients lost to follow-up or still alive at study end) and compares survival curves between different groups using the log-rank test."),
      tags$h4("When to use survival analysis"),
      tags$p("Use this analysis when you have time-to-event data such as patient survival time, time to disease recurrence, or time to treatment response. It's particularly valuable when some observations are censored (incomplete follow-up). This helps identify whether different patient groups or biomarker levels are associated with different survival outcomes."),
      tags$h4("How to interpret results"),
      tags$p("The Kaplan-Meier plot shows survival probability (y-axis) over time (x-axis) for different groups. Each step down represents an event (e.g., death). Censored observations are marked with tick marks. Wider separation between curves indicates greater differences in survival. The log-rank p-value tests whether survival curves are significantly different - values less than 0.05 indicate significant differences between groups. Median survival time is where the curve crosses the 50% survival line.")
    ),
    uiOutput(ns("censor_output")),
    plotOutput(ns("kmplot"))
  )

}

#' @description Utility survival KM module
#' @param id module id
#' @param reactive_resp reponse value, in the format 1, 345, 345+, 23, 45, 355+
#' @param reactive_strata strata variable
#' @param reactive_checkpoint checkpoint
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}); when given, the censor-time
#'   slider registers under this view (the sample_general owner passes its
#'   own child so the key lands at
#'   \code{resultspace.sample_general.survival_censor}). The push is gated
#'   on the checkpoint because the slider only exists in the survival view.
#' @importFrom survminer ggsurvplot surv_pvalue
#' @importFrom survival survfit Surv
#' @examples
#' #' # library(shiny)
#' # library(survminer)
#' # library(survival)
#' #
#' # ui <- fluidPage(
#' #   survival_ui("surv")
#' # )
#' #
#' # server <- function(input, output, session) {
#   survival_module('surv', reactive_resp = reactive(t), reactive_strata = reactive#' (strata))
#' # }
#' #
#' # # significant
#' # v <- runif(n = 1000, 1, 1000)
#' # t <- paste0(v, sample(c("", "+"), replace = TRUE, size = 1000))
#' # strata <- c("a", "b")[as.integer(v < 500)+1]
#' # shinyApp(ui, server)
#' #
#' # # insignificant
#' # v <- runif(n = 1000, 1, 1000)
#' # t <- paste0(v, sample(c("", "+"), replace = TRUE, size = 1000))
#' # strata <- sample(c("a", "b"), replace = TRUE, size = 1000)
#' # shinyApp(ui, server)
#'
survival_module <- function(
  id, reactive_resp, reactive_strata, reactive_checkpoint = reactive(TRUE),
  store = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # The censor slider is the module's only user-editable widget; its
  # bounds derive from the data, so the binding keeps them open and the
  # push clamps to the rendered range. Observer retention mandatory.
  # ------------------------------------------------------------------
  .sv_store_observers <- list()
  .sv_keep <- function(obs) {
    .sv_store_observers[[length(.sv_store_observers) + 1L]] <<- obs
    invisible(obs)
  }
  if (!is.null(store)) {
    .sv_range <- function() {
      d <- tryCatch(dat(), shiny.silent.error = function(e) NULL,
                    error = function(e) NULL)
      if (is.null(d) || is.null(d$time)) return(NULL)
      rg <- suppressWarnings(range(d$time, na.rm = TRUE))
      if (any(!is.finite(rg))) return(NULL)
      rg
    }
    store_register(
      store,
      widget_binding("survival_censor", "numeric", label = "Censor time",
        help = paste("Time at which the Kaplan-Meier analysis",
                     "right-censors observations"),
        choices_provider = NULL)
    )
    .sv_root_store <- if (is.null(store$parent)) store else store$parent
    .sv_keep(observeEvent(input$censor, {
      if (!is.null(input$censor))
        store_sync_from_ui(store, "survival_censor", input$censor)
    }, ignoreInit = TRUE))
    # seed once the slider exists (restore-first-wins)
    .sv_seeded <- FALSE
    .sv_keep(observe({
      if (.sv_seeded) return(NULL)
      if (is.null(input$censor)) return(NULL)
      .sv_seeded <<- TRUE
      store_seed(store, list(survival_censor = input$censor))
    }))
    # store -> UI push, gated on the survival view (the slider is
    # renderUI'd only there) and clamped to the rendered range
    .sv_epoch <- store_epoch(store)
    .sv_keep(observe({
      .sv_epoch()
      cv <- store_read(store, "survival_censor")[[1]]
      if (is.null(cv) ||
          is.null(.sv_root_store$pending[[paste0(store$prefix, ".survival_censor")]]))
        return(NULL)
      ck <- tryCatch(reactive_checkpoint(), shiny.silent.error = function(e) FALSE,
                     error = function(e) FALSE)
      if (!isTRUE(ck)) return(NULL)
      rg <- .sv_range()
      if (!is.null(rg))
        cv <- min(max(cv, rg[1]), rg[2])
      updateSliderInput(session, "censor", value = cv)
    }))
  }

  dat <- reactive({
    req(reactive_checkpoint())
    y <- reactive_resp()
    data.frame(
      time = as.numeric(sub("\\+$", "", y)),
      event = as.integer(grepl("\\+$", y)), 
      strata = reactive_strata(),
      stringsAsFactors = FALSE
    )
  })
  
  output$censor_output <- renderUI({
    nm <- max(dat()$time, na.rm = TRUE)
    fluidRow(
      column(12, offset = 0, style='padding-left:5px; padding-right:5px; padding-top:0px; padding-bottom:0px',
             div(style="display: inline-block;vertical-align:top;", h5("Censor at")),
             div(style="padding-left:25px; display: inline-block;vertical-align:top; width:65%;", 
                 sliderInput(ns("censor"), label = NULL, min = min(dat()$time, na.rm = TRUE), max = nm, value = nm)
             )
      ))
  })
  
  output$kmplot <- renderPlot({
    req(input$censor)
    df <- dat()
    i <- which(df$time > input$censor)
    df$time[i] <- input$censor
    df$event[i] <- 0
    fit <- survfit(Surv(time, event) ~ strata, data = df)
    lab <- ""
    if (length(df$strata > 1)) {
      r <- surv_pvalue(fit, data = df, method = "survdiff")
      lab <- paste(r$method, r$pval.txt)
    }
    
    ggsurvplot(fit, data = df, risk.table = TRUE, conf.int = TRUE, pval = lab, surv.median.line = "hv")

  })

  }) # end moduleServer
}