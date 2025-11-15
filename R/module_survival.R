
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
  id, reactive_resp, reactive_strata, reactive_checkpoint = reactive(TRUE)
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

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