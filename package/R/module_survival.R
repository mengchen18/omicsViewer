
survival_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("censor_output")),
    plotOutput(ns("kmplot"))
  )
  
}

#' Utility survival KM module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_resp reponse value, in the format 1, 345, 345+, 23, 45, 355+
#' @param reactive_strata strata variable
#' @param reactive_checkpoint checkpoint
#' @importFrom survminer ggsurvplot surv_pvalue
#' @importFrom survival survfit
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
#   callModule(survival_module, id = 'surv', reactive_resp = reactive(t), reactive_strata = reactive#' (strata))
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
  input, output, session, reactive_resp, reactive_strata, reactive_checkpoint = reactive(TRUE)
) {
  
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
    suppressWarnings(
      ggsurvplot(fit, data = df, risk.table = TRUE, conf.int = TRUE, pval = lab, surv.median.line = "hv")
    )
  })
}



