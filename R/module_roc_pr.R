
plot_roc_pr_ui <- function(id) {
  ns <- NS(id)
  shinycssloaders::withSpinner(
    plotOutput(ns('roc_pr'))
  )
}

#' Shiny module for boxplot using plotly - Module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_param_plotly_boxplot reactive value; argument passed to plotly_boxplot
#' @param reactive_checkpoint reactive_value; check this value before render any plot/executing any calculation
#' @examples 
#' if (interactive()) {
#'   library(shiny)
#'   
#'   ui <- fluidPage(
#'     sliderInput("ngrp", label = "Number of groups", min = 2, max = 5, value = 2),
#'     plot_roc_pr_ui("testplot")
#'   )
#'   
#'   server <- function(input, output, session) {
#'     ng <- reactive(
#'       sample(letters[1:input$ngrp], size = 100, replace = TRUE)
#'     )
#'     callModule(
#'       plot_roc_pr_module, id = "testplot",
#'       reactive_param = reactive(list(
#'         x = ng(),
#'         y = rnorm(100)
#'       ))
#'     )
#'   }
#'  shinyApp(ui, server)
#' }
#' @return do not return any values


plot_roc_pr_module <- function(input, output, session, reactive_param, reactive_checkpoint = reactive(TRUE)) {
  
  output$roc_pr <- renderPlot({
    req(reactive_checkpoint())
    draw_roc_pr(reactive_param()$y, reactive_param()$x)
  })
}



