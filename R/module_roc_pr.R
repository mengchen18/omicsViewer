
plot_roc_pr_ui <- function(id) {
  ns <- NS(id)
  shinycssloaders::withSpinner(
    plotOutput(ns('roc_pr'))
  )
}

#' Shiny module for ROC/PR plot - Module
#' @param id module id
#' @param reactive_param reactive value; argument pass to draw_roc_pr
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
#'     plot_roc_pr_module("testplot",
#'       reactive_param = reactive(list(
#'         x = ng(),
#'         y = rnorm(100)
#'       ))
#'     )
#'   }
#'  shinyApp(ui, server)
#' }
#' @return do not return any values


plot_roc_pr_module <- function(id, reactive_param, reactive_checkpoint = reactive(TRUE)) {

  moduleServer(id, function(input, output, session) {

  output$roc_pr <- renderPlot({
    req(reactive_checkpoint())
    draw_roc_pr(reactive_param()$y, reactive_param()$x)
  })

  }) # end moduleServer
}



