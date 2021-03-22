#' Independency test for two categorical vectors
#' @param x a vector of character values
#' @param y the second vector of character values, the same length as x
#' @examples 
#' 1
#' # # examples
#' # t1 <- sample(c('A', "B"), size = 30, replace = TRUE)
#' # t2 <- sample(c('a', "b"), size = 30, replace = TRUE)
#' # factorIndependency(t1, t2)
#' # t3 <- sample(c('a', "b", "c", 'd'), size = 30, replace = TRUE)
#' # factorIndependency(t1, t3)

factorIndependency <- function(x, y) {
  
  tab <- table(x, y)
  suppressWarnings( r1 <- chisq.test(tab) )
  r2 <- fisher.test(tab)
  
  df <- data.frame(
    method = c(r1$method, r2$method),
    pvalue = signif(c(r1$p.value, r2$p.value), digits = 3),
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
  )
  
  toMat <- function(x) apply(x, 2, function(i) i)
  
  if (nrow(tab) == 2 && ncol(tab) == 2){
    odds_wald <- (tab[1, 1]/tab[1, 2])/(tab[2, 1]/tab[2, 2]) # conditional maximum likelihood
    odds_fisher <- r2$estimate[[1]] # conditional maximum likelihood
    df$OR <- signif(c(odds_wald, odds_fisher))
    df$OR.method <- c("Wald", "Fisher")
  }
  list(cont.table = data.frame(toMat(tab)),
       residual.ratio = signif(toMat(r1$observed)/toMat(r1$expected), digits = 3),
       p.table = df
  )
}

#' utility - factorIndependency shiny UI
#' @param id id
#' 
factorIndependency_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    wellPanel(
      uiOutput(ns("error")),
      DT::dataTableOutput(ns("count.table.output")),
      DT::dataTableOutput(ns("residual.ratio.output")),
      DT::dataTableOutput(ns("p.table.output")))
  )
}

#' utility - factorIndependency shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param x x, char vector
#' @param y y, char vector, same length as x
#' @param reactive_checkpoint checkpoint
#' @examples
#' 1 #
#' # library(shiny)
#' # 
#' # ui <- fluidPage(
#' #   factorIndependency_ui("cti")
#' # )
#' # 
#' # server <- function(input, output, session) {
#' #   x0 <- reactive( sample(c('A', "B"), size = 30, replace = TRUE) )
#' #   y0 <- reactive( sample(c('a', "b", "c"), size = 30, replace = TRUE) )
#' #   callModule(factorIndependency_module, id = "cti", x = x0, y = y0)
#' # }
#' # 
#' # shinyApp(ui, server)
#' # 
#' # 
#' # t1 <- sample(c('A', "B"), size = 30, replace = TRUE)
#' # t2 <- sample(c('a', "b", "c"), size = 30, replace = TRUE)
#' # u <- factorIndependency(t1, t2)
#' # u <- u$cont.table

factorIndependency_module <- function(
  input, output, session, x, y, reactive_checkpoint = reactive(TRUE)
) {
  
  ns <- session$ns
  
  stats <- reactive({
    req(reactive_checkpoint())
    tx <- table(x())
    ty <- table(y())
    if (length(tx) < 2 || length(ty) < 2 || max(tx) < 2 || max(ty) < 2 || length(tx) > 12 || length(ty) > 12)
      return(NULL)
    factorIndependency(x = x(), y = y()) 
  })
  
  output$error <- renderUI(
    verbatimTextOutput(ns("errorMsg"))
  )
  
  output$errorMsg <- renderText({
    req(is.null(stats()))
    "The selected variable is not suitable for independence test, too many distinct values or no duplicate value!"
  })
  
  output$count.table.output <- DT::renderDataTable({
    req(stats())
    DT::datatable(stats()$cont.table,
                  options = list(searching = FALSE, lengthChange = FALSE, dom = 't', scrollX = TRUE), 
                  rownames = TRUE, class = "compact",
                  caption = "Contingency table")
  })
  
  output$residual.ratio.output <- DT::renderDataTable({
    req(stats())
    DT::datatable(stats()$residual.ratio, 
                  options = list(searching = FALSE, lengthChange = FALSE, dom = 't', scrollX = TRUE), 
                  rownames = TRUE, class = "compact",
                  caption = "Contingency table fold change (observed/expect)"
    )
  })
  
  output$p.table.output <- DT::renderDataTable({
    req(stats())
    DT::datatable(stats()$p.table,
                  options = list(searching = FALSE, lengthChange = FALSE, dom = 't'), 
                  rownames = FALSE, class = "compact",
                  caption = "Significance test of the independence"
    )
  })
}



