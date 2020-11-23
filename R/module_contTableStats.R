# # examples
# t1 <- sample(c('A', "B"), size = 30, replace = TRUE)
# t2 <- sample(c('a', "b"), size = 30, replace = TRUE)
# factorIndependency(t1, t2)
# t3 <- sample(c('a', "b", "c", 'd'), size = 30, replace = TRUE)
# factorIndependency(t1, t3)

# x - always two levels
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

factorIndependency_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    wellPanel(
      DT::dataTableOutput(ns("count.table.output")),
      DT::dataTableOutput(ns("residual.ratio.output")),
      DT::dataTableOutput(ns("p.table.output")))
  )
}

factorIndependency_module <- function(
  input, output, session, x, y, reactive_checkpoint = reactive(TRUE)
) {
  
  stats <- reactive({
    req(reactive_checkpoint())
    tx <- table(x())
    ty <- table(y())
    req(length(tx) > 1 && length(ty) > 1)
    factorIndependency(x = x(), y = y()) 
  })
  
  output$count.table.output <- DT::renderDataTable(
    DT::datatable(stats()$cont.table,
                  options = list(searching = FALSE, lengthChange = FALSE, dom = 't'), 
                  rownames = TRUE, class = "compact",
                  caption = "Contingency table")
  )
  
  output$residual.ratio.output <- DT::renderDataTable(
    DT::datatable(stats()$residual.ratio, 
                  options = list(searching = FALSE, lengthChange = FALSE, dom = 't'), 
                  rownames = TRUE, class = "compact",
                  caption = "Contingency table fold change (observed/expect)"
                  )
  )
  
  output$p.table.output <- DT::renderDataTable(
    DT::datatable(stats()$p.table,
                  options = list(searching = FALSE, lengthChange = FALSE, dom = 't'), 
                  rownames = FALSE, class = "compact",
                  caption = "Significance test of the independence"
                  )
  )
}

# library(shiny)
# 
# ui <- fluidPage(
#   factorIndependency_ui("cti")
# )
# 
# server <- function(input, output, session) {
#   x0 <- reactive( sample(c('A', "B"), size = 30, replace = TRUE) )
#   y0 <- reactive( sample(c('a', "b", "c"), size = 30, replace = TRUE) )
#   callModule(factorIndependency_module, id = "cti", x = x0, y = y0)
# }
# 
# shinyApp(ui, server)
# 
# 
# t1 <- sample(c('A', "B"), size = 30, replace = TRUE)
# t2 <- sample(c('a', "b", "c"), size = 30, replace = TRUE)
# u <- factorIndependency(t1, t2)
# u <- u$cont.table

