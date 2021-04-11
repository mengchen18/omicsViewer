#' plotly boxplot wrapper
#' @param x a matrix object
#' @param i index of which row in x should be highlighted
#' @param highlight index of which column in x to be highlighted
#' @param ylab y-axis label
#' @param extvar external variable, put on top of the boxplot
#' @param ylab.extvar y-axis label of external variable
#' @importFrom stringr str_pad
#' @importFrom reshape2 melt
#' @examples
#' 1 # not for uesrs
#' # ## examples
#' # library(plotly)
#' # x <- cbind(matrix(rnorm(10000, mean = 3), 1000, 10), matrix(rnorm(20000), 1000, 20))
#' # x[sample(1:length(x), size = 0.3*length(x))] <- NA
#' # rownames(x) <- paste("R", 1:nrow(x), sep = "")
#' # colnames(x) <- paste("C", 1:ncol(x), sep = "")
#' # 
#' # plotly_boxplot(x)
#' # plotly_boxplot(x = x, i = 1)
#' # plotly_boxplot(x = x, i  = c(4, 20, 80))
#' # plotly_boxplot(x, i = 1:10)
#' # plotly_boxplot(x, i = 1:1000)
#' # 
#' # plotly_boxplot(x = x, highlight = c(1, 4, 5, 20))
#' # plotly_boxplot(x = x, i = 1, highlight = c(1, 4, 5, 20))
#' # plotly_boxplot(x = x, i  = c(4, 20, 80), highlight = c(1, 4, 5, 20))
#' # plotly_boxplot(x, i = 1:10, highlight = c(1, 4, 5, 20))
#' # plotly_boxplot(x, i = 1:1000, highlight = c(1, 4, 5, 20))
#' # 
#' # plotly_boxplot(x, i = c(4, 20, 80), extvar = 1:30)

plotly_boxplot <- function(x, i = NULL, highlight = NULL, ylab = "ylab", extvar = NULL, ylab.extvar = "ylab.extvar") {
  
  o.name <- colnames(x)
  colnames(x) <- paste0("B", str_pad(1:ncol(x), pad = "0", nchar(ncol(x))))
  mat_i <- x[i, , drop = FALSE]
  tlpmap <- structure(o.name, names = colnames(x))
  if (length(highlight) == 0)
    highlight <- NULL
  if (length(i) == 0)
    i <- NULL
  
  convertDF <- function(m, c, maxr = 200) {
    if (nrow(m) > maxr)
      m <- apply(m, 2, quantile, probs = seq(0, 1, by = 0.02), na.rm = TRUE)
    df <- melt(m)
    hp <- colnames(m)[c]
    df$Cat <- c("n", "c")[as.integer(df$Var2 %in% hp)+1]
    na.omit(df)
  }
  
  df <- convertDF(x, c = highlight)
  
  fig <- plot_ly(showlegend = FALSE)
  fig <- add_boxplot(
    fig, data = df, x = ~ Var2, y = ~ value, color = ~ Cat, 
    colors = c( "#35c636", "#8a8588"), boxpoints = is.null(i), 
    opacity=ifelse(is.null(i), 1, 0.1), 
    text = paste0("<b>Sample: </b>", tlpmap[df$Var2]), hoverinfo = "text")
  
  if (nrow(mat_i) > 20) { # boxplot
    
    dfi <- convertDF(mat_i, c = highlight, maxr = 100)
    fig <- add_boxplot(fig, data = dfi, x = ~ Var2, y = ~ value, boxpoints = FALSE, color = ~ Cat)
    
  } else if (nrow(mat_i) > 0 ) { # point plot
    
    if (nrow(mat_i) <=3) { #  add line
      for (ii in 1:nrow(mat_i)) {
        fig <- add_lines(fig, x = colnames(mat_i), y = mat_i[ii, ], name = rownames(mat_i)[ii],
                         line = list(color = "rgb(150, 150, 150)", width = 1), showlegend = FALSE)# TRUE
        
      }
    }
    dfi <- dfi <- convertDF(mat_i, c = highlight)
    fig <- add_markers(fig, data = dfi, x = ~ Var2, y = ~ value, color = ~ Cat, 
                       text = paste0("<b>Feature: </b>", dfi$Var1, "<br>",
                                     "<b>Sample: </b>", tlpmap[dfi$Var2]), 
                       hoverinfo = "text")
  }
  fig <- plotly::layout(
    fig, 
    xaxis = list(
      categoryarray = colnames(x), 
      categoryorder = "array", 
      ticktext = o.name,
      tickvals = as.list(1:ncol(x) - 1), 
      title = ""))
  
  if (is.null(extvar))
    return( plotly::layout(fig, yaxis = list(title = ylab)) )
  
    figext <- plot_ly(
      x = colnames(x), y = extvar, type = "bar", showlegend = FALSE
    )
  ff <- subplot(figext, fig, shareX = TRUE, nrows = 2, heights = c(0.2, 0.8), margin = 0, titleY = TRUE)
  plotly::layout(ff, yaxis = list(title = ylab.extvar), yaxis2 = list(title = ylab))
  
}


###
#' utility - boxplot shiny UI
#' @param id id
#' 
plotly_boxplot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    plotlyOutput(ns('boxplotly'))
  )
}

#' utility - boxplot shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_param_plotly_boxplot reactive value; argument passed to plotly_boxplot
#' @param reactive_checkpoint reactive_value; check this value before render any plot/executing any calculation
#' @examples 
#' # # # ### examples
#' # library(shiny)
#' # 
#' # ui <- fluidPage(
#' #   plotly_boxplot_ui("testplotly")
#' # )
#' # 
#' # server <- function(input, output, session) {
#' # 
#' #   x <- cbind(matrix(rnorm(10000, mean = 3), 1000, 10), matrix(rnorm(20000), 1000, 20))
#' #   x[sample(1:length(x), size = 0.3*length(x))] <- NA
#' #   rownames(x) <- paste("R", 1:nrow(x), sep = "")
#' #   colnames(x) <- paste("C", 1:ncol(x), sep = "")
#' #   callModule(plotly_boxplot_module, id = "testplotly",
#' #              reactive_param_plotly_boxplot = reactive(list(
#' #                x = x# , i  = c(4, 20, 80)# , highlight = c(1, 4, 5, 20), extvar = 1:30
#' #                ))
#' #              )
#' # }
#' # 
#' # shinyApp(ui, server)
plotly_boxplot_module <- function(input, output, session, reactive_param_plotly_boxplot, reactive_checkpoint = reactive(TRUE)) {
  
  output$boxplotly <- renderPlotly({
    req(reactive_checkpoint())
    do.call(plotly_boxplot, args = reactive_param_plotly_boxplot()) 
  })
}



