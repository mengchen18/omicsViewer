# barplot utility
# 
# # ## examples
# set.seed(100)
# n <- 100
# s_x <- rnorm(n)
# s_n <- paste(sample(letters, replace = TRUE, size = n), 1:n)
# plotly_barplot(x = s_x, names = s_n)
# plotly_barplot(x = s_x, names = s_n, sort = "inc")
# plotly_barplot(x = s_x, names = s_n, sort = "dec")
# 
# plotly_barplot(x = s_x, names = s_n, highlight = c(5, 6, 9, 10, NA), background = c(15, 18, 30))
# plotly_barplot(x = s_x, names = s_n, highlight = c(5, 6, 9, 10, NA), background = c(15, 18, 30), sort = "inc")
# # 
# dat <- readRDS("Dat/exampleEset.RDS")
# fd <- fData(dat)
# plotly_barplot(fd$`t-test|OV_BR|md`, rownames(fd))
# plotly_barplot(fd$`t-test|OV_BR|md`, rownames(fd), 
#                highlight = c(1, 5, 10), highlight_width = 6,
#                sort = "inc")
# 
# plotly_barplot(fd$`t-test|OV_BR|md`, rownames(fd), 
#                highlight = c(1, 5, 10), highlight_width = 6,
#                background = c(20, 50, 70), background_width = 6,
#                sort = "inc")

plotly_barplot <- function(
  x, names, 
  highlight = NULL, highlight_color = "red", highlight_width = 1, highlight_legend = "highlighted",
  background = NULL, background_color = "gray", background_width = 1, background_legend = "background", 
  ylab = "ylab", xlab = 'xlab', sort = c("none", "increasing", "decreasing")[1], tooltips = NULL,
  source = "plotlybarchart"
  ) {
  
  background <- na.omit(background)
  highlight <- na.omit(highlight)
  data <- data.frame(x = x, names = names, stringsAsFactors = FALSE)
  sort <- match.arg(sort, choices = c("none", "increasing", "decreasing"))
  set.seed(100)
  if (sort == "increasing") {
    data$xpos <- rank(data$x, ties.method = "random")
  } else if (sort == "decreasing") {
    rk <- rank(data$x, ties.method = "random")
    data$xpos <- (max(rk)+1)-rk
  } else
    data$xpos <- 1:nrow(data)
  
  fig <- plot_ly(data, source = source)
  fig <- add_lines(fig, x = data$xpos, y = data$x, type = "lines", line = list(color = "lightgray"), name = "Ranking stats")
  
  st <- (max(data$x, na.rm = TRUE) - min(data$x, na.rm = TRUE))/20
  
  if (!is.null(tooltips))
    txt <- tooltips else
      txt <- names
  
  if (!is.null(highlight)) {
    y <- data$x[highlight]
    y <- sign(y)*pmax(abs(y), st)
    fig <- add_bars(fig, x = data$xpos[highlight], y = y, type = "bar", name = highlight_legend,
                    marker = list(color = highlight_color), opacity = 0.8, width = highlight_width,
                    text = txt[highlight], hoverinfo = 'text')
    fig <- plotly::layout(fig, xaxis= list(
      showticklabels = TRUE, tickvals = data$xpos[highlight], ticktext = data$names[highlight])
      )
  }
  
  if (!is.null(background)) {
    y <- data$x[background]
    y <- sign(y)*pmax(abs(y), st)
    fig <- add_bars(fig, x = data$xpos[background], y = y, type = "bar", name = background_legend,
                    marker = list(color = background_color), opacity = 0.8, width = background_width,
                    text = txt[background], hoverinfo = 'text')
  }
  
  plotly::layout(fig, yaxis = list(title = ylab), xaxis = list(title = xlab), legend = list(
    orientation = 'h', xanchor = "center", x = 0.5, yancher = "center", y = 1.1)
    )
}


###
plotly_barplot_ui <- function(id) {
  ns <- NS(id)
  plotlyOutput(ns("plot"))
}

plotly_barplot_module <- function(input, output, session, ...) {
  
  ns <- session$ns
  output$plot <- renderPlotly(
    do.call(plotly_barplot, args = list(..., source = ns("barplotly")))
  )
  
  reactive(
    event_data("plotly_click", source = ns('barplotly') )
  )
}

# library(shiny)
# dat <- readRDS("Dat/exampleEset.RDS")
# fd <- fData(dat)
# 
# # # examples
# ui <- fluidPage(
#   plotly_barplot_ui("tx")
# )
# server <- function(input, output, session) {
#   set.seed(100)
#   # n <- 100
#   s_x <- fd$`t-test|OV_BR|md`
#   s_n <- rownames(fd)
#   v <- callModule(plotly_barplot_module, id = "tx",
#              x = s_x, names = s_n, highlight <- c(5, 6, 9, 10, NA) , sort = "inc")
#   observe(print(v()))
# }
# shinyApp(ui, server)



