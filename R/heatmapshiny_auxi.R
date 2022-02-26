#' @description convert value to colors
#' @param x a vector of numeric, character or factor
#' @param n number of distinct values, only used when x is a numeric vector
#' @importFrom RColorBrewer brewer.pal
#' 
value2color <- function(x, n=10) {
  if (is.numeric(x)) {
    nl <- as.factor(x)
    if (length(unique(x)) > n)
      nl <- cut(x,n)
    cp <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(n)
    names(cp) <- levels(nl)
    cc <- cp[nl]
  } else if (is.factor(x) || is.character(x) || is.logical(x)) {
    x <- as.character(x)
    nl <- sort(unique(x))
    if (length(nl) > 40) {
      cp <- nColors(k = 40)
      cp <- rep(cp, ceiling(length(nl)/40))[seq_along(nl)]
    } else {
      cp <- nColors(k = length(nl))
    }
    names(cp) <- nl
    cp <- cp[sort(names(cp))]
    cc <- cp[x]
  } else 
    stop("x is not one of: 1) numeric vectoer 2) char vector 3) factor 4) logical!")
  cc[is.na(cc)] <- "white"
  list(color = cc, key = cp)
}

#' @description plot heatmap annotations
#' @param object an object returned by \code{addHeatmapAnnotation}
#' @param ... other parameters passed to plot, usually xlim or ylim

addHeatmapAnnotation_plot <- function( object, ... ) {
  
  x <- object$x
  type <- object$type
  column <- object$column
  var.name <- object$var.name
  ll <- object$key
  mm <- object$mm
  cmat <- object$cmat
  
  minColShow <- 50
  minRowShow <- 100
  borderColor <- NA
  
  if ( object$type == 2 ) {
    if (length(object$key$color) < ifelse(column, minColShow, minRowShow))
      borderColor <- "gray25"
  } else if ( object$type == 3 ) {
    if (ncol(cmat) < ifelse(column, minColShow, minRowShow)) 
      borderColor <- "gray25"
  }
  lim <- list(...)
  if (!is.null(lim$xlim)) {
    if (lim$xlim[2] - lim$xlim[1] < minColShow)
      borderColor <- "gray25"
  }
  if (!is.null(lim$ylim)) {
    if (lim$ylim[2] - lim$ylim[1] < minRowShow)
      borderColor <- "gray25"
  }
  
  if (type == 1) { # single categorical side color
    if (!column)
      x <- -x
    barplot(x, horiz = !column, xaxs = "i", yaxs = "i", space = 0, xpd = FALSE, axes = FALSE, ...)
    if (!column) {
      pp <- par('xaxp')
      at <- seq(pp[1], pp[2], length.out = pp[3]+1)
      lb <- round(rev(seq(-pp[2], -pp[1], length.out = pp[3]+1)), digits = 3)
      axis(side = 1, at = at, labels = lb, las = 2)
    } else {
      axis(side = 4, las = 2)
    }
  } else if (type == 2) { # barplot
    barplot(rep(1, length(ll$color)), col = ll$color, space = 0, border = borderColor,
            xaxs = "i", yaxs = "i", axes = FALSE, horiz = !column, xpd = FALSE, ...)
    mtext(side = ifelse(column, 4, 1), at = 0.5, text = var.name, las = 2, line = 0.5)
  } else if (type == 3) { # multiple categorical
    for (i in nrow(mm):1) 
      barplot(mm[i, ], space = 0, col = cmat[i, ], xaxs = "i", yaxs = "i", axes = FALSE, border = borderColor,
              add = i != nrow(mm), horiz = !column, xpd = FALSE, ...)
    mtext(side = ifelse(column, 4, 1), at = mm[, 1]-0.5, text = rownames(mm), las = 2, line = 0.5)
  } else 
    stop("Unreal type!")
}

#' @description Prepare data for heatmap annotations
#' @param x a vector or matrix used for annotation. nrow of matrix should equal ncol/nrow of heatmap matrix
#' @param column whether the annotations are for columns
#' @param var.name variable name
#' @importFrom matrixStats colCumsums
#' 
addHeatmapAnnotation <- function(x, column = TRUE, var.name="") {
  
  object <- list()
  object$column <- column
  object$var.name <- var.name

  if (is.numeric(x) && is.vector(x) && !is.matrix(x)) {
    object$type <- 1
    object$x <- x
  } else if ((is.character(x) && is.vector(x)) || is.factor(x) || is.logical(x)) {
    ll <- value2color(x)
    object$key <- ll
    object$type <- 2
  } else if (is.matrix(x) || is.data.frame(x)) {
    ll <- lapply(seq_len(ncol(x)), function(i) value2color(x[, i]))
    # cmat <- t(sapply(ll, "[[", "color"))
    cmat <- t(
      vapply( ll, function(x) x[["color"]], character(length(ll[[1]][["color"]])) )
      )
    keylist <- lapply(ll, "[[", "key")
    
    mm <- matrix(1, nrow(cmat), ncol(cmat))
    mm <- colCumsums(mm)
    rownames(mm) <- rownames(cmat) <- names(keylist) <- colnames(x) 
    object$key <- ll
    object$mm <- mm
    object$cmat <- cmat
    object$type <- 3
  }
  object
}

#' @description Prepare data for heatmap annotations
#' @param x a matrix, the distance between rows are calculated
#' @param method method used for distance calculation
#' @importFrom stats dist cor
adist <- function(x, method="pearson") {
  arg.dist <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  arg.cor <- c("spearman", "pearson")
  method <- match.arg(method, choices = c(arg.dist, arg.cor))
  
  if (method %in% arg.dist) 
    dd <- dist(x, method = method) else
      dd <- as.dist(1-cor(t(x), use = "pairwise"))
  dd
}

#' @description Plot heatmap key
#' @param range range of values 
#' @param colors colors 
#' @importFrom grid grid.pretty
#' 
heatmapKey <- function(range, colors) {
  par(mar = c(2, 1, 0, 1))
  bp <- barplot(rep(1, length(colors)), col = colors, space = 0, axes = FALSE, xaxs = "i", yaxs = "i", border = colors)
  lab <- grid.pretty(range)
  at <- (lab - min(range))/(max(range) - min(range))*max(bp)
  mtext(text = lab, side = 1, at = at)
}

#' @description Plot side color keys
#' @param x object
#' @param label labels
#' 
sideCorKey <- function(x, label) {
  par(mar = c(0.1, 0.1, 2, 0.1))
  bp <- barplot(rep(1, length(x$key)), col = x$key, space = 0, xaxs = 'i', yaxs = 'i', 
                horiz = TRUE, border = x$key, axes = FALSE, main = label)
  text(0.5, y = bp, names(x$key))
}

############## tooltips ################
#' @description shiny UI for tooltips
#' @param id id
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' 
shinyPlotTooltipsUI <- function(id) {
  
  ns <- NS(id)
  
  htmlCode <- sprintf(
    "// Get mouse coordinates
    var mouseX, mouseY;
    $(document).mousemove(function(e) {
    mouseX = e.clientX;
    mouseY = e.clientY;
    }).mouseover();
    
    // Function to possition draggable, place on current mouse coordinates
    Shiny.addCustomMessageHandler ('placeDraggable',function (message) {
    var element = $('#%s').parent();
    element.css({'top': mouseY + 'px', 'left' : mouseX + 'px'})
    });", 
    # ns("hover_info"),
    # ns("hover_info"),
    ns("hover_info")
  )
  
  tagList(
    
    tags$head(
      tags$script(
        HTML(htmlCode)
      )
    ),
    absolutePanel(fixed=TRUE, draggable = TRUE, htmlOutput(ns("hover_info"), style = "z-index: 9999;")
    )
  )
}

#' @description Tooltips for shiny plot, Module
#' @param input input
#' @param output output
#' @param session session
#' @param points legend
#' @details modified from https://stackoverflow.com/questions/39346273/how-to-place-shiny-panel-near-the-pointer-mouse/50623683
#' @name tooltipsShinyPlot
#' @examples
#' 1
#' # require(shiny)
#' # require(ggplot2)
#' # ui <- shinyUI(fluidPage(
#' # shinyPlotTooltipsUI("stp"),
#' # plotOutput( "myplot", hover = hoverOpts(id ="myplot_hover") )
#' # ))
#' 
#' # server <- shinyServer(function(input, output, session) {
#'   
#' # output$myplot <- renderPlot({
#' # ggplot(mtcars) + geom_point(aes(mpg,cyl))
#' # })
#' 
#' ## Create reactive variable
#' #points <- reactive({
#' # req(input$myplot_hover)
#' # unlist(nearPoints(mtcars, input$myplot_hover, maxpoints=1))
#' # })
#'   
#' # callModule(shinyPlotTooltips, id = "stp", points = points)
#' # })
#' 
#' # shinyApp(ui, server)
#' 
shinyPlotTooltips <- function(input, output, session, points) {
  
  output$hover_info <- renderText({ 
    req(points())
    session$sendCustomMessage(type = 'placeDraggable', message = list()) 
    paste0(
      '<p style="background-color:rgb(250,250,250); padding:5px; border-radius:3px; border: 2px solid #CCC">', 
      points(), '</p>'
      )
    })
}
