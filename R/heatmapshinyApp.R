#' interactive heatmap
#' @param x a matrix object or ExpressionSet
#' @param fData feature data, ignored if x is an ExpressionSet
#' @param pData phenotype data, ignored if x is an ExpressionSet
#' @export
#' @import Biobase
#' @importFrom matrixStats rowSums2
#' 
iheatmap <- function(x, fData = NULL, pData = NULL, impute = FALSE) {
  
  if (inherits(x, "ExpressionSet") || inherits(x, "xcmsFeatureSet")) {
    fData <- Biobase::fData(x)
    pData <- Biobase::pData(x)
    x <- Biobase::exprs(x)
  }
  
  ir <- unique(
    c(which(matrixStats::rowSums2(!is.na(x)) == 0), 
      which(matrixStats::rowVars(x) == 0))
  )
  if (length(ir) > 0) {
    fData <- fData[-ir, ]
    x <- x[-ir, ]
  }
  
  if (impute) {
    x <- apply(x, 1, function(xx) {
      xx[is.na(xx)] <- min(xx, na.rm = TRUE)*0.9
      xx
    })
    x <- t(x)
  }
  
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        tabsetPanel(
          tabPanel("Parameters", iheatmapInput(id = "test")),
          tabPanel("Legend", iheatmapLegend(id = "test"))
        )
      ),
      mainPanel = mainPanel(
        iheatmapOutput(id = "test")
      )
    )
  )
  
  server <- function(input, output) {
    callModule(
      iheatmapModule, 'test', 
      mat = reactive(x), 
      pd = reactive(pData), 
      fd = reactive(fData))
  }
  
  shinyApp(ui, server)
}



#' iheatmap input
#' @param id id
#' @describeIn iheatmap
#' 
iheatmapInput <- function(id) {
  
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        6, 
        # clustering row/distance/linkage
        selectInput(ns("rowSortBy"), "Sorting rows by", choices = NULL, selectize = TRUE, multiple = FALSE),
        conditionalPanel(
          sprintf("input['%s'] == 'hierarchical cluster'", ns("rowSortBy")), 
          selectInput(ns("clusterRowDist"), "Distance", 
                      choices = c("Pearson correlation", "Euclidean", "Maximum", "Manhattan", "Canberra", "Binary", 
                                  "Minkowski", "Spearman correlation"), 
                      selectize = TRUE),
          selectInput(ns("clusterRowLink"), "Linkage", 
                      choices = c("ward.D", "ward.D2", "single", "complete", "average", 
                                  "mcquitty", "median", "centroid"), 
                      selectize = TRUE)
        ),
        hr(),
        # clustering column/distance/linkage
        selectInput(ns("colSortBy"), "Sorting columns by", choices = NULL, selectize = TRUE),
        conditionalPanel(
          sprintf("input['%s'] == 'hierarchical cluster'", ns("colSortBy")), 
          selectInput(ns("clusterColDist"), "Distance", 
                      choices = c("Pearson correlation", "Euclidean", "Maximum", "Manhattan", "Canberra", "Binary", 
                                  "Minkowski", "Spearman correlation"), 
                      selectize = TRUE),
          selectInput(ns("clusterColLink"), "Linkage", 
                      choices = c("ward.D", "ward.D2", "single", "complete", "average", 
                                  "mcquitty", "median", "centroid"), 
                      selectize = TRUE)
        ),
        hr(),
        # scale none/row/col
        selectInput(ns("scale"), "Scale on", 
                    choices = c("row", "none", "column"), 
                    selectize = TRUE)
      ),
      column(
        6, 
        # row annotations (1 vs 1+; numeric vs categorith)
        selectInput(ns("annotCol"), label = "Column annotations", choices = NULL, multiple = TRUE),
        
        # column annotations
        selectInput(ns("annotRow"), label = "Row annotations", choices = NULL, multiple = TRUE),
        
        # column annotations
        selectInput(ns("tooltipInfo"), label = "Tooltips", choices = NULL, multiple = TRUE),
        hr(),
        
        # margin
        sliderInput(ns("marginBottom"), "Bottom margin", min = 1, max = 20, value = 4),
        sliderInput(ns("marginRight"), "Right margin", min = 1, max = 20, value = 4),
        
        # color pallete
        selectInput(ns("heatmapColors"), label = "Heatmap color panel", choices = c(
          "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn"
        ), selected = "RdYlBu", selectize = TRUE, multiple = FALSE)   
      )
    )
  )
}

#' iheatmap Output
#' @param id id
#' @describeIn iheatmap
#' 
iheatmapOutput <- function(id) {
  ns <- NS(id) 
  tagList(
    fluidRow(
      uiOutput(ns("topleft")),
      uiOutput(ns("topmiddle")),
      uiOutput(ns("column_dend_ui")),
      uiOutput(ns("middleleft")),
      uiOutput(ns("middlemiddle")),
      uiOutput(ns("column_sideColor_ui")),
      uiOutput(ns("row_dend_ui")),
      uiOutput(ns("row_sideColor_ui")),
      uiOutput(ns("heatmap_ui"))
    ),
    shinyPlotTooltipsUI(ns("tlp_heatmap"))
  )
}

#' iheatmap legend
#' @param id id
#' @describeIn iheatmap
#' 
iheatmapLegend <- function(id) {
  ns <- NS(id) 
  tagList(
    # verbatimTextOutput(ns("info")),
    uiOutput(ns("key_heatmap_ui")),
    uiOutput(ns("key_colSideCor_ui")),
    uiOutput(ns("key_rowSideCor_ui"))
  )
}

#' iheatmap module
#' @param input input
#' @param output output
#' @param session session
#' @name iheatmap
#' 
iheatmapModule <- function(input, output, session, mat, pd, fd) {
  
  ns <- session$ns
  
  clsRow <- reactive({
    if (nrow(mat()) > 750)
      return("none")
    "hierarchical cluste"
  })
  ######## update selectize input ########
  observe( updateSelectInput(session, "annotCol", choices = colnames(pd())) )
  observe( updateSelectInput(session, "annotRow", choices = colnames(fd())) )
  observe( updateSelectInput(session, "colSortBy", choices = c("hierarchical cluster", "none", colnames(pd()))) )
  observe( updateSelectInput(session, "rowSortBy", choices = c("hierarchical cluster", "none", colnames(fd())), selected = clsRow()) )
  observe( updateSelectInput(session, "tooltipInfo", choices = c(colnames(fd()), colnames(pd())) ) )
  
  
  ######## prepare heatmap data ########
  mm <- reactive({
    
    if (input$scale == "row") {
      mm <- t(scale(t(mat()))) 
      brk <- c(min(mm, na.rm = TRUE), seq(-2, 2, length.out = 99), max(mm, na.rm = TRUE))
    } else if (input$scale == "column") {
      mm <- scale(mat())
      brk <- c(min(mm, na.rm = TRUE), seq(-2, 2, length.out = 99), max(mm, na.rm = TRUE))
    } else {
      mm <- mat()
      brk <- seq(min(mm, na.rm = TRUE), max(mm, na.rm = TRUE), length.out = 101)
    }
    
    list(mat = mm, breaks = brk)
  })
  
  rowSB <- eventReactive(list(input$rowSortBy, mm()$mat), {
    req(input$rowSortBy)
    hcl_r <- NULL
    ord_r <- 1:nrow(mm()$mat)
    if (!input$rowSortBy %in% c("", "none", "hierarchical cluster")) {
      ord_r <- order(fd()[, input$rowSortBy])
    } else if (input$rowSortBy == "hierarchical cluster") { #clusterRow
      dd <- tolower(strsplit(input$clusterRowDist, " ")[[1]][1])
      hcl_r <- hclust(adist(mm()$mat, method = dd), method = input$clusterRowLink)
      ord_r <- hcl_r$order 
      hcl_r <- as.dendrogram(hcl_r)
    }
    list(ord = ord_r, hcl = hcl_r)
  })
  
  colSB <- eventReactive(list(input$colSortBy, mm()$mat), {
    input$colSortBy
    hcl_c <- NULL
    ord_c <- 1:ncol(mm()$mat)
    if (!input$colSortBy %in% c("", "none", "hierarchical cluster")) {
      ord_c <- order(pd()[, input$colSortBy])
    } else if (input$colSortBy == "hierarchical cluster") {
      dd <- tolower(strsplit(input$clusterColDist, " ")[[1]][1])
      hcl_c <- hclust(adist(t(mm()$mat), method = dd), method = input$clusterColLink)
      ord_c <- hcl_c$order
      hcl_c <- as.dendrogram(hcl_c)
    }
    list(ord = ord_c, hcl = hcl_c)
  })
  
  hm <- reactive({
    mmo <- mm()$mat[rowSB()$ord, colSB()$ord]
    list(
      dend_c = colSB()$hcl,
      dend_r = rowSB()$hcl,
      mat = t(mmo),
      ord_c = colSB()$ord,
      ord_r = rowSB()$ord,
      brk =  mm()$breaks
    )
  })
  ######## render plot ########
  dat_colSideCol <- reactive({
    req(input$annotCol)
    addHeatmapAnnotation(pd()[hm()$ord_c, input$annotCol],  var.name = input$annotCol)
  })
  output$colSideCol <- renderPlot({
    par(mar = c(0, 0, 0, input$marginRight))
    addHeatmapAnnotation_plot( dat_colSideCol(), xlim = ranges$x-0.5)
  })
  
  dat_rowSideCol <- reactive({
    req(input$annotRow)
    addHeatmapAnnotation(fd()[hm()$ord_r, input$annotRow], column = FALSE, var.name = input$annotRow)
  })
  output$rowSideCol <- renderPlot({
    par(mar= c(input$marginBottom, 0, 0, 0))
    addHeatmapAnnotation_plot(dat_rowSideCol(), ylim = ranges$y-0.5)
  })
  
  output$heatmap <- renderPlot({
    par(mar = c(input$marginBottom, 0, 0, input$marginRight))
    image(hm()$mat, x = 1:nrow(hm()$mat), y = 1:ncol(hm()$mat), xlim = ranges$x, ylim = ranges$y,
          col = heatColor(), axes = FALSE, xlab = "", ylab = "", breaks = hm()$brk)
    
    if (ranges$x[2] - ranges$x[1] <= 60) {
      irn <- .rg(ranges$x, rownames(hm()$mat))
      mtext(side = 1, at = irn$at, text = irn$lab, las = 2, line = 0.5)
      abline(v = seq(ranges$x[1], ranges$x[2], by = 1), col = "white")
    }
    
    if (ranges$y[2] - ranges$y[1] <= 100) {
      rrn <- .rg(ranges$y, colnames(hm()$mat))
      abline(h = seq(ranges$y[1], ranges$y[2], by = 1), col = "white")
      mtext(side = 4, at = rrn$at, text = rrn$lab, las = 2, line = 0.5)
    }
    
  })
  
  output$dendCol <- renderPlot({
    req(hm()$dend_c)
    par(mar = c(0, 0, 1, input$marginRight))
    plot(hm()$dend_c, xaxs="i", yaxs = "i", axes = FALSE, xlim = ranges$x, center = TRUE)
    axis(side = 4)
  })
  
  output$dendRow <- renderPlot({
    req(hm()$dend_r)
    par(mar = c(input$marginBottom, 1, 0, 0))
    plot(hm()$dend_r, horiz = TRUE, yaxs="i", xaxs = "i", axes = FALSE, ylim = ranges$y)
    axis(side = 1)
  })
  
  ######## update range - heatmap ########
  ranges <- reactiveValues(x = NULL, y = NULL)
  observe({
    ranges$x <- c(0, nrow(hm()$mat))+0.5
    ranges$y <- c(0, ncol(hm()$mat))+0.5
  })
  
  .rg <- function(x, tx) {
    v1 <- ceiling(x[1])
    v2 <- floor(x[2])
    list(at = v1:v2, lab = tx[v1:v2])
  }
  heatColor <- reactive({
    colorRampPalette(
      rev(RColorBrewer::brewer.pal(n = 7, name = input$heatmapColors))
    )(100)
  })
  
  observeEvent(input$heatmap_dblclick, {
    brush <- input$heatmap_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin)-0.5, 0.5), min(round(brush$xmax)+0.5, nrow(hm()$mat) + 0.5))
      ranges$y <- c(max(floor(brush$ymin)-0.5, 0.5), min(round(brush$ymax)+0.5, ncol(hm()$mat) + 0.5))
      # ranges$x <- c(max(round(brush$xmin)-0.5, 0.5), min(round(brush$xmax)+0.5, nrow(hm()$mat)))
      # ranges$y <- c(max(round(brush$ymin)-0.5, 0.5), min(round(brush$ymax)+0.5, ncol(hm()$mat)))
    } else {
      ranges$x <- c(0, nrow(hm()$mat))+0.5
      ranges$y <- c(0, ncol(hm()$mat))+0.5
    }
  })
  
  ######## update range - dendrogram ########
  observeEvent(input$dendRow_dblclick, {
    brush <- input$dendRow_brush
    if (!is.null(brush)) {
      ranges$y <- c(max(floor(brush$ymin)-0.5, 0.5), min(ceiling(brush$ymax)+0.5, ncol(hm()$mat) + 0.5))
    } else {
      ranges$y <- c(0, ncol(hm()$mat))+0.5
    }
  })
  observeEvent(input$dendCol_dblclick, {
    brush <- input$dendCol_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin)-0.5, 0.5), min(ceiling(brush$xmax)+0.5, nrow(hm()$mat) + 0.5))
    } else {
      ranges$x <- c(0, nrow(hm()$mat))+0.5
    }
  })
  
  ######## update range - sidebar ########
  observeEvent(input$rowSideCol_dblclick, {
    brush <- input$rowSideCol_brush
    if (!is.null(brush)) {
      ranges$y <- c(max(floor(brush$ymin)-0.5, 0.5), min(ceiling(brush$ymax)+0.5, ncol(hm()$mat) + 0.5))
    } else {
      ranges$y <- c(0, ncol(hm()$mat))+0.5
    }
  })
  observeEvent(input$colSideCol_dblclick, {
    brush <- input$colSideCol_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin)-0.5, 0.5), min(ceiling(brush$xmax)+0.5, nrow(hm()$mat) + 0.5))
    } else {
      ranges$x <- c(0, nrow(hm()$mat))+0.5
    }
  })
  
  ######### render keys ##########
  output$key_heatmap <- renderPlot(
    heatmapKey(range(hm()$mat, na.rm = TRUE), heatColor())
  )
  output$key_heatmap_ui <- renderUI({ 
    plotOutput(ns("key_heatmap"), height = "45px")
  })
  
  output$key_colSideCor <- renderPlot({
    lg <- dat_colSideCol()
    req(lg$key)
    if (lg$type == 2) {
      sideCorKey(x = lg$key, label = lg$var.name)
    } else {
      layout(matrix(1:length(lg$key), 1))
      for (i in 1:length(lg$key))
        sideCorKey(x = lg$key[[i]], label = lg$var.name[i])
    }
  })
  output$key_colSideCor_ui <- renderUI({
    if (is.null(dat_colSideCol()$key))
      return()
    tagList(
      h4("Column bars"),
      plotOutput(ns("key_colSideCor"), height = 266)
    )
  })
  
  output$key_rowSideCor <- renderPlot({
    lg <- dat_rowSideCol()
    req(lg$key)
    if (lg$type == 2) {
      sideCorKey(x = lg$key, label = lg$var.name)
    } else {
      layout(matrix(1:length(lg$key), 1))
      for (i in 1:length(lg$key))
        sideCorKey(x = lg$key[[i]], label = lg$var.name[i])
    }
  })
  output$key_rowSideCor_ui <- renderUI({
    if (is.null(dat_rowSideCol()$key))
      return()
    tagList(
      h4("Row bars"),
      plotOutput(ns("key_rowSideCor"), height = 266)
    )
  })
  
  
  
  ######### place holder plots ##########
  output$empty1 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty2 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty3 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty4 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  
  ######### dynamic UI render; dynamic layout ##########
  show_col_dend <- reactiveVal(TRUE)
  observe( show_col_dend(input$colSortBy == "hierarchical cluster") )
  
  show_col_sideColor <- reactiveVal(FALSE)
  observe( show_col_sideColor(length(input$annotCol) != 0) )
  
  show_row_dend <- reactiveVal(TRUE)
  observe( show_row_dend(input$rowSortBy == "hierarchical cluster") ) #clusterRow
  
  show_row_sideColor <- reactiveVal(FALSE)
  observe( show_row_sideColor(length(input$annotRow) != 0) )
  
  width_left <- reactive( 2 )
  width_mid <- reactive ({ 
    if (length(input$annotRow) <= 3)
      r <- 1 else if (length(input$annotRow) < 6)
        r <- 2 else
          r <- 3
  })
  width_right <- reactive({ 
    if (show_row_sideColor() & show_row_dend()) {
      wid <- 12 - width_mid() - width_left()
    } else if (!show_row_sideColor() & show_row_dend()) {
      wid <- 12 - width_left()
    } else if (show_row_sideColor() & !show_row_dend()) {
      wid <- 12 - width_mid()
    } else
      wid <- 12
    wid
  })
  
  # top left 
  height_colSideCol <- reactive({
    if ( dat_colSideCol()$type == 1 )
      return("50px")
    paste0(20*length(input$annotCol), "px")
  })
  
  output$topleft <- renderUI({
    if (!show_col_dend() || !show_row_dend())
      return(NULL)
    column(width_left(), plotOutput(ns("empty1"), height = "100px" ), offset = 0, style='padding:0px;')
  })
  # top middle
  output$topmiddle <- renderUI({
    if (!show_col_dend() || !show_row_sideColor())
      return(NULL)
    column(width_mid(), plotOutput(ns("empty2"), height = "100px" ), offset = 0, style='padding:0px;')
  })
  # top right
  output$column_dend_ui <- renderUI({
    if (!show_col_dend())
      return(NULL)
    column(width_right(), plotOutput(ns("dendCol"),
                                     dblclick = ns("dendCol_dblclick"),
                                     brush = brushOpts(id = ns("dendCol_brush"), resetOnNew = TRUE),
                                     height = "100px" ), offset = 0, style='padding:0px;')
  })
  
  # middle left
  output$middleleft <- renderUI({ 
    if (!show_col_sideColor() || !show_row_dend() )
      return(NULL)
    column(width_left(), plotOutput(ns("empty3"), height = height_colSideCol() ), offset = 0, style='padding:0px;')
  })
  # middle middle
  output$middlemiddle <- renderUI({
    if (!show_col_sideColor() || !show_row_sideColor() )
      return(NULL)
    column(width_mid(), plotOutput(ns("empty4"), height = height_colSideCol() ), offset = 0, style='padding:0px;')
  })
  # middle right
  output$column_sideColor_ui <- renderUI({
    if ( !show_col_sideColor() )
      return(NULL)
    column(width_right(), plotOutput(ns("colSideCol"), 
                                     click = ns("colSideCol_click"), 
                                     dblclick = ns("colSideCol_dblclick"), 
                                     hover = hoverOpts(id = ns("colSideCol_hover"), delay = 150), 
                                     brush = brushOpts(id = ns("colSideCol_brush"), resetOnNew = TRUE),
                                     height = height_colSideCol() ), offset = 0, style='padding:0px;')
  })
  # bottom left
  output$row_dend_ui <- renderUI({
    if ( ! show_row_dend() )
      return(NULL)
    column(width_left(), plotOutput(ns("dendRow"), 
                                    dblclick = ns("dendRow_dblclick"), 
                                    brush = brushOpts(id = ns("dendRow_brush"), resetOnNew = TRUE),
                                    height = "800px" ), offset = 0, style='padding:0px;')
  })
  # bottom middle
  output$row_sideColor_ui <- renderUI({
    if ( !show_row_sideColor() )
      return(NULL)
    column(width_mid(), plotOutput(ns("rowSideCol"), 
                                   click = ns("rowSideCol_click"), 
                                   dblclick = ns("rowSideCol_dblclick"), 
                                   hover = hoverOpts(id = ns("rowSideCol_hover"), delay = 150), 
                                   brush = brushOpts(id = ns("rowSideCol_brush"), resetOnNew = TRUE),
                                   height = "800px" ), offset = 0, style='padding:0px;')
  })
  # bottom right
  output$heatmap_ui <- renderUI({
    column(width_right(),  
           plotOutput(ns("heatmap"), click = ns("heatmap_click"), 
                      dblclick = ns("heatmap_dblclick"), 
                      hover = hoverOpts(id = ns("heatmap_hover"), delay = 150), 
                      brush = brushOpts(id = ns("heatmap_brush"), resetOnNew = TRUE),
                      height = "800px"),
           offset = 0, style='padding:0px;')
    
  })
  
  ############### tooltips ################
  
  dat_tooltip <- reactive({
    list(
      pd = pd()[hm()$ord_c, colnames(pd()) %in% input$tooltipInfo, drop = FALSE],
      fd = fd()[hm()$ord_r, colnames(fd()) %in% input$tooltipInfo, drop = FALSE]
    )
  })
  tooltip_mouse <- reactiveVal(NULL)
  # col side color tooltips
  observeEvent(list(input$annotCol, input$annotRow), { 
    tooltip_mouse(NULL)
  })
  observe({
    res <- NULL
    req(input$colSideCol_hover)
    req(input$colSideCol_click)
    hover_x <- ceiling(input$colSideCol_hover$x)
    hover_y <- ceiling(input$colSideCol_hover$y)
    x <- ceiling(input$colSideCol_click$x)
    y <- ceiling(input$colSideCol_click$y)
    if (hover_x == x && hover_y == y && dat_colSideCol()$type %in% 2:3) {
      if (dat_colSideCol()$type == 2) {
        varname <- dat_colSideCol()$var.name 
        key <- dat_colSideCol()$key
        leg <- names(dat_colSideCol()$key$color)[x]
      } else {
        varname <- rownames(dat_colSideCol()$cmat)[y]
        key <- dat_colSideCol()$key[[match(varname, dat_colSideCol()$var.name)]]
        color <- match(dat_colSideCol()$cmat[y, x], dat_colSideCol()$cmat[y, ])
        leg <- names(key$color)[color]
      }
      res <- sprintf("<b>%s:</b> %s", varname, leg)
    }
    tooltip_mouse(res)
  })
  
  # row side color tooltips
  observe({
    res <- NULL
    req(input$rowSideCol_hover)
    req(input$rowSideCol_click)
    hover_x <- ceiling(input$rowSideCol_hover$x)
    hover_y <- ceiling(input$rowSideCol_hover$y)
    x <- ceiling(input$rowSideCol_click$x)
    y <- ceiling(input$rowSideCol_click$y)
    if (hover_x == x && hover_y == y && dat_rowSideCol()$type %in% 2:3) {
      if (dat_rowSideCol()$type == 2) {
        varname <- dat_rowSideCol()$var.name 
        key <- dat_rowSideCol()$key
        leg <- names(dat_rowSideCol()$key$color)[y]
      } else {
        varname <- rownames(dat_rowSideCol()$cmat)[x]
        key <- dat_rowSideCol()$key[[match(varname, dat_rowSideCol()$var.name)]]
        color <- match(dat_rowSideCol()$cmat[x, y], dat_rowSideCol()$cmat[x, ])
        leg <- names(key$color)[color]
      }
      res <- sprintf("<b>%s:</b> %s", varname, leg)
    }
    tooltip_mouse(res)
  })
  
  clickedName <- reactiveVal(NULL)
  # heatmap tooltips
  observe({
    res <- NULL
    req(input$heatmap_hover)
    req(input$heatmap_click)
    hover_x <- round(input$heatmap_hover$x)
    hover_y <- round(input$heatmap_hover$y)
    x <- round(input$heatmap_click$x)
    y <- round(input$heatmap_click$y)
    
    if (hover_x == x && hover_y == y ) {
      
      l <- c(dat_tooltip()$pd[x, , drop = FALSE], dat_tooltip()$fd[y, , drop = FALSE])
      tl <- sapply(names(l), function(x) {
        if (is.numeric(cv <- l[[x]]))
          cv <- round(cv, digits = 2)
        sprintf("<b>%s:</b> %s", x, cv)
      })
      tl <- paste(tl, collapse = "<br>")
      res <- sprintf(
        "<b>Sample:</b> %s <br> <b>Feature:</b> %s",
        rownames(hm()$mat)[x], colnames(hm()$mat)[y]
      )
      res <- paste(res, tl, sep = "<br>")
    }
    
    clickedName( c(col = rownames(hm()$mat)[x], row = colnames(hm()$mat)[y] ) )
    tooltip_mouse(res)
  })
  
  callModule(shinyPlotTooltips, id = "tlp_heatmap", points = tooltip_mouse) 
  
  ############### return ##############
  
  brushedValues <- reactive({
    brush <- input$heatmap_brush
    l <- list(col = NULL, row = NULL)
    if (!is.null(brush)) {
      x1 <- ceiling(max(round(brush$xmin), 0.5))
      x2 <- floor(min(round(brush$xmax), nrow(hm()$mat)))
      y1 <- ceiling(max(round(brush$ymin), 0.5))
      y2 <- floor(min(round(brush$ymax), ncol(hm()$mat)))
      l <- list(
        row = colnames(hm()$mat)[y1:y2],
        col = rownames(hm()$mat)[x1:x2]
      )
    }
    l
  })
  
  reactive({
    list(clicked = clickedName(),
         brushed = brushedValues())
  })
}

