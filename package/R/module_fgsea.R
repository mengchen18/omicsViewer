#' Utility - fgsea shiny ui
#' @param id id
enrichment_fgsea_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # select variable
    triselector_ui(ns("tris_fgsea")),
    # plotly barplot
    plotlyOutput(ns("bplot")),
    # table
    # DT::dataTableOutput(ns("stab"))
    dataTableDownload_ui(ns("stab"))
  )
}

#' Utility - fgsea shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_featureData reactive feature data
#' @importFrom stringr str_split_fixed
#' @importFrom DT renderDataTable datatable
#' @importFrom fastmatch fmatch
#' @examples 
#' # library(shiny)
#' # source("Git/R/module_triselector.R")
#' # source("Git/R/auxi_fgsea.R")
#' # source("Git/R/module_barplotGsea.R")
#' # 
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # fd <- Biobase::fData(dat)
#' # 
#' # ui <- fluidPage(
#' #   enrichment_fgsea_ui("fgsea")
#' # )
#' # 
#' # server <- function(input, output, session) {
#' #   callModule(enrichment_fgsea_module, id = "fgsea", reactive_featureData = reactive(fd) )
#' # }
#' # 
#' # shinyApp(ui, server)
enrichment_fgsea_module <- function(input, output, session, reactive_featureData) {
  
  ns <- session$ns
  
  # select stats
  triset <- reactive({
    fd <- reactive_featureData()
    cn <- colnames(fd)[sapply(fd, is.numeric) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })
  v1 <- callModule(triselector_module, id = "tris_fgsea", reactive_x = triset, label = "Value")
  
  # run fgsea
  tab <- reactive({
    req( ! v1()$variable %in% c("Select a variable!", ""))
    fdgs <- reactive_featureData()[, grepl("^GS\\|", colnames(reactive_featureData()))]
    stats <- reactive_featureData()[, paste(v1(), collapse = "|")]
    res <- fgsea0(fdgs, stats = stats, nperm = 1000, minSize = 5, 
                  maxSize = 500, gs_desc = sub("GS\\|All\\|", "", colnames(fdgs)))
    list(
      pathway_mat = fdgs,
      table = res,
      stats = stats,
      statsNames = rownames(fdgs)
    )
  })
  
  vi <- callModule(
    dataTableDownload_module, id = "stab", 
    reactive_table = reactive(tab()$table), 
    reactive_cols = reactive(setdiff(colnames(tab()$table), "leadingEdge")), 
    prefix = "fgsea_"
    )
  # output$stab <- DT::renderDataTable({
  #   t0 <- tab()$table[, setdiff(colnames(tab()$table), "leadingEdge")]
  #   ic <- which(sapply(t0, function(x) is.numeric(x) & !is.integer(x)))
  #   t0[, ic] <- lapply(t0[, ic], signif, digits = 3)
  #   DT::datatable(t0, options = list(scrollX = TRUE), rownames = FALSE, selection = "single")
  # })
  
  output$bplot <- renderPlotly({
    
    hid <- bid <- NULL
    
    if (!is.null(i <- vi() )) {
      i <- tab()$table[i, ]
      hid <- fmatch(i$leadingEdge[[1]], tab()$statsNames)
      bid <- setdiff(which(tab()$pathway_mat[, i$pathway] != 0), hid)
    }
    
    plotly_barplot(
      x = tab()$stats, names = tab()$statsNames, 
      highlight = hid, highlight_color = "red", highlight_width = 2, highlight_legend = "Leading edges",
      background = bid, background_color = "gray", background_width = 2, background_legend = "background", 
      ylab = "ylab", xlab = '', sort = "dec", source = ns("plotlybarchart")
    )
    
  })
}

# 
