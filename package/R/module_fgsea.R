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
    scc <- paste(v1(), collapse = "|")
    req(scc %in% colnames(reactive_featureData()))
    fdgs <- reactive_featureData()[, grepl("^GS\\|", colnames(reactive_featureData()))]
    stats <- reactive_featureData()[, scc]

    ii <- !is.na(stats)
    stats <- stats[ii]
    fdgs <- fdgs[ii, ]

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
  
  output$bplot <- renderPlotly({
    
    hid <- bid <- NULL
    
    if (!is.null(i <- vi() ) && length(vi()) > 0) {
      i <- tab()$table[i, ]
      hid <- fmatch(i$leadingEdge[[1]], tab()$statsNames)
      bid <- setdiff(which(tab()$pathway_mat[, i$pathway] != 0), hid)
      if (length(bid) == 0)
        bid <- NULL
    }
    
    plotly_barplot(
      x = tab()$stats, names = tab()$statsNames, 
      highlight = hid, highlight_color = "red", highlight_width = 2, highlight_legend = "Leading edges",
      background = bid, background_color = "gray", background_width = 2, background_legend = "background", 
      ylab = "Rankding stats", xlab = '', sort = "dec", source = ns("plotlybarchart")
    )
    
  })
}

# 
