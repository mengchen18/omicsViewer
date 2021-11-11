#' @description Utility - fgsea shiny ui
#' @param id id
enrichment_fgsea_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # select variable
    triselector_ui(ns("tris_fgsea")),
    # plotly barplot
    plotlyOutput(ns("bplot")),
    # table
    dataTableDownload_ui(ns("stab"))
  )
}

#' @description Utility - fgsea shiny module
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
#' 
# dat <- readRDS("inst/extdata/demo.RDS")
# obj <- tallGS(dat)
# fd <- fData(obj)
# 
# ui <- fluidPage(
#   enrichment_fgsea_ui("fgsea")
# )
# server <- function(input, output, session) {
#   callModule(enrichment_fgsea_module, id = "fgsea", reactive_featureData = reactive(fd) )
# }
# shinyApp(ui, server)

enrichment_fgsea_module <- function(input, output, session, reactive_featureData, reactive_status = reactive(NULL)) {
  
  ns <- session$ns
  
  # select stats
  triset <- reactive({
    fd <- reactive_featureData()
    cn <- colnames(fd)[sapply(fd, is.numeric) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })

  xax <- reactiveVal()
  v1 <- callModule(
    triselector_module, id = "tris_fgsea", reactive_x = triset, label = "Value",
    reactive_selector1 = reactive(xax()$v1), 
    reactive_selector2 = reactive(xax()$v2), 
    reactive_selector3 = reactive(xax()$v3)
    )
  
  gsInfo <- reactive({
    fdgs <- attr(reactive_featureData(), "GS")
    uniqueGs <- unique(fdgs$gsId)
    names(uniqueGs) <- uniqueGs
    list(gs = fdgs, desc = uniqueGs)
    })
  
  # run fgsea
  tab <- reactive({
    req( ! v1()$variable %in% c("Select a variable!", ""))
    scc <- paste(v1(), collapse = "|")
    req(scc %in% colnames(reactive_featureData()))
    stats <- reactive_featureData()[, scc]
    names(stats) <- rownames(reactive_featureData())
    stats <- na.omit(stats)
    fdgs <- gsInfo()$gs[gsInfo()$gs$featureId %fin% names(stats), ]
    if (nrow(fdgs) < 3) {
      message("Perhaps a problem ... enrichment_fgsea_module")
      return(NULL)
    }
    
    res <- fgsea1(
      fdgs, stats = stats, nperm = 1000, minSize = 3, maxSize = 500, 
      gs_desc = gsInfo()$desc)
    
    cn <- colnames(res)
    cn[cn == "ES"] <- "enrichment score (ES)"
    cn[cn == "NES"] <- "normalized ES"
    cn[cn == "desc"] <- "description"
    colnames(res) <- cn
    list(
      pathway_mat = fdgs,
      table = res[order(abs(res$"normalized ES"), decreasing = TRUE), , drop = FALSE],
      stats = stats,
      statsNames = names(stats)
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
      # hid <- fmatch(i$leadingEdge[[1]], tab()$statsNames)
      hid <- i$leadingEdge[[1]]
      bid <- setdiff(tab()$pathway_mat$featureId[tab()$pathway_mat$gsId == i$pathway], hid)
      if (length(bid) == 0)
        bid <- NULL
      hid <- fmatch(hid, tab()$statsNames)
      bid <- fmatch(bid, tab()$statsNames)
    }
    
    plotly_barplot(
      x = tab()$stats, names = tab()$statsNames, 
      highlight = hid, highlight_color = "red", highlight_width = 2, highlight_legend = "Leading edges",
      background = bid, background_color = "gray", background_width = 2, background_legend = "background", 
      ylab = "Rankding stats", xlab = '', sort = "dec", source = ns("plotlybarchart")
    )
    
  })

  # save status and restore
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    })

  reactive(list(xax = v1()))
}

# 
