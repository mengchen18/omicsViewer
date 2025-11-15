#' @description Utility - fgsea shiny ui
#' @param id id
enrichment_fgsea_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Fast Gene Set Enrichment Analysis (fGSEA)"),
      tags$p("fGSEA is a computational method that determines whether predefined gene sets show statistically significant, concordant differences between two biological states. Unlike ORA which only uses a cutoff, fGSEA uses the entire ranked list of genes, making it more sensitive to subtle but coordinated changes in gene expression across a pathway."),
      tags$h4("When to use fGSEA"),
      tags$p("Use fGSEA when you have a ranked list of all genes/proteins (e.g., ranked by fold change or t-statistic) and want to identify which gene sets are enriched at the top or bottom of your ranking. This is especially useful when there are no clear cutoffs or when you want to detect pathway-level changes that might not be apparent from individual gene changes."),
      tags$h4("How to interpret results"),
      tags$p("The Normalized Enrichment Score (NES) indicates the degree of enrichment: positive NES means the gene set is enriched at the top of your ranked list, negative NES means enrichment at the bottom. The p-value and FDR indicate statistical significance. The leading edge genes are the core subset of genes driving the enrichment signal. The bar plot visualizes NES values, with longer bars indicating stronger enrichment.")
    ),
    tags$h3("Ranking Variable Selection", class = "sr-only", `aria-label` = "Select the ranking statistic for gene set enrichment analysis such as t-statistic, log fold change, or p-value"),
    triselector_ui(ns("tris_fgsea"), right_margin = "5"),
    tags$h3("Enrichment Score Visualization", class = "sr-only", `aria-label` = "Bar plot showing normalized enrichment scores with leading edge genes highlighted when a pathway is selected"),
    shinycssloaders::withSpinner(
      plotlyOutput(ns("bplot")),
      type = 8, color = "green"
    ),
    tags$h3("Detailed Results Table", class = "sr-only", `aria-label` = "Table with enrichment scores, p-values, FDR, pathway size, and leading edge genes for all gene sets"),
    dataTableDownload_ui(ns("stab")),
    # Hidden text summary for AI browsers
    div(class = "sr-only", `aria-live` = "polite", `aria-atomic` = "true",
        uiOutput(ns("plotSummary")))
  )
}

#' @description Utility - fgsea shiny module
#' @param id module id
#' @param reactive_featureData reactive feature data
#' @param reactive_status reactive status for restoring saved sessions
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
#   enrichment_fgsea_module("fgsea", reactive_featureData = reactive(fd) )
# }
# shinyApp(ui, server)

enrichment_fgsea_module <- function(id, reactive_featureData, reactive_status = reactive(NULL)) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  triset <- reactive({
    fd <- reactive_featureData()
    cn <- colnames(fd)[vapply(fd, is.numeric, logical(1)) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })

  xax <- reactiveVal()
  v1 <- triselector_module(
    "tris_fgsea", reactive_x = triset, label = "Input variable",
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
      fdgs, stats = stats, minSize = 3, maxSize = 500, 
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
  
  vi <- dataTableDownload_module(
    "stab",
    reactive_table = reactive(tab()$table),
    reactive_cols = reactive(setdiff(colnames(tab()$table), "leadingEdge")),
    prefix = "fgsea_"
  )
  
  output$bplot <- renderPlotly({
    
    hid <- bid <- NULL


    if (!is.null(i <- vi() ) && length(vi()) > 0) {
      i <- tab()$table[i, ]
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

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    })

  reactive(list(xax = v1()))

  }) # end moduleServer
}
