#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
enrichment_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # table
    uiOutput(ns("error")),
    # DT::dataTableOutput(ns("stab")),
    dataTableDownload_ui(ns("stab")),
    dataTableDownload_ui(ns("overlapTab"))
    # plotly barplot
    # plotlyOutput(ns("bplot"))
  )
}

#' @description Utility enrichment analysis shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_i reactive index of rows to be selected (for ORA)
#' @param reactive_featureData reactive feature data
#' @importFrom fastmatch fmatch
#' @examples 
#' #' # source("Git/R/auxi_fgsea.R")
#' # source("Git/R/auxi_vectORA.R")
#' # source("Git/R/module_barplotGsea.R")
# dat <- readRDS("inst/extdata/demo.RDS")
# obj <- tallGS(dat)
# fd <- Biobase::fData(obj)
# fdgs <- attr(fd, "GS")
# selected_ids <- rownames(fd)[fd$`PCA|All|PC1(10.1%)` > 0.02]
# ui <- fluidPage(
#   enrichment_analysis_ui("ea")
# )
# server <- function(input, output, session) {
#   callModule(
#     enrichment_analysis_module, id = "ea",
#     reactive_featureData = reactive(fd), reactive_i = reactive(selected_ids)
#   )
# }
# shinyApp(ui, server)


enrichment_analysis_module <- function(
  input, output, session, reactive_featureData, reactive_i
) {
  
  ns <- session$ns
  
  reactive_pathway <- reactive({
    attr(reactive_featureData(), "GS")
  })
  
  rii <- reactive({
    if (length(reactive_i()) <=  3)
      return("notest")
    if (is.integer(reactive_i() ))
      return(reactive_i())
    # fmatch(reactive_i(), rownames(reactive_pathway_mat()))
    reactive_i()
    })
  
  oraTab <- reactive({
    # print(rii()[1])
    notest <- "No geneset has been tested, please try to include more input feature IDs!" 
    if (rii()[1] == "notest")
      return(notest)
    tab <- vectORATall(reactive_pathway(), i = rii(), background = nrow(reactive_featureData()))
    # print(head(tab))
    if (is.null(tab))
      return(notest)
    ic <- which(sapply(tab, function(x) is.numeric(x) & !is.integer(x)))
    tab[, ic] <- lapply(tab[, ic], signif, digits = 3)
    tab <- tab[which(tab$p.adjusted < 0.1 | tab$p.value < 0.05 | tab$OR >= 3), ]
    tab
  })
  
  output$errorMsg <- renderText({
    req(is.character(oraTab()))
    oraTab()
  })
  
  output$error <- renderUI(
    verbatimTextOutput(ns("errorMsg"))
  )
  
  vi <- callModule(
    dataTableDownload_module, id = "stab", 
    reactive_table = reactive({
      req(is.data.frame(oraTab()))
      oraTab()
    }), 
    reactive_cols = reactive( setdiff(colnames(oraTab()), "overlap_ids") ), 
    prefix = "ORA_", sortBy = "p.value", decreasing = FALSE
  )

  hd <- reactive({
    req(is.data.frame(oraTab()))
    req( i <- vi() )

    ii <- grep("^General", colnames(reactive_featureData()), ignore.case = TRUE)
    if (length(ii) == 0)
      ii <- 1:min(3, ncol(reactive_featureData()))

    i <- oraTab()[i, ]        
    # positive_ids <- i$overlap_ids[[1]]
    # hid <- fmatch(positive_ids, rownames(reactive_pathway_mat()))
    hid <- i$overlap_ids[[1]]
    req(hid)
    df1 <- reactive_featureData()[hid, ii, drop = FALSE]
    df1 <- cbind(Overlap = "+", df1)
    
    apath <- reactive_pathway()[reactive_pathway()$gsId == i$pathway, ]
    aid <- setdiff(apath$featureId, hid)
    # pathway_name <- as.character(i$pathway)
    # aid <- which(reactive_pathway_mat()[, pathway_name] != 0)
    # aid <- setdiff(aid, hid)
    if (length(aid) > 0) {
      df2 <- reactive_featureData()[aid, ii, drop = FALSE]
      df2 <- cbind(Overlap = "", df2)
      df1 <- rbind(df1, df2)
    }    
    df1
    })

  vi2 <- callModule(
    dataTableDownload_module, id = "overlapTab", 
    reactive_table = hd,
    # reactive_cols = reactive( setdiff(colnames(oraTab()), "overlap_ids") ), 
    prefix = "ORA_overlapGenes_"
  )
}
