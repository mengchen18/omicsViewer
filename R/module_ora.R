#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
enrichment_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("error")),
    triselector_ui(ns("tris_ora"), right_margin = "5"),
    dataTableDownload_ui(ns("stab")),
    dataTableDownload_ui(ns("overlapTab"))
  )
}

#' @description Utility enrichment analysis shiny module
#' @param id module id
#' @param reactive_featureData reactive feature data
#' @param reactive_i reactive index of rows to be selected (for ORA)
#' @importFrom fastmatch fmatch
#' @importFrom stats cutree
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
#   enrichment_analysis_module("ea",
#     reactive_featureData = reactive(fd), reactive_i = reactive(selected_ids)
#   )
# }
# shinyApp(ui, server)


enrichment_analysis_module <- function(
  id, reactive_featureData, reactive_i
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  reactive_pathway <- reactive({
    attr(reactive_featureData(), "GS")
  })
  
  triset <- reactive({
    trisetter(meta = reactive_featureData(), combine = "none")
  })

  v1 <- triselector_module(
    "tris_ora", reactive_x = triset, label = "Collapse features on"
    )

  size_bg <- reactiveVal()
  rii <- reactiveVal()
  reactive_pathway_collapsed <- reactiveVal( NULL )
  col_key <- reactiveVal( NULL )

  observeEvent(list(
    reactive_i(),
    reactive_featureData(),
    reactive_pathway(),
    v1()
    ), {

    req(reactive_i())
    req(reactive_featureData())
    req(rp <- reactive_pathway())
    req(v1()$variable)

    if (v1()$variable %in% c("", "Select a variable!")) {
      size_bg( nrow(reactive_featureData()) )
      reactive_pathway_collapsed(NULL)
      rii(reactive_i())
      col_key( NULL )
      if ( length(rii()) <= 1 )
        rii(NULL) else if ( length(rii()) <= 3 )
          rii("notest")
      return()
    }
    
    cs <- do.call(paste, list(v1(), collapse = "|"))
    if (!cs %in% colnames(reactive_featureData()))
      return(NULL)
    val <- reactive_featureData()[, cs]
    names(val) <- rownames(reactive_featureData())
    ck <- val[reactive_i()]
    col_key( ck )

    rii(unique(ck))
    if ( length(rii()) <= 1)
      rii(NULL) else if (length(rii()) <=  3)
        rii("notest")

    rp$featureId <- as.factor( val[ as.character( rp$featureId ) ] )
    reactive_pathway_collapsed( unique(rp) )

    size_bg( length(unique(val)) )
    })
  
  OT <- reactive({
  
    req(size_bg())
    req(rii())
  
    notest <- "No geneset has been tested, please try to include more input feature IDs!" 

    if (rii()[1] == "notest")
      return(notest)
  
    if (is.null(reactive_pathway_collapsed()))
      rp <- reactive_pathway() else
        rp <- reactive_pathway_collapsed()
  
    tab <- vectORATall(rp, i = rii(), background = size_bg())
  
    if (is.null(tab))
      return(notest)
  
    ic <- which(vapply(tab, function(x) is.numeric(x) & !is.integer(x), logical(1)))  
    tab[, ic] <- lapply(tab[, ic], signif, digits = 3)  
    tab <- tab[which(tab$p.adjusted < 0.1 | tab$p.value < 0.05 | tab$OR >= 3), ]    
  
    if (nrow(tab) > 3) {    
      hcl <- hclust(jaccardList(tab$overlap_ids))
      cls <- cutree(hcl, h = 0.5)
      tab$desc <- paste("cluster", cls, tab$desc, sep = "_")    
    }    
  
    tab
  })
  
  oraTab <- reactiveVal( NULL )
  observe({
    if (!is.null(rii()))
      oraTab( OT() )
  })
  
  
  output$errorMsg <- renderText({
    req(is.character(oraTab()))
    oraTab()
  })
  
  output$error <- renderUI(
    verbatimTextOutput(ns("errorMsg"))
  )
  
  vi <- dataTableDownload_module(
    "stab",
    reactive_table = reactive({
      req(is.data.frame(oraTab()))
      oraTab()
    }),
    reactive_cols = reactive( setdiff(colnames(oraTab()), "overlap_ids") ),
    prefix = "ORA_", sortBy = "p.value", decreasing = FALSE, pageLength = ENRICHMENT_TABLE_PAGE_LENGTH
  )

  hd <- reactive({
    req(is.data.frame(oraTab()))
    req( i <- vi() )

    ii <- grep("^General", colnames(reactive_featureData()), ignore.case = TRUE)
    if (length(ii) == 0)
      ii <- seq_len( min(3, ncol(reactive_featureData())) )
    i <- oraTab()[i, ]
    hid <- i$overlap_ids[[1]]
    req(hid)
    if (!is.null(ck <- col_key()))
      hid <- names(ck)[ck %in% hid]
    df1 <- reactive_featureData()[hid, ii, drop = FALSE]
    df1 <- cbind(Overlap = "+", df1)
    apath <- reactive_pathway()[reactive_pathway()$gsId == i$pathway, ]
    aid <- setdiff(apath$featureId, hid)
    if (length(aid) > 0) {
      df2 <- reactive_featureData()[aid, ii, drop = FALSE]
      df2 <- cbind(Overlap = "", df2)
      df1 <- rbind(df1, df2)
    }    
    df1
    })

  vi2 <- dataTableDownload_module(
    "overlapTab",
    reactive_table = hd,
    prefix = "ORA_overlapGenes_", pageLength = ENRICHMENT_TABLE_PAGE_LENGTH
  )

  }) # end moduleServer
}
