enrichment_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # table
    uiOutput(ns("error")),
    DT::dataTableOutput(ns("stab")),
    # plotly barplot
    plotlyOutput(ns("bplot"))
  )
}

enrichment_analysis_module <- function(
  input, output, session, reactive_pathway_mat, reactive_i
) {
  
  ns <- session$ns
  
  rii <- reactive({
    req(length(reactive_i()) > 3)
    if (is.integer(reactive_i() ))
      return(reactive_i())
    fastmatch::fmatch(reactive_i(), rownames(reactive_pathway_mat()))
    })
  
  oraTab <- reactive({
    tab <- vectORA(reactive_pathway_mat(), i = rii())
    if (is.null(tab))
      return("No geneset has been tested, please try to include more input feature IDs!")
    ic <- which(sapply(tab, function(x) is.numeric(x) & !is.integer(x)))
    tab[, ic] <- lapply(tab[, ic], signif, digits = 3)
    tab
  })
  
  output$errorMsg <- renderText({
    req(is.character(oraTab()))
    oraTab()
  })
  
  output$error <- renderUI(
    verbatimTextOutput(ns("errorMsg"))
  )

  output$stab <- DT::renderDataTable({
    req(is.data.frame(oraTab()))
    DT::datatable(oraTab()[, setdiff(colnames(oraTab()), "overlap_ids")],
                  options = list(scrollX = TRUE), rownames = FALSE, selection = "single")
  })

  hd <- reactive({
    req(i <- input$stab_rows_selected )
    i <- oraTab()[i, ]
    pathway_name <- i$pathway
    aid <- which(reactive_pathway_mat()[, pathway_name] != 0)
    stats <- rep(-1, nrow(reactive_pathway_mat()))
    stats[aid] <- 1
    positive_ids <- i$overlap_ids[[1]]
    hid <- fastmatch::fmatch(positive_ids, rownames(reactive_pathway_mat()))
    bid <- setdiff(rii(), hid)
    list(hid = hid, bid = bid, stats = stats, names = rownames(reactive_pathway_mat()))
  })

  output$bplot <- renderPlotly(
    plotly_barplot(
      x = hd()$stats, names = hd()$names,
      highlight = hd()$hid, highlight_color = "red", highlight_width = 1, highlight_legend = "highlighted",
      background = hd()$bid, background_color = "gray", background_width = 1, background_legend = "background",
      ylab = "ylab", xlab = '', sort = "dec", source = "plotlybarchart"
    )
  )
}

# # #######
# source("Git/R/auxi_fgsea.R")
# source("Git/R/auxi_vectORA.R")
# source("Git/R/module_barplotGsea.R")
# 
# dat <- readRDS("Dat/exampleEset.RDS")
# fd <- fData(dat)
# fdgs <- fd[, grep("^GS\\|", colnames(fd))]
# selected_ids <- which(fd$`PCA|All|PC1(9.2%)` > 0.02 )
# 
# ui <- fluidPage(
#   enrichment_analysis_ui("ea")
# )
# 
# server <- function(input, output, session) {
#   callModule(enrichment_analysis_module, id = "ea",
#              reactive_pathway_mat = reactive(fdgs), reactive_i = reactive(selected_ids)
#   )
# }
# 
# shinyApp(ui, server)
