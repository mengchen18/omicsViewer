enrichment_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # table
    DT::dataTableOutput(ns("stab")),
    # plotly barplot
    plotlyOutput(ns("bplot"))
  )
}

enrichment_analysis_module <- function(
  input, output, session, reactive_pathway_mat, reactive_i
) {
  
  # observe({
  #   print(length(reactive_i()))
  # })
  
  rii <- reactive({
    req(length(reactive_i()) > 3)
    if (is.integer(reactive_i() ))
      return(reactive_i())
    fastmatch::fmatch(reactive_i(), rownames(reactive_pathway_mat()))
    })
  
  oraTab <- reactive({
    tab <- vectORA(reactive_pathway_mat(), i = rii())
    req(tab)
      # return(data.frame("No signicant result" = "Too little input, overlap < minOverlap!"))
    ic <- which(sapply(tab, function(x) is.numeric(x) & !is.integer(x)))
    tab[, ic] <- lapply(tab[, ic], signif, digits = 3)
    tab
  })
  # observe(print(head(oraTab())))

  output$stab <- DT::renderDataTable({
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

# #######
# source("shiny/auxi_fgsea.R")
# source("shiny/auxi_vectORA.R")
# source("shiny/module_barplotGsea.R")
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
