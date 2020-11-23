sample_general_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(1, style = "margin-top: 0px;", attr4selector_ui(ns("a4_gp"), circle = FALSE)),
      column(11, triselector_ui(ns("tris_sample_general")))
    ),
    uiOutput(ns("sample_general_plot")),
    DT::dataTableOutput(ns('mtab'))
  )
}

sample_general_module <- function(input, output, session, reactive_phenoData, reactive_j = reactive(NULL)) {
  
  ns <- session$ns
  
  triset <- reactive({
    trisetter(meta = reactive_phenoData(), combine = "none")
  })
  v1 <- callModule(triselector_module, id = "tris_sample_general", reactive_x = triset, label = "Compare")
  attr4select <- callModule(
    attr4selector_module, id = "a4_gp", reactive_meta = reactive_phenoData, reactive_triset = triset
  )
  
  pheno <- reactive({
    req(!v1()$variable %in% c("", "Select a variable!"))
    cs <- do.call(paste, list(v1(), collapse = "|"))
    if (!cs %in% colnames(reactive_phenoData()))
      return(NULL)
    val <- reactive_phenoData()[, cs]
    
    if (v1()$analysis == "Surv") {
      type <- "surv"
    } else if (is.numeric(val)) {
      type <- "beeswarm"  
    } else if (is.character(val) || is.facet(val)) {
      type <- "table"
    } else {
      warnings("Unknown type of val: sample_general_module, return NULL!")
      return(NULL)
    }
    list(value = val, type = type)
  })
  
  select <- reactive({
    select <- rep("Unselected", nrow(reactive_phenoData()))
    if (!is.null(reactive_j()))
      select[rownames(reactive_phenoData()) %in% reactive_j()] <- "selected"
    select
  })
  
  output$sample_general_plot <- renderUI({
    if (pheno()$type == "beeswarm")
      return( plotly_scatter_ui(ns("sample_general_beeswarm")) )
    if (pheno()$type == "table")
      return( factorIndependency_ui(ns("sample_general_contab")) )
    if (pheno()$type == "surv")
      return( survival_ui(ns("sample_general_surv")) )
  })
  
  ## beeswarm
  showRegLine <- reactiveVal(FALSE)
  vs_scatter <- callModule(
    plotly_scatter_module, id = "sample_general_beeswarm", 
    reactive_param_plotly_scatter = reactive({
      req(reactive_j())
      req(pheno()$value)
      tooltips <- attr4select$tooltips
      if (is.null(tooltips))
        tooltips <- rownames(reactive_phenoData())
      l <- list(
        x = select(), 
        y = pheno()$value,
        xlab = "", 
        ylab = do.call(paste, list(v1(), collapse = "|")),
        tooltips = tooltips
      )
      l$color <- attr4select$color
      l$shape = attr4select$shape
      l$size = attr4select$size
      l$highlight = attr4select$highlight
      l$lighlightName = attr4select$highlightName
      l
    }), 
    reactive_regLine = showRegLine,
    reactive_checkpoint = reactive(pheno()$type == "beeswarm"))
  observe({
    req(pheno()$type == "beeswarm")
    showRegLine(vs_scatter()$regline)
  })
  
  # cont table stats
  callModule(factorIndependency_module, id = "sample_general_contab", 
             x = select, y = reactive(pheno()$value),
             reactive_checkpoint = reactive(pheno()$type == "table")
  )
  
  ## survival
  callModule(survival_module, id = 'sample_general_surv', 
             reactive_resp = reactive(pheno()$value), reactive_strata = select,
             reactive_checkpoint = reactive(pheno()$type == "surv")
  )
  
  ## table
  
  metatab <- reactive({
    req(reactive_j())
    tab <- reactive_phenoData()
    tab <- tab[, grep("^General\\|", colnames(tab)), drop = FALSE]
    tab <- tab[reactive_j(), , drop = FALSE]
    ic <- sapply(tab, is.numeric) & sapply(tab, is.integer)
    tab[ic] <- lapply(tab[ic], signif, digits = 2)
    colnames(tab) <- sub('General\\|All\\|', "", colnames(tab))
    tab
  })
  
  output$mtab <- DT::renderDataTable(
    DT::datatable(metatab(), options = list(scrollX = TRUE), rownames = FALSE, selection = "single")
  )
  
  
}

# ####
# library(shiny)
# 
# source("Git/R/module_triselector.R")
# source("Git/R/module_scatter.R")
# source("Git/R/module_contTableStats.R")
# source("Git/R/module_survival.R")
# 
# dat <- readRDS("Dat/exampleEset.RDS")
# pd <- pData(dat)
# 
# ui <- fluidPage(
#   sample_general_ui("sample_general")
# )
# 
# server <- function(input, output, session) {
#   callModule(sample_general_module, id = "sample_general", reactive_phenoData = reactive(pd), reactive_j = reactive(sample(rownames(pd), size = 20)) )
# }
# 
# shinyApp(ui, server)
