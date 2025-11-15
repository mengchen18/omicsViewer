#' The three-step selector - the ui function
#' @description Function should only be used for the developers
#' @param id id
#' @param right_margin margin on the right side, in px. For example, "20" translates to "20px".
#' @importFrom shinyWidgets pickerInput updatePickerInput
#' @export
#' @examples
#' if (interactive()) {
#'   library(shiny)
#'   library(Biobase)
#'   
#'   file <- system.file("extdata/demo.RDS", package = "omicsViewer")
#'   dat <- readRDS(file)
#'   fData <- fData(dat)
#'   triset <- stringr::str_split_fixed(colnames(fData), '\\|', n= 3)
#'   
#'   ui <- fluidPage(
#'     triselector_ui("tres"),
#'     triselector_ui("tres2")
#'   )
#'   server <- function(input, output, session) {
#'     v1 <- triselector_module("tres", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("mean.diff")
#'     )
#'     v2 <- triselector_module("tres2", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("log.fdr"))
#'     observe({
#'       print("/////////////////////////")
#'       print(v1())
#'     })
#'   }
#'   
#'   shinyApp(ui, server)
#' }
#' @return a tagList of UI components

triselector_ui <- function(id, right_margin = "20") {
  rmar <- sprintf("padding-left:2px; padding-right:%spx; padding-top:2px; padding-bottom:2px", right_margin)
  ns <- NS(id)
  tagList(
    fluidRow(
      column(2, offset = 0, align = "right",
             style='padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px',
             uiOutput(ns("groupLabel"))
      ),
      column(3, offset = 0, style='padding:2px;',
        selectInput(inputId = ns("analysis"), label = NULL, choices = NULL, selectize = TRUE, width = "100%") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-analysis-selector"))),
      column(4, offset = 0, style='padding:2px;',
        selectInput(inputId = ns("subset"), label = NULL, choices = NULL, selectize = TRUE, width = "100%") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-subset-selector"))),
      column(3, offset = 0, style=rmar, #'padding:2px;',
        selectInput(inputId = ns("variable"), label = NULL, choices = NULL, selectize = TRUE, width = "100%") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-variable-selector")))
      # column(3, offset = 0, style='padding:2px;'
        # pickerInput(inputId = ns("analysis"), label = NULL, choices = NULL, options = list(`live-search` = TRUE))),
      # column(4, offset = 0, style='padding:2px;',
        # pickerInput(inputId = ns("subset"), label = NULL, choices = NULL, options = list(`live-search` = TRUE))),
      # column(3, offset = 0, style="padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px", #'padding:2px;',
        # pickerInput(inputId = ns("variable"), label = NULL, choices = NULL, options = list(`live-search` = TRUE)))
    )
  )
}

#' The three-step selector - the module function
#' @description The selector is used to select columns of phenotype and feature data.
#' Function should only be used for the developers.
#' @param id module id
#' @param reactive_x an nx3 matrix
#' @param reactive_selector1 default value for selector 1
#' @param reactive_selector2 default value for selector 2
#' @param reactive_selector3 default value for selector 3
#' @param label of the triselector
#' @export
#' @examples
#' if (interactive()) {
#'   library(shiny)
#'   library(Biobase)
#'
#'   file <- system.file("extdata/demo.RDS", package = "omicsViewer")
#'   dat <- readRDS(file)
#'   fData <- fData(dat)
#'   triset <- stringr::str_split_fixed(colnames(fData), '\\|', n= 3)
#'
#'   ui <- fluidPage(
#'     triselector_ui("tres"),
#'     triselector_ui("tres2")
#'   )
#'   server <- function(input, output, session) {
#'     v1 <- triselector_module("tres", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("mean.diff")
#'     )
#'     v2 <- triselector_module("tres2", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("log.fdr"))
#'     observe({
#'       print("/////////////////////////")
#'       print(v1())
#'     })
#'   }
#'
#'   shinyApp(ui, server)
#' }
#' @return an reactive object containing the selected values

triselector_module <- function(id,
                               reactive_x,
                               reactive_selector1 = reactive(NULL),
                               reactive_selector2 = reactive(NULL),
                               reactive_selector3 = reactive(NULL),
                               label = "Group Label:") {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # Validate reactive_x input - silently handle NULL during startup
  validated_x <- reactive({
    # Allow NULL during initial app startup
    if (is.null(reactive_x())) {
      return(NULL)
    }

    # Perform validation
    validation <- validate_triselector_input(
      reactive_x(),
      name = "triselector input",
      allow_null = TRUE
    )

    # Only show error notifications after app has initialized
    # (prevents error spam during startup)
    if (!validation$valid && length(names(input)) > 0) {
      showNotification(
        validation$message,
        type = "error",
        duration = 5
      )
      return(NULL)
    }

    # Return validated data or NULL
    if (validation$valid) {
      reactive_x()
    } else {
      NULL
    }
  })

  output$groupLabel <- renderUI({
    h5(HTML(sprintf("<b>%s</b>", label)))
  })

  inte <- reactive({
    aa <- grepl("tris_feature_general", ns("x"))
    length(aa) > 0 && aa
  })
  
  observeEvent(list(reactive_selector1()), {
    req(vx <- validated_x())
    if (length(names(input)) == 0)
      return(NULL)
    cc <- unique(vx[, 1])
    if (!is.null(reactive_selector1()))
      ss <- reactive_selector1() else
        ss <- cc[1]
    updateSelectInput(session, inputId = "analysis", choices = cc, selected = ss)
  })
  
  observeEvent(list(names(input), validated_x()), {
    req(vx <- validated_x())
    if (length(names(input)) == 0)
      return(NULL)
    cc <- unique(vx[, 1])
    if (input$analysis %in% cc)
      ss <- input$analysis else if (!is.null(reactive_selector1()))
        ss <- reactive_selector1() else
          ss <- cc[1]
    updateSelectInput(session, inputId = "analysis", choices = cc, selected = ss)
    # updatePickerInput(session, inputId = "analysis", choices = cc, selected = ss)
  })
  
  # bug fix
  observeEvent(input$analysis, {
    req(vx <- validated_x())
    req (input$analysis == "")
    cc <- unique(vx[, 1])
    if (is.null(reactive_selector1()))
      ss <- cc[1] else
        ss <- reactive_selector1()    
      updateSelectInput(session, inputId = "analysis", choices = cc, selected = ss)
      # updatePickerInput(session, inputId = "analysis", choices = cc, selected = ss)
    })
  
  # updat selectize input when reactive_x is given
  observe({
    req(vx <- validated_x())
    input$analysis
    req(input$analysis)
    cc <- unique(vx[vx[, 1] == input$analysis, 2])
    updateSelectInput(session, inputId = "subset", choices = cc, selected = reactive_selector2())
    # updatePickerInput(session, inputId = "subset", choices = cc, selected = reactive_selector2())
  })
  
  observe({
    input$analysis
    input$subset
    req(input$analysis)
    req(input$subset)
    req(vx <- validated_x())

    cc <- vx[, 3][vx[, 1] == input$analysis & vx[, 2] == input$subset]
    cc <- c("--select--", cc)
    preselected <- try(match.arg(reactive_selector3(), cc), silent = TRUE)
      if (inherits(preselected, "try-error"))
        preselected <- NULL
    updateSelectizeInput(session, inputId = "variable", choices = cc, selected = preselected, server = TRUE)
    # updatePickerInput(session, inputId = "variable", choices = cc, selected = preselected)
  })
  
  reactive({
    req(input$variable)
    list(analysis = input$analysis, subset = input$subset, variable = input$variable)
  })

  }) # end moduleServer
}



