#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
ptmotif_ui <- function(id) {
  
  ns <- NS(id)
  tagList(
    uiOutput(ns("msg_ui")),
    fluidRow(
      column(
        width = 3, textInput(ns("cent.res"), label = "Center residue", width = "100%", value = "STY")
      ),
      column(
        width = 3, sliderInput(ns("min.seqs"), label = "Minimum seqs", min = 3, max = 50, value = 10, width = "100%")
      ),
      column(
        width = 3, textInput(ns("pval.cut"), label = "P-value cutoff", value = 1e-6, width = "100%")
      ),
      column(
        width = 3, offset = 0, style='padding-left:5px; padding-right:5px; padding-top:25px; padding-bottom:5px',
        actionButton(ns("submit"), label = "Run", width = "100%")
      )
    ),
    dataTableDownload_ui(ns("tbl"))
  )
}

ptmotif_module <- function(
  input, output, session, pdata, fdata, expr, feature_selected, sample_selected, background
) {
  ns <- session$ns
  
  bg.seqs <- reactiveVal(NULL)
  errText <- reactiveVal(NULL)
  observe({
    req(fdata())
    req(ic <- grep("^PTMSeq\\|", colnames(fdata())))
    x <- background()
    if (!is.null(x)) {
      errText(NULL)
      bg.seqs(x)
    } else {
      errText(
        "The length of sequences is different. The input of PTM motif analysis requires 
the sequences have the same length, the middle of the sequence should be the modified AA."
      )
      bg.seqs(NULL)
    }
  })

  output$errorMsg <- renderText({
    if (is.null(errText()))
      return(NULL)
    errText() 
    })

  output$msg_ui <- renderUI({
    if (!is.null( errText() ))
      verbatimTextOutput(ns("errorMsg")) else
        plotOutput(ns("plt"))
    })
  
  foregroundSeqs <- reactive({
    req(fdata())
    req(ic <- grep("^PTMSeq\\|", colnames(fdata())))
    fg.seq <- fdata()[feature_selected(), ic]
    unique(unlist(strsplit(fg.seq, ";")))
  })
  
  logo <- reactive({
    req(foregroundSeqs() >= input$min.seqs )
    req(r0 <- bg.seqs())
    motifRF(foregroundSeqs(), bg.pfm = r0)
  })
  
  mot <- eventReactive(input$submit, {
    fg.seq <- foregroundSeqs()
    midpos <- (nchar(fg.seq[1])+1)/2
    fg.seq <- fg.seq[substr(fg.seq, midpos, midpos) %in% strsplit(input$cent.res, split = "|")[[1]]]
    req( length(fg.seq) >= input$min.seqs )
    pc <- try(as.numeric(input$pval.cut))
    req(is.numeric(pc))
    show_modal_spinner(text = "Calculating ...")
    tab <- rmotifx::motifx(
      fg.seqs = fg.seq, bg.seqs = bg.seqs(), central.res = input$cent.res, min.seqs = input$min.seqs, pval.cutoff = pc
      )
    if (is.null(tab))
      tab <- data.frame(
        motif = character(0), score = numeric(0), fg.matches = numeric(0), fg.size = numeric(0),
        bg.matches = numeric(0), bg.size = numeric(0), fold.increase = numeric(0),
        stringsAsFactors = FALSE
      )
    remove_modal_spinner()
    tab
  })
  
  ##
  output$plt <- renderPlot({
    req( logo() )
    ggseqlogo::ggseqlogo( data = logo() )
  })
  
  motifTab <- callModule(
    dataTableDownload_module,    
    id = "tbl", reactive_table = reactive({
      req(mot())
      mot()
    }), prefix = "motif", pageLength = 10)
}