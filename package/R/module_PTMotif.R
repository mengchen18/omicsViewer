#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
ptmotif_ui <- function(id) {
  
  ns <- NS(id)
  tagList(
    plotOutput(ns("plt")),
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
  input, output, session, pdata, fdata, expr, feature_selected, sample_selected, object
) {
  ns <- session$ns
  
  bg.seqs <- reactiveVal(NULL)
  observe({
    req(ic <- grep("^PTMSeq\\|", colnames(fdata())))
    l <- unique(unlist(strsplit(fdata()[, ic], ";")))
    if (length(unique(nchar(l))) > 1) {
      showModal(modalDialog(
        title = "Oops!",
        "The length of sequences is different. The input of PTM motif analysis requires the sequences 
        have the same length, the middle of the sequence should be the modified AA."
      ))
      return()
    }
    bg.seqs(l)
  })
  observeEvent(fdata(), {
    if (exists("__pfm.bg__", envir = .GlobalEnv))
      rm(list = "__pfm.bg__", envir = .GlobalEnv)
  })
  observe({
    req(bg.seqs())
    if (exists("__pfm.bg__", envir = .GlobalEnv))
      return()
    assign("__pfm.bg__", aaFreq( bg.seqs() ), envir = .GlobalEnv)
  })
  
  foregroundSeqs <- reactive({
    req(ic <- grep("^PTMSeq\\|", colnames(fdata())))
    fg.seq <- fdata()[feature_selected(), ic]
    unique(unlist(strsplit(fg.seq, ";")))
  })
  
  logo <- reactive({
    req(foregroundSeqs() >= input$min.seqs )
    req(r0 <- get("__pfm.bg__", envir = .GlobalEnv))
    motifRF(foregroundSeqs(), bg.pfm = r0)
  })
  
  mot <- eventReactive(input$submit, {
    fg.seq <- foregroundSeqs()
    midpos <- (nchar(fg.seq[1])+1)/2
    fg.seq <- fg.seq[substr(fg.seq, midpos, midpos) %in% strsplit(input$cent.res, split = "|")[[1]]]
    print(length(fg.seq))
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