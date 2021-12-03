#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
ptmotif_ui <- function(id) {
  
  ns <- NS(id)
  tagList(
    triselector_ui(ns("tris_seqlogo")),
    uiOutput(ns("msg_ui")),
    tabsetPanel(
      tabPanel("Selected seqs", dataTableDownload_ui(ns("seqtable"))),
      tabPanel("PFM selected", dataTableDownload_ui(ns("seqtable_fg"))),
      tabPanel("PFM all", dataTableDownload_ui(ns("seqtable_bg"))),
      tabPanel("PFM selected/all", dataTableDownload_ui(ns("seqtable_rat")))#,
      # tabPanel("Motifx", fluidRow(
      #   column(
      #     width = 3, textInput(ns("cent.res"), label = "Center residue", width = "100%", value = "all AA", placeholder = "e.g. all AA, STY, ST")
      #   ),
      #   column(
      #     width = 3, sliderInput(ns("min.seqs"), label = "Minimum seqs", min = 3, max = 50, value = 10, width = "100%")
      #   ),
      #   column(
      #     width = 3, textInput(ns("pval.cut"), label = "P-value cutoff", value = 1e-6, width = "100%")
      #   ),
      #   column(
      #     width = 3, offset = 0, style='padding-left:5px; padding-right:5px; padding-top:25px; padding-bottom:5px',
      #     actionButton(ns("submit"), label = "Run", width = "100%")
      #   )),
      #   dataTableDownload_ui(ns("tbl"))
      #   )
      )
    )
}

ptmotif_module <- function(
  input, output, session, pdata, fdata, expr, feature_selected, sample_selected, background
) {
  ns <- session$ns

  triset <- reactive({    
    req(fdata())    
    i <- grep("^SeqLogo\\|", colnames(fdata()), value = TRUE)    
    req( length(i) > 0 )    
    str_split_fixed(i, "\\|", n = 3)
    })

  xax <- reactiveVal()
  observe({
    xax(list(
      v1 = triset()[1, 1],
      v2 = triset()[1, 2],
      v3 = triset()[1, 3]
      ))
    })

  v1 <- callModule(
    triselector_module, id = "tris_seqlogo", reactive_x = triset, label = 'Sequence',
    reactive_selector1 = reactive(xax()$v1), 
    reactive_selector2 = reactive(xax()$v2), 
    reactive_selector3 = reactive(xax()$v3)
    )

  scc <- reactive({
    req(v1())
    cs <- do.call(paste, list(v1(), collapse = "|"))
    req(cs %in% colnames(fdata()))
    cs
    })

  cleanSeqs <- function(x) {
    x <- unique(unlist(strsplit(x, ";")))
    x[which(nchar(x) > 0)]
  }

  bg.seqs <- reactiveVal(NULL)
  errText <- reactiveVal(NULL)

  observe({
    req(fdata())
    req(scc())
    bg.seqs( cleanSeqs( fdata()[, scc()] ) )
  }) 

  observe({
    req(bg.seqs())
    if (length(unique(nchar(bg.seqs()))) > 1)
      errText(
        "The length of sequences is different. The input of seqLogo analysis requires the sequences have the same length"
      )
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
    req(scc())
    req(feature_selected())    
    req(is.null(errText()))
    cleanSeqs( fdata()[feature_selected(), scc()] )
  })
  
  bg.pfm <- reactive({
    req(bg.seqs())
    aaFreq(bg.seqs())
    })

  fg.pfm <- reactive({
    req(foregroundSeqs())
    aaFreq(foregroundSeqs())
    })

  logo <- reactive({    
    req(bg.pfm())    
    req(fg.pfm())
    motifRF(fg.pfm = fg.pfm(), bg.pfm = bg.pfm())
  })
  
  # mot <- eventReactive(input$submit, {
  #   fg.seq <- foregroundSeqs()
  #   midpos <- (nchar(fg.seq[1])+1)/2
    
  #   ctres <- trimws(input$cent.res)
  #   if (ctres == "all AA")
  #     ctres <- LETTERS
  #   fg.seq <- fg.seq[substr(fg.seq, midpos, midpos) %in% strsplit(ctres, split = "|")[[1]]]
  #   req( length(fg.seq) >= input$min.seqs )
  #   pc <- try(as.numeric(input$pval.cut))
  #   req(is.numeric(pc))
  #   show_modal_spinner(text = "Calculating ...")
  #   tab <- motifx(
  #     fg.seqs = fg.seq, bg.seqs = bg.seqs(), central.res = input$cent.res, min.seqs = input$min.seqs, pval.cutoff = pc
  #     )
  #   if (is.null(tab))
  #     tab <- data.frame(
  #       motif = character(0), score = numeric(0), fg.matches = numeric(0), fg.size = numeric(0),
  #       bg.matches = numeric(0), bg.size = numeric(0), fold.increase = numeric(0),
  #       stringsAsFactors = FALSE
  #     )
  #   remove_modal_spinner()
  #   tab
  # })
  
  ##
  output$plt <- renderPlot({
    req( logo() )
    ggseqlogo::ggseqlogo( data = logo() )
  })

  mat2df <- function(x) {
    data.frame(Name = rownames(x), x, stringsAsFactors = FALSE)
  }

  callModule(
    dataTableDownload_module,    
    id = "seqtable", reactive_table = reactive({
      fg <- foregroundSeqs()
      fg <- fg[which(nchar(fg)>0)]
      do.call(rbind, strsplit(fg, "|"))
    }), prefix = "motif", pageLength = 10)

  callModule(
    dataTableDownload_module,    
    id = "seqtable_fg", reactive_table = reactive(mat2df(fg.pfm())), prefix = "seqLogoPFM_foreground", pageLength = 10)

  callModule(
    dataTableDownload_module,    
    id = "seqtable_bg", reactive_table = reactive(mat2df(bg.pfm())), prefix = "seqLogoPFM_background", pageLength = 10)

  callModule(
    dataTableDownload_module,    
    id = "seqtable_rat", reactive_table = reactive(mat2df(logo())), prefix = "seqLogoPFM_ratio", pageLength = 10)

  
  # motifTab <- callModule(
  #   dataTableDownload_module,    
  #   id = "tbl", reactive_table = reactive({
  #     req(mot())
  #     mot()
  #   }), prefix = "motif", pageLength = 10)
}