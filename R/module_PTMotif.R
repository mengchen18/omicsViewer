#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' @importFrom ggplot2 geom_vline
ptmotif_ui <- function(id) {
  
  ns <- NS(id)
  tagList(
    triselector_ui(ns("tris_seqlogo")),
    uiOutput(ns("msg_ui"))
  )
}

ptmotif_module <- function(
  id, pdata, fdata, expr, feature_selected, sample_selected, background
) {

  moduleServer(id, function(input, output, session) {

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

  v1 <- triselector_module(
    "tris_seqlogo", reactive_x = triset, label = 'Sequence',
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
      tabsetPanel(
        tabPanel("Ratio selected/all",
          plotOutput(ns("plt"), height = PTMOTIF_PLOT_HEIGHT),
          tabsetPanel(
            tabPanel("Selected seqs", dataTableDownload_ui(ns("seqtable"))),
            tabPanel("Position weighted matrix", dataTableDownload_ui(ns("seqtable_rat")))            
            )
          ),
        tabPanel("Selected", 
          plotOutput(ns("plt.fg")),
          dataTableDownload_ui(ns("seqtable_fg"))
          ),
        tabPanel("All", 
          plotOutput(ns("plt.bg")),
          dataTableDownload_ui(ns("seqtable_bg"))
          )
        )
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

  output$plt <- renderPlot({
    req( logo() )
    ggseqlogo::ggseqlogo( data = logo() ) + geom_vline(
      xintercept = (ncol(logo())+1)/2, linetype="dashed", color = "orange", size=1.5
      )
  })

  output$plt.fg <- renderPlot({
    req( d <- fg.pfm() )
    ggseqlogo::ggseqlogo( data = d ) + geom_vline(
      xintercept = (ncol(d)+1)/2, linetype="dashed", color = "orange", size=1.5
      )
    })

  output$plt.bg <- renderPlot({
    req( d <- bg.pfm() )
    ggseqlogo::ggseqlogo( data = d ) + geom_vline(
      xintercept = (ncol(d)+1)/2, linetype="dashed", color = "orange", size=1.5
      )
    })

  mat2df <- function(x) {
    data.frame(Name = rownames(x), x, stringsAsFactors = FALSE)
  }

  dataTableDownload_module(
    "seqtable", reactive_table = reactive({
      fg <- foregroundSeqs()
      fg <- fg[which(nchar(fg)>0)]
      do.call(rbind, strsplit(fg, "|"))
    }), prefix = "motif", pageLength = DEFAULT_TABLE_PAGE_LENGTH)

  dataTableDownload_module(
    "seqtable_fg", reactive_table = reactive(mat2df(fg.pfm())), prefix = "seqLogoPFM_foreground", pageLength = DEFAULT_TABLE_PAGE_LENGTH)

  dataTableDownload_module(
    "seqtable_bg", reactive_table = reactive(mat2df(bg.pfm())), prefix = "seqLogoPFM_background", pageLength = DEFAULT_TABLE_PAGE_LENGTH)

  dataTableDownload_module(
    "seqtable_rat", reactive_table = reactive(mat2df(logo())), prefix = "seqLogoPFM_ratio", pageLength = DEFAULT_TABLE_PAGE_LENGTH)

  }) # end moduleServer
}