#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' @importFrom ggplot2 geom_vline
ptmotif_ui <- function(id) {

  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Post-Translational Modification (PTM) Motif Analysis"),
      tags$p("PTM motif analysis identifies enriched sequence patterns surrounding post-translational modification sites (e.g., phosphorylation, acetylation, ubiquitination) in your selected proteins. This analysis uses sequence logo visualization to show the consensus amino acid sequence motif and performs statistical testing to determine if certain motifs are over-represented in your dataset compared to background."),
      tags$h4("When to use PTM motif analysis"),
      tags$p("Use this analysis when you have identified PTM sites (typically from mass spectrometry) and want to understand which kinases, enzymes, or regulatory mechanisms might be responsible for those modifications. This helps identify active signaling pathways and regulatory networks. It's particularly valuable for phosphoproteomics data to predict kinase-substrate relationships."),
      tags$h4("How to interpret results"),
      tags$p("The sequence logo shows the consensus motif with letter height indicating amino acid frequency at each position relative to the PTM site (position 0). Larger letters indicate more conserved positions. The enrichment statistics show whether your selected PTM sites have a different motif pattern compared to all PTM sites in the dataset. Significant enrichment suggests specific kinases or regulatory proteins are active in your experimental condition.")
    ),
    triselector_ui(ns("tris_seqlogo")),
    uiOutput(ns("msg_ui"))
  )
}

ptmotif_module <- function(
  id, pdata, fdata, expr, feature_selected, sample_selected, background,
  reactive_status = reactive(NULL), store = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  triset <- reactive({
    req(fdata())
    i <- grep("^SeqLogo\\|", colnames(fdata()), value = TRUE)
    req( length(i) > 0 )
    str_split_fixed(i, "\\|", n = 3)
    })

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # The sequence triselector is driven by store_watch selectors; the
  # first-choice auto-default seeds the store once (unset keys only,
  # restore/user-first-wins) instead of resetting on every fdata change.
  # ------------------------------------------------------------------
  .ptm_store_observers <- list()
  .ptm_keep <- function(obs) {
    .ptm_store_observers[[length(.ptm_store_observers) + 1L]] <<- obs
    invisible(obs)
  }
  xax <- reactiveVal()
  if (!is.null(store)) {
    .ptm_ts <- function() {
      fd <- tryCatch(fdata(), shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.null(fd)) return(NULL)
      i <- grep("^SeqLogo\\|", colnames(fd), value = TRUE)
      if (length(i) == 0L) return(NULL)
      str_split_fixed(i, "\\|", n = 3)
    }
    .ptm_ts1 <- function() {
      ts <- .ptm_ts()
      if (is.null(ts)) character(0) else unique(ts[, 1])
    }
    .ptm_ts2 <- function(a) {
      ts <- .ptm_ts()
      if (is.null(ts) || is.null(a) || !nzchar(a)) character(0)
      else unique(ts[ts[, 1] %in% a, 2])
    }
    .ptm_ts3 <- function(a, b) {
      ts <- .ptm_ts()
      if (is.null(ts)) character(0)
      else {
        i <- rep(TRUE, nrow(ts))
        if (!is.null(a) && nzchar(a)) i <- i & ts[, 1] %in% a
        if (!is.null(b) && nzchar(b)) i <- i & ts[, 2] %in% b
        unique(ts[i, 3])
      }
    }
    kp1 <- paste0(store$prefix, ".xax_analysis")
    kp2 <- paste0(store$prefix, ".xax_subset")
    store_register(
      store,
      widget_binding("xax_analysis", "select", label = "Sequence category",
        help = "Annotation category of the sequence-window column",
        choices_provider = function(v) .ptm_ts1()),
      widget_binding("xax_subset", "select", label = "Sequence subcategory",
        help = "Subcategory within the sequence category",
        depends_on = "xax_analysis",
        choices_provider = function(v) .ptm_ts2(v[[kp1]])),
      widget_binding("xax_variable", "select_cascaded", label = "Sequence column",
        help = paste("Feature annotation column holding the sequence windows",
                     "around modification sites"),
        depends_on = c("xax_analysis", "xax_subset"),
        choices_provider = function(v) .ptm_ts3(v[[kp1]], v[[kp2]]))
    )
    v1 <- triselector_module(
      "tris_seqlogo", reactive_x = triset, label = 'Sequence',
      reactive_selector1 = store_watch(store, "xax_analysis"),
      reactive_selector2 = store_watch(store, "xax_subset"),
      reactive_selector3 = store_watch(store, "xax_variable"),
      reactive_axis_request = store_epoch(store))
    # UI -> store: settled triples only (WP2 store_bind_triselector)
    store_bind_triselector(store,
      keys = c(analysis = "xax_analysis", subset = "xax_subset",
               variable = "xax_variable"),
      sel = v1, keep = .ptm_keep)
    # auto-default to the first available sequence column, but only while
    # the keys are unset (a restore or user pick wins and sticks)
    .ptm_keep(observe({
      ts <- tryCatch(triset(), shiny.silent.error = function(e) NULL,
                    error = function(e) NULL)
      if (is.null(ts) || nrow(ts) == 0L) return(NULL)
      held <- store_read(store, c("xax_analysis", "xax_subset", "xax_variable"))
      patch <- list()
      if (is.null(held[[paste0(store$prefix, ".xax_analysis")]]))
        patch$xax_analysis <- ts[1, 1]
      if (is.null(held[[paste0(store$prefix, ".xax_subset")]]))
        patch$xax_subset <- ts[1, 2]
      if (is.null(held[[paste0(store$prefix, ".xax_variable")]]))
        patch$xax_variable <- ts[1, 3]
      if (length(patch))
        tryCatch(store_apply(store, patch, origin = "system", strict = FALSE),
                 error = function(e) NULL)
    }))
  } else {
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
  }

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

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    if (!is.null(store)) {
      # single transactional path (per-key resilient, meta_scatter style)
      if (identical(length(s$xax), 3L)) {
        tr <- lapply(s$xax, function(x)
          if (is.null(x) || !nzchar(x) || identical(x, "--select--")) NULL else x)
        patch <- list()
        if (!is.null(tr[[1]])) patch$xax_analysis <- tr[[1]]
        if (!is.null(tr[[2]])) patch$xax_subset <- tr[[2]]
        if (!is.null(tr[[3]])) patch$xax_variable <- tr[[3]]
        if (length(patch))
          tryCatch(store_apply(store, patch, origin = "restore", strict = FALSE),
                   error = function(e) NULL)
      }
    } else {
      xax(NULL)
      xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    }
  })

  reactive(list(xax = v1()))

  }) # end moduleServer
}