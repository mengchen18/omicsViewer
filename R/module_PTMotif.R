#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' @importFrom ggplot2 geom_vline

# ggseqlogo 0.2.2 uses ggplot2 idioms that are deprecated upstream:
# aes_string() (deprecated in ggplot2 3.0.0) and guides(scale = FALSE)
# (deprecated in ggplot2 3.5.0). Until that is fixed upstream, mute
# exactly those lifecycle warnings at our call sites; any other warning
# still surfaces.
.ggseqlogo_quiet <- function(data) {
  withCallingHandlers(
    ggseqlogo::ggseqlogo(data = data),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("aes_string", msg, fixed = TRUE) ||
          grepl("argument of `guides()`", msg, fixed = TRUE))
        invokeRestart("muffleWarning")
    }
  )
}

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
  reactive_status = reactive(NULL), store = NULL, context = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # ------------------------------------------------------------------
  # WP-5 (panel debt, D13): the panel context, built at the L1 call
  # site (WP-11 collapses the per-panel builds into the registry loop).
  # Context wins when present; the legacy formals remain the
  # direct-instantiation path (the string/dose precedent). This module
  # consumes store + keep from the context; fdata/feature_selected/
  # reactive_status keep arriving through their formals (nothing moved
  # in this WP).
  # ------------------------------------------------------------------
  if (!is.null(context)) {
    store <- context$store %||% store
  }

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
  #
  # WP-5 (panel debt): retention is the shared per-instance panel_keep()
  # factory (D3 - the context's closure when context-built, else a fresh
  # local one for direct instantiation), and the cascade boilerplate
  # (guarded SeqLogo triset source, ts1/ts2/ts3 choices providers, the
  # three xax_* registrations, triselector wiring through
  # store_watch/store_epoch and the UI->store sync) is the shared
  # panel_triset_bind helper (D4). Labels/help stay panel-owned; the
  # registered keys/kinds/labels/help/depends_on are byte-identical to
  # the hand-rolled block (frozen pre-migration golden,
  # tests/test_panelRegistry.R).
  # ------------------------------------------------------------------
  keep <- if (is.null(context)) panel_keep() else context$keep
  xax <- reactiveVal()
  if (!is.null(store)) {
    v1 <- panel_triset_bind(
      store, keep,
      tri_id = "tris_seqlogo", tri_label = "Sequence",
      triset_source = function() {
        i <- grep("^SeqLogo\\|", colnames(fdata()), value = TRUE)
        if (length(i) == 0L) return(NULL)
        str_split_fixed(i, "\\|", n = 3)
      },
      bindings = list(
        analysis = list(
          label = "Sequence category",
          help = "Annotation category of the sequence-window column"),
        subset = list(
          label = "Sequence subcategory",
          help = "Subcategory within the sequence category"),
        variable = list(
          label = "Sequence column",
          help = paste("Feature annotation column holding the sequence windows",
                       "around modification sites"))))
    # auto-default to the first available sequence column, but only while
    # the keys are unset (a restore or user pick wins and sticks)
    keep(observe({
      ts <- tryCatch(triset(), shiny.silent.error = function(e) NULL,
                    error = function(e) NULL)
      if (is.null(ts) || nrow(ts) == 0L) return(NULL)
      store_seed(store, list(xax_analysis = ts[1, 1],
                             xax_subset = ts[1, 2],
                             xax_variable = ts[1, 3]))
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
    # R-M10: clear the message on a good column - it used to stick
    # forever after one bad SeqLogo column
    errText(
      if (length(unique(nchar(bg.seqs()))) > 1)
        "The length of sequences is different. The input of seqLogo analysis requires the sequences have the same length" else
        NULL
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
    .ggseqlogo_quiet( data = logo() ) + geom_vline(
      xintercept = (ncol(logo())+1)/2, linetype="dashed", color = "orange", linewidth=1.5
      )
  })

  output$plt.fg <- renderPlot({
    req( d <- fg.pfm() )
    .ggseqlogo_quiet( data = d ) + geom_vline(
      xintercept = (ncol(d)+1)/2, linetype="dashed", color = "orange", linewidth=1.5
      )
    })

  output$plt.bg <- renderPlot({
    req( d <- bg.pfm() )
    .ggseqlogo_quiet( data = d ) + geom_vline(
      xintercept = (ncol(d)+1)/2, linetype="dashed", color = "orange", linewidth=1.5
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