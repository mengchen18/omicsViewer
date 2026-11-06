#' @name panelRegistryHelpers
#' @title Result-space panel registry and shared panel context (WP-0)
#' @description Internal registry + socket for the built-in result-space
#'   panels (PANEL_DEBT_PLAN.md sections 3-4, WP-0). The registry is ONE
#'   ordered list of validated records describing every built-in panel
#'   (module pair, store prefix, tab value, status key, availability
#'   predicate); \code{\link{rs_tab_set}} is the ONE availability function
#'   derived from it. \code{\link{panel_context}} bundles everything a
#'   panel receives (data reactives, selections, child store, sliced
#'   status, retention helper, triset/table helpers, logger, session,
#'   active flag) so later work packages can converge the modules on one
#'   shape without touching the widget-store control plane or the
#'   triselector cascade (AGENTS.md blast-radius rule: this file only
#'   CALLS \code{widget_store_child}, \code{store_register},
#'   \code{store_watch}, \code{store_epoch}, \code{store_bind_triselector},
#'   \code{triselector_module} - wrap-only).
#'
#'   Adoption is deliberately staged: WP-0 ships the registry and helpers
#'   with ZERO call sites inside the app (call-site adoption is WP-1+),
#'   so nothing here may change any user-visible behavior yet.
#' @keywords internal
NULL

############################################################################
### [1] registry records (plan section 3.1)
############################################################################

# Composition pins these two tab values last, exactly as the L1 renderUI
# does today (additionalTabs compose between the core panels and this
# pinned pair). rs_tab_set restores that order (plan section 3.1: registry
# order EQUALS today's literal composition order).
.panel_pinned_tabs <- c("Geneshot", "Sample")

.panel_id_pattern <- "^[a-zA-Z][a-zA-Z0-9_]*$"

# dotted, namespaced (at least one dot), mirroring the widget-store id
# shape; kept local so this file has no load-order dependency on
# auxi_widgetStore.R (R sources R/ in collate order: auxi_* first)
.panel_dotted_pattern <- "^[a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z][a-zA-Z0-9_]*)+$"

#' Describe one built-in result-space panel
#'
#' Plain validated list; the registry (\code{\link{panel_records}}) is an
#' ordered list of these. Field semantics (plan section 3.1):
#' \itemize{
#'   \item \code{module_id} - the namespace id L1 passes to
#'     \code{moduleServer} TODAY, verbatim (renaming changes DOM ids -
#'     forbidden).
#'   \item \code{store_id} - the canonical widget-store prefix; verbatim
#'     for the 7 existing children, NEW \code{"resultspace.dose_response"}
#'     for dose (D2/D18; the module id \code{"rescurve"} is historical).
#'   \item \code{tab_value} - the navbar tab value (drives
#'     \code{input$analyst} and the \code{analyst_tab} store choices).
#'   \item \code{label}/\code{icon} - the tab label / FontAwesome icon
#'     name used by the navbar composition (WP-12 consumes these).
#'   \item \code{tab_title}/\code{sr_title}/\code{aria_label} - the
#'     tab tooltip (the title span's \code{title} attribute), the
#'     screen-reader-only h2 text and its \code{aria-label}; the navbar
#'     composition renders all three (WP-12). NULL (absent) means the
#'     generic tab builder omits that attribute - the seed records all
#'     carry them.
#'   \item \code{order} - registry list order equals today's literal
#'     composition order: Feature, ORA, fGSEA, Response, StringDB,
#'     SeqLogo, Geneshot, Sample.
#'   \item \code{requires(fd, pd)} - availability predicate over the
#'     feature data (and, for future use, phenotype data); must be
#'     defensive: \code{fd} may be NULL (hidden tab) or make reads error
#'     (\code{rs_tab_set} guards evaluation).
#'   \item \code{ui}/\code{server} - the module UI/server pair.
#'   \item \code{stateful} - whether the legacy \code{analyst_*} status
#'     transport applies to the module's own state.
#'   \item \code{status_key} - the \code{analyst_*} snapshot key (note
#'     \code{analyst_gene_shot} differs from the module id - historical,
#'     kept).
#' }
#'
#' @param module_id Single-segment namespace id (letters, digits,
#'   underscore), e.g. \code{"feature_general"}.
#' @param store_id Dotted canonical store prefix, e.g.
#'   \code{"resultspace.ora"}.
#' @param tab_value Navbar tab value, e.g. \code{"ORA"}.
#' @param label Tab label, e.g. \code{"ORA"}.
#' @param icon FontAwesome icon name, e.g. \code{"object-group"}.
#' @param order Scalar integer-ish composition rank (1-based).
#' @param requires Function of \code{(fd, pd)} returning TRUE when the
#'   tab is available.
#' @param ui Module UI function.
#' @param server Module server function.
#' @param stateful Scalar logical; whether legacy status transport
#'   applies.
#' @param status_key NULL or the \code{analyst_*} snapshot key.
#' @param tab_title NULL or the tab tooltip string (title attribute of
#'   the tab label span).
#' @param sr_title NULL or the screen-reader-only h2 text.
#' @param aria_label NULL or the h2's aria-label string.
#' @return A validated record (named list).
#' @keywords internal
#' @rdname panelRegistryHelpers
panel_record <- function(module_id, store_id, tab_value, label, icon,
                         order, requires, ui, server, stateful = TRUE,
                         status_key = NULL, tab_title = NULL,
                         sr_title = NULL, aria_label = NULL) {
  .panel_chr <- function(x, what) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
      stop("panel_record: ", what, " must be a single non-empty character string.",
           call. = FALSE)
    x
  }
  module_id <- .panel_chr(module_id, "module_id")
  if (!grepl(.panel_id_pattern, module_id))
    stop("panel_record: module_id must be a single-segment identifier ",
         "(letters, digits, underscore), not a dotted path: ", module_id,
         call. = FALSE)
  store_id <- .panel_chr(store_id, "store_id")
  if (!grepl(.panel_dotted_pattern, store_id))
    stop("panel_record: store_id must be a dotted namespaced identifier ",
         "such as 'resultspace.ora': ", store_id, call. = FALSE)
  tab_value <- .panel_chr(tab_value, "tab_value")
  label <- .panel_chr(label, "label")
  icon <- .panel_chr(icon, "icon")
  if (!is.numeric(order) || length(order) != 1L || is.na(order) ||
      order != floor(order))
    stop("panel_record: order must be a single finite integer-like number.",
         call. = FALSE)
  if (!is.function(requires))
    stop("panel_record: requires must be a function of (fd, pd).", call. = FALSE)
  if (!is.function(ui))
    stop("panel_record: ui must be a function.", call. = FALSE)
  if (!is.function(server))
    stop("panel_record: server must be a function.", call. = FALSE)
  if (!is.logical(stateful) || length(stateful) != 1L || is.na(stateful))
    stop("panel_record: stateful must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(status_key))
    status_key <- .panel_chr(status_key, "status_key")
  for (.panel_fld in c("tab_title", "sr_title", "aria_label"))
    if (!is.null(get(.panel_fld)))
      assign(.panel_fld, .panel_chr(get(.panel_fld), .panel_fld))
  list(
    module_id = module_id, store_id = store_id, tab_value = tab_value,
    label = label, icon = icon, order = as.integer(order),
    requires = requires, ui = ui, server = server,
    stateful = isTRUE(stateful), status_key = status_key,
    tab_title = tab_title, sr_title = sr_title,
    aria_label = aria_label
  )
}

#' The ordered registry of built-in result-space panels
#'
#' Seed content, WP-0 (plan section 3.1 + appendix cheat sheet; names are
#' verbatim and must never drift). WP-11 record-order data fix: the LIST
#' order is the L1 SERVER composition order - the order L1's registry
#' loop invokes the module servers, which is the binding registration
#' order and therefore the store-snapshot byte order (call order =
#' snapshot/enumeration order; the pre-WP-11 literal call sites composed
#' feature_general, fgsea, ora, stringdb, sample_general, geneshotTab,
#' ptm, rescurve). The TAB composition order (Feature, ORA, fGSEA,
#' Response, StringDB, SeqLogo, [additionalTabs], Geneshot, Sample) is
#' the records' \code{order} field; \code{\link{rs_tab_set}} sorts by
#' it, so its output is byte-identical to the pre-WP-11 registry list
#' order.
#'
#' Dose/STRING records carry their current semantics faithfully: dose
#' gained the store in WP-3 (\code{resultspace.dose_response}; the
#' \code{analyst_dose_response} status entry stays the state_support
#' marker), STRING's status emission was removed in WP-4 (widget state
#' rides \code{resultspace.stringdb.*}; the record keeps
#' \code{status_key = "analyst_stringdb"} for the context's status
#' slice, which is dormant).
#'
#' @format A named list of \code{\link{panel_record}} lists keyed by
#'   module id (a lazily materialized binding - see the note below).
#' @keywords internal
#' @rdname panelRegistryHelpers
# The registry references the module UI/server FUNCTIONS, but collate
# order sources this file (auxi_*) before the module_* definitions - so
# the list is a DELAYED binding: it materializes on first access, when
# the namespace is complete. Data-shaped for consumers (plan section
# 3.1/3.2 sketch uses it as a plain list); tests force it explicitly.
delayedAssign(
  "panel_records",
  value = stats::setNames(
  list(
  panel_record(
    module_id = "feature_general",
    store_id = "resultspace.feature_general",
    tab_value = "Feature", label = "Feature", icon = "chart-column",
    tab_title = "Selected-feature analysis.\nCompares the abundance of the currently selected feature(s) across sample groups:\nbeeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,\nANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary\noutcome is chosen.\nInterpretation: the box spans the interquartile range with the median as center line.\nA small p-value (e.g. < 0.05) indicates a statistically significant group difference -\nalways judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0\nmeans the feature separates the two classes well; 0.5 is no better than chance.",
    sr_title = "Feature Analysis",
    aria_label = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes",
    order = 1L,
    requires = function(fd, pd) TRUE,
    ui = feature_general_ui, server = feature_general_module,
    stateful = TRUE, status_key = "analyst_feature_general"
  ),
  panel_record(
    module_id = "fgsea",
    store_id = "resultspace.fgsea",
    tab_value = "fGSEA", label = "fGSEA", icon = "arrow-down-short-wide",
    tab_title = "Fast gene set enrichment analysis (GSEA on a ranked list).\nUnlike ORA, fGSEA uses ALL features, ranked by a statistic you choose (e.g. fold\nchange or signed p-value), and asks whether a gene set's members cluster at the\ntop or bottom of that ranking - no hard significance cutoff needed.\nInterpretation: NES (normalized enrichment score) > 0 means the set is enriched\namong up-ranked features, NES < 0 depleted/enriched at the other end; padj < 0.05\nis a common significance threshold; the leading edge lists the members driving\nthe enrichment signal.",
    sr_title = "Fast Gene Set Enrichment Analysis",
    aria_label = "Ranked gene set enrichment analysis computing normalized enrichment scores and identifying leading edge genes",
    order = 3L,
    requires = function(fd, pd) !is.null(fd) && !is.null(attr(fd, "GS")),
    ui = enrichment_fgsea_ui, server = enrichment_fgsea_module,
    stateful = TRUE, status_key = "analyst_fgsea"
  ),
  panel_record(
    module_id = "ora",
    store_id = "resultspace.ora",
    tab_value = "ORA", label = "ORA", icon = "object-group",
    tab_title = "Over-representation analysis (hypergeometric test).\nTests whether your currently selected features belong to a gene set (pathway, GO\ncategory, custom set) more often than expected by chance. Needs a well-defined\nselection, e.g. significant hits from the Feature tab.\nInterpretation: the p-value is the chance probability of the observed overlap; FDR\ncorrects for testing many gene sets (FDR < 0.05 is a common significance cutoff);\nthe odds ratio / gene ratio measures enrichment strength. The second table lists\nthe overlapping features per gene set. Results depend on selection size and the\nbackground gene set collection.",
    sr_title = "Over-Representation Analysis",
    aria_label = "Gene set over-representation analysis using hypergeometric test to identify enriched pathways and functional categories",
    order = 2L,
    requires = function(fd, pd) !is.null(fd) && !is.null(attr(fd, "GS")),
    ui = enrichment_analysis_ui, server = enrichment_analysis_module,
    stateful = TRUE, status_key = "analyst_ora"
  ),
  panel_record(
    module_id = "stringdb",
    store_id = "resultspace.stringdb",
    tab_value = "StringDB", label = "StringDB", icon = "share-nodes",
    tab_title = "Protein-protein interaction network (STRING database).\nQueries STRING for known and predicted interactions among the currently selected\nfeatures (max ~300; set the taxonomy code, e.g. 9606 = human) and draws the\nnetwork - nodes are proteins, edges are interactions. Requires internet access.\nInterpretation: edge confidence reflects combined evidence (experiments,\ndatabases, text mining); dense clusters often correspond to functional modules or\ncomplexes; highly connected hub proteins may be key regulators.",
    sr_title = "STRING Protein Interaction Network",
    aria_label = "Protein-protein interaction network from STRING database showing physical and functional associations",
    order = 5L,
    requires = function(fd, pd)
      !is.null(fd) && any(grepl("^StringDB\\|", colnames(fd))),
    ui = string_ui, server = string_module,
    # stateful = FALSE: WP-4 removed the analyst_stringdb status emission;
    # widget state rides the store section, and the context's status slice
    # (analyst_stringdb) is dormant - always NULL
    stateful = FALSE, status_key = "analyst_stringdb"
  ),
  panel_record(
    module_id = "sample_general",
    store_id = "resultspace.sample_general",
    tab_value = "Sample", label = "Sample", icon = "arrows-down-to-people",
    tab_title = "Selected-sample analysis.\nAnalyses the currently selected samples against the rest of the cohort:\nbeeswarm plots for continuous variables, contingency tables of categorical\nannotations (chi-square / Fisher's exact test), and Kaplan-Meier survival\ncurves with the log-rank test when time-to-event data are available.\nInterpretation: a small contingency p-value means the selected samples differ in\nthat annotation more than expected by chance. Kaplan-Meier curves show the\nevent-free proportion over time per group; the log-rank p-value tests whether the\ncurves differ. Check group sizes (n) - small groups make tests unreliable.",
    sr_title = "Sample Analysis",
    aria_label = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves",
    order = 8L,
    requires = function(fd, pd) TRUE,
    ui = sample_general_ui, server = sample_general_module,
    stateful = TRUE, status_key = "analyst_sample_general"
  ),
  panel_record(
    module_id = "geneshotTab",
    store_id = "resultspace.geneshot",
    tab_value = "Geneshot", label = "Geneshot", icon = "book-open",
    tab_title = "Literature-based gene discovery.\nSearches PubMed-associated literature (AutoRIF) for genes co-mentioned with your\nsearch term(s); separate multiple terms with ';'. Useful for hypothesis generation.\nInterpretation: the scatter plot shows publication count (x) versus the share of a\ngene's publications mentioning your term (y). Genes in the upper right are both\nwell studied and specifically associated with your topic. A '+' marks genes that\nare present in your current feature selection. Ranking reflects literature bias,\nnot statistical evidence from your data.",
    sr_title = "Geneshot Literature Search",
    aria_label = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms",
    order = 7L,
    requires = function(fd, pd) TRUE,
    ui = geneshot_ui, server = geneshot_module,
    stateful = TRUE, status_key = "analyst_gene_shot"
  ),
  panel_record(
    module_id = "ptm",
    store_id = "resultspace.ptm",
    tab_value = "SeqLogo", label = "SeqLogo", icon = "font",
    tab_title = "PTM motif enrichment (sequence logo).\nTakes the amino-acid sequence windows around the selected modification sites and\ntests which residues are over-represented at each position, shown as a sequence\nlogo around the modification site (position 0).\nInterpretation: taller letters mean the amino acid (or property class) is more\nconserved/enriched at that position; significant enrichment versus all sites in\nthe dataset hints at the kinase or recognition machinery involved.",
    sr_title = "Post-Translational Modification Motifs",
    aria_label = "PTM motif enrichment analysis identifying overrepresented amino acid sequence patterns around modification sites",
    order = 6L,
    requires = function(fd, pd)
      !is.null(fd) && any(grepl("^SeqLogo\\|", colnames(fd))),
    ui = ptmotif_ui, server = ptmotif_module,
    stateful = TRUE, status_key = "analyst_ptm"
  ),
  panel_record(
    module_id = "rescurve",
    store_id = "resultspace.dose_response",  # NEW in WP-3 (D2/D18)
    tab_value = "Response", label = "Response", icon = "capsules",
    tab_title = "Dose-response curves (4-parameter logistic fit).\nFits abundance vs concentration for the selected features and estimates EC50/IC50\n(concentration giving half-maximal effect), Hill slope (curve steepness) and the\ntop/bottom response plateaus.\nInterpretation: a lower EC50/IC50 means higher potency; the Hill slope shows how\nabruptly the response switches; the parameter table gives fitted values with\nconfidence intervals. Points are observed data, the line the model fit.",
    sr_title = "Dose-Response Curves",
    aria_label = "Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model",
    order = 4L,
    requires = function(fd, pd)
      !is.null(fd) && any(grepl("^ResponseCurve\\|", colnames(fd))),
    ui = dose_response_ui, server = dose_response_module,
    # WP-3 joined the store; the legacy status channel carries only the
    # state_support marker (the plot is a computed payload)
    stateful = FALSE, status_key = "analyst_dose_response"
  )
  ),
  c("feature_general", "fgsea", "ora", "stringdb", "sample_general",
    "geneshotTab", "ptm", "rescurve")
  ),
  eval.env = environment(),
  assign.env = environment()
)

############################################################################
### [2] the ONE availability function (plan section 3.2)
############################################################################

#' Compute the available analysis (result-space) tab values
#'
#' Registry-backed replacement for the triplicated availability logic
#' (L1 \code{.rs_tab_choices}, the L1 renderUI conditions, L0
#' \code{agent_analysis_tabs}). Since the WP-11 record-order data fix the
#' registry list is in SERVER-composition order, so records are evaluated
#' in \code{order}-field (tab-composition) order here - the returned
#' vector is byte-identical to the pre-WP-11 output. Predicates are
#' NULL-safe and their
#' evaluation is guarded: a failed read degrades to
#' Feature/Geneshot/Sample exactly as today's tryCatch-wrapped
#' \code{.rs_tab_choices}. Additional tab names compose between the core
#' panels and the pinned Geneshot/Sample pair (today's literal order), and
#' the result is uniqued.
#'
#' Call sites adopt this in WP-1 (L1 \code{.rs_tab_choices} and L0
#' \code{agent_analysis_tabs}); the renderUI conditions convert in WP-12.
#'
#' @param fd Feature data (may be NULL when the tab/data is unavailable).
#' @param additional Character vector of additional tab names
#'   (\code{additionalTabs$tabName}); \code{character(0)} by default.
#' @return Unique character vector of available tab values in composition
#'   order.
#' @keywords internal
#' @rdname panelRegistryHelpers
rs_tab_set <- function(fd, additional = character(0)) {
  if (is.null(additional))
    additional <- character(0)
  additional <- as.character(additional)
  # WP-11: evaluate in TAB-composition order (the order field), not
  # registry list order (the server-composition order since the WP-11
  # record-order data fix) - the output vector stays byte-identical
  recs <- panel_records[order(vapply(panel_records,
                                     function(r) r$order, integer(1)))]
  ok <- vapply(
    recs,
    function(r)
      isTRUE(tryCatch(
        r$requires(fd, NULL),
        shiny.silent.error = function(e) FALSE,
        error = function(e) FALSE)),
    logical(1))
  vals <- vapply(recs[ok], function(r) r$tab_value, character(1))
  # pinned pair composes last; their records (always available) are held
  # back from the leading segment so additional tabs sit between the core
  # panels and Geneshot/Sample, byte-identical to today's vector
  unique(c(vals[!vals %in% .panel_pinned_tabs], additional,
           .panel_pinned_tabs))
}

############################################################################
### [3] shared panel context + helpers (plan section 3.3/3.4)
############################################################################

#' Observer-retention factory (the \code{.xx_keep} idiom, owned once)
#'
#' Shiny holds observer dependencies weakly: store-glue observers that are
#' not referenced anywhere are garbage collected between flushes and later
#' store-to-UI pushes silently never fire (AGENTS.md observer-GC lore).
#' Every module hand-rolled this closure; this is the one shared
#' definition. PER-INSTANCE by construction (each call returns a fresh
#' closure over its own list) - never share one factory result across
#' sessions or modules (D3). Contexts hold one instance; modules running
#' without a context (\code{context = NULL}) fall back to a local one.
#'
#' @return A closure \code{function(obs)} that appends and returns
#'   \code{obs}, keeping it referenced for the lifetime of the closure's
#'   owning frame.
#' @keywords internal
#' @rdname panelRegistryHelpers
panel_keep <- function() {
  obs <- list()
  function(o) {
    obs[[length(obs) + 1L]] <<- o
    invisible(o)
  }
}

#' Bind a triselector cascade to the widget store (the boilerplate, once)
#'
#' Absorbs ONLY the per-panel cascade boilerplate (D4): the guarded triset
#' source, the ts1/ts2/ts3 choices providers, the
#' xax_analysis/xax_subset/xax_variable binding registration, the
#' \code{\link{triselector_module}} wiring through
#' \code{\link{store_watch}}/\code{\link{store_epoch}}, and
#' \code{\link{store_bind_triselector}}. Label/help strings, sync/seed/
#' store-to-UI push observers, extra bindings and actions stay in the
#' module files - they are panel semantics, not boilerplate.
#'
#' The registered bindings are exactly today's shape: xax_analysis
#' (\code{"select"}), xax_subset (\code{"select"}, depends on
#' xax_analysis), xax_variable (\code{"select_cascaded"}, depends on
#' both), with the caller's labels/help and choices providers that mirror
#' the triset source (\code{unique} per cascade level, empty when the
#' source is cold - the exact semantics of the current \code{.xx_ts1/
#' .xx_ts2/.xx_ts3} helpers).
#'
#' @param store Child view of the canonical widget store for the panel
#'   (e.g. \code{resultspace.ora}).
#' @param keep Retention closure from \code{\link{panel_keep}}; keeps the
#'   UI-to-store sync observer referenced (observer-GC rule).
#' @param tri_id Module id of the triselector (e.g. \code{"tris_ora"}).
#' @param tri_label Group label shown by the triselector (e.g.
#'   \code{"Collapse features on"}).
#' @param triset_source Panel-owned function returning the nx3 character
#'   cascade matrix (e.g. \code{function() trisetter(meta = fd(),
#'   combine = "none")}). Evaluated defensively: NULL/errors/non-matrix
#'   results degrade to empty choices.
#' @param bindings List with elements \code{analysis}, \code{subset},
#'   \code{variable}, each a list with \code{label} and \code{help}
#'   strings for that binding (panel-owned strings; they do not move into
#'   this file).
#' @param allow_unset Passed to \code{\link{triselector_module}}; TRUE for
#'   panels whose variable may be unset (\code{"--select--"} placeholder).
#' @return The triselector module's returned reactive (settled triples,
#'   or NULL while the first selection has not settled).
#' @keywords internal
#' @rdname panelRegistryHelpers
panel_triset_bind <- function(store, keep, tri_id, tri_label,
                              triset_source, bindings,
                              allow_unset = FALSE) {
  if (is.null(store) || !is.environment(store))
    stop("panel_triset_bind: store must be a widget store/child view ",
         "(e.g. from widget_store_child).", call. = FALSE)
  if (!is.function(keep))
    stop("panel_triset_bind: keep must be a panel_keep() closure.", call. = FALSE)
  if (!is.character(tri_id) || length(tri_id) != 1L || !nzchar(tri_id))
    stop("panel_triset_bind: tri_id must be a non-empty module id.",
         call. = FALSE)
  if (!is.character(tri_label) || length(tri_label) != 1L ||
      !nzchar(tri_label))
    stop("panel_triset_bind: tri_label must be a non-empty group label.",
         call. = FALSE)
  if (!is.function(triset_source))
    stop("panel_triset_bind: triset_source must be a function returning ",
         "an nx3 character matrix.", call. = FALSE)
  if (!is.list(bindings) ||
      !setequal(names(bindings), c("analysis", "subset", "variable")))
    stop("panel_triset_bind: bindings must be a list with elements ",
         "analysis, subset, variable.", call. = FALSE)
  for (nm in names(bindings)) {
    b <- bindings[[nm]]
    if (!is.list(b) || !all(c("label", "help") %in% names(b)) ||
        !is.character(b$label) || length(b$label) != 1L || !nzchar(b$label) ||
        !is.character(b$help) || length(b$help) != 1L || !nzchar(b$help))
      stop("panel_triset_bind: bindings$", nm,
           " must be a list with non-empty label and help strings.",
           call. = FALSE)
  }
  if (!is.logical(allow_unset) || length(allow_unset) != 1L || is.na(allow_unset))
    stop("panel_triset_bind: allow_unset must be TRUE or FALSE.", call. = FALSE)

  # guarded triset source (the .xx_ts shape shared by every panel today)
  .ts <- function() {
    ts <- tryCatch(triset_source(),
                   shiny.silent.error = function(e) NULL,
                   error = function(e) NULL)
    if (is.null(ts) || !is.matrix(ts) || nrow(ts) == 0L) NULL else ts
  }
  .ts1 <- function() {
    ts <- .ts()
    if (is.null(ts)) character(0) else unique(ts[, 1])
  }
  .ts2 <- function(a) {
    ts <- .ts()
    if (is.null(ts) || is.null(a) || !nzchar(a)) character(0)
    else unique(ts[ts[, 1] %in% a, 2])
  }
  .ts3 <- function(a, b) {
    ts <- .ts()
    if (is.null(ts)) character(0)
    else {
      i <- rep(TRUE, nrow(ts))
      if (!is.null(a) && nzchar(a)) i <- i & ts[, 1] %in% a
      if (!is.null(b) && nzchar(b)) i <- i & ts[, 2] %in% b
      unique(ts[i, 3])
    }
  }
  k1 <- .widget_store_key(store, "xax_analysis")
  k2 <- .widget_store_key(store, "xax_subset")
  store_register(
    store,
    widget_binding(
      "xax_analysis", "select",
      label = bindings$analysis$label, help = bindings$analysis$help,
      choices_provider = function(v) .ts1()),
    widget_binding(
      "xax_subset", "select",
      label = bindings$subset$label, help = bindings$subset$help,
      depends_on = "xax_analysis",
      choices_provider = function(v) .ts2(v[[k1]])),
    widget_binding(
      "xax_variable", "select_cascaded",
      label = bindings$variable$label, help = bindings$variable$help,
      depends_on = c("xax_analysis", "xax_subset"),
      choices_provider = function(v) .ts3(v[[k1]], v[[k2]]))
  )
  v1 <- triselector_module(
    tri_id,
    reactive_x = reactive(.ts()), label = tri_label,
    reactive_selector1 = store_watch(store, "xax_analysis"),
    reactive_selector2 = store_watch(store, "xax_subset"),
    reactive_selector3 = store_watch(store, "xax_variable"),
    reactive_axis_request = store_epoch(store),
    allow_unset = allow_unset)
  store_bind_triselector(
    store,
    keys = c(analysis = "xax_analysis", subset = "xax_subset",
             variable = "xax_variable"),
    sel = v1, keep = keep)
  v1
}

# D14: dormant logger. The assistant owns real logging
# (agent_logger_new); L0 does not construct a logger in this pass. The
# context field exists so later threading is mechanical, and its default
# accepts any event call and does nothing.
.panel_logger_noop <- function(...) invisible(NULL)

#' Construct the per-panel context bundle
#'
#' The bundle every result-space panel receives (plan section 3.3). Built
#' once per panel (registry loop in WP-11) and passed as the optional
#' \code{context} formal each module server gains; modules coalesce
#' (\code{store <- context$store \%||\% store}) so all legacy formals keep
#' working. Fields:
#' \describe{
#'   \item{id, tab_value}{Panel identity (module id; navbar tab value).}
#'   \item{data}{Named list of data reactives; subset of
#'     \code{expr/pdata/fdata/object}.}
#'   \item{selection}{Named list of selection reactives; subset of
#'     \code{feature/sample}.}
#'   \item{store}{The panel's child store view.}
#'   \item{status}{Reactive carrying the panel's sliced legacy status
#'     (\code{status()$analyst_<key>}); \code{reactive(NULL)} by default.
#'     Status OUT stays the module's return value this pass (the
#'     convention is documented, not re-plumbed).}
#'   \item{keep}{Retention closure (a fresh \code{\link{panel_keep()}}
#'     instance per context - D3, never package-level).}
#'   \item{triset}{Partial application of \code{\link{panel_triset_bind}}
#'     over this context's store/keep.}
#'   \item{table}{Passthrough to \code{\link{dataTableDownload_module}}
#'     with this context's store forwarded (the sanctioned table pattern;
#'     dormant until adopted).}
#'   \item{logger}{Dormant no-op by default (D14); accepts any event call.}
#'   \item{session}{The caller's Shiny session (dormant field); panels
#'     may ignore it.}
#'   \item{active}{When \code{active_tab} is given (which requires
#'     \code{tab_value}): \code{reactive(active_tab() == tab_value)} -
#'     the string module's \code{active} wiring re-expressed;
#'     \code{reactive(NULL)} otherwise.}
#' }
#'
#' Deliberately NOT in the context (plan section 3.3): agent tool
#' registration, agent excerpt providers, prose/labels, snapshot policy,
#' additionalTabs handling.
#'
#' @param id Module id (namespace id L1 passes to moduleServer).
#' @param tab_value Navbar tab value; needed for \code{active}.
#' @param store The panel's child store view, or NULL (legacy/no-store
#'   use; \code{triset}/\code{table} then refuse with an actionable
#'   error).
#' @param data Named list of data reactives (\code{expr}, \code{pdata},
#'   \code{fdata}, \code{object} - any subset).
#' @param selection Named list of selection reactives (\code{feature},
#'   \code{sample} - any subset).
#' @param status Reactive carrying the panel's sliced status; defaults to
#'   \code{reactive(NULL)}.
#' @param active_tab Optional reactive yielding the session-wide active
#'   analysis tab name (L1's \code{status()$analyst_active_tab}); drives
#'   \code{active}. Requires \code{tab_value} (error otherwise - the
#'   comparison needs it).
#' @param keep Retention closure; defaults to a FRESH
#'   \code{\link{panel_keep}()} instance (per-context, D3).
#' @param logger Optional logger closure; defaults to the dormant no-op
#'   (D14).
#' @param session Optional Shiny session (dormant field).
#' @return A named list (the context bundle).
#' @keywords internal
#' @rdname panelRegistryHelpers
panel_context <- function(id, tab_value = NULL, store = NULL,
                          data = list(), selection = list(),
                          status = shiny::reactive(NULL),
                          active_tab = NULL,
                          keep = panel_keep(),
                          logger = NULL,
                          session = NULL) {
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id))
    stop("panel_context: id must be a single non-empty character string.",
         call. = FALSE)
  if (!is.null(tab_value) &&
      (!is.character(tab_value) || length(tab_value) != 1L ||
       is.na(tab_value) || !nzchar(tab_value)))
    stop("panel_context: tab_value must be NULL or a non-empty string.",
         call. = FALSE)
  if (!is.null(store) && !is.environment(store))
    stop("panel_context: store must be a widget store/child view.",
         call. = FALSE)
  .panel_reactives <- function(x, what, allowed) {
    if (!is.list(x))
      stop("panel_context: ", what, " must be a named list of reactives.",
           call. = FALSE)
    bad <- setdiff(names(x), allowed)
    if (length(bad))
      stop("panel_context: ", what, " has unknown element(s): ",
           paste(bad, collapse = ", "), "; allowed: ",
           paste(allowed, collapse = ", "), ".", call. = FALSE)
    for (nm in names(x))
      if (!is.function(x[[nm]]))
        stop("panel_context: ", what, "$", nm, " must be a reactive ",
             "(function).", call. = FALSE)
    invisible(NULL)
  }
  .panel_reactives(data, "data", c("expr", "pdata", "fdata", "object"))
  .panel_reactives(selection, "selection", c("feature", "sample"))
  if (!is.function(status))
    stop("panel_context: status must be a reactive (function).", call. = FALSE)
  if (!is.null(active_tab)) {
    if (!is.function(active_tab))
      stop("panel_context: active_tab must be NULL or a reactive (function).",
           call. = FALSE)
    if (is.null(tab_value))
      stop("panel_context: active_tab requires tab_value (the panel's ",
           "navbar tab value that 'active' compares against); without it ",
           "active() would silently yield logical(0).", call. = FALSE)
  }
  if (!is.function(keep))
    stop("panel_context: keep must be a panel_keep() closure.", call. = FALSE)
  if (is.null(logger))
    logger <- .panel_logger_noop
  else if (!is.function(logger))
    stop("panel_context: logger must be NULL or a closure.", call. = FALSE)
  ctx <- list(
    id = id,
    tab_value = tab_value,
    data = data,
    selection = selection,
    store = store,
    status = status,
    keep = keep,
    triset = function(tri_id, tri_label, triset_source, bindings,
                      allow_unset = FALSE) {
      if (is.null(store))
        stop("panel_context: panel '", id, "' has no store; the triset ",
             "cascade helper requires a store-backed panel.", call. = FALSE)
      panel_triset_bind(store, keep, tri_id, tri_label, triset_source,
                        bindings, allow_unset)
    },
    table = function(id, ..., store_key = NULL) {
      if (is.null(store))
        stop("panel_context: panel '", id, "' has no store; the table ",
             "helper requires a store-backed panel.", call. = FALSE)
      dataTableDownload_module(id, ..., store = store, store_key = store_key)
    },
    logger = logger,
    session = session,
    active = if (is.null(active_tab)) shiny::reactive(NULL)
             else shiny::reactive(active_tab() == tab_value)
  )
  ctx
}

############################################################################
### [4] additionalTabs entry validation (P8, WP-2)
############################################################################

# describes an offending value for the actionable-error convention
# ("a numeric of length 1", "a function", "a character of length 2", ...)
.panel_val_desc <- function(v) {
  if (is.function(v))
    "a function"
  else
    paste0("a ", paste(class(v), collapse = "/"), " of length ", length(v))
}

#' Validate the \code{additionalTabs} argument (fail fast; P8, WP-2)
#'
#' Entry guard for the legacy 5-field additionalTabs convention, which
#' STAYS the convention this pass (plan sections 4 WP-2, 7; D8/D15: no
#' deprecation shims, no registry absorption). Validates exactly the
#' fields the consumers read - the two L1 loops (server instantiation,
#' tab-panel composition, status builder) and the L0 availability read
#' use \code{tabName} (navbar title; \code{vapply} \code{character(1)}
#' in L0), \code{moduleName} (moduleServer id and status key),
#' \code{moduleUi}/\code{moduleServer} (functions) and the optional
#' \code{stateful} (only ever fed to \code{isTRUE()}) - nothing
#' stricter: extra fields, unnamed top-level lists, NA/absent/NULL
#' \code{stateful} and duplicate names all pass because the loops
#' tolerate them.
#'
#' Called at BOTH entries - \code{\link{omicsViewer}} (first statement,
#' before the app object is built) and \code{\link{app_module}} (before
#' \code{moduleServer}) - so malformed input fails fast with a message
#' naming the offending entry index, field and expected type, instead
#' of the deep-dive explosions an unnamed list element used to produce
#' inside reactive contexts (\code{do.call}/\code{vapply}). NULL (the
#' default) is the fast path (one \code{is.null} check); the empty
#' list passes (the loops' \code{length(x) > 0} guards).
#'
#' @param x The \code{additionalTabs} argument: NULL or a list of tab
#'   records.
#' @return \code{invisible(NULL)}; throws an actionable error
#'   (\code{call. = FALSE}) naming the first malformed entry.
#' @keywords internal
#' @rdname panelRegistryHelpers
panel_validate_additional_tabs <- function(x) {
  if (is.null(x))
    return(invisible(NULL))
  if (!is.list(x))
    stop("additionalTabs must be NULL or a list of tab records; got ",
         .panel_val_desc(x), ".", call. = FALSE)
  for (i in seq_along(x)) {
    lo <- x[[i]]
    if (!is.list(lo))
      stop("additionalTabs[[", i, "]] must be a list (a tab record); got ",
           .panel_val_desc(lo), ".", call. = FALSE)
    for (fld in c("tabName", "moduleName")) {
      v <- lo[[fld]]
      if (is.null(v))
        stop("additionalTabs[[", i, "]] is missing required field '", fld,
             "' (a non-empty character string of length 1).", call. = FALSE)
      if (!is.character(v) || length(v) != 1L || is.na(v) || !nzchar(v))
        stop("additionalTabs[[", i, "]]$", fld,
             " must be a non-empty character string of length 1; got ",
             .panel_val_desc(v), ".", call. = FALSE)
    }
    for (fld in c("moduleUi", "moduleServer")) {
      v <- lo[[fld]]
      if (is.null(v))
        stop("additionalTabs[[", i, "]] is missing required field '", fld,
             "' (a function).", call. = FALSE)
      if (!is.function(v))
        stop("additionalTabs[[", i, "]]$", fld, " must be a function; got ",
             .panel_val_desc(v), ".", call. = FALSE)
    }
    st <- lo$stateful
    if (!is.null(st) && !(is.logical(st) && length(st) == 1L))
      stop("additionalTabs[[", i, "]]$stateful must be TRUE or FALSE ",
           "(a logical of length 1), or absent/NULL; got ",
           .panel_val_desc(st), ".", call. = FALSE)
  }
  invisible(NULL)
}
