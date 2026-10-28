#' Application Constants
#'
#' @description
#' This file contains named constants used throughout the omicsViewer application.
#' Using named constants instead of "magic numbers" improves code readability,
#' maintainability, and makes it easier to update values in a single location.
#'
#' @return Documentation-only topic; no return value.
#' @name constants
#' @keywords internal
NULL

# =============================================================================
# External API Limits
# =============================================================================

#' @description Maximum number of genes allowed in STRING database query.
#' This is a limit imposed by the STRING API.
STRING_MAX_GENES <- 300

#' @description Maximum number of network edges to display from STRING database.
#' Limited to prevent performance issues in network visualization.
STRING_MAX_NETWORK_EDGES <- 999

#' @description Required score threshold for STRING database queries (0-1000).
#' Higher values mean more confident interactions.
STRING_REQUIRED_SCORE <- 500

# =============================================================================
# Display and Table Settings
# =============================================================================

#' @description Number of top results to display in Geneshot analysis.
#' Shows the most relevant genes based on publication ranking.
GENESHOT_TOP_RESULTS <- 20

#' @description Default number of rows per page in data tables.
#' Used across multiple table displays in the application.
DEFAULT_TABLE_PAGE_LENGTH <- 10

#' @description Page length for large data tables (feature/sample tables).
DEFAULT_TABLE_PAGE_LENGTH_LARGE <- 25

#' @description Page length for enrichment analysis result tables.
ENRICHMENT_TABLE_PAGE_LENGTH <- 8

#' @description Maximum number of unique labels for ROC/PR curve analysis.
#' Prevents cluttered visualizations with too many categories.
MAX_ROC_PR_LABELS <- 12

#' @description Minimum number of unique labels required for ROC/PR analysis.
MIN_ROC_PR_LABELS <- 2

# =============================================================================
# Statistical Analysis Defaults
# =============================================================================

#' @description Minimum number of samples required for correlation analysis.
#' Ensures statistical reliability of correlation coefficients.
MIN_SAMPLES_CORRELATION <- 12

#' @description Analysis categories whose `Category|Subcategory|Variable` columns
#'   form a volcano view (`<category>|<contrast>|mean.diff` paired with
#'   `<category>|<contrast>|log.fdr` or `log.pvalue`). Drives the default
#'   feature-space axes, the auto-detected volcano quick views and the volcano
#'   corner auto-selection. Override per dataset through
#'   `prepOmicsViewer(volcano.categories = ...)` (stored as the
#'   `volcanoCategories` feature-data attribute).
VOLCANO_DEFAULT_ANALYSES <- c("ttest", "DE")

#' @description Quantile threshold for NA imputation (default 15th percentile).
#' Used when filling missing values in expression data.
IMPUTATION_QUANTILE_THRESHOLD <- 0.15

#' @description Minimum filtering threshold for AutoRIF gene filtering.
#' Genes must appear in at least this fraction of publications.
AUTORIF_MIN_THRESHOLD <- 3

#' @description Divisor for AutoRIF dynamic threshold calculation.
#' Used to calculate: ceiling(nrow(df)/AUTORIF_THRESHOLD_DIVISOR)
AUTORIF_THRESHOLD_DIVISOR <- 200

# =============================================================================
# UI Layout Dimensions
# =============================================================================

#' @description Case-insensitive grep pattern matched against annotation
#'   variable names (the third triselector level) when a Tooltips
#'   attribute cascade starts unset: the first matching column is
#'   preselected so hover labels show gene symbols / gene names out of
#'   the box (see \code{attr4selector_module}).
ATTR4_TOOLTIP_DEFAULT_PATTERN <- "symbol|name"

#' @description Default height for scatter plots in pixels.
DEFAULT_SCATTER_PLOT_HEIGHT <- "400px"

#' @description Height for main scatter plot in meta analysis view.
META_SCATTER_PLOT_HEIGHT <- "666px"

#' @description Height for dose response plots in pixels.
DOSE_RESPONSE_PLOT_HEIGHT <- "550px"

#' @description Height for PTMotif sequence logo plots in pixels.
PTMOTIF_PLOT_HEIGHT <- "300px"

#' @description Height for heatmap key/legend in pixels.
HEATMAP_KEY_HEIGHT <- "45px"

#' @description Default height for empty placeholder plots in pixels.
EMPTY_PLOT_HEIGHT <- "100px"

#' @description Default height for main heatmap display in pixels.
HEATMAP_MAIN_HEIGHT <- "800px"

#' @description Height for scrollable data table views in pixels.
TABLE_SCROLL_HEIGHT <- "450px"

#' @description Top margin for dose response parameter table in pixels.
DOSE_RESPONSE_TABLE_MARGIN_TOP <- "30px"

#' @description Standard margin for UI elements in pixels.
STANDARD_MARGIN <- "10px"

#' @description Margin for figure attribute controls in pixels.
FIGURE_ATTR_MARGIN <- "25px"

#' @description Width for figure attribute control panel in pixels.
FIGURE_ATTR_PANEL_WIDTH <- "788px"

#' @description Width for file selection dropdown in pixels.
FILE_SELECT_WIDTH <- "500px"

#' @description Width for file download button in pixels.
DOWNLOAD_BUTTON_WIDTH <- "115px"

#' @description Default height for annotation row in heatmap (per annotation).
HEATMAP_ANNOT_ROW_HEIGHT_PX <- 20

#' @description Minimum height when no annotations in heatmap.
HEATMAP_MIN_HEIGHT <- "50px"

# =============================================================================
# Bootstrap Grid System
# =============================================================================

#' @description Total columns in Bootstrap grid system.
#' Used for dynamic column width calculations in responsive layouts.
BOOTSTRAP_GRID_COLUMNS <- 12

# =============================================================================
# Margin and Slider Limits
# =============================================================================

#' @description Maximum value for heatmap margin sliders.
HEATMAP_MARGIN_MAX <- 20

#' @description Minimum value for heatmap margin sliders.
HEATMAP_MARGIN_MIN <- 1

#' @description Default margin value for heatmap bottom.
HEATMAP_MARGIN_DEFAULT <- 4

# =============================================================================
# Network Visualization
# =============================================================================

#' @description Opacity for nodes in STRING network visualization (0-1).
STRING_NETWORK_OPACITY <- 0.7

#' @description Font size for labels in STRING network visualization.
STRING_NETWORK_FONT_SIZE <- 12

# =============================================================================
# Agent Tool Boundaries
# =============================================================================

#' @description Literal sentinel strings some LLM providers serialize for
#' omitted optional tool arguments instead of JSON null (observed with glm
#' flash models: "null"; also "{}"/"[]" for object/array-shaped optionals).
#' Every agent-facing boundary treats these exactly like omitted values.
#' One shared definition so the canonical widget store, the semantic apply
#' validators, the generic widget tier, and the figure grammar cannot drift
#' apart.
AGENT_SENTINEL_STRINGS <- c("null", "NULL", "{}", "[]")

#' @description Single-string sentinel test used at tool and validation
#' boundaries.
#' @keywords internal
#' @rdname constants
agent_sentinel_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) &&
    x %in% AGENT_SENTINEL_STRINGS
}

#' @description Full-detail sections of \code{get_omics_viewer_state} that a
#' model can request beyond the compact overview (WP1 progressive
#' disclosure). The overview itself always carries dataset, active and
#' available tabs, selection counts with bounded example IDs, quick-view
#' id+label lists, the widget-store scatter view, and this menu.
#' \code{panel_results} (WP-7, F4) returns the head of the ORA/fGSEA
#' enrichment result tables after a run.
AGENT_STATE_SECTIONS <- c("annotations", "quick_views", "panels", "panel_results", "figure_grammar")

# R-H3 (todo 4.3): ORA result rows fed into the jaccard/hclust
# description clustering - the pairwise pass is quadratic in row count
ORA_CLUSTER_MAX_ROWS <- 200

# R-H6 (todo 4.3): hard row cap for adist()-based clustering. The dense
# n x n distance matrix costs ~n^2 * 8 bytes (~20 GB at the 50k cap plus
# copies) - large but bounded, and clustering beyond it is not informative
# anyway. Above the cap adist refuses with an actionable message.
ADIST_MAX_ROWS <- 50000

# WP-G1: canonical owning-tab label per namespace prefix (longest prefix
# wins; the bare dataspace/resultspace entries are the fallback for the
# tab-select widgets themselves). One static table - the single source of
# the tab field in capability records.
AGENT_TAB_LABELS <- c(
  "dataspace.feature_space" = "Feature",
  "dataspace.sample_space" = "Sample",
  "dataspace.tab_pheno" = "Sample table",
  "dataspace.tab_feature" = "Feature table",
  "dataspace.tab_expr" = "Expression",
  "dataspace.cor_heatmap" = "Cor",
  "dataspace.expr_heatmap" = "Heatmap",
  "dataspace.dyn_heatmap" = "Dynamic heatmap",
  "dataspace" = "Data space",
  "resultspace.feature_general" = "Feature",
  "resultspace.ora" = "ORA",
  "resultspace.fgsea" = "fGSEA",
  "resultspace.stringdb" = "StringDB",
  "resultspace.geneshot" = "Geneshot",
  "resultspace.sample_general" = "Sample",
  "resultspace.ptm" = "SeqLogo",
  "resultspace" = "Analysis space"
)

# WP-G1: panel capability records sourced from the modules' sr-only
# screen-reader prose (verbatim first paragraph per panel; Feature/Sample
# Analysis use their sr-only aria-labels). Surfaced so find_controls can
# return real explanations of what each panel does.
AGENT_PANEL_PROSE <- list(
  list(space = "Data space", tab = "Feature", title = "Feature Scatter Plot",
    prose = "Feature scatter plots visualize relationships between two feature-level variables (e.g., fold change vs. p-value for volcano plots, or PCA coordinates). Each point represents one gene/protein. This view helps identify interesting features, assess data quality, and select features for downstream analysis."),
  list(space = "Data space", tab = "Feature table", title = "Feature Metadata Table",
    prose = "The feature table displays all annotations for genes/proteins in your dataset. Each row represents one feature with columns for identifiers, statistical results, functional annotations, and gene set memberships. This is your primary reference for feature information."),
  list(space = "Data space", tab = "Sample", title = "Sample Scatter Plot",
    prose = "Sample scatter plots visualize relationships between sample-level variables or show sample distributions in reduced dimensions (PCA, t-SNE, UMAP). Each point represents one sample/patient. This view reveals sample clustering, batch effects, and relationships between clinical variables."),
  list(space = "Data space", tab = "Sample table", title = "Sample Metadata Table",
    prose = "The sample table displays all clinical and experimental annotations for samples in your dataset. Each row represents one sample with columns for identifiers, experimental conditions, phenotypes, and other metadata. This is your primary reference for sample information."),
  list(space = "Data space", tab = "Cor", title = "Correlation Heatmap",
    prose = "The correlation heatmap displays pairwise correlations between all samples in your dataset. Hierarchical clustering organizes samples with similar expression profiles together, revealing sample relationships, experimental groupings, and potential batch effects. This is a key quality control and exploratory visualization."),
  list(space = "Data space", tab = "Heatmap", title = "Expression Heatmap",
    prose = "The expression heatmap displays abundance values for features (rows) across samples (columns). Hierarchical clustering groups similar features and samples together, revealing expression patterns, co-regulated genes, and sample groupings. This is essential for exploring large-scale expression patterns."),
  list(space = "Data space", tab = "Dynamic heatmap", title = "Dynamic Heatmap with Selection",
    prose = "The dynamic heatmap allows interactive row selection to subset features based on expression patterns. Click individual rows or drag to select multiple features showing interesting patterns. Selected features are automatically updated in the analysis panel for downstream investigation."),
  list(space = "Data space", tab = "Expression", title = "Expression Matrix Table",
    prose = "The expression matrix table displays raw quantitative values for all features across all samples. Each row represents a feature (gene/protein) and each column represents a sample. This is the numerical data underlying all visualizations and analyses in the application."),
  list(space = "Data space", tab = "GSList", title = "Gene Set Membership Table",
    prose = "The gene set list shows which features belong to which functional categories, pathways, or custom gene sets. This table is used by enrichment analysis tools (ORA, fGSEA) to test for over-representation or enrichment. Each row links a feature to its gene set memberships."),
  list(space = "Analysis space", tab = "Feature", title = "Feature Analysis",
    prose = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"),
  list(space = "Analysis space", tab = "Sample", title = "Sample Analysis",
    prose = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"),
  list(space = "Analysis space", tab = "Geneshot", title = "Geneshot Literature Search",
    prose = "Geneshot is a search engine that finds genes and proteins associated with arbitrary text queries by mining the scientific literature. It ranks genes by how frequently they co-occur with your search terms in PubMed abstracts and other biomedical texts."),
  list(space = "Analysis space", tab = "ORA", title = "Over-Representation Analysis",
    prose = "Over-Representation Analysis identifies gene sets, pathways, and functional categories that are statistically enriched in your selected features. It uses the hypergeometric test to determine if your selected genes/proteins overlap with known biological gene sets more often than would be expected by chance. After a run, the top rows of the results table are readable by the AI assistant through the panel_results state section; column sorting is click-only."),
  list(space = "Analysis space", tab = "fGSEA", title = "Fast Gene Set Enrichment Analysis",
    prose = "fGSEA is a computational method that determines whether predefined gene sets show statistically significant, concordant differences between two biological states. Unlike ORA which only uses a cutoff, fGSEA uses the entire ranked list of genes, making it more sensitive to subtle but coordinated changes in gene expression across a pathway. After a run, the top rows of the results table are readable by the AI assistant through the panel_results state section; column sorting is click-only."),
  list(space = "Analysis space", tab = "Response", title = "Dose-Response Curves",
    prose = "Dose-response analysis models the relationship between drug/compound concentration and biological response using a 4-parameter logistic curve. This analysis fits sigmoid curves to your data and estimates key pharmacological parameters including EC50, IC50, Hill slope, and maximum/minimum response values."),
  list(space = "Analysis space", tab = "StringDB", title = "STRING Protein Interaction Network",
    prose = "STRING (Search Tool for the Retrieval of Interacting Genes/Proteins) is a database of known and predicted protein-protein interactions. This module queries the STRING database to build an interaction network for your selected proteins, showing how they might work together in biological processes."),
  list(space = "Analysis space", tab = "SeqLogo", title = "Post-Translational Modification Motifs",
    prose = "PTM motif analysis identifies enriched sequence patterns surrounding post-translational modification sites (e.g., phosphorylation, acetylation, ubiquitination) in your selected proteins. This analysis uses sequence logo visualization to show the consensus amino acid sequence motif.")
)
