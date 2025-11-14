#' Application Constants
#'
#' @description
#' This file contains named constants used throughout the omicsViewer application.
#' Using named constants instead of "magic numbers" improves code readability,
#' maintainability, and makes it easier to update values in a single location.
#'
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
