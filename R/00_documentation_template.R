#' Documentation Template for omicsViewer Modules
#'
#' @description
#' This file provides the standard documentation template for all Shiny modules
#' in the omicsViewer package. All module functions should follow this structure
#' to ensure consistency and completeness.
#'
#' @name documentation_template
#' @keywords internal
NULL

# =============================================================================
# UI FUNCTION DOCUMENTATION TEMPLATE
# =============================================================================

#' [Module Name] UI Function
#'
#' @description
#' Creates the user interface for the [module name] module. [Brief description
#' of what this module displays/does - 1-2 sentences].
#'
#' [Optional: More detailed description of the UI layout, key components,
#' and user interaction elements - 2-4 sentences]
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in the corresponding module server function.
#' @param [other_param] [Type]. [Description]. Default: [value].
#'
#' @return
#' A \code{tagList} (or \code{fluidRow}) containing the UI elements:
#' \itemize{
#'   \item [Component 1]: [Description]
#'   \item [Component 2]: [Description]
#'   \item ...
#' }
#'
#' @family [module_family] modules
#' @seealso
#' \code{\link{[module_name]_module}} for the corresponding server logic.
#'
#' @examples
#' if (interactive()) {
#'   ui <- fluidPage(
#'     [module_name]_ui("module_id")
#'   )
#'   server <- function(input, output, session) {
#'     [module_name]_module("module_id", ...)
#'   }
#'   shinyApp(ui, server)
#' }
#'
#' @keywords internal
#' @importFrom shiny NS tagList fluidRow column

# =============================================================================
# MODULE SERVER FUNCTION DOCUMENTATION TEMPLATE
# =============================================================================

#' [Module Name] Server Function
#'
#' @description
#' Server logic for the [module name] module. [Brief description of what
#' this module does - 1-2 sentences].
#'
#' [Optional: More detailed description of the processing logic, data flow,
#' and key functionality - 2-4 sentences]
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in the corresponding UI function.
#'
#' @param reactive_param1 Reactive expression. [Description of what this
#'   reactive returns, expected format/structure]. Example: reactive expression
#'   returning a numeric matrix with features as rows and samples as columns.
#'
#' @param reactive_param2 Reactive expression. [Description]. Optional.
#'   Default: \code{reactive(NULL)}.
#'
#' @param non_reactive_param [Type]. [Description of non-reactive parameter].
#'   Default: [value].
#'
#' @details
#' [Optional section for additional implementation details]
#'
#' ## Input Requirements
#' \itemize{
#'   \item [Required input 1]: [Specification]
#'   \item [Required input 2]: [Specification]
#' }
#'
#' ## Processing Steps
#' \enumerate{
#'   \item [Step 1 description]
#'   \item [Step 2 description]
#'   \item ...
#' }
#'
#' ## Output Format
#' [Description of what the module returns and how it can be used]
#'
#' @return
#' A reactive expression that returns [description of return value and
#' structure]. Returns NULL if [error/edge case conditions].
#'
#' Alternatively for modules with no return value:
#' NULL (invisibly). The module updates reactive values internally and
#' communicates with child modules via reactive expressions.
#'
#' @family [module_family] modules
#' @seealso
#' \code{\link{[module_name]_ui}} for the corresponding UI function.
#' \code{\link{[related_module]}} for related functionality.
#'
#' @examples
#' if (interactive()) {
#'   # Example with sample data
#'   data <- matrix(rnorm(1000), nrow = 100)
#'   rownames(data) <- paste0("Feature", 1:100)
#'   colnames(data) <- paste0("Sample", 1:10)
#'
#'   ui <- fluidPage(
#'     [module_name]_ui("module_id")
#'   )
#'
#'   server <- function(input, output, session) {
#'     result <- [module_name]_module(
#'       "module_id",
#'       reactive_param1 = reactive(data),
#'       reactive_param2 = reactive(NULL)
#'     )
#'   }
#'
#'   shinyApp(ui, server)
#' }
#'
#' @keywords internal
#' @importFrom shiny moduleServer reactive observeEvent req

# =============================================================================
# DOCUMENTATION CHECKLIST
# =============================================================================

# For each module, ensure documentation includes:
#
# UI Functions:
#   [ ] Brief description (1-2 sentences)
#   [ ] @param id with standard text
#   [ ] @param for all other parameters with type and description
#   [ ] @return describing UI structure
#   [ ] @family tag for grouping related modules
#   [ ] @seealso linking to server function
#   [ ] @examples with runnable code (if applicable)
#   [ ] @keywords internal
#   [ ] @importFrom for all used packages
#
# Server Functions:
#   [ ] Brief description (1-2 sentences)
#   [ ] @param id with standard text
#   [ ] @param for all reactive parameters with "Reactive expression" prefix
#   [ ] @param for all non-reactive parameters
#   [ ] Distinguish between required and optional parameters
#   [ ] @details section for complex logic (optional but recommended)
#   [ ] @return describing return value and structure
#   [ ] @family tag matching UI function
#   [ ] @seealso linking to UI function and related modules
#   [ ] @examples with runnable code (if applicable)
#   [ ] @keywords internal
#   [ ] @importFrom for all used packages
#
# Additional Best Practices:
#   [ ] Use consistent terminology (e.g., "feature" not "gene" in general docs)
#   [ ] Specify data structure expectations (data.frame, matrix, vector, etc.)
#   [ ] Document default values explicitly
#   [ ] Note NULL handling behavior
#   [ ] Describe error/edge case behavior
#   [ ] Link to related modules and helper functions
#   [ ] Include practical examples where helpful
