
#' Launch the omicsViewer Shiny Application
#'
#' @description
#' Starts an interactive Shiny application for exploring omics data, including visualization
#' of expression matrices, feature and sample metadata, statistical analyses, and functional
#' enrichment results. The viewer supports both \code{ExpressionSet} and
#' \code{SummarizedExperiment} objects.
#'
#' @param dir Character. Path to directory containing the \code{ExpressionSet} or
#'   \code{SummarizedExperiment} object saved as .RDS file. Provide only the directory path,
#'   not the full file path. The viewer will list all compatible files in this directory.
#' @param additionalTabs List. Optional custom Shiny modules to add as tabs in the "Analyst" panel.
#'   Each element should be a list with: \code{tabName} (character), \code{moduleName} (character),
#'   \code{moduleUi} (UI function), and \code{moduleServer} (server function).
#' @param filePattern Character. Regular expression pattern to filter files displayed in the
#'   directory. Default: \code{".(RDS|DB|SQLITE|SQLITE3)$"} (case-insensitive).
#' @param ESVObj ExpressionSet or SummarizedExperiment. Optional pre-loaded object to view
#'   directly without file selection. If provided, the file dropdown will show "ESVObj.RDS".
#' @param esetLoader Function. Custom loader for reading saved objects. Default: \code{readESVObj}.
#'   Should accept a file path and return an ExpressionSet or SummarizedExperiment object.
#' @param exprsGetter Function. Extracts expression matrix from loaded object.
#'   Default: \code{getExprs}. Should return a numeric matrix.
#' @param pDataGetter Function. Extracts phenotype/sample metadata. Default: \code{getPData}.
#'   Should return a data.frame with rownames matching sample names.
#' @param fDataGetter Function. Extracts feature metadata. Default: \code{getFData}.
#'   Should return a data.frame with rownames matching feature names.
#' @param defaultAxisGetter Function. Determines default axes for plots. Takes two arguments:
#'   \code{x} (the loaded object) and \code{what} (one of "sx", "sy", "fx", "fy" for
#'   sample/feature space x/y axes). Should return column name from metadata.
#' @param appName Character. Application title displayed in the UI. Default: "omicsViewer".
#' @param appVersion Character or package_version. Version number displayed in UI.
#'   Default: current package version.
#' @param log_file Character, FALSE or NULL. Where to write the app console
#'   log (see Details). Default NULL: a per-run file named
#'   \code{app-<timestamp>-<pid>.log} inside \code{tempdir()/omicsviewer-logs}.
#'   Pass an explicit file path to change it, or \code{FALSE} to disable console
#'   logging. By default logging is always on; it can also be disabled through
#'   the environment (see Details).
#'
#' @export
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @examples
#' if (interactive()) {
#'   # Basic usage with example data
#'   omicsViewer(system.file("extdata", package = "omicsViewer"))
#'
#'   # With pre-loaded object
#'   packdir <- system.file("extdata", package = "omicsViewer")
#'   eset <- readRDS(file.path(packdir, "exampleEset.RDS"))
#'   omicsViewer(packdir, ESVObj = eset)
#' }
#'
#' @return NULL (invisibly). Launches the Shiny application. The app runs until stopped by the user.
#'
#' @details
#' The app includes an optional AI assistant. Enable it by installing the
#' suggested packages \code{ellmer (>= 0.5.0)} and \code{shinychat (>= 0.5.0)}.
#' A user can supply a provider, model, and session-only API key through the
#' assistant settings dialog, or an administrator can configure credentials
#' before launching the app.
#'
#' Supported environment variables are:
#' \itemize{
#'   \item \code{OMICSVIEWER_LLM_PROVIDER}: \code{"openai"} or \code{"anthropic"}
#'   \item \code{OMICSVIEWER_LLM_API_KEY}, or the provider's usual
#'     \code{OPENAI_API_KEY} / \code{ANTHROPIC_API_KEY}
#'   \item \code{OMICSVIEWER_LLM_MODEL}, or
#'     \code{OMICSVIEWER_OPENAI_MODEL} / \code{OMICSVIEWER_ANTHROPIC_MODEL}
#'   \item \code{OMICSVIEWER_LLM_BASE_URL}, or a provider-specific
#'     \code{OMICSVIEWER_OPENAI_BASE_URL} / \code{OMICSVIEWER_ANTHROPIC_BASE_URL}
#'   \item \code{OMICSVIEWER_LLM_MAX_REQUESTS}: integer from 1--200; default 40
#'   \item \code{OMICSVIEWER_LLM_LOG}: opt in to local assistant diagnostics;
#'     default false
#'   \item \code{OMICSVIEWER_LLM_LOG_DIR}: writable directory for JSONL logs
#'   \item \code{OMICSVIEWER_LLM_LOG_MAX_BYTES}: log size limit from
#'     1 MB--100 MB; default 10 MB
#' }
#'
#' The assistant is session-local and does not store credentials or conversations
#' in omicsViewer snapshots. Its tools expose only bounded state and annotation
#' summaries and can change only explicitly validated interface views. Figures
#' use an allowlisted declarative ggplot2 grammar: the model never supplies R
#' code, and the chat shows a low-resolution preview with a high-resolution PNG
#' download.
#'
#' Optional diagnostic logging records user prompts, assistant responses, tool
#' requests/results, provider request boundaries, configuration changes, and
#' failures to a local JSONL file. Logs may contain dataset-derived summaries;
#' they never contain API keys. Logging can also be toggled in the assistant
#' drawer for the current session.
#'
#' In addition, the app console is logged by default for debugging: everything
#' the R process prints (\code{print}/\code{cat} output) plus all
#' \code{message}/\code{warning}/\code{error} conditions raised during the
#' session are teed to a local file, while the interactive console display is
#' unchanged. The default location is one file per run,
#' \code{tempdir()/omicsviewer-logs/app-<timestamp>-<pid>.log}, following the
#' same session-scoped convention as the assistant diagnostics; supply
#' \code{log_file} to choose another destination. Two environment variables
#' adjust the default behaviour:
#' \itemize{
#'   \item \code{OMICSVIEWER_LOG_DIR}: writable directory for console logs
#'   \item \code{OMICSVIEWER_LOG}: \code{off}/\code{false}/\code{0} opts out
#'     of console logging entirely
#' }
#'
#' @seealso
#' \code{\link{prepOmicsViewer}} for preparing data objects for visualization.
#' \code{\link{app_module}} for the main application module (developers only).
#' 
omicsViewer <- function(
  dir, additionalTabs = NULL, filePattern = ".(RDS|DB|SQLITE|SQLITE3)$", ESVObj = NULL,
  esetLoader = readESVObj, 
  exprsGetter = getExprs, pDataGetter = getPData, fDataGetter = getFData, 
  defaultAxisGetter = getAx,
  appName = "omicsViewer", appVersion = packageVersion("omicsViewer"),
  log_file = NULL
  ) {
  
  app <- list(
    ui = fluidPage(
      app_ui("app")
    ),
    server = function(input, output, session, aTabs = additionalTabs,
                      f_eset = esetLoader, f_exprs = exprsGetter, f_pd = pDataGetter, f_fd = fDataGetter,
                      axg = defaultAxisGetter) {
      app_module(
        "app", .dir = reactive(dir), additionalTabs = aTabs, filePattern = filePattern,
        esetLoader = f_eset, exprsGetter = f_exprs, pDataGetter = f_pd, fDataGetter = f_fd,
        defaultAxisGetter = axg, appName = appName, appVersion = appVersion, ESVObj = reactive(ESVObj)
        )
    }
  )
  # console log (default on; see ?omicsViewer) ------------------------
  applog <- applog_begin(log_file)
  on.exit(applog_end(applog), add = TRUE)
  if (!is.null(applog))
    message("omicsViewer console log: ", applog$path)

  withCallingHandlers(
    runApp(app),
    message = function(c) applog_write(applog, "MESSAGE", c),
    warning = function(c) applog_write(applog, "WARNING", c),
    error   = function(c) applog_write(applog, "ERROR",   c)
  )
}


