#' Batch Comparison UI Function
#'
#' @description
#' Creates the user interface for the batch comparison module that automatically
#' compares either all phenotype variables or all features between selected
#' and unselected sample groups.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{batch_comparison_module}}.
#'
#' @return
#' A \code{tagList} containing:
#' \itemize{
#'   \item Toggle buttons to select comparison mode (phenotype vs features)
#'   \item Results table with p-values and FDR correction
#'   \item Download functionality for results
#'   \item Selected row information panel
#' }
#'
#' @family analysis modules
#' @seealso
#' \code{\link{batch_comparison_module}} for the corresponding server logic.
#'
#' @keywords internal
#' @importFrom shinyWidgets radioGroupButtons
#' @importFrom DT dataTableOutput
#'
batch_comparison_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Batch Comparison Analysis"),
      tags$p("Batch comparison analysis automatically tests all phenotype variables or all molecular features for differences between your selected and unselected sample groups. This comprehensive approach applies appropriate statistical tests based on variable types (numeric, categorical, or survival data) and corrects for multiple testing using FDR. It provides a systematic way to identify which clinical variables or molecular features significantly differ between your sample groups."),
      tags$h4("When to use batch comparison"),
      tags$p("Use this analysis when you have selected a group of samples (e.g., based on clustering, clinical criteria, or molecular signature) and want to systematically identify: (1) which clinical/phenotype variables characterize your selected group, or (2) which genes/proteins are differentially expressed in your selected group. This is essential for biomarker discovery, patient stratification validation, and comprehensive group characterization."),
      tags$h4("How to interpret results"),
      tags$p("The results table shows one row per variable/feature tested. P-values indicate statistical significance of group differences. FDR-corrected p-values (q-values) control for false discoveries when testing many variables - use FDR < 0.05 as a stringent threshold. For phenotype comparisons, summary statistics show distributions in each group. For feature comparisons, mean differences indicate effect sizes. Sort by FDR to identify the most significant differences.")
    ),
    fluidRow(
      column(
        6,
        checkboxInput(ns("show_phenotype"), "Compare All Phenotype Variables", value = TRUE)
      ),
      column(
        6,
        checkboxInput(ns("show_features"), "Compare All Features/Expression", value = FALSE)
      ),
      tabsetPanel(
        tabPanel(
          title = "Phenotype variable",
          uiOutput(ns("phenotype_results"))
        ),
        tabPanel(
          title = "Feature differential expression",
          uiOutput(ns("features_results"))
        )
      )
    )
  )
}

#' Batch Comparison Server Function
#'
#' @description
#' Server logic for the batch comparison module. Automatically tests all
#' phenotype variables or all features for differences between selected
#' and unselected samples with appropriate statistical tests and FDR correction.
#'
#' @param id Character. Namespace ID for the Shiny module.
#' @param reactive_expr Reactive expression. Returns numeric matrix with
#'   features as rows and samples as columns.
#' @param reactive_phenoData Reactive expression. Returns data.frame of
#'   sample metadata.
#' @param reactive_featureData Reactive expression. Returns data.frame of
#'   feature metadata.
#' @param reactive_i_samples Reactive expression. Returns integer vector of
#'   selected sample indices.
#'
#' @return
#' Reactive value containing information about the selected row from the
#' results table.
#'
#' @family analysis modules
#' @seealso
#' \code{\link{batch_comparison_ui}} for the corresponding UI function.
#'
#' @keywords internal
#' @importFrom stats wilcox.test fisher.test p.adjust
#' @importFrom survival survdiff survfit coxph Surv
#'
batch_comparison_module <- function(
    id,
    reactive_expr,
    reactive_phenoData,
    reactive_featureData,
    reactive_i_samples
) {

  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # Helper function: Detect variable type
    detectVarType <- function(col_name, values, pd) {
      # Priority 1: Check if survival variable
      if (grepl("^Surv\\|", col_name)) {
        # Validate it's actually a Surv object
        if (inherits(values, "Surv")) {
          return("survival")
        } else {
          return("exclude")
        }
      }

      # Priority 2: Check data type
      if (is.numeric(values)) {
        # Need at least 2 non-NA values
        if (sum(!is.na(values)) >= 2) {
          return("numeric")
        } else {
          return("exclude")
        }
      } else if (is.character(values) || is.factor(values)) {
        # Need 2-12 levels with sufficient observations
        vals <- values[!is.na(values)]
        n_levels <- length(unique(vals))
        if (n_levels >= 2 && n_levels <= 12) {
          # Check that we have at least 2 observations
          if (length(vals) >= 2) {
            return("categorical")
          }
        }
        return("exclude")
      }

      return("exclude")
    }

    # Helper function: Compare numeric variable
    compareNumericVar <- function(values_selected, values_unselected) {
      # Remove NAs
      v1 <- values_selected[!is.na(values_selected)]
      v2 <- values_unselected[!is.na(values_unselected)]

      # Need at least 2 values in each group
      if (length(v1) < 2 || length(v2) < 2) {
        return(list(
          test_method = "Wilcoxon rank-sum",
          mean_diff = NA,
          p_value = NA
        ))
      }

      # Wilcoxon rank-sum test (non-parametric)
      test_result <- tryCatch(
        wilcox.test(v1, v2, exact = FALSE),
        error = function(e) NULL
      )

      if (is.null(test_result)) {
        p_val <- NA
      } else {
        p_val <- test_result$p.value
      }

      # Mean difference (effect size)
      mean_diff <- mean(v1) - mean(v2)

      list(
        test_method = "Wilcoxon rank-sum",
        mean_diff = mean_diff,
        p_value = p_val
      )
    }

    # Helper function: Compare categorical variable
    compareCategoricalVar <- function(values_selected, values_unselected) {
      # Create grouping vector
      groups <- c(rep("selected", length(values_selected)),
                  rep("unselected", length(values_unselected)))
      all_values <- c(values_selected, values_unselected)

      # Remove NAs
      keep_idx <- !is.na(all_values)
      groups <- groups[keep_idx]
      all_values <- all_values[keep_idx]

      # Need at least some observations
      if (length(all_values) < 2) {
        return(list(
          test_method = "Fisher's exact",
          odds_ratio = NA,
          p_value = NA
        ))
      }

      # Create contingency table
      tab <- table(all_values, groups)

      # Fisher's exact test
      test_result <- tryCatch({
        fisher.test(tab)
      }, error = function(e) {
        # If table too large, use simulation
        tryCatch(
          fisher.test(tab, simulate.p.value = TRUE, B = 1e5),
          error = function(e2) NULL
        )
      })

      if (is.null(test_result)) {
        p_val <- NA
        odds_ratio <- NA
      } else {
        p_val <- test_result$p.value
        # Extract odds ratio if available (2x2 table)
        if (!is.null(test_result$estimate)) {
          odds_ratio <- as.numeric(test_result$estimate)
        } else {
          odds_ratio <- NA
        }
      }

      list(
        test_method = "Fisher's exact",
        odds_ratio = odds_ratio,
        p_value = p_val
      )
    }

    # Helper function: Compare survival variable
    compareSurvivalVar <- function(surv_selected, surv_unselected) {
      # Combine survival objects
      surv_combined <- c(surv_selected, surv_unselected)
      groups <- factor(c(rep("selected", length(surv_selected)),
                         rep("unselected", length(surv_unselected))),
                       levels = c("unselected", "selected"))

      # Log-rank test
      test_result <- tryCatch(
        survival::survdiff(surv_combined ~ groups),
        error = function(e) NULL
      )

      if (is.null(test_result)) {
        p_val <- NA
      } else {
        p_val <- 1 - pchisq(test_result$chisq, df = 1)
      }

      # Hazard ratio using Cox proportional hazards
      hazard_ratio <- tryCatch({
        cox_model <- survival::coxph(surv_combined ~ groups)
        exp(coef(cox_model)[1])
      }, error = function(e) NA)

      list(
        test_method = "Log-rank",
        hazard_ratio = hazard_ratio,
        p_value = p_val
      )
    }

    # Master function: Compare all phenotype variables
    comparePhenotypeVars <- function(pd, selected_idx, unselected_idx) {
      results_list <- list()

      for (col_name in colnames(pd)) {
        # Get values for both groups
        values_all <- pd[[col_name]]
        values_selected <- values_all[selected_idx]
        values_unselected <- values_all[unselected_idx]

        # Detect variable type
        var_type <- detectVarType(col_name, values_all, pd)

        if (var_type == "exclude") {
          next  # Skip this variable
        }

        # Perform appropriate test
        if (var_type == "numeric") {
          result <- compareNumericVar(values_selected, values_unselected)
        } else if (var_type == "categorical") {
          result <- compareCategoricalVar(values_selected, values_unselected)
        } else if (var_type == "survival") {
          result <- compareSurvivalVar(values_selected, values_unselected)
        } else {
          next
        }

        # Store result with effect size (mean_diff, odds_ratio, or hazard_ratio)
        row_data <- data.frame(
          variable_name = col_name,
          variable_type = var_type,
          test_method = result$test_method,
          p_value = result$p_value,
          stringsAsFactors = FALSE
        )

        # Add effect size column based on variable type
        if (var_type == "numeric") {
          row_data$effect_size = result$mean_diff
        } else if (var_type == "categorical") {
          row_data$effect_size = result$odds_ratio
        } else if (var_type == "survival") {
          row_data$effect_size = result$hazard_ratio
        }

        results_list[[col_name]] <- row_data
      }

      # Combine all results
      if (length(results_list) == 0) {
        return(data.frame(
          variable_name = character(),
          variable_type = character(),
          test_method = character(),
          p_value = numeric(),
          effect_size = numeric(),
          fdr_p_value = numeric(),
          stringsAsFactors = FALSE
        ))
      }

      results_df <- do.call(rbind, results_list)
      rownames(results_df) <- NULL

      # FDR correction
      results_df$fdr_p_value <- p.adjust(results_df$p_value, method = "fdr")

      # Sort by FDR
      results_df <- results_df[order(results_df$fdr_p_value, decreasing = FALSE), ]

      return(results_df)
    }

    # Function: Compare all features using multi.t.test
    compareFeatures <- function(expr, pd, selected_idx, unselected_idx) {
      # Create grouping vector
      groups <- rep("unselected", ncol(expr))
      groups[selected_idx] <- "selected"

      # Create temporary phenoData with group column
      temp_pd <- data.frame(group = factor(groups, levels = c("selected", "unselected")))
      rownames(temp_pd) <- colnames(expr)

      # Define comparison
      compare_spec <- c("group", "selected", "unselected")

      # Use multi.t.test
      test_results <- tryCatch(
        multi.t.test(x = expr, pheno = temp_pd, compare = compare_spec, fillNA = TRUE),
        error = function(e) NULL
      )

      if (is.null(test_results)) {
        return(data.frame(
          feature_name = character(),
          mean_selected = numeric(),
          mean_unselected = numeric(),
          mean_diff = numeric(),
          n_selected = integer(),
          n_unselected = integer(),
          p_value = numeric(),
          fdr_p_value = numeric(),
          stringsAsFactors = FALSE
        ))
      }

      # Extract relevant columns
      results_df <- data.frame(
        feature_name = rownames(test_results),
        mean_selected = test_results[["mean|group|selected"]],
        mean_unselected = test_results[["mean|group|unselected"]],
        mean_diff = test_results[["ttest|selected_vs_unselected|mean.diff"]],
        n_selected = test_results[["n value|group|selected"]],
        n_unselected = test_results[["n value|group|unselected"]],
        p_value = test_results[["ttest|selected_vs_unselected|pvalue"]],
        fdr_p_value = test_results[["ttest|selected_vs_unselected|fdr"]],
        stringsAsFactors = FALSE
      )

      # Sort by FDR
      results_df <- results_df[order(results_df$fdr_p_value, decreasing = FALSE), ]
      rownames(results_df) <- NULL

      return(results_df)
    }

    # Reactive: Get selected and unselected sample indices
    sample_groups <- reactive({
      req(reactive_i_samples())
      req(reactive_phenoData())

      all_idx <- seq_len(nrow(reactive_phenoData()))
      selected_idx <- reactive_i_samples()
      unselected_idx <- setdiff(all_idx, selected_idx)

      list(
        selected = selected_idx,
        unselected = unselected_idx
      )
    })

    # Reactive: Phenotype comparison results
    phenotype_results <- reactive({
      req(input$show_phenotype)
      req(sample_groups())
      req(reactive_phenoData())

      comparePhenotypeVars(
        pd = reactive_phenoData(),
        selected_idx = sample_groups()$selected,
        unselected_idx = sample_groups()$unselected
      )
    })

    # Reactive: Feature comparison results
    features_results <- reactive({
      req(input$show_features)
      req(sample_groups())
      req(reactive_expr())

      compareFeatures(
        expr = reactive_expr(),
        pd = reactive_phenoData(),
        selected_idx = sample_groups()$selected,
        unselected_idx = sample_groups()$unselected
      )
    })

    # Output: Phenotype results UI
    output$phenotype_results <- renderUI({
      if (!input$show_phenotype) return(NULL)

      tagList(
        h3("Phenotype Variables Comparison"),
        dataTableDownload_ui(ns("phenotype_table"))
      )
    })

    # Output: Features results UI
    output$features_results <- renderUI({
      if (!input$show_features) return(NULL)

      tagList(
        h3("Features/Expression Comparison"),
        dataTableDownload_ui(ns("features_table"))
      )
    })

    # Render phenotype table and capture selection
    phenotype_selection <- dataTableDownload_module(
      "phenotype_table",
      reactive_table = phenotype_results,
      prefix = "batch_comparison_phenotype_"
    )

    # Render features table and capture selection
    features_selection <- dataTableDownload_module(
      "features_table",
      reactive_table = features_results,
      prefix = "batch_comparison_features_"
    )

    # Track the last selected row from either table
    last_selected <- reactiveVal(NULL)

    # Update when phenotype table selection changes
    observe({
      req(input$show_phenotype)
      sel <- phenotype_selection()
      if (!is.null(sel) && length(sel) > 0) {
        # Get the actual row data
        row_idx <- sel[1]  # Get first selected row
        row_data <- phenotype_results()[row_idx, ]

        last_selected(list(
          source = "phenotype",
          row_index = row_idx,
          data = row_data
        ))
      }
    })

    # Update when features table selection changes
    observe({
      req(input$show_features)
      sel <- features_selection()
      if (!is.null(sel) && length(sel) > 0) {
        # Get the actual row data
        row_idx <- sel[1]  # Get first selected row
        row_data <- features_results()[row_idx, ]

        last_selected(list(
          source = "features",
          row_index = row_idx,
          data = row_data
        ))
      }
    })
    
    observe(last_selected())

    # Return reactive containing last selected row info
    return(reactive({
      last_selected()
    }))

  }) # end moduleServer
}
