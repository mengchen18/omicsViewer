#' Input Validation Utilities
#'
#' @description
#' Provides comprehensive validation functions for Shiny module inputs to ensure
#' data integrity and provide user-friendly error messages.
#'
#' @keywords internal

#' Validate reactive data frame input
#'
#' @param data Reactive expression or data.frame to validate
#' @param name Character. Name of the data for error messages
#' @param required Logical. Whether NULL values are allowed. Default: TRUE
#' @param min_rows Integer. Minimum number of rows required. Default: 1
#' @param min_cols Integer. Minimum number of columns required. Default: 1
#'
#' @return List with $valid (logical) and $message (character) fields
#'
#' @examples
#' \dontrun{
#' validation <- validate_dataframe(
#'   pData(eset),
#'   name = "phenotype data",
#'   min_rows = 3,
#'   min_cols = 2
#' )
#' if (!validation$valid) {
#'   showNotification(validation$message, type = "error")
#'   return(NULL)
#' }
#' }
#'
#' @export
validate_dataframe <- function(data, name = "data", required = TRUE,
                                min_rows = 1, min_cols = 1) {

  # Handle reactive expressions
  if (is.function(data)) {
    data <- try(data(), silent = TRUE)
    if (inherits(data, "try-error")) {
      return(list(
        valid = FALSE,
        message = sprintf("Error evaluating reactive %s", name)
      ))
    }
  }

  # Check for NULL
  if (is.null(data)) {
    if (required) {
      return(list(
        valid = FALSE,
        message = sprintf("%s is required but is NULL", name)
      ))
    }
    return(list(valid = TRUE, message = ""))
  }

  # Check if data.frame
  if (!is.data.frame(data)) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must be a data.frame, got %s", name, class(data)[1])
    ))
  }

  # Check dimensions
  if (nrow(data) < min_rows) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must have at least %d rows, has %d",
                        name, min_rows, nrow(data))
    ))
  }

  if (ncol(data) < min_cols) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must have at least %d columns, has %d",
                        name, min_cols, ncol(data))
    ))
  }

  return(list(valid = TRUE, message = ""))
}


#' Validate numeric matrix input
#'
#' @param mat Reactive expression or matrix to validate
#' @param name Character. Name of the matrix for error messages
#' @param required Logical. Whether NULL values are allowed. Default: TRUE
#' @param min_rows Integer. Minimum number of rows required. Default: 1
#' @param min_cols Integer. Minimum number of columns required. Default: 1
#' @param allow_na Logical. Whether NA values are allowed. Default: TRUE
#'
#' @return List with $valid (logical) and $message (character) fields
#'
#' @export
validate_numeric_matrix <- function(mat, name = "matrix", required = TRUE,
                                     min_rows = 1, min_cols = 1, allow_na = TRUE) {

  # Handle reactive expressions
  if (is.function(mat)) {
    mat <- try(mat(), silent = TRUE)
    if (inherits(mat, "try-error")) {
      return(list(
        valid = FALSE,
        message = sprintf("Error evaluating reactive %s", name)
      ))
    }
  }

  # Check for NULL
  if (is.null(mat)) {
    if (required) {
      return(list(
        valid = FALSE,
        message = sprintf("%s is required but is NULL", name)
      ))
    }
    return(list(valid = TRUE, message = ""))
  }

  # Check if matrix or data.frame with numeric values
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must be a matrix or data.frame, got %s",
                        name, class(mat)[1])
    ))
  }

  # Check if numeric
  if (!is.numeric(mat)) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must be numeric", name)
    ))
  }

  # Check dimensions
  if (nrow(mat) < min_rows) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must have at least %d rows, has %d",
                        name, min_rows, nrow(mat))
    ))
  }

  if (ncol(mat) < min_cols) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must have at least %d columns, has %d",
                        name, min_cols, ncol(mat))
    ))
  }

  # Check for NA values
  if (!allow_na && any(is.na(mat))) {
    return(list(
      valid = FALSE,
      message = sprintf("%s contains NA values which are not allowed", name)
    ))
  }

  return(list(valid = TRUE, message = ""))
}


#' Validate character vector input
#'
#' @param vec Reactive expression or vector to validate
#' @param name Character. Name of the vector for error messages
#' @param required Logical. Whether NULL values are allowed. Default: TRUE
#' @param min_length Integer. Minimum length required. Default: 1
#' @param max_length Integer. Maximum length allowed. Default: Inf
#' @param allow_empty Logical. Whether empty strings are allowed. Default: FALSE
#'
#' @return List with $valid (logical) and $message (character) fields
#'
#' @export
validate_character_vector <- function(vec, name = "vector", required = TRUE,
                                       min_length = 1, max_length = Inf,
                                       allow_empty = FALSE) {

  # Handle reactive expressions
  if (is.function(vec)) {
    vec <- try(vec(), silent = TRUE)
    if (inherits(vec, "try-error")) {
      return(list(
        valid = FALSE,
        message = sprintf("Error evaluating reactive %s", name)
      ))
    }
  }

  # Check for NULL
  if (is.null(vec)) {
    if (required) {
      return(list(
        valid = FALSE,
        message = sprintf("%s is required but is NULL", name)
      ))
    }
    return(list(valid = TRUE, message = ""))
  }

  # Check if character
  if (!is.character(vec)) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must be character, got %s", name, class(vec)[1])
    ))
  }

  # Check length
  if (length(vec) < min_length) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must have at least %d elements, has %d",
                        name, min_length, length(vec))
    ))
  }

  if (length(vec) > max_length) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must have at most %d elements, has %d",
                        name, max_length, length(vec))
    ))
  }

  # Check for empty strings
  if (!allow_empty && any(nchar(vec) == 0)) {
    return(list(
      valid = FALSE,
      message = sprintf("%s contains empty strings which are not allowed", name)
    ))
  }

  return(list(valid = TRUE, message = ""))
}


#' Validate triselector input matrix
#'
#' @param triset Matrix or reactive expression. Should be nx3 matrix from str_split_fixed
#' @param name Character. Name for error messages
#' @param allow_null Logical. Whether NULL is allowed (e.g., during app startup). Default: FALSE
#'
#' @return List with $valid (logical) and $message (character) fields
#'
#' @export
validate_triselector_input <- function(triset, name = "triselector data", allow_null = FALSE) {

  # Handle reactive expressions
  if (is.function(triset)) {
    triset <- try(triset(), silent = TRUE)
    if (inherits(triset, "try-error")) {
      return(list(
        valid = FALSE,
        message = sprintf("Error evaluating reactive %s", name)
      ))
    }
  }

  # Check for NULL
  if (is.null(triset)) {
    if (allow_null) {
      return(list(valid = TRUE, message = ""))
    }
    return(list(
      valid = FALSE,
      message = sprintf("%s is NULL - cannot populate selectors", name)
    ))
  }

  # Check if matrix
  if (!is.matrix(triset)) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must be a matrix, got %s", name, class(triset)[1])
    ))
  }

  # Check dimensions
  if (ncol(triset) != 3) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must have exactly 3 columns, has %d",
                        name, ncol(triset))
    ))
  }

  if (nrow(triset) == 0) {
    return(list(
      valid = FALSE,
      message = sprintf("%s has no rows - no options available", name)
    ))
  }

  # Check for empty values in first column (analysis type)
  if (any(is.na(triset[, 1]) | nchar(triset[, 1]) == 0)) {
    return(list(
      valid = FALSE,
      message = sprintf("%s has empty values in first column (analysis type)", name)
    ))
  }

  return(list(valid = TRUE, message = ""))
}


#' Validate numeric range
#'
#' @param value Numeric value to validate
#' @param name Character. Name for error messages
#' @param min Numeric. Minimum allowed value (inclusive). Default: -Inf
#' @param max Numeric. Maximum allowed value (inclusive). Default: Inf
#' @param allow_na Logical. Whether NA is allowed. Default: FALSE
#'
#' @return List with $valid (logical) and $message (character) fields
#'
#' @export
validate_numeric_range <- function(value, name = "value",
                                    min = -Inf, max = Inf, allow_na = FALSE) {

  # Check for NA
  if (is.na(value)) {
    if (!allow_na) {
      return(list(
        valid = FALSE,
        message = sprintf("%s is NA which is not allowed", name)
      ))
    }
    return(list(valid = TRUE, message = ""))
  }

  # Check if numeric
  if (!is.numeric(value)) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must be numeric, got %s", name, class(value)[1])
    ))
  }

  # Check range
  if (value < min) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must be >= %g, got %g", name, min, value)
    ))
  }

  if (value > max) {
    return(list(
      valid = FALSE,
      message = sprintf("%s must be <= %g, got %g", name, max, value)
    ))
  }

  return(list(valid = TRUE, message = ""))
}


#' Validate matching dimensions
#'
#' @param obj1 First object (matrix or data.frame)
#' @param obj2 Second object (matrix or data.frame)
#' @param dimension Character. Either "rows" or "cols" to check
#' @param name1 Character. Name of first object for error messages
#' @param name2 Character. Name of second object for error messages
#'
#' @return List with $valid (logical) and $message (character) fields
#'
#' @export
validate_matching_dimensions <- function(obj1, obj2, dimension = c("rows", "cols"),
                                          name1 = "object1", name2 = "object2") {

  dimension <- match.arg(dimension)

  # Get dimensions
  if (dimension == "rows") {
    dim1 <- nrow(obj1)
    dim2 <- nrow(obj2)
    dim_name <- "rows"
  } else {
    dim1 <- ncol(obj1)
    dim2 <- ncol(obj2)
    dim_name <- "columns"
  }

  # Check if matching
  if (dim1 != dim2) {
    return(list(
      valid = FALSE,
      message = sprintf("%s and %s must have matching %s: %d vs %d",
                        name1, name2, dim_name, dim1, dim2)
    ))
  }

  return(list(valid = TRUE, message = ""))
}
