#' Error Handling Utilities for API Calls
#'
#' @description
#' Utility functions for robust error handling of external API calls.
#' Provides consistent error messages and graceful degradation when
#' network requests fail.
#'
#' @name error_handling
#' @keywords internal
NULL

#' Safe HTTP GET wrapper with error handling
#'
#' @description
#' Wraps httr::GET with comprehensive error handling including:
#' - Network connectivity errors
#' - HTTP status code errors (4xx, 5xx)
#' - Timeout handling
#' - SSL/TLS errors
#'
#' @param url The URL to fetch
#' @param query Query parameters (optional)
#' @param timeout Timeout in seconds (default: 30)
#' @param api_name Name of the API for error messages (default: "API")
#' @param ... Additional parameters passed to httr::GET
#'
#' @return A list with:
#'   - success: logical, TRUE if request succeeded
#'   - data: response content if successful, NULL otherwise
#'   - error: error message if failed, NULL otherwise
#'   - status_code: HTTP status code if available
#'
#' @importFrom httr GET content status_code http_error http_status timeout
#' @keywords internal
#'
safe_GET <- function(url, query = NULL, timeout = 30, api_name = "API", ...) {

  result <- list(
    success = FALSE,
    data = NULL,
    error = NULL,
    status_code = NULL
  )

  tryCatch({
    # Perform GET request with timeout
    response <- httr::GET(
      url = url,
      query = query,
      httr::timeout(timeout),
      ...
    )

    result$status_code <- httr::status_code(response)

    # Check for HTTP errors
    if (httr::http_error(response)) {
      status_info <- httr::http_status(response)
      result$error <- sprintf(
        "%s request failed with HTTP %d: %s",
        api_name,
        result$status_code,
        status_info$message
      )
      return(result)
    }

    # Extract content
    content_data <- httr::content(response)

    result$success <- TRUE
    result$data <- content_data

  }, error = function(e) {
    error_msg <- conditionMessage(e)

    # Provide more specific error messages
    if (grepl("timeout|timed out", error_msg, ignore.case = TRUE)) {
      result$error <- sprintf(
        "%s request timed out after %d seconds. Please check your internet connection and try again.",
        api_name, timeout
      )
    } else if (grepl("could not resolve host|name or service not known", error_msg, ignore.case = TRUE)) {
      result$error <- sprintf(
        "Cannot connect to %s. Please check your internet connection.",
        api_name
      )
    } else if (grepl("SSL|certificate", error_msg, ignore.case = TRUE)) {
      result$error <- sprintf(
        "%s SSL/certificate error: %s",
        api_name, error_msg
      )
    } else {
      result$error <- sprintf(
        "%s request failed: %s",
        api_name, error_msg
      )
    }
  })

  return(result)
}

#' Safe HTTP POST wrapper with error handling
#'
#' @description
#' Wraps httr::POST with comprehensive error handling including:
#' - Network connectivity errors
#' - HTTP status code errors (4xx, 5xx)
#' - Timeout handling
#' - SSL/TLS errors
#'
#' @param url The URL to post to
#' @param body Request body
#' @param encode Encoding method (default: "json")
#' @param timeout Timeout in seconds (default: 30)
#' @param api_name Name of the API for error messages (default: "API")
#' @param ... Additional parameters passed to httr::POST
#'
#' @return A list with:
#'   - success: logical, TRUE if request succeeded
#'   - data: response content if successful, NULL otherwise
#'   - error: error message if failed, NULL otherwise
#'   - status_code: HTTP status code if available
#'
#' @importFrom httr POST content status_code http_error http_status timeout
#' @keywords internal
#'
safe_POST <- function(url, body = NULL, encode = "json", timeout = 30, api_name = "API", ...) {

  result <- list(
    success = FALSE,
    data = NULL,
    error = NULL,
    status_code = NULL
  )

  tryCatch({
    # Perform POST request with timeout
    response <- httr::POST(
      url = url,
      body = body,
      encode = encode,
      httr::timeout(timeout),
      ...
    )

    result$status_code <- httr::status_code(response)

    # Check for HTTP errors
    if (httr::http_error(response)) {
      status_info <- httr::http_status(response)
      result$error <- sprintf(
        "%s request failed with HTTP %d: %s",
        api_name,
        result$status_code,
        status_info$message
      )
      return(result)
    }

    # Extract content
    content_data <- httr::content(response)

    result$success <- TRUE
    result$data <- content_data

  }, error = function(e) {
    error_msg <- conditionMessage(e)

    # Provide more specific error messages
    if (grepl("timeout|timed out", error_msg, ignore.case = TRUE)) {
      result$error <- sprintf(
        "%s request timed out after %d seconds. Please check your internet connection and try again.",
        api_name, timeout
      )
    } else if (grepl("could not resolve host|name or service not known", error_msg, ignore.case = TRUE)) {
      result$error <- sprintf(
        "Cannot connect to %s. Please check your internet connection.",
        api_name
      )
    } else if (grepl("SSL|certificate", error_msg, ignore.case = TRUE)) {
      result$error <- sprintf(
        "%s SSL/certificate error: %s",
        api_name, error_msg
      )
    } else {
      result$error <- sprintf(
        "%s request failed: %s",
        api_name, error_msg
      )
    }
  })

  return(result)
}

#' Validate API response data
#'
#' @description
#' Validates that API response data is in expected format and not empty.
#'
#' @param data The response data to validate
#' @param expected_type Expected type ("data.frame", "list", "character", etc.)
#' @param api_name Name of the API for error messages
#'
#' @return A list with:
#'   - valid: logical, TRUE if data is valid
#'   - error: error message if invalid, NULL otherwise
#'
#' @keywords internal
#'
validate_api_response <- function(data, expected_type = "data.frame", api_name = "API") {

  result <- list(
    valid = FALSE,
    error = NULL
  )

  # Check if data is NULL
  if (is.null(data)) {
    result$error <- sprintf("%s returned no data.", api_name)
    return(result)
  }

  # Check expected type
  if (expected_type == "data.frame" && !is.data.frame(data)) {
    result$error <- sprintf(
      "%s returned unexpected data format (expected data.frame, got %s).",
      api_name, class(data)[1]
    )
    return(result)
  }

  if (expected_type == "list" && !is.list(data)) {
    result$error <- sprintf(
      "%s returned unexpected data format (expected list, got %s).",
      api_name, class(data)[1]
    )
    return(result)
  }

  # Check if data is empty
  if (is.data.frame(data) && nrow(data) == 0) {
    result$error <- sprintf("%s returned empty results.", api_name)
    return(result)
  }

  if (is.list(data) && length(data) == 0) {
    result$error <- sprintf("%s returned empty results.", api_name)
    return(result)
  }

  result$valid <- TRUE
  return(result)
}
