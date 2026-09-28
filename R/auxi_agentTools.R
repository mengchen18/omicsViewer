#' @name agentToolHelpers
#' @title Canonical agent tool-argument transport (the boundary seam)
#' @description One anti-corruption layer between ellmer and every
#'   assistant tool handler. Tools registered through
#'   \code{\link{agent_tool}} receive their arguments in the canonical
#'   JSON-document shape (named lists, atomic vectors, NULL for absent),
#'   regardless of how a provider or ellmer version materializes the
#'   wire payload. See section WP15 of AGENT_ACCURACY_PLAN.md.
NULL

############################################################################
### [WP15] canonical tool-argument transport (the agent boundary seam)
###
### One anti-corruption layer between ellmer and every tool handler.
###
### Problem this solves (AGENT_ACCURACY_PLAN.md section WP15): ellmer's
### default `convert = TRUE` materializes JSON tool arguments into R shapes
### that no validator was designed for -- `type_array(type_object)`
### arguments become column-wise tibbles (a layer list whose length()
### counts COLUMNS, an absent nested object becomes a 1-row all-NA
### tibble, a nested combinator array becomes list(<tibble>)). Unit tests
### feed plain nested lists, production received tibbles: the tested shape
### was never the shipped shape. Some providers additionally serialize
### omitted optionals as literal "null"/"{}"/"[]" strings
### (AGENT_SENTINEL_STRINGS).
###
### The seam installs ONE canonical representation -- the JSON document
### shape (named lists, atomic scalars and vectors, NULL for absent) --
### and deletes the per-consumer shape knowledge instead of extending it:
###
###   * tools are registered with `ellmer::tool(convert = FALSE)`, so
###     handlers receive the raw parsed JSON (jsonlite simplifyVector =
###     FALSE) rather than ellmer's tibble/factor conversion;
###   * `agent_args_sanitize()` runs once at tool entry and restores the
###     canonical shape the raw JSON almost has: homogeneous scalar
###     arrays become atomic vectors, empty and all-null arrays become
###     NULL, and provider sentinel strings become NULL;
###   * validators keep domain logic only (allowlists, ranges, operators)
###     and never see a tibble, a factor, or a sentinel string again.
###
### `agent_args_sanitize` is idempotent: canonicalizing an already
### canonical document returns it unchanged, so direct programmatic
### callers and the unit tests construct the very shape production sees.

#' Canonicalize one JSON tool-argument value
#'
#' Recursively maps a parsed JSON value (or any near-canonical R value)
#' onto the canonical document shape used by every agent tool handler:
#' named lists for objects, atomic vectors for homogeneous scalar
#' arrays, and \code{NULL} for absent values. Provider sentinel strings
#' (\code{AGENT_SENTINEL_STRINGS} -- some models serialize omitted
#' optionals as literal \code{"null"}/\code{"{}"}/\code{"[]"} strings)
#' become \code{NULL} at any depth. Empty arrays and arrays whose every
#' element is absent stay explicit empty lists -- an empty array is a
#' real value (\dQuote{clear the selection}, \dQuote{reject via the
#' min rule}), distinct from an omitted optional.
#'
#' The function is idempotent: a canonical document is returned
#' unchanged, so unit tests and programmatic callers can feed documents
#' straight through it before validators run.
#'
#' @param x A parsed JSON value, a named list of tool arguments, or any
#'   already-canonical R value.
#' @param depth Recursion guard (internal).
#' @return The canonical representation of \code{x}.
#' @keywords internal
#' @rdname agentToolHelpers
agent_args_sanitize <- function(x, depth = 0L) {
  if (depth > .agent_args_max_depth)
    stop("Tool arguments nest at most ", .agent_args_max_depth, " levels deep.")
  if (is.null(x))
    return(NULL)
  if (!is.list(x))
    return(if (agent_sentinel_string(x)) NULL else x)

  nms <- names(x)
  if (is.null(nms) || any(!nzchar(nms))) {
    # JSON array: drop absent elements; a homogeneous scalar array
    # becomes an atomic vector (the shape every array-argument consumer
    # already expects), anything else stays a list of elements. An
    # empty (or all-absent) array stays an explicit empty list -- it is
    # a real value (e.g. "clear the selection" / "reject via the min
    # rule"), distinct from an omitted optional.
    items <- lapply(x, agent_args_sanitize, depth = depth + 1L)
    items <- items[!vapply(items, is.null, logical(1))]
    if (!length(items))
      return(list())
    scalar <- vapply(items, function(v) is.atomic(v) && length(v) == 1L,
                     logical(1))
    if (all(scalar))
      return(unlist(items, use.names = FALSE))
    items
  } else {
    # JSON object: recurse per property, keeping explicit nulls as
    # named NULL entries (consumers test absence with is.null()).
    lapply(x, agent_args_sanitize, depth = depth + 1L)
  }
}

#' Depth cap for sanitized tool arguments
#'
#' Real argument documents nest at most a handful of levels (figure
#' specs: layers -> filter -> combinators -> leaves). The cap turns a
#' pathological payload into a self-correctable tool error instead of
#' unbounded recursion.
#' @keywords internal
#' @rdname agentToolHelpers
.agent_args_max_depth <- 64L

#' Register one agent tool on the canonical argument transport
#'
#' Drop-in replacement for \code{ellmer::tool()} for every assistant
#' tool. The handler is wrapped so that exactly one
#' \code{\link{agent_args_sanitize}} pass runs at tool entry, and the
#' tool is registered with \code{convert = FALSE}: ellmer then hands the
#' raw parsed JSON to the wrapper instead of converting
#' \code{type_array(type_object)} arguments into tibbles and enum
#' arrays into factors. Handlers, validators, unit tests, and the WP3
#' echo therefore all see the same canonical shape.
#'
#' Direct calls on the returned ToolDef (as the unit tests make) run
#' through the identical wrapper, so the tested entry is the shipped
#' entry.
#'
#' @param fun Tool handler; formals must match \code{names(arguments)}
#'   exactly (ellmer enforces this).
#' @param description Tool description sent to the model.
#' @param arguments Named list of ellmer type specifications.
#' @param name Tool name (letters, numbers, - and _ only).
#' @param annotations ellmer tool annotations.
#' @return An \code{ellmer::ToolDef} with \code{convert = FALSE}.
#' @keywords internal
#' @rdname agentToolHelpers
agent_tool <- function(fun, description, arguments = list(), name = NULL,
                       annotations = list()) {
  arg_names <- names(arguments)
  # Embed the handler and the sanitizer by value so the wrapper needs no
  # name resolution beyond base R at call time.
  sanitize <- agent_args_sanitize
  handler <- eval(bquote(function() {
    args <- .(sanitize)(as.list(environment()))
    # Drop absent top-level arguments before dispatch, exactly as
    # ellmer's convert = TRUE path did: a handler's own formal default
    # (e.g. max_results = 20L) then applies instead of an explicit NULL.
    do.call(.(fun), args[!vapply(args, is.null, logical(1))])
  }))
  # ellmer::tool() requires the formals to match the declared argument
  # names exactly; every formal defaults to NULL so omitted optionals and
  # explicit JSON nulls reach the dispatch above as absent values.
  formals(handler) <- stats::setNames(
    rep(list(NULL), length(arg_names)), arg_names)
  ellmer::tool(
    handler,
    name = name,
    description = description,
    arguments = arguments,
    convert = FALSE,
    annotations = annotations
  )
}
