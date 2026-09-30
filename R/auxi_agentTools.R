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
#' become \code{NULL} -- but only in \emph{object field} positions,
#' never inside data arrays: a feature list like
#' \code{["TP53", "NULL"]} is data whose entries may legitimately equal
#' those strings, and silently dropping them corrupts the selection
#' (todo 3.3). Empty arrays and arrays whose every element is absent
#' stay explicit empty lists -- an empty array is a real value
#' (\dQuote{clear the selection}, \dQuote{reject via the min rule}),
#' distinct from an omitted optional.
#'
#' The function is idempotent: a canonical document is returned
#' unchanged, so unit tests and programmatic callers can feed documents
#' straight through it before validators run.
#'
#' @param x A parsed JSON value, a named list of tool arguments, or any
#'   already-canonical R value.
#' @param depth Recursion guard (internal).
#' @param field_position Whether \code{x} sits in an object-field
#'   position (sentinel strings map to \code{NULL}); array elements and
#'   the top-level document keep literal strings (internal).
#' @return The canonical representation of \code{x}.
#' @keywords internal
#' @rdname agentToolHelpers
agent_args_sanitize <- function(x, depth = 0L, field_position = TRUE) {
  if (depth > .agent_args_max_depth)
    stop("Tool arguments nest at most ", .agent_args_max_depth, " levels deep.")
  if (is.null(x))
    return(NULL)
  if (!is.list(x))
    return(if (field_position && agent_sentinel_string(x)) NULL else x)

  nms <- names(x)
  if (is.null(nms) || any(!nzchar(nms))) {
    # JSON array: drop absent elements; a homogeneous scalar array
    # becomes an atomic vector (the shape every array-argument consumer
    # already expects), anything else stays a list of elements. An
    # empty (or all-absent) array stays an explicit empty list -- it is
    # a real value (e.g. "clear the selection" / "reject via the min
    # rule"), distinct from an omitted optional. Array ELEMENTS are
    # data (ids, values): the sentinel-string rule never applies to
    # them, only to the fields of objects nested inside the array.
    items <- lapply(x, function(el) {
      if (is.list(el) && !is.null(names(el)) && all(nzchar(names(el))))
        agent_args_sanitize(el, depth + 1L, field_position = TRUE)
      else
        agent_args_sanitize(el, depth + 1L, field_position = FALSE)
    })
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
    lapply(x, agent_args_sanitize, depth = depth + 1L, field_position = TRUE)
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

############################################################################
### [Stage 2 / todo 3.2] output codec at the same seam
###
### WP15 fixed the INPUT seam; this is the OUTPUT seam. Tool values were
### handed to ellmer as raw R objects, and the provider body serializer
### (jsonlite::toJSON(auto_unbox = TRUE)) mangled them in flight:
###   * named vectors lost their names ("dimensions":[2702,60]);
###   * NULL fields serialized as {} (a main source of the literal-"{}"
###     sentinel strings the input seam then has to strip);
###   * nothing bounded the size -- one oversized tool result inflates
###     the context (volcano echoes grew with the dataset, todo 3.1).
### The codec converts every non-string, non-error tool value to a
### jsonlite `json` string ONCE, at the seam, so what the model receives
### is exactly what we serialized: names kept, nulls null, size capped.
### `tool_string()` returns a json-class value verbatim and the request
### body embeds it as the tool-output string -- ellmer's own documented
### "return toJSON(...) from a tool" path. Internal consumers (context
### stubber, compaction summariser, logging) read the structured copy
### kept in `extra$data` (todo 3.4), which also survives the snapshot
### slim path.

#' Tool-result output byte cap
#'
#' Env knob \code{OMICSVIEWER_LLM_TOOL_OUTPUT_BYTES} (default 8192,
#' clamped to [1024, 65536]): every model-facing tool result value is at
#' most this many bytes; oversized results are replaced by a truncation
#' marker that tells the model how to narrow the request.
#' @keywords internal
#' @rdname agentToolHelpers
.agent_tool_output_max_bytes <- function() {
  raw <- suppressWarnings(as.integer(Sys.getenv(
    "OMICSVIEWER_LLM_TOOL_OUTPUT_BYTES", "")))
  if (length(raw) != 1L || is.na(raw))
    return(8192L)
  max(1024L, min(65536L, raw))
}

#' Prepare an R value for faithful JSON serialization
#'
#' jsonlite's \code{auto_unbox} drops the names of atomic vectors
#' (\code{c(a = 1)} becomes the scalar \code{1}); converting named
#' vectors to named lists restores the object shape the value intended.
#' Everything else passes through recursively (explicit \code{NULL}s are
#' kept and serialize as \code{null} with \code{null = "null"}).
#'
#' @param x An R value from a tool handler.
#' @return A JSON-faithful R representation.
#' @keywords internal
#' @rdname agentToolHelpers
agent_json_prep <- function(x) {
  if (is.list(x)) {
    out <- lapply(x, agent_json_prep)
    if (!is.null(names(x)))
      names(out) <- names(x)
    return(out)
  }
  if (is.atomic(x) && length(x) >= 1L && !is.null(names(x)))
    return(as.list(x))
  x
}

#' Apply the output codec to one tool result
#'
#' Runs inside the \code{\link{agent_tool}} wrapper on every handler
#' return value: non-string, non-error \code{ContentToolResult} values
#' become \code{json}-class strings (names preserved, \code{NULL} as
#' \code{null}), bounded by \code{max_bytes} (default from
#' \code{OMICSVIEWER_LLM_TOOL_OUTPUT_BYTES}); the structured original is
#' kept in \code{extra$data} for internal consumers. Errored results,
#' plain strings, and values that are already \code{json} pass through
#' unchanged. Results that are not \code{ContentToolResult} objects are
#' returned untouched (ellmer wraps plain strings itself).
#'
#' @param result Handler return value.
#' @param max_bytes Model-facing byte cap for the serialized value.
#' @return The codec-processed result.
#' @keywords internal
#' @rdname agentToolHelpers
agent_result_codec <- function(result, max_bytes = .agent_tool_output_max_bytes()) {
  if (!inherits(result, "ellmer::ContentToolResult"))
    return(result)
  if (!is.null(result@error))
    return(result)
  value <- result@value
  if (is.character(value) && length(value) == 1L && !inherits(value, "json"))
    return(result)  # plain string output: already the wire shape

  extra <- result@extra
  if (is.null(extra)) extra <- list()
  if (is.null(extra$data) && !is.null(value))
    extra$data <- value

  encoded <- tryCatch(
    jsonlite::toJSON(
      agent_json_prep(value), auto_unbox = TRUE, null = "null", na = "null"),
    error = function(e) NULL)
  if (is.null(encoded))
    return(result)  # unserializable value: let ellmer's own path handle it

  if (nchar(encoded, type = "bytes") <= max_bytes) {
    out_value <- encoded
  } else {
    head_chars <- max(200L, floor(max_bytes / 4L))
    out_value <- paste0(
      "Tool result truncated: ", nchar(encoded, type = "bytes"),
      " bytes exceeded the ", max_bytes,
      "-byte tool-output limit. Narrow the request (fewer ids, sections,",
      " rows, or columns) and call the tool again. First ", head_chars,
      " characters:\n", substr(encoded, 1L, head_chars)
    )
  }
  ellmer::ContentToolResult(
    value = out_value,
    error = result@error,
    extra = extra,
    request = result@request
  )
}

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
  codec <- agent_result_codec
  handler <- eval(bquote(function() {
    args <- .(sanitize)(as.list(environment()))
    # Drop absent top-level arguments before dispatch, exactly as
    # ellmer's convert = TRUE path did: a handler's own formal default
    # (e.g. max_results = 20L) then applies instead of an explicit NULL.
    out <- do.call(.(fun), args[!vapply(args, is.null, logical(1))])
    # Output seam (todo 3.2): one canonical serialization of the value
    # the model receives -- names preserved, nulls null, size capped.
    .(codec)(out)
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
