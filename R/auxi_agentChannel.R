#' @name agentChannel
#' @title Transport-neutral agent channel (Stage 0)
#' @description Stage 0 of the agent-channel plan (AGENT_CHANNEL_PLAN.md):
#'   a transport-neutral dispatcher, a JSON-Schema emitter, and a
#'   session-token registry that together make the assistant tool surface
#'   a supported extension point. Everything here is a \emph{new caller}
#'   of the existing validated plane: ToolDefs built by
#'   \code{\link{agent_tool}} (sanitize-in / codec-out) are invoked
#'   directly, so in-process agents (Stage 1) and the external channel
#'   (Stages 2-3, JSON-RPC over WS) drive exactly the pathway a model
#'   drives. No reactive or store internals are touched; handlers keep
#'   their own \code{isolate}/\code{withReactiveDomain} discipline.
#'
#'   Plan decisions honoured here: numeric ranges stay description-only
#'   in emitted schemas (D15); nothing in this file edits the widget
#'   store or the sanitizer seam (cluster rule); serialization to JSON
#'   text happens only at the transport edge - the returned specs are
#'   JSON-ready R lists (jsonlite::toJSON with \code{auto_unbox = TRUE},
#'   \code{null = "null"} projects them directly).
NULL

############################################################################
### [Stage 0] JSON-Schema emitter
###
### Converts ellmer 0.5.0 type objects (S7 classes; introspected slots)
### into plain named lists shaped as JSON Schema fragments:
###
###   TypeBasic  -> {type}
###   TypeEnum   -> {type: "string", enum: values}
###   TypeArray  -> {type: "array", items: <spec>}
###   TypeObject -> {type: "object", properties: {...},
###                  required: [...], additionalProperties: false}
###
### A property's own `required` slot is argument-position metadata: it
### is expressed in the PARENT's required array, never inside the
### emitted property spec. Empty/NULL descriptions are omitted. Numeric
### ranges live in descriptions only (D15) - no minimum/maximum keys.
###
### Edge note for Stage 2: length-1 character vectors (a single enum
### value, a single required name) unbox to JSON scalars under
### auto_unbox = TRUE; the WS edge must wrap those with I() when it
### serializes tools/list.

#' Convert one ellmer type object to a JSON-Schema fragment
#'
#' Recursive converter behind \code{\link{agent_channel_toolspecs}}.
#' Branches on the ellmer 0.5.0 S7 class strings via base
#' \code{inherits()} (S7 objects carry their class attributes, so no
#' runtime S7 dependency is needed). Unknown type classes are a loud
#' error, not a silently degraded schema: if a future ellmer type (or a
#' new constructor in the tool registry) reaches this converter, the
#' caller fails fast instead of advertising a wrong schema. The guard
#' suite pins this with a static constructor-coverage sweep over
#' auxi_agentToolRegistry.R / auxi_agentFigureGrammar.R.
#'
#' @param type An ellmer type object (TypeBasic, TypeEnum, TypeArray,
#'   or TypeObject).
#' @return A JSON-ready named list (JSON-Schema fragment); the
#'   \code{required} flag of \code{type} itself is not part of the
#'   fragment (the parent carries it).
#' @keywords internal
#' @rdname agentChannel
.agent_channel_type_spec <- function(type) {
  if (inherits(type, "ellmer::TypeEnum")) {
    spec <- list(type = "string", enum = type@values)
  } else if (inherits(type, "ellmer::TypeArray")) {
    spec <- list(type = "array",
                 items = .agent_channel_type_spec(type@items))
  } else if (inherits(type, "ellmer::TypeObject")) {
    properties <- lapply(type@properties, .agent_channel_type_spec)
    required <- names(type@properties)[vapply(
      type@properties,
      function(property) isTRUE(property@required),
      logical(1))]
    spec <- list(type = "object")
    if (length(properties))
      spec$properties <- properties
    if (length(required))
      spec$required <- required
    if (identical(type@additional_properties, FALSE))
      spec$additionalProperties <- FALSE
  } else if (inherits(type, "ellmer::TypeBasic")) {
    spec <- list(type = type@type)
  } else {
    stop("Unsupported ellmer type for channel JSON-Schema emission: ",
         paste(class(type), collapse = "/"),
         ". Add a converter branch before exposing it on the channel.")
  }
  description <- type@description
  if (length(description) == 1L && !is.na(description) && nzchar(description))
    spec$description <- description
  spec
}

#' Emit JSON-Schema tool specs for a list of ToolDefs
#'
#' Takes the named list of ellmer ToolDefs that
#' \code{agent_tool_registry(...)$tools} returns and projects every tool
#' onto the wire shape a channel client needs:
#' \code{list(name, description, inputSchema)} per tool, where
#' \code{inputSchema} is the JSON-Schema object for the tool's arguments
#' (properties, required array, enum values, nested objects/arrays,
#' \code{additionalProperties} reflecting ellmer's
#' \code{additional_properties} slot). The result is keyed by each
#' ToolDef's wire name (\code{@name}, falling back to the list key), so
#' \code{names(agent_channel_toolspecs(tools))} lines up with the names
#' \code{\link{agent_channel_dispatch}} accepts. The return value is
#' JSON-ready R lists, not strings; serialize with
#' \code{jsonlite::toJSON(..., auto_unbox = TRUE, null = "null")} only
#' at the transport edge.
#'
#' @param tools Named list of ellmer ToolDefs (the
#'   \code{$tools} element of \code{agent_tool_registry(...)}).
#' @return Named list of JSON-Schema tool specs keyed by wire tool
#'   name; \code{description} is omitted when the ToolDef carries none.
#' @keywords internal
#' @rdname agentChannel
agent_channel_toolspecs <- function(tools) {
  if (!is.list(tools))
    stop("agent_channel_toolspecs: tools must be a named list of ToolDefs.")
  if (!length(tools))
    return(list())
  specs <- lapply(seq_along(tools), function(i) {
    tool_def <- tools[[i]]
    key <- names(tools)[[i]]
    if (!inherits(tool_def, "ellmer::ToolDef"))
      stop("agent_channel_toolspecs: element '", key,
           "' is not an ellmer::ToolDef.")
    name <- tool_def@name
    if (length(name) != 1L || is.na(name) || !nzchar(name))
      name <- key
    spec <- list(name = name)
    description <- tool_def@description
    if (length(description) == 1L && !is.na(description) &&
        nzchar(description))
      spec$description <- description
    spec$inputSchema <- .agent_channel_type_spec(tool_def@arguments)
    spec
  })
  names(specs) <- vapply(specs, function(spec) spec$name, character(1))
  specs
}

############################################################################
### [Stage 0] transport-neutral dispatcher
###
### The verified transport-clean entry into the tool surface: call the
### ToolDef directly as a function (module_agentTestHooks.R reuses the
### same surface; the ToolDef wraps sanitize-in / codec-out, so every
### caller sees the canonical argument shape and the bounded output).
### Direct ToolDef calls let handler conditions ESCAPE (unlike ellmer's
### invoke_tool), so the dispatcher catches them and always answers
### with a result record - a channel transport must never propagate an
### R condition across the wire.

#' Dispatch one tool call through the channel plane
#'
#' Calls the requested ToolDef from \code{tools} directly with the
#' named arguments \code{args} (the JSON-decoded object shape; explicit
#' \code{NULL} entries are dropped before the call, mirroring
#' \code{\link{agent_tool}}'s absent-argument semantics so handler
#' defaults apply). The ToolDef's own wrapper runs the argument
#' sanitizer and the output codec, so sentinel strings
#' (\code{AGENT_SENTINEL_STRINGS}) are treated as omitted and oversized
#' values come back as bounded truncation markers.
#'
#' Conditions never escape: handler errors, errored tool results, and
#' malformed requests all return \code{ok = FALSE} with the condition
#' message in \code{error}. An unknown \code{name} fails with the list
#' of available tools so an agent can self-correct. Tools are matched
#' by wire name first (the ToolDef's \code{@name}), then by their list
#' key, so both the channel client and internal callers can dispatch.
#'
#' @param tools Named list of ellmer ToolDefs (the
#'   \code{$tools} element of \code{agent_tool_registry(...)}).
#' @param name Tool name (single non-empty string).
#' @param args Named list of tool arguments in the JSON-decoded shape
#'   (object fields, atomic scalars/vectors, \code{NULL} for absent);
#'   \code{NULL} is treated as an empty call.
#' @return \code{list(ok = TRUE, value = <ToolDef result>, error =
#'   NULL)} on success - \code{value} is the full
#'   \code{ellmer::ContentToolResult} (model-facing \code{@value} plus
#'   the structured original in \code{extra$data}) - or
#'   \code{list(ok = FALSE, value = NULL, error = <message>)}.
#' @keywords internal
#' @rdname agentChannel
agent_channel_dispatch <- function(tools, name, args) {
  fail <- function(message) list(ok = FALSE, value = NULL, error = message)
  if (!is.list(tools))
    return(fail("tools must be a named list of ToolDefs."))
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(name))
    return(fail("Tool name must be a single non-empty string."))
  if (is.null(args))
    args <- list()
  if (!is.list(args))
    return(fail(
      "Tool arguments must be a named list (the JSON object shape) or NULL."))
  arg_names <- names(args)
  if (length(args) && (is.null(arg_names) || any(!nzchar(arg_names))))
    return(fail(
      "Tool arguments must be a fully named list; unnamed entries cannot be dispatched."))

  wire_names <- vapply(tools, function(tool_def) {
    if (inherits(tool_def, "ellmer::ToolDef") && length(tool_def@name) == 1L &&
        !is.na(tool_def@name) && nzchar(tool_def@name))
      tool_def@name
    else
      NA_character_
  }, character(1))
  index <- match(name, wire_names)
  if (is.na(index) && !is.null(names(tools)))
    index <- match(name, names(tools))
  if (is.na(index))
    return(fail(paste0(
      "Unknown tool '", name, "'. Available tools: ",
      paste(sort(unique(c(wire_names[!is.na(wire_names)], names(tools)))),
            collapse = ", "), ".")))

  # Absent-argument semantics (mirrors agent_tool): drop explicit NULLs
  # so the handler's own formal defaults apply.
  args <- args[!vapply(args, is.null, logical(1))]
  result <- tryCatch(
    do.call(tools[[index]], args),
    error = function(e) e)
  if (inherits(result, "condition"))
    return(fail(conditionMessage(result)))
  if (inherits(result, "ellmer::ContentToolResult") &&
      !is.null(result@error))
    return(fail(conditionMessage(result@error)))
  list(ok = TRUE, value = result, error = NULL)
}

############################################################################
### [Stage 0] session-token registry
###
### Package-env registry keyed by session token; the caller (assistant
### module init, D14) stores whatever the channel endpoint needs
### (registry/tools/session_domain/logger) and wires session$onEnded
### eviction through `on_ended`. The registry itself stays dumb: no
### reactivity, no policy, no logging - multi-tab refusal (D4) and rate
### caps (D6) live at the endpoint, not here.

#' The channel registry environment (token -> entry)
#' @keywords internal
#' @rdname agentChannel
.agent_channel_registry <- new.env(parent = emptyenv())

#' Bind a channel entry under a session token
#'
#' Stores \code{entry} (an arbitrary list - callers keep their own
#' registry/tools/session_domain/logger) under \code{token},
#' overwriting any previous entry for the same token. When
#' \code{on_ended} is non-NULL it is called exactly once, immediately
#' after storing; the caller passes its session-lifetime wiring (e.g.
#' \code{function() session$onEnded(function() unbind(token))}). The
#' registry keeps no policy and never calls back on its own.
#'
#' @param token Session token (single non-empty string).
#' @param entry Arbitrary list to store.
#' @param on_ended Optional callback invoked right after storing
#'   (caller-supplied session$onEnded wiring); \code{NULL} (default)
#'   binds without lifetime wiring.
#' @return \code{token}, invisibly.
#' @keywords internal
#' @rdname agentChannel
agent_channel_registry_bind <- function(token, entry, on_ended = NULL) {
  if (!is.character(token) || length(token) != 1L || is.na(token) ||
      !nzchar(token))
    stop("agent_channel_registry_bind: token must be a single non-empty string.")
  if (!is.list(entry))
    stop("agent_channel_registry_bind: entry must be a list.")
  assign(token, entry, envir = .agent_channel_registry)
  if (!is.null(on_ended))
    on_ended()
  invisible(token)
}

#' Remove a channel registry entry
#'
#' Unbinds \code{token}; absent tokens are removed silently (idempotent
#' - a session that already ended never fails its own eviction).
#'
#' @param token Session token to unbind.
#' @return \code{NULL}, invisibly.
#' @keywords internal
#' @rdname agentChannel
agent_channel_registry_unbind <- function(token) {
  if (exists(token, envir = .agent_channel_registry, inherits = FALSE))
    rm(list = token, envir = .agent_channel_registry)
  invisible(NULL)
}

#' List bound channel tokens
#'
#' @return Sorted character vector of currently bound session tokens
#'   (empty character vector when none).
#' @keywords internal
#' @rdname agentChannel
agent_channel_registry_list <- function() {
  ls(.agent_channel_registry, sorted = TRUE)
}
