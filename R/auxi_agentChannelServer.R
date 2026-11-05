#' @name agentChannelServer
#' @title Localhost JSON-RPC agent channel endpoint (Stage 2)
#' @description Stage 2 of the agent-channel plan (AGENT_CHANNEL_PLAN.md):
#'   a JSON-RPC 2.0-over-WebSocket endpoint on 127.0.0.1 through which an
#'   EXTERNAL agent (editor bridge, Python driver, ...) can drive the
#'   currently-bound assistant session's read-tier tools. The endpoint is a
#'   thin transport shell around the Stage 0 plane: every request is answered
#'   by the pure handler \code{\link{agent_channel_handle}}, which dispatches
#'   through \code{\link{agent_channel_dispatch}} against the session-bound
#'   registry (exactly the pathway a model drives). Nothing here touches the
#'   widget store, the sanitizer seam, or any triselector internals.
#'
#'   Security posture (plan decisions D3/D10): the server binds 127.0.0.1
#'   ONLY - there is no host argument and none may be added without
#'   revisiting this posture; it is compiled in but starts only when
#'   \code{OMICSVIEWER_AGENT_CHANNEL} is truthy; connections must present
#'   the per-process random pairing token as the WebSocket path
#'   (\code{/<token>}) - wrong paths are closed immediately. The FULL token
#'   is printed to the console at server start (that is the pairing
#'   mechanism) and never written to a log or response; log events carry
#'   only a short non-reversible tag (\code{\link{.agent_channel_token_tag}}).
#'   Stage 2 exposes the READ tier only (decision: enabled tier "read");
#'   write tools arrive in Stage 3 behind their own env flag.
NULL

############################################################################
### [Stage 2] env knobs, tokens, tags
###
### Truthiness: a value that is set and not one of the recognized falsy
### words (0/false/f/no/n/off, case-insensitive). Port default 7790 sits
### OUTSIDE the 7775-7786 e2e range (AGENTS.md). The per-connection rate
### cap (D6) defaults to 600 tool calls, 0 = unlimited.

.agent_channel_falsy <- c("", "0", "false", "f", "no", "n", "off")

.agent_channel_enabled <- function() {
  value <- tolower(trimws(Sys.getenv("OMICSVIEWER_AGENT_CHANNEL", unset = "")))
  !(value %in% .agent_channel_falsy)
}

.agent_channel_port_default <- function() {
  raw <- suppressWarnings(as.integer(
    trimws(Sys.getenv("OMICSVIEWER_AGENT_CHANNEL_PORT", unset = ""))))
  if (length(raw) != 1L || is.na(raw) || raw < 1L || raw > 65535L)
    return(7790L)
  raw
}

.agent_channel_max_calls <- function() {
  raw <- suppressWarnings(as.integer(
    trimws(Sys.getenv("OMICSVIEWER_AGENT_CHANNEL_MAX_CALLS", unset = ""))))
  if (length(raw) != 1L || is.na(raw) || raw < 0L)
    return(600L)
  raw
}

#' Channel schema version (D8)
#'
#' Bump on any breaking change to the wire shapes (\code{initialize}
#' result, tool specs, error codes). Clients re-fetch per connection.
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_schema_version <- 1L

#' Enabled tool tiers
#'
#' Stage 2 serves the \code{read} tier only. The \code{write} tier is
#' listed for Stage 3 readiness and is NOT reachable through this
#' server until Stage 3 enables it behind its own env flag
#' (\code{OMICSVIEWER_AGENT_CHANNEL_WRITES}).
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_tiers <- list(
  read = c("get_omics_viewer_state", "search_annotations",
           "summarize_annotation", "find_controls", "describe_control",
           "get_figure"),
  write = c("get_omics_viewer_state", "search_annotations",
            "summarize_annotation", "set_omics_viewer_state",
            "set_scatter_view", "find_controls", "describe_control",
            "set_widgets", "run_control", "set_enrichment_parameters",
            "set_table_view", "create_figure", "get_figure",
            "update_figure"))

#' Generate a fresh pairing token
#'
#' Cryptographically-flavored randomness is not required for a token
#' that is only ever compared inside one process over loopback; a
#' uniform 48-hex-char draw from \code{stats::runif} exceeds the >=32
#' hex-char requirement (D10) by a wide margin and adds no dependency.
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_new_token <- function() {
  paste(sprintf("%04x", as.integer(floor(stats::runif(12L, 0, 65536)))),
        collapse = "")
}

#' Short non-reversible tag for log events
#'
#' Log events must never carry a full token (D5/D10); this rolling
#' polynomial hash mod 2^32 (exact in doubles: every intermediate is
#' < 2^52) yields a stable 16-hex-char correlation tag. It is a log
#' label, not a security primitive.
#' @param token Token (or session token) to tag.
#' @return 16-character hex string.
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_token_tag <- function(token) {
  ints <- utf8ToInt(paste0(as.character(token), collapse = ""))
  if (!length(ints))
    return("0000000000000000")
  h1 <- 0
  h2 <- 0
  for (i in seq_along(ints)) {
    h1 <- (h1 * 1000003 + ints[[i]] + i) %% 4294967296
    h2 <- (h2 * 33 + ints[[i]] * i) %% 4294967296
  }
  # doubles > .Machine$integer.max: hex-format via two exact 16-bit
  # halves (%% and %/% are exact below 2^32 for doubles)
  half16 <- function(v) sprintf("%04x%04x", v %/% 65536, v %% 65536)
  paste0(half16(h1), half16(h2))
}

#' Resolve the package version for `initialize`
#'
#' \code{utils::packageVersion} fails when the package runs from a
#' working tree via pkgload without an installed copy; the loaded
#' namespace records the version it was built from, so fall back to
#' that before giving up with "unknown".
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_package_version <- function() {
  version <- tryCatch(
    as.character(utils::packageVersion("omicsViewer")),
    error = function(e) NULL)
  if (is.null(version) || !length(version) || is.na(version) ||
      !nzchar(version))
    version <- tryCatch({
      meta <- get(".__NAMESPACE__.", envir = getNamespace("omicsViewer"),
                  inherits = TRUE)
      as.character(meta$spec[["version"]])
    }, error = function(e) NULL)
  if (is.null(version) || !length(version) || is.na(version) ||
      !nzchar(version))
    "unknown"
  else
    version
}

############################################################################
### [Stage 2] pure JSON-RPC handler
###
### One function, no sockets, no package-env reads: everything it needs
### rides in `state` (an environment built by agent_channel_state_new /
### the connection layer), and the only mutation it makes is the
### per-connection rate-cap counter (D6) - so the guard suite drives the
### exact wire behavior without a server. The handler NEVER throws: every
### branch answers with a JSON-RPC response string (internal faults map
### to -32603). Wire projection obligations (Stage-0 review F3): tool
### results cross the wire as `value@value` ONLY - `extra$data` and
### `extra$display` are in-process payload that never serializes here.
### Logging (D5) rides the bound session's logger; agent_logger_event is
### failure-silent by design.

#' Construct per-connection handler state
#'
#' The connection layer (and the guard suite) builds one state per
#' WebSocket connection: the bound session's tools, its logger, the
#' served tier, and the per-connection rate cap. \code{tools = NULL}
#' means "no session bound" - \code{tools/list} and \code{tools/call}
#' answer -32000 while \code{initialize} keeps working. The state is an
#' environment because \code{\link{agent_channel_handle}} increments
#' \code{state$calls} for every executed \code{tools/call} (the only
#' mutation); \code{max_calls = 0} disables the cap.
#'
#' @param tools Named list of ellmer ToolDefs, or NULL when no session
#'   is bound.
#' @param logger Session logger (\code{\link{agent_logger_new}}) or
#'   NULL.
#' @param session_token Bound session token (used only for log tags),
#'   or NULL.
#' @param tier Enabled tool tier ("read" in Stage 2).
#' @param max_calls Per-connection tools/call cap; NULL reads
#'   \code{OMICSVIEWER_AGENT_CHANNEL_MAX_CALLS} (default 600, 0 =
#'   unlimited).
#' @return An environment (handler state).
#' @keywords internal
#' @rdname agentChannelServer
agent_channel_state_new <- function(tools = NULL, logger = NULL,
                                    session_token = NULL, tier = "read",
                                    max_calls = NULL) {
  tier <- match.arg(tier, names(.agent_channel_tiers))
  if (is.null(max_calls))
    max_calls <- .agent_channel_max_calls()
  state <- new.env(parent = emptyenv())
  state$tools <- tools
  state$logger <- logger
  state$session_token <- session_token
  state$tier <- tier
  state$max_calls <- max_calls
  state$calls <- 0L
  state
}

#' JSON-RPC 2.0 response serialization (the wire edge)
#'
#' The single place response objects become text: auto_unbox for
#' scalars, explicit nulls. Length-1 arrays that must stay arrays
#' (tool-schema \code{required}/\code{enum}) are protected upstream by
#' wrapping them in \code{I()} (Stage-0 review obligation) - see
#' \code{\link{.agent_channel_wire_spec}}.
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_res <- function(response) {
  as.character(jsonlite::toJSON(response, auto_unbox = TRUE,
                                null = "null", na = "null"))
}

.agent_channel_ok <- function(id, result)
  .agent_channel_res(list(jsonrpc = "2.0", id = id, result = result))

.agent_channel_err <- function(id, code, message)
  .agent_channel_res(list(jsonrpc = "2.0", id = id,
                          error = list(code = code, message = message)))

#' Parse-error response
#'
#' @param message Human-readable detail appended to the standard
#'   parse-error text.
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_parse_error_response <- function(message = NULL) {
  text <- "Parse error: the message is not valid JSON."
  if (!is.null(message) && length(message) == 1L && !is.na(message) &&
      nzchar(message))
    text <- paste(text, message)
  .agent_channel_err(NULL, -32700L, text)
}

#' Keep length-1 schema arrays as arrays on the wire
#'
#' jsonlite's \code{auto_unbox = TRUE} collapses a length-1 atomic
#' vector to a scalar, but JSON Schema requires \code{required} and
#' \code{enum} to be arrays regardless of length (a single required
#' argument must serialize as \code{["name"]}, not \code{"name"}).
#' This recursive pre-pass wraps exactly those fields in \code{I()} so
#' the edge serializer keeps them as arrays (Stage-0 review
#' obligation); every other field passes through untouched.
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_wire_spec <- function(spec) {
  if (!is.list(spec))
    return(spec)
  out <- lapply(names(spec), function(nm) {
    value <- spec[[nm]]
    if (nm %in% c("required", "enum") && is.atomic(value) &&
        !is.null(value) && length(value) >= 1L)
      return(I(value))
    .agent_channel_wire_spec(value)
  })
  names(out) <- names(spec)
  out
}

.agent_channel_log <- function(state, event, details) {
  if (is.null(state$logger))
    return(invisible(FALSE))
  try(agent_logger_event(state$logger, event, details), silent = TRUE)
  invisible(TRUE)
}

############################################################################
### [Stage 2] the handler itself

#' Answer one JSON-RPC message (pure handler)
#'
#' Parses \code{msg_json} (a single JSON-RPC 2.0 request object) and
#' returns the response as a JSON string - no sockets, no package-env
#' access, never throws. Methods:
#' \itemize{
#'   \item \code{initialize} - protocol identity: name
#'     \code{"omicsviewer-agent-channel"}, jsonrpc \code{"2.0"},
#'     \code{package_version}, \code{channel_schema_version} (D8), and
#'     the enabled \code{tool_tier} ("read" in Stage 2). No bound
#'     session required.
#'   \item \code{tools/list} - the enabled tier's tool specs
#'     (\code{\link{agent_channel_toolspecs}} filtered to the tier,
#'     \code{required}/\code{enum} kept as arrays) plus the version
#'     echo (D8). Requires a bound session.
#'   \item \code{tools/call} (params \code{name}, optional \code{args}
#'     object) - dispatch through \code{\link{agent_channel_dispatch}}
#'     against the bound session's tools. The wire result carries
#'     \code{value@value} ONLY, wrapped as
#'     \code{content = [{type: "text", text}]} (Stage-0 review F3:
#'     \code{extra$data}/\code{extra$display} never cross the wire).
#'     Dispatcher \code{ok = FALSE} maps to error -32602 (arg/tool
#'     errors) and the connection stays open so an agent can
#'     self-correct.
#' }
#' Error codes: -32700 parse error; -32600 invalid request (batch
#' arrays, non-object messages, missing id - notifications are not
#' supported in v1; wrong \code{jsonrpc} field); -32601 unknown method;
#' -32602 tool/argument errors (including tools outside the enabled
#' tier); -32000 channel-state errors (no session bound; rate cap
#' exceeded - the cap is per connection,
#' \code{OMICSVIEWER_AGENT_CHANNEL_MAX_CALLS}, default 600, 0 =
#' unlimited, and an exceeded cap does NOT close the connection);
#' -32603 internal error (never expected).
#'
#' Rate-capped and executed \code{tools/call} requests emit
#' \code{channel_request_start} / \code{channel_request_end} (method,
#' tool, ok, duration) and \code{channel_rate_cap_hit} events through
#' the state's logger when one is present (D5).
#'
#' @param state Handler state from \code{\link{agent_channel_state_new}}
#'   (an environment; \code{state$calls} is incremented once per
#'   executed tools/call).
#' @param msg_json Single character string (one JSON-RPC request
#'   object, text frame).
#' @return JSON-RPC response string.
#' @keywords internal
#' @rdname agentChannelServer
agent_channel_handle <- function(state, msg_json) {
  tryCatch(
    .agent_channel_handle_impl(state, msg_json),
    error = function(e)
      .agent_channel_err(NULL, -32603L,
                         paste0("Internal error: ", conditionMessage(e))))
}

.agent_channel_handle_impl <- function(state, msg_json) {
  if (!is.character(msg_json) || length(msg_json) != 1L || is.na(msg_json))
    return(.agent_channel_parse_error_response(
      "The channel speaks JSON-RPC 2.0 over text frames."))
  parsed <- tryCatch(
    jsonlite::fromJSON(msg_json, simplifyVector = FALSE),
    error = function(e) NULL)
  if (is.null(parsed))
    return(.agent_channel_parse_error_response())
  if (!is.list(parsed))
    return(.agent_channel_err(NULL, -32600L,
      "Invalid Request: the message must be a JSON-RPC 2.0 request object."))
  if (is.null(names(parsed)) || any(!nzchar(names(parsed))))
    return(.agent_channel_err(NULL, -32600L,
      "Invalid Request: batch arrays are not supported; send one request object per text frame."))
  if (!is.null(parsed$jsonrpc) && !identical(parsed$jsonrpc, "2.0"))
    return(.agent_channel_err(NULL, -32600L,
      "Invalid Request: jsonrpc must be exactly \"2.0\"."))
  if (is.null(parsed$id) || length(parsed$id) != 1L ||
      (!is.character(parsed$id) && !is.numeric(parsed$id)))
    return(.agent_channel_err(NULL, -32600L,
      "Invalid Request: a single string or number id is required (notifications are not supported in v1)."))
  id <- parsed$id
  method <- parsed$method
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method))
    return(.agent_channel_err(id, -32600L,
                              "Invalid Request: method must be a single non-empty string."))
  started_at <- proc.time()[["elapsed"]]

  finish <- function(response, ok, tool = NULL) {
    .agent_channel_log(state, "channel_request_end", list(
      method = method, tool = tool, ok = ok,
      duration_ms = round((proc.time()[["elapsed"]] - started_at) * 1000)))
    response
  }
  .agent_channel_log(state, "channel_request_start", list(method = method))

  if (identical(method, "initialize")) {
    result <- list(
      protocol = "omicsviewer-agent-channel",
      jsonrpc = "2.0",
      package_version = .agent_channel_package_version(),
      channel_schema_version = .agent_channel_schema_version,
      tool_tier = state$tier)
    return(finish(.agent_channel_ok(id, result), TRUE))
  }

  if (identical(method, "tools/list")) {
    if (is.null(state$tools))
      return(finish(.agent_channel_err(
        id, -32000L,
        paste("No session bound: open the omicsViewer assistant in a",
              "browser (with OMICSVIEWER_AGENT_CHANNEL enabled) so a",
              "session can pair with the channel.")), FALSE))
    specs <- agent_channel_toolspecs(state$tools)
    tier_names <- .agent_channel_tiers[[state$tier]]
    tier_names <- tier_names[tier_names %in% names(specs)]
    tools_array <- if (length(tier_names))
      unname(lapply(specs[tier_names], .agent_channel_wire_spec))
    else
      list()
    result <- list(
      tools = tools_array,
      package_version = .agent_channel_package_version(),
      channel_schema_version = .agent_channel_schema_version)
    return(finish(.agent_channel_ok(id, result), TRUE))
  }

  if (identical(method, "tools/call")) {
    params <- parsed$params
    if (is.null(params))
      params <- list()
    if (!is.list(params) || is.null(names(params)) ||
        any(!nzchar(names(params))))
      return(finish(.agent_channel_err(
        id, -32602L,
        "Invalid params: tools/call needs a params object with `name` and optional `args`."),
        FALSE))
    name <- params$name
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !nzchar(name))
      return(finish(.agent_channel_err(
        id, -32602L,
        "Invalid params: `name` must be a single non-empty tool name."), FALSE))
    args <- params$args
    if (is.null(args))
      args <- list()
    if (!is.list(args) || is.null(names(args)) ||
        any(!nzchar(names(args))))
      return(finish(.agent_channel_err(
        id, -32602L,
        "Invalid params: `args` must be a JSON object of tool arguments or omitted."),
        FALSE))
    if (is.null(state$tools))
      return(finish(.agent_channel_err(
        id, -32000L,
        paste("No session bound: open the omicsViewer assistant in a",
              "browser (with OMICSVIEWER_AGENT_CHANNEL enabled) so a",
              "session can pair with the channel.")), FALSE, name))
    if (!(name %in% .agent_channel_tiers[[state$tier]]))
      return(finish(.agent_channel_err(
        id, -32602L,
        paste0("Tool '", name, "' is not available on the '", state$tier,
               "' tier. Tier tools: ",
               paste(.agent_channel_tiers[[state$tier]], collapse = ", "),
               ".")), FALSE, name))
    if (state$max_calls > 0L && state$calls >= state$max_calls) {
      .agent_channel_log(state, "channel_rate_cap_hit", list(
        tool = name, calls = state$calls, max_calls = state$max_calls))
      # Documented behavior: the connection stays open; only tools/call
      # keeps answering -32000 (D6).
      return(finish(.agent_channel_err(
        id, -32000L,
        paste0("Rate cap exceeded: this connection has used its ",
               state$max_calls,
               " tool calls (OMICSVIEWER_AGENT_CHANNEL_MAX_CALLS).",
               " The connection stays open; reconnect for a fresh cap.")),
        FALSE, name))
    }
    state$calls <- state$calls + 1L
    result <- agent_channel_dispatch(state$tools, name, args)
    if (!isTRUE(result$ok))
      return(finish(.agent_channel_err(id, -32602L, result$error),
                    FALSE, name))
    # F3 wire projection: value@value ONLY (bounded codec JSON text or
    # plain string); extra$data / extra$display stay in-process.
    text <- result$value@value
    if (!is.character(text) || length(text) != 1L || is.na(text))
      text <- tryCatch(
        as.character(jsonlite::toJSON(text, auto_unbox = TRUE,
                                      null = "null")),
        error = function(e) "")
    return(finish(.agent_channel_ok(
      id, list(content = list(list(type = "text", text = text)))),
      TRUE, name))
  }

  finish(.agent_channel_err(id, -32601L,
    paste0("Method not found: '", method,
           "'. Supported: initialize, tools/list, tools/call.")), FALSE)
}

############################################################################
### [Stage 2] httpuv server (127.0.0.1, token-paired WS path)

#' Channel server state (package env)
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_server <- local({
  env <- new.env(parent = emptyenv())
  env$running <- FALSE
  env$server <- NULL
  env$token <- NULL
  env$port <- NULL
  env$max_calls <- 600L
  env$next_id <- 1L
  env
})

.agent_channel_server_logger <- function() {
  tokens <- agent_channel_registry_list()
  if (!length(tokens))
    return(NULL)
  entry <- tryCatch(get(tokens[[1]], envir = .agent_channel_registry),
                    error = function(e) NULL)
  if (is.null(entry) || is.null(entry$logger))
    NULL
  else
    entry$logger
}

#' Build one connection's handler state from the bound session
#'
#' Reads the (single) bound registry entry at connection open: the
#' first-bound session's tools and logger become this connection's
#' dispatch target (first-bind-wins, plan D4 interpretation). Falls
#' back to a no-session state when nothing is bound.
#' @keywords internal
#' @rdname agentChannelServer
.agent_channel_conn_state <- function(server) {
  tokens <- agent_channel_registry_list()
  entry <- NULL
  if (length(tokens))
    entry <- tryCatch(get(tokens[[1]], envir = .agent_channel_registry),
                      error = function(e) NULL)
  agent_channel_state_new(
    tools = if (!is.null(entry)) entry$tools else NULL,
    logger = if (!is.null(entry) && !is.null(entry$logger))
      entry$logger else .agent_channel_server_logger(),
    session_token = if (length(tokens)) tokens[[1]] else NULL,
    tier = "read",
    max_calls = server$max_calls)
}

.agent_channel_request_path <- function(req) {
  path <- NULL
  if (!is.null(req)) {
    path <- req$PATH_INFO
    if (is.null(path) || length(path) != 1L || is.na(path))
      path <- req$path
  }
  if (is.null(path) || length(path) != 1L || is.na(path))
    "/"
  else
    path
}

#' Start the agent channel WebSocket server
#'
#' Starts (or reports) the JSON-RPC 2.0 endpoint on 127.0.0.1. The
#' host is fixed: the server never binds a non-loopback interface, and
#' a \code{host} argument is deliberately absent - asking for another
#' interface is impossible by construction (D10). Does nothing (and
#' says why, returning \code{FALSE}) unless
#' \code{OMICSVIEWER_AGENT_CHANNEL} is truthy; a second start while
#' running reports the existing server and returns \code{TRUE}.
#'
#' The pairing token is generated fresh at start (>= 32 hex chars), so
#' a restart invalidates prior clients; the FULL token is printed to
#' the console (the pairing surface) while log events carry only
#' \code{\link{.agent_channel_token_tag}}. Connections must hit
#' \code{ws://127.0.0.1:<port>/<token>}; anything else is closed
#' immediately (nothing about the token is echoed back).
#'
#' @param port Port to listen on; NULL reads
#'   \code{OMICSVIEWER_AGENT_CHANNEL_PORT} (default 7790, outside the
#'   7775-7786 e2e range).
#' @return TRUE when a server is listening after the call, FALSE
#'   otherwise (invisibly).
#' @keywords internal
#' @rdname agentChannelServer
agent_channel_server_start <- function(port = NULL) {
  if (!.agent_channel_enabled()) {
    message(paste(
      "omicsViewer agent channel: disabled.",
      "Set OMICSVIEWER_AGENT_CHANNEL=1 (and optionally",
      "OMICSVIEWER_AGENT_CHANNEL_PORT) to expose the assistant tools",
      "on a localhost JSON-RPC channel."))
    return(invisible(FALSE))
  }
  server_env <- .agent_channel_server
  if (isTRUE(server_env$running)) {
    message("omicsViewer agent channel: already running on port ",
            server_env$port, " (token tag ",
            .agent_channel_token_tag(server_env$token), ").")
    return(invisible(TRUE))
  }
  if (is.null(port))
    port <- .agent_channel_port_default()
  port <- suppressWarnings(as.integer(port[1]))
  if (length(port) != 1L || is.na(port) || port < 1L || port > 65535L)
    stop("agent_channel_server_start: port must be an integer in [1, 65535].")
  token <- .agent_channel_new_token()
  server_env$token <- token
  server_env$port <- port
  server_env$max_calls <- .agent_channel_max_calls()

  handle <- NULL
  ok_start <- tryCatch({
    handle <- httpuv::startServer("127.0.0.1", port, list(
      onWSOpen = function(ws) {
        path <- .agent_channel_request_path(ws$request)
        expected <- paste0("/", server_env$token)
        tag <- .agent_channel_token_tag(server_env$token)
        if (!identical(path, expected)) {
          # Wrong path: refuse and say NOTHING that could be used to
          # probe for the right one (D5: no secrets in events).
          .agent_channel_log(
            agent_channel_state_new(logger = .agent_channel_server_logger()),
            "channel_connection_refused",
            list(path = "<redacted>", token_tag = tag))
          try(ws$close(), silent = TRUE)
          return(invisible(NULL))
        }
        conn_id <- server_env$next_id
        server_env$next_id <- conn_id + 1L
        conn_state <- .agent_channel_conn_state(server_env)
        session_tag <- if (is.null(conn_state$session_token)) NULL else
          .agent_channel_token_tag(conn_state$session_token)
        .agent_channel_log(conn_state, "channel_connection_open", list(
          connection = conn_id, token_tag = tag, session = session_tag))
        ws$onMessage(function(binary, message) {
          response <- tryCatch(
            if (isTRUE(binary))
              .agent_channel_parse_error_response(
                "Binary frames are not supported; the channel speaks JSON-RPC 2.0 over text frames.")
            else
              agent_channel_handle(conn_state, message),
            error = function(e)
              .agent_channel_err(NULL, -32603L,
                                 paste0("Internal error: ",
                                        conditionMessage(e))))
          try(ws$send(response), silent = TRUE)
        })
        ws$onClose(function() {
          .agent_channel_log(conn_state, "channel_connection_close",
                             list(connection = conn_id))
        })
      },
      onRequest = function(req) {
        # Plain HTTP requests carry no endpoint in Stage 2 (figures
        # arrive in Stage 3); answer generically without echoing any
        # path content.
        list(status = 404L,
             body = "Not Found",
             headers = list("Content-Type" = "text/plain"))
      }
    ))
    TRUE
  }, error = function(e) {
    message("omicsViewer agent channel: failed to start on port ", port,
            ": ", conditionMessage(e))
    FALSE
  })
  if (!isTRUE(ok_start))
    return(invisible(FALSE))
  server_env$server <- handle
  server_env$running <- TRUE
  # The pairing surface: full token to the console ONLY. Log events
  # below carry the tag, never the token (D5/D10).
  message("omicsViewer agent channel: listening on ws://127.0.0.1:",
          port, "/", token)
  logger <- .agent_channel_server_logger()
  if (!is.null(logger))
    try(agent_logger_event(logger, "channel_server_start", list(
      port = port, token_tag = .agent_channel_token_tag(token),
      tier = "read", max_calls = server_env$max_calls)), silent = TRUE)
  invisible(TRUE)
}

#' Stop the agent channel WebSocket server
#'
#' Idempotent: stopping an already-stopped channel is a silent no-op.
#' The pairing token is discarded (a later start generates a fresh
#' one).
#' @return TRUE if a running server was stopped (invisibly).
#' @keywords internal
#' @rdname agentChannelServer
agent_channel_server_stop <- function() {
  server_env <- .agent_channel_server
  if (!isTRUE(server_env$running))
    return(invisible(FALSE))
  logger <- .agent_channel_server_logger()
  port <- server_env$port
  if (!is.null(server_env$server))
    try(httpuv::stopServer(server_env$server), silent = TRUE)
  server_env$server <- NULL
  server_env$token <- NULL
  server_env$port <- NULL
  server_env$running <- FALSE
  message("omicsViewer agent channel: stopped (was on port ", port, ").")
  if (!is.null(logger))
    try(agent_logger_event(logger, "channel_server_stop",
                           list(port = port)), silent = TRUE)
  invisible(TRUE)
}

#' Introspect the channel server (internal)
#'
#' Test/introspection helper. The full token is included so a local
#' client can pair; it is never on the wire or in logs.
#' @return List: running, port, token, token_tag, tier, max_calls.
#' @keywords internal
#' @rdname agentChannelServer
agent_channel_server_info <- function() {
  server_env <- .agent_channel_server
  running <- isTRUE(server_env$running)
  list(
    running = running,
    port = if (running) server_env$port else NULL,
    token = if (running) server_env$token else NULL,
    token_tag = if (running)
      .agent_channel_token_tag(server_env$token) else NULL,
    tier = "read",
    max_calls = if (running) server_env$max_calls else NULL)
}

############################################################################
### [Stage 2] session binding (D14): the assistant-module hook
###
### One defaulted call from module_aiAssistant.R (right after the tool
### registry is built): bind this session into the Stage-0 registry and
### start the server when the channel is enabled. First-bind-wins (the
### plan's D4 "refuse when >1 session" narrowed for Stage 2): a second
### session logs a console+logger warning and does NOT rebind - the
### channel stays attached to the first session until it ends, and the
### next session to attach after that wins. Everything here is
### failure-tolerant: a channel problem must never break app boot.

#' Bind a session to the agent channel (assistant-module hook)
#'
#' Called once per session from \code{module_aiAssistant.R} with the
#' freshly built tool registry. Inert (returns \code{FALSE}) unless
#' \code{OMICSVIEWER_AGENT_CHANNEL} is truthy. When enabled and no
#' session is bound, stores the session's tools/domain/logger in the
#' Stage-0 registry under \code{session$token}, wires
#' \code{session$onEnded} eviction, and starts the server. When a
#' session is already bound, logs a console+logger warning and does
#' not rebind (first-bind-wins). Never throws - a channel failure is
#' reported to the console and swallowed.
#'
#' @param session Shiny session (needs \code{$token} and
#'   \code{$onEnded}).
#' @param tools Named list of ellmer ToolDefs (the registry's
#'   \code{$tools}).
#' @param session_domain The session domain the tool handlers use.
#' @param logger Session logger (\code{\link{agent_logger_new}}).
#' @return TRUE when this session bound the channel (invisibly).
#' @keywords internal
#' @rdname agentChannelServer
agent_channel_session_attach <- function(session, tools, session_domain,
                                         logger) {
  tryCatch(
    {
      if (!.agent_channel_enabled())
        return(invisible(FALSE))
      existing <- agent_channel_registry_list()
      if (length(existing) >= 1L) {
        tags <- vapply(existing, .agent_channel_token_tag, character(1))
        message(paste(
          "omicsViewer agent channel: a session is already bound (tag",
          paste(tags, collapse = ", "),
          "); this session will NOT be exposed on the channel",
          "(first-bind-wins)."))
        agent_logger_event(logger, "channel_bind_refused",
                           list(bound_sessions = tags))
        return(invisible(FALSE))
      }
      token <- session$token
      agent_channel_registry_bind(
        token,
        list(tools = tools, session_domain = session_domain,
             logger = logger,
             bound_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z",
                               tz = "UTC")),
        on_ended = function()
          session$onEnded(function() agent_channel_session_detach(token)))
      agent_logger_event(
        logger, "channel_bind",
        list(session = .agent_channel_token_tag(token), tier = "read"))
      agent_channel_server_start()
      invisible(TRUE)
    },
    error = function(e) {
      message("omicsViewer agent channel: session attach failed: ",
              conditionMessage(e))
      invisible(FALSE)
    })
}

#' Unbind a session and stop the server when it was the last one
#'
#' Wired through \code{session$onEnded} by
#' \code{\link{agent_channel_session_attach}}; safe to call manually
#' (tests) and idempotent for absent tokens.
#' @param token Session token to unbind.
#' @return NULL (invisibly); never throws.
#' @keywords internal
#' @rdname agentChannelServer
agent_channel_session_detach <- function(token) {
  tryCatch(
    {
      if (!is.character(token) || length(token) != 1L || is.na(token) ||
          !nzchar(token))
        return(invisible(NULL))
      agent_channel_registry_unbind(token)
      if (!length(agent_channel_registry_list()))
        agent_channel_server_stop()
      invisible(NULL)
    },
    error = function(e) {
      message("omicsViewer agent channel: session detach failed: ",
              conditionMessage(e))
      invisible(NULL)
    })
}
