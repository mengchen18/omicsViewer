# Regression tests for the triselector cascade derivation.
#
# History: the S2 control-plane migration made cascaded choice sets derive
# unconditionally from reactive_selector1/2 (the canonical store values).
# Modules that drive triselectors WITHOUT a store (feature_general, fgsea,
# geneshot, tables, attr4) pass restore-only reactive_selectors backed by an
# empty reactiveVal - those cascades never populated, and the analysis panel
# stayed blank when features were selected. The cascade must work in BOTH
# regimes: store-less (inputs drive; selector1/2 NULL until a restore) and
# store-driven (requested state wins over in-flight inputs).
#
# Server-initiated updates are captured by spying on the mock session's
# sendInputMessage, which is exactly the channel the regression broke.

library(shiny)
library(omicsViewer)
library(unittest, quietly = TRUE)

ts_matrix <- matrix(c(
  "General", "All", "x",
  "General", "All", "y",
  "PCA",     "All", "PC1",
  "PCA",     "All", "PC2"
), ncol = 3, byrow = TRUE)

.spy_input_messages <- function(session, sink) {
  # sink must be an environment: list element assignment inside the spy
  # closure would modify a local copy
  orig <- session$sendInputMessage
  session$sendInputMessage <- function(id, msg) {
    sink$msgs[[length(sink$msgs) + 1L]] <- list(id = id, msg = msg)
    orig(id, msg)
  }
}

.sent_choices <- function(sink, id) {
  # updateSelectInput messages carry an options HTML string, not a choices
  # vector; extract the encoded choice values
  hits <- Filter(function(m) grepl(paste0("(^|\\.)", id, "$"), m$id), sink$msgs)
  if (!length(hits)) return(NULL)
  opts <- utils::tail(hits, 1)[[1]]$msg$options
  if (is.null(opts)) return(NULL)
  m <- regmatches(opts, gregexpr('value="[^"]*"', opts))[[1]]
  sub('^value="', "", sub('"$', "", m))
}

.sent_msg <- function(sink, id) {
  hits <- Filter(function(m) grepl(paste0("(^|\\.)", id, "$"), m$id), sink$msgs)
  if (!length(hits)) return(NULL)
  utils::tail(hits, 1)[[1]]$msg
}
.sent_value <- function(sink, id) {
  msg <- .sent_msg(sink, id)
  if (is.null(msg)) return(NULL)
  msg$value %||% msg$selected
}

## ---------------------------------------------------------- store-less ----
# feature_general / fgsea pattern: reactive_selectors read an empty
# reactiveVal until a status restore writes it.
sent1 <- new.env(); sent1$msgs <- list()
app1 <- function(input, output, session) {
  xax <- reactiveVal()
  res <- omicsViewer:::triselector_module(
    "t",
    reactive_x = reactive(ts_matrix),
    reactive_selector1 = reactive(xax()$v1),
    reactive_selector2 = reactive(xax()$v2),
    reactive_selector3 = reactive(xax()$v3)
  )
  exported <<- list(result = res, xax = xax)
}

testServer(app1, {
  .spy_input_messages(session, sent1)
  # simulate the user choosing the analysis: the subset cascade must react
  session$setInputs(`t-analysis` = "General")
  session$flushReact()
  ok(
    ut_cmp_identical(
      .sent_choices(sent1, "subset"), "All"),
    "store-less: subset choices populate after the user picks an analysis"
  )
  session$setInputs(`t-subset` = "All")
  session$flushReact()
  # server-side selectize messages carry no options HTML; the regression
  # signature was the variable observer never firing at all, so assert the
  # update message exists
  has_variable_update <- any(vapply(sent1$msgs, function(m)
    grepl("(^|\\.)variable$", m$id), logical(1)))
  ok(
    ut_cmp_identical(has_variable_update, TRUE),
    "store-less: variable cascade fires after the user picks a subset"
  )
  session$setInputs(`t-variable` = "y")
  session$flushReact()
  ok(
    ut_cmp_identical(exported$result()$variable, "y"),
    "store-less: manual selections are reported by the module"
  )
})

## ---------------------------------------------------------- store-driven ----
# meta_scatter pattern: selector values come from live reactiveVals and a
# version bump (reactive_axis_request) re-asserts a requested triple.
sent2 <- new.env(); sent2$msgs <- list()
app2 <- function(input, output, session) {
  sel1 <- reactiveVal("General"); sel2 <- reactiveVal("All"); sel3 <- reactiveVal("x")
  ver <- reactiveVal(0L)
  res <- omicsViewer:::triselector_module(
    "t",
    reactive_x = reactive(ts_matrix),
    reactive_selector1 = reactive(sel1()),
    reactive_selector2 = reactive(sel2()),
    reactive_selector3 = reactive(sel3()),
    reactive_axis_request = reactive(ver())
  )
  exported <<- list(
    result = res,
    request = function(v1, v2, v3) {
      sel1(v1); sel2(v2); sel3(v3); ver(isolate(ver()) + 1L)
    }
  )
}

testServer(app2, {
  .spy_input_messages(session, sent2)
  session$setInputs(`t-analysis` = "General", `t-subset` = "All", `t-variable` = "x")
  session$flushReact()
  # request a different triple (as a store transaction would): the choice
  # sets and selections must follow the REQUESTED state, not the in-flight
  # input values
  exported$request("PCA", "All", "PC2")
  session$flushReact()
  session$flushReact()
  ok(
    ut_cmp_identical(.sent_value(sent2, "analysis"), "PCA"),
    "store-driven: analysis update carries the requested selection"
  )
  ok(
    ut_cmp_identical(.sent_choices(sent2, "subset"), "All") &&
      ut_cmp_identical(.sent_value(sent2, "subset"), "All"),
    "store-driven: subset update carries requested choices and selection"
  )
  ok(
    ut_cmp_identical(.sent_value(sent2, "variable"), "PC2"),
    "store-driven: variable update carries the requested selection"
  )
})
