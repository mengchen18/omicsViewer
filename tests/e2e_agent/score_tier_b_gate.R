# WP7 gate scorer — machine-score the full Tier B benchmark runs.
#
# Reads artifacts/<label>/run<K>/task*/record.json (+ log.jsonl) produced by
# tier_b_full.mjs and prints a per-task / per-run score table with the gate
# metric: tasks 10-11 FIRST-ATTEMPT success (plan §4 WP7, decision 5).
#
# Definitions:
#   * prompt span = log events after a user_message up to the next one.
#   * probe = the LAST prompt of a task (setup prompts exist for 8/10/11).
#   * first-attempt success = the probe's outcome was achieved AND no
#     tool_result error occurred anywhere in the probe span (a retry that
#     eventually succeeds downgrades the task to outcome-ok / not-first).
#
# Usage: Rscript tests/e2e_agent/score_tier_b_gate.R [label]     (default: gate)

args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args)) args[1] else "gate"
root <- file.path("tests", "e2e_agent", "artifacts", label)
stopifnot(dir.exists(root))

runs <- list.dirs(root, recursive = FALSE)
runs <- runs[grepl("run[0-9]+$", basename(runs))]
stopifnot(length(runs) > 0)

read_events <- function(log_path) {
  if (!file.exists(log_path)) return(list())
  lines <- readLines(log_path, warn = FALSE)
  events <- list()
  for (l in lines) {
    if (!nzchar(trimws(l))) next
    ev <- tryCatch(jsonlite::fromJSON(l, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(ev)) events[[length(events) + 1L]] <- ev
  }
  events
}

# split events into spans keyed by user_message; return list of spans,
# each = list(prompt, events, tool_calls) where tool_calls carry args+error
span_events <- function(events) {
  spans <- list()
  cur <- NULL
  flush <- function() {
    if (!is.null(cur) && length(cur$events)) spans[[length(spans) + 1L]] <<- cur
  }
  for (ev in events) {
    d <- ev$details
    if (identical(ev$event, "user_message")) {
      flush()
      cur <- list(prompt = d$message %||% "",
                  events = list())
    } else if (!is.null(cur)) {
      cur$events[[length(cur$events) + 1L]] <- ev
    }
  }
  flush()
  spans
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# first figure-operation tool call within a span's events; returns
# list(name, args_json, error) of the FIRST update/create call and flags
# for any tool errors in the span
probe_figure_call <- function(spans) {
  if (!length(spans)) return(NULL)
  sp <- spans[[length(spans)]] # probe = last prompt
  first_fig <- NULL
  any_tool_error <- FALSE
  fig_ops_after_error <- FALSE
  for (ev in sp$events) {
    d <- ev$details
    if (identical(ev$event, "tool_result")) {
      if (!is.null(d$error)) any_tool_error <- TRUE
      next
    }
    if (identical(ev$event, "tool_request") &&
        !is.null(d$tool_name) &&
        d$tool_name %in% c("update_figure", "create_figure") && is.null(first_fig)) {
      first_fig <- d
      if (any_tool_error) fig_ops_after_error <- TRUE
    }
  }
  list(span = sp, first_fig = first_fig,
       any_tool_error = any_tool_error, fig_after_error = fig_ops_after_error)
}

flat_json <- function(x) {
  if (is.null(x)) return("")
  paste(utils::capture.output(str(x, max.level = 6)), collapse = " ")
}

score_task <- function(task_id, record, events) {
  spans <- span_events(events)
  probe <- probe_figure_call(spans)
  obs <- if (length(record$observations)) record$observations[[length(record$observations)]] else list()
  err_in_probe <- !is.null(probe) && probe$any_tool_error

  first <- FALSE # first-attempt success
  outcome <- FALSE
  note <- ""

  switch(as.character(task_id),
    "1" = {
      txt <- obs$chat_tail %||% ""
      outcome <- grepl("2702|2,702|NCI|60 samples|60, |proteom", txt, ignore.case = TRUE)
      first <- outcome && !err_in_probe
      note <- "answer names the dataset"
    },
    "2" = {
      txt <- obs$chat_tail %||% ""
      sel <- record$baseline$feature_selected %||% 0
      pat <- if (sel == 0) "no genes|none are|0 genes|zero genes|no features" else as.character(sel)
      outcome <- grepl(pat, txt, ignore.case = TRUE)
      first <- outcome && !err_in_probe
      note <- sprintf("answer states selection count (truth: %s)", sel)
    },
    "3" = {
      outcome <- identical(obs$tab, "Sample")
      first <- outcome && !err_in_probe
      note <- sprintf("data-space tab = %s", obs$tab %||% "NULL")
    },
    "4" = {
      ax <- obs$axis_titles %||% c("", "")
      outcome <- any(grepl("log\\.fdr|log\\.pvalue", ax, ignore.case = TRUE)) &&
        any(grepl("mean\\.diff", ax, ignore.case = TRUE))
      first <- outcome && !err_in_probe
      note <- paste(ax, collapse = " / ")
    },
    "5" = {
      ax <- obs$axis_titles %||% c("", "")
      outcome <- any(grepl("RE_vs_ME\\|mean\\.diff|mean\\.diff.*RE_vs_ME", ax)) &&
        any(grepl("log\\.fdr", ax, ignore.case = TRUE))
      first <- outcome && !err_in_probe
      note <- paste(ax, collapse = " / ")
    },
    "106" = , # dataset-faithful variant of 6 (see tier_b_full.mjs)
    "6" = {
      outcome <- identical(obs$feature_selected, 5L) || identical(obs$feature_selected, 5)
      first <- outcome && !err_in_probe
      note <- sprintf("selected features = %s", obs$feature_selected %||% "NULL")
    },
    "7" = {
      ok_summary <- FALSE
      if (!is.null(probe)) {
        for (ev in probe$span$events) {
          if (identical(ev$event, "tool_result") && identical(ev$details$tool_name, "summarize_annotation") &&
              is.null(ev$details$error)) ok_summary <- TRUE
        }
      }
      txt <- obs$chat_tail %||% ""
      outcome <- ok_summary || grepl("Origin|group|TP53|cell line", txt, ignore.case = TRUE)
      first <- ok_summary && !err_in_probe
      note <- if (ok_summary) "summarize_annotation ok" else "text-only answer"
    },
    "108" = , # dataset-faithful variant of 8
    "8" = {
      outcome <- (obs$figures %||% 0) >= 1
      first <- outcome && !err_in_probe
      note <- sprintf("figures = %s", obs$figures %||% 0)
    },
    "9" = {
      outcome <- (obs$figures %||% 0) >= 1
      first <- outcome && !err_in_probe
      note <- sprintf("figures = %s", obs$figures %||% 0)
    },
    "10" = {
      fig <- probe$first_fig
      args_txt <- if (!is.null(fig)) flat_json(fig$arguments) else ""
      title_ok <- grepl("expression-wide t-test significance", args_txt, ignore.case = TRUE)
      outcome <- (obs$figures %||% 0) >= 2 && !is.null(fig) && title_ok
      first <- outcome && !err_in_probe
      note <- sprintf("fig tool = %s, title in args = %s, errors in probe = %s",
                      if (is.null(fig)) "none" else fig$tool_name, title_ok, err_in_probe)
    },
    "11" = {
      fig <- probe$first_fig
      args_txt <- if (!is.null(fig)) flat_json(fig$arguments) else ""
      color_ok <- grepl("Intensity", args_txt)
      outcome <- (obs$figures %||% 0) >= 2 && !is.null(fig) && color_ok
      first <- outcome && !err_in_probe
      note <- sprintf("fig tool = %s, Intensity color in args = %s, errors in probe = %s",
                      if (is.null(fig)) "none" else fig$tool_name, color_ok, err_in_probe)
    },
    "12" = {
      txt <- obs$chat_tail %||% ""
      outcome <- grepl("scatter|heatmap|table|enrich|figure|tab", txt, ignore.case = TRUE)
      first <- outcome && !err_in_probe
      note <- "capability disclosure text"
    },
    {
      outcome <- FALSE; first <- FALSE
      note <- "unknown task id (skipped)"
    }
  )

  list(task_id = task_id, status = record$status, outcome = outcome,
       first = first, note = note)
}

rows <- list()
for (rd in runs) {
  run <- basename(rd)
  tdirs <- list.dirs(rd, recursive = FALSE)
  tdirs <- tdirs[grepl("^task[0-9]+$", basename(tdirs))]
  for (td in sort(tdirs)) {
    tid <- as.integer(sub("task", "", basename(td)))
    rec <- tryCatch(jsonlite::fromJSON(file.path(td, "record.json"), simplifyVector = FALSE),
                    error = function(e) NULL)
    if (is.null(rec)) { rows[[length(rows)+1]] <- list(run=run, task_id=tid, status="missing", outcome=FALSE, first=FALSE, note="no record.json"); next }
    events <- read_events(file.path(td, "log.jsonl"))
    s <- score_task(tid, rec, events)
    s$run <- run
    rows[[length(rows) + 1L]] <- s
  }
}

df <- do.call(rbind, lapply(rows, function(r)
  data.frame(run = r$run, task = r$task_id, status = r$status,
             outcome = r$outcome, first_attempt = r$first, note = r$note,
             stringsAsFactors = FALSE)))

cat("\n=== WP7 gate — Tier B full benchmark score (label:", label, ") ===\n\n")
print(df, row.names = FALSE)

cat("\n--- per-run outcome / first-attempt totals ---\n")
agg <- aggregate(cbind(outcome, first_attempt) ~ run, data = df, FUN = sum)
agg$n <- as.integer(table(df$run)[agg$run])
print(agg, row.names = FALSE)

cat("\n--- GATE METRIC: tasks 10-11 first-attempt success ---\n")
gate <- df[df$task %in% c(10, 11), ]
for (t in c(10, 11)) {
  g <- gate[gate$task == t, ]
  n_first <- sum(g$first_attempt)
  n_out <- sum(g$outcome)
  cat(sprintf("task %d: first-attempt %d/%d, outcome-ok %d/%d  -> %s\n",
              t, n_first, nrow(g), n_out, nrow(g),
              if (n_first / max(nrow(g), 1) >= 2/3) "GATE: no patch mode needed" else "GATE: below 2/3"))
}
n_first_all <- sum(gate$first_attempt)
cat(sprintf("tasks 10+11 combined first-attempt: %d/%d (%.0f%%)\n",
            n_first_all, nrow(gate), 100 * n_first_all / max(nrow(gate), 1)))
if (nrow(gate) < 6) {
  cat("\n(incomplete gate data: need 2 tasks x 3 runs; no verdict yet)\n")
} else {
cat(sprintf("\nDecision rule (plan decision 5): build WP7 patch-mode update_figure iff\ncombined first-attempt success < 2/3. Verdict: %s\n",
            if (n_first_all / max(nrow(gate), 1) < 2/3) "BUILD patch mode" else "SKIP patch mode, proceed to WP8"))
}

out <- file.path(root, "gate_scores.csv")
utils::write.csv(df, out, row.names = FALSE)
cat("scores written to", out, "\n")
