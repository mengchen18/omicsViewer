# WP0 Tier A - UI-effect regression for the agent state bridge.
#
# Runner: spawns the Playwright harness (tests/e2e_agent/tier_a.mjs) against
# a locally spawned omicsViewer app with demo.RDS preloaded, then reports the
# harness assertions through unittest. No LLM provider is required.
#
# Run:  Rscript tests/test_agentUiEffects.R
# Needs: node + playwright module (tests/e2e_agent), system Chrome.

library(omicsViewer)
library(unittest, quietly = TRUE)

e2e_dir <- file.path("tests", "e2e_agent")
if (!dir.exists(e2e_dir))
  e2e_dir <- "e2e_agent" # allow running from tests/ directly

skip <- function(reason) {
  ok(TRUE, paste("agent UI-effect tests skipped:", reason))
  quit(save = "no", status = 0)
}

if (!nzchar(Sys.which("node")))
  skip("node unavailable")
if (!file.exists(file.path(e2e_dir, "node_modules", "playwright")))
  skip("playwright module unavailable (run npm install in tests/e2e_agent)")
chrome <- Sys.which(c("google-chrome", "google-chrome-stable", "chromium", "chromium-browser"))
if (!any(nzchar(chrome)))
  skip("no Chrome/Chromium found")

# Stale app processes from crashed runs cause port conflicts; the harness
# kills its own child on exit, but sweep the port defensively first.
port_pids <- tryCatch(
  {
    out <- system2("ss", c("-tlnp"), stdout = TRUE, stderr = TRUE)
    pids <- regmatches(out, gregexpr("pid=\\K[0-9]+", out, perl = TRUE))
    unique(unlist(pids[grepl(":7778", out)]))
  },
  error = function(e) character()
)
if (length(port_pids)) {
  for (pid in port_pids) try(system2("kill", c("-9", pid)), silent = TRUE)
  Sys.sleep(1)
}

status <- suppressWarnings(system2(
  "node", c(file.path(e2e_dir, "tier_a.mjs")),
  stdout = TRUE, stderr = TRUE
))

results_file <- file.path(e2e_dir, "tier_a_results.json")
if (file.exists(results_file)) {
  results <- jsonlite::fromJSON(results_file, simplifyVector = FALSE)
  for (entry in results$results) {
    ok(isTRUE(entry$pass), paste0("tier-a: ", entry$name,
                                  if (!isTRUE(entry$pass) && nzchar(entry$detail)) paste0(" [", entry$detail, "]")))
  }
  ok(isTRUE(results$passed == results$total),
     sprintf("tier-a suite complete (%s/%s)", results$passed, results$total))
} else {
  ok(FALSE, "tier-a harness produced no results file")
  cat(paste(status, collapse = "\n"), "\n")
}

if (!isTRUE(attr(status, "status") == 0) && file.exists(results_file)) {
  cat(paste(status, collapse = "\n"), "\n")
}
