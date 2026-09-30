library(omicsViewer)
library(unittest, quietly = TRUE)

## todo 1.4: a warning emitted while loading a dataset must NOT abort the
## load. The old tryCatch(warning = ...) handler swallowed the warning and
## returned NULL, so any coercion/reshape warning surfaced as
## "Failed to load data - file may be corrupted".

demo_path <- system.file("extdata", "demo.RDS", package = "omicsViewer")
ok(ut_cmp_identical(nzchar(demo_path), TRUE), "demo dataset is available")

ldir <- file.path(tempdir(), paste0("loader-warn-", Sys.getpid()))
if (dir.exists(ldir)) unlink(ldir, recursive = TRUE)
dir.create(ldir, recursive = TRUE, showWarnings = FALSE)
file.copy(demo_path, file.path(ldir, "warn.RDS"))
file.copy(demo_path, file.path(ldir, "err.RDS"))

app_loader <- function(input, output, session) {
  omicsViewer:::app_module(
    "app",
    .dir = shiny::reactive(ldir),
    esetLoader = function(path) {
      if (grepl("warn", path, fixed = TRUE))
        warning("test coercion warning")
      if (grepl("err", path, fixed = TRUE))
        stop("test load error")
      readRDS(path)
    }
  )
}

shiny::testServer(app_loader, {
  # the data-space navbar input must exist before the first flush, otherwise
  # the status assembly errors under the mock session (see test_appState.R)
  session$setInputs(`app-dataspace-eset` = "Feature")

  loading_text <- function() {
    session$flushReact()
    session$flushOutput()
    as.character(output[["app-loadingStatus"]])
  }

  # warnings must pass through without aborting the load
  session$setInputs(`app-selectFile` = "warn.RDS")
  ok(
    ut_cmp_identical(
      any(grepl("Dataset loaded successfully", loading_text(), fixed = TRUE)),
      TRUE
    ),
    "a warning during loading does not abort the load"
  )

  # genuine errors still fail cleanly (NULL dataset, session alive)
  session$setInputs(`app-selectFile` = "err.RDS")
  ok(
    ut_cmp_identical(
      any(grepl("Loading dataset", loading_text(), fixed = TRUE)),
      TRUE
    ),
    "an error during loading keeps the dataset empty"
  )
  ok(
    ut_cmp_identical(isFALSE(session$isClosed()), TRUE),
    "the session survives both loading outcomes"
  )
})
