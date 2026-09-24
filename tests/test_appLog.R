library(omicsViewer)
library(unittest, quietly = TRUE)

# tests for the omicsViewer() console log (R/auxi_appLog.R)
applog_dir <- omicsViewer:::applog_dir
applog_new_file <- omicsViewer:::applog_new_file
applog_begin <- omicsViewer:::applog_begin
applog_end <- omicsViewer:::applog_end
applog_write <- omicsViewer:::applog_write

ok(sink.number() == 0L, "clean sink stack at start")

# ---- 1. directory resolution ------------------------------------------
Sys.setenv(OMICSVIEWER_LOG_DIR = "")
ok(grepl("omicsviewer-logs$", applog_dir()),
   "default log dir is tempdir()/omicsviewer-logs")
Sys.setenv(OMICSVIEWER_LOG_DIR = file.path(tempdir(), "custom-logs"))
ok(grepl("custom-logs$", applog_dir()), "OMICSVIEWER_LOG_DIR overrides default")
Sys.setenv(OMICSVIEWER_LOG_DIR = "")

# ---- 2. per-run file naming -------------------------------------------
d <- file.path(tempdir(), "applog-unit")
unlink(d, recursive = TRUE)
f <- applog_new_file(d)
ok(dir.exists(d), "log directory is created")
ok(grepl("^app-\\d{8}-\\d{6}-p\\d+\\.log$", basename(f)),
   "file name carries timestamp and pid")

# ---- 3. round trip: output sink + condition handlers ------------------
tf <- tempfile(fileext = ".log")
env <- applog_begin(tf)
ok(!is.null(env), "applog_begin returns a logger")
ok(identical(env$path, tf), "explicit path is honoured")

cat("cat-output-line\n")
print("print-output-line")
withCallingHandlers({
  message("a-message")
  warning("a-warning")
  # uncaught propagation: the handler must wrap the stop() closer than
  # any tryCatch, exactly like the runApp() wrapper in omicsViewer()
  tryCatch(
    withCallingHandlers(stop("an-error"),
      error = function(c) applog_write(env, "ERROR", c)),
    error = function(e) NULL
  )
},
  message = function(c) applog_write(env, "MESSAGE", c),
  warning = function(c) applog_write(env, "WARNING", c)
)
applog_end(env)

txt <- readLines(tf)
ok(any(grepl("cat-output-line", txt)), "cat output captured")
ok(any(grepl("print-output-line", txt)), "print output captured")
ok(any(grepl("MESSAGE.*a-message", txt)), "message condition captured with level")
ok(any(grepl("WARNING.*a-warning", txt)), "warning condition captured with level")
ok(any(grepl("ERROR.*an-error", txt)), "error condition captured with level")
ok(any(grepl("^== omicsViewer console log ==", txt)), "header written")
ok(any(grepl("^== ended", txt)), "footer written")
ok(sink.number() == 0L, "output sink fully restored")

# timestamped entries carry an ISO-looking prefix
i <- grep("MESSAGE.*a-message", txt)
ok(grepl("^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}", txt[i]),
   "condition entries are timestamped")

# ---- 4. warning carries the call site ---------------------------------
tf2 <- tempfile(fileext = ".log")
env2 <- applog_begin(tf2)
withCallingHandlers(
  warning("with-call"),
  warning = function(c) applog_write(env2, "WARNING", c)
)
applog_end(env2)
ok(any(grepl("with-call.*\\[call:", readLines(tf2))),
   "condition call site is appended")

# ---- 5. nested sinks are unwound safely ------------------------------
tf3 <- tempfile(fileext = ".log")
env3 <- applog_begin(tf3)
sink(tempfile())  # stray inner sink, e.g. from app code
applog_end(env3)
ok(sink.number() == 0L, "nested sink unwound by restore")
# clean up any sink left over by the stray one
while (sink.number() > 0L) sink()

# ---- 6. opt-out paths --------------------------------------------------
ok(is.null(applog_begin(FALSE)), "log_file = FALSE disables logging")
Sys.setenv(OMICSVIEWER_LOG = "off")
ok(is.null(applog_begin(NULL)), "OMICSVIEWER_LOG=off disables logging")
Sys.setenv(OMICSVIEWER_LOG = "")
dflt <- applog_begin(NULL)
ok(!is.null(dflt), "default-on when env unset")
ok(identical(dirname(dflt$path), applog_dir()), "default path in default dir")
ok(file.exists(dflt$path), "default log file reserved on disk")
applog_end(dflt)

# applog_end / applog_write tolerate NULL
ok(!applog_end(NULL), "applog_end(NULL) is a no-op")
ok(!applog_write(NULL, "MESSAGE", simpleMessage("x")),
   "applog_write(NULL, ..) is a no-op")

ok(sink.number() == 0L, "clean sink stack at end")
