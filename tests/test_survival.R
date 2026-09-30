## Survival module convention tests (todo item 1.1, Option A)
##
## Pins the author-intended data convention: a trailing "+" marks an EVENT
## at that time; values without "+" are censored at the recorded time.
## Administrative censoring (slider) truncates and censors, never creates
## events. Also covers the prepEsetViewer surv validation, which was a
## no-op before (wrong variable + match-everything regex).

library(unittest)

parse_surv <- omicsViewer:::.survival_parse_response
apply_censor <- omicsViewer:::.survival_apply_censor

## ---- convention: trailing '+' marks an event ------------------------------
p <- parse_surv(c("969+", "227", "120+", "45"))
ok(ut_cmp_equal(p$time, c(969, 227, 120, 45)), "parse - times")
ok(ut_cmp_equal(p$event, c(1, 0, 1, 0)), "parse - trailing + marks event")

## numeric-only column: no '+' means no events (all censored)
ok(ut_cmp_equal(parse_surv(c("100", "200"))$event, c(0, 0)),
   "parse - no + means censored")

## fractional times with '+'
ok(ut_cmp_equal(parse_surv(c("45.5+", "3.2"))$event, c(1, 0)),
   "parse - decimal times keep convention")

## ---- demo data: the OS column has 17 '+'-suffixed of 47 entries -----------
svf <- system.file("extdata", "sampleSurv.tsv", package = "omicsViewer")
sv <- read.delim(svf)
p2 <- parse_surv(as.character(sv$OS))
ok(ut_cmp_equal(sum(p2$event), sum(grepl("\\+$", as.character(sv$OS)))),
   "sampleSurv - event count equals '+' entries")
ok(ut_cmp_equal(sum(p2$event), 17L), "sampleSurv - demo OS has 17 events")
ok(!any(is.na(p2$time)), "sampleSurv - all times parse")

## ---- administrative censoring ----------------------------------------------
df <- data.frame(time = c(100, 200, 300), event = c(1, 0, 1),
                 strata = c("a", "b", "a"))
d2 <- apply_censor(df, 200)
ok(ut_cmp_equal(d2$time, c(100, 200, 200)), "censor - truncates time")
ok(ut_cmp_equal(d2$event, c(1, 0, 0)), "censor - beyond-censor becomes event 0")
ok(ut_cmp_equal(apply_censor(df, 300)$event, df$event),
   "censor - censor at max time is a no-op")

## ---- strata gate (regression: was length(df$strata > 1), always TRUE) -----
ok(!isTRUE(length(unique(rep("selected", 3))) > 1),
   "strata gate - single stratum is FALSE")
ok(isTRUE(length(unique(c("a", "b", "a"))) > 1),
   "strata gate - two strata is TRUE")

## ---- prepEsetViewer surv validation ----------------------------------------
expr <- matrix(rnorm(15), 5, 3,
               dimnames = list(paste0("g", 1:5), paste0("s", 1:3)))
pd <- data.frame(row.names = paste0("s", 1:3), group = c("a", "a", "b"))
fd <- data.frame(row.names = paste0("g", 1:5), symbol = paste0("g", 1:5))

res <- try(
  omicsViewer::prepOmicsViewer(expr, pd, fd, surv = c("120+", "45", "969+"),
                              PCA = FALSE, correlation = FALSE),
  silent = TRUE)
ok(!inherits(res, "try-error"), "prep - accepts numeric times with trailing +")

res <- try(
  omicsViewer::prepOmicsViewer(expr, pd, fd, surv = c("120+", "45 months", "969+"),
                              PCA = FALSE, correlation = FALSE),
  silent = TRUE)
ok(inherits(res, "try-error"), "prep - rejects non-numeric surv values")

## NA entries are tolerated (missing survival for a sample)
res <- try(
  omicsViewer::prepOmicsViewer(expr, pd, fd, surv = c("120+", NA, "969+"),
                              PCA = FALSE, correlation = FALSE),
  silent = TRUE)
ok(!inherits(res, "try-error"), "prep - accepts NA survival entries")

# ---------------- R-M7: censor slider survives data changes ----------------
# The slider used to be rebuilt AT MAX on every dat() change (a
# sample-selection change silently reset the user's censor time), and
# store pushes landing before the renderUI slider existed were lost.
# The renderUI now sets the value at render time from the module's
# desired-value reactiveVal.
sv_resp_r7 <- shiny::reactiveVal(c("100", "200+", "300", "150+"))
sv_ui_r7 <- NULL
shiny::testServer(function(input, output, session) {
  omicsViewer:::survival_module("sv", reactive_resp = sv_resp_r7,
    reactive_strata = reactive(NULL), reactive_checkpoint = reactive(TRUE))
  shiny::outputOptions(output, "sv-censor_output", suspendWhenHidden = FALSE)
}, {
  for (i in 1:5) session$flushReact()
  session$setInputs(`sv-censor` = 300)   # warm ignoreInit observer (slider init)
  for (i in 1:3) session$flushReact()
  session$setInputs(`sv-censor` = 150)   # user picks a censor time
  for (i in 1:3) session$flushReact()
  sv_resp_r7(c("100", "200+", "300", "400"))  # sample-selection change
  for (i in 1:4) session$flushReact()
  sv_ui_r7 <<- tryCatch(paste(output[["sv-censor_output"]], collapse = ""),
                        error = function(e) conditionMessage(e))
})
ok(grepl('data-from="150"', sv_ui_r7, fixed = TRUE) &&
     grepl('data-max="400"', sv_ui_r7, fixed = TRUE),
   "survival: censor slider keeps the user value across data changes (R-M7)")

# R-L4: a 2-column gs data.frame (featureId, gsId; no weight) is
# documented input - as.integer(NULL) used to error "replacement has 0
# rows". The fix defaults weight to 1.
gsdf <- data.frame(featureId = c(1, 2), gsId = c("gs1", "gs1"))
res <- try(
  omicsViewer::prepOmicsViewer(expr, pd, fd, gs = gsdf,
                               PCA = FALSE, correlation = FALSE,
                               SummarizedExperiment = FALSE),
  silent = TRUE)
ok(!inherits(res, "try-error"), "prep - 2-column gs accepted (R-L4)")
ok(!inherits(res, "try-error") &&
     identical(as.integer(attr(Biobase::fData(res), "GS")$weight), c(1L, 1L)),
   "prep - missing gs weight defaults to 1 (R-L4)")
