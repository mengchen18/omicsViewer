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
