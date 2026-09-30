library(omicsViewer)
library(unittest, quietly = TRUE)

# todo 4.3 (R-H3/R-H6/M7): performance-gate regressions

# ---------------- R-H6: adist size guard ----------------
xm <- matrix(rnorm(120 * 40), nrow = 120)
ok(inherits(tryCatch(omicsViewer:::adist(xm, "pearson"), error = function(e) e), "dist"),
   "adist: within the cap returns a dist object (R-H6)")
big <- matrix(rnorm((omicsViewer:::ADIST_MAX_ROWS + 1) * 4),
              nrow = omicsViewer:::ADIST_MAX_ROWS + 1)
err <- tryCatch(omicsViewer:::adist(big, "euclidean"), error = conditionMessage)
ok(grepl("limited to", err, fixed = TRUE),
   "adist: refuses above ADIST_MAX_ROWS with an actionable message (R-H6)")

# ---------------- output_visible helper ----------------
ok(ut_cmp_identical(
  omicsViewer:::output_visible(session = NULL, "any"), TRUE),
   "output_visible: no session (non-reactive) is visible (4.3)")

# ---------------- M7/R-M11: boxplot background cache ----------------
bx <- matrix(rnorm(250 * 12), nrow = 250,
             dimnames = list(paste0("F", 1:250), paste0("S", 1:12)))
bgq <- apply(bx, 2, quantile, probs = seq(0, 1, by = 0.02), na.rm = TRUE)
fig <- tryCatch(omicsViewer:::plotly_boxplot(bx, i = 1:5, bg_quantiles = bgq),
                error = function(e) e)
ok(!inherits(fig, "error") && inherits(fig, "plotly"),
   "plotly_boxplot: pre-quantiled background renders (R-M11)")
fig2 <- tryCatch(omicsViewer:::plotly_boxplot(bx, i = 1:5), error = function(e) e)
ok(!inherits(fig2, "error") && inherits(fig2, "plotly"),
   "plotly_boxplot: no background (uncached path) still renders (control)")
