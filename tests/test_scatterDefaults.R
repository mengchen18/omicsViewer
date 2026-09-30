## Scatter default-axis seeding tests (todo follow-up: demo.RDS truncated sx/sy)
##
## 1. .scatter_resolve_axis: exact match, truncated/stale-label fallback, NULL cases
## 2. demo.RDS: sx/sy attributes resolve exactly against the live sample-space triset
## 3. demoSurv.RDS: mock survival endpoints exist and keep the '+' convention
## 4. testServer probe: booting the sample-space scatter with the OLD truncated
##    defaults ("PCA|All|PC1(") seeds the live PC1/PC2 axes, not a dead triple

local({
  library(unittest, quietly = TRUE)
  pkgload::load_all(".", quiet = TRUE, attach_testthat = FALSE)
})
dat <- readRDS("inst/extdata/demo.RDS")
pd <- pData(dat); ex <- exprs(dat)

## ---- 1. resolve helper ------------------------------------------------------
ts <- matrix(c("PCA", "All", "PC1(10.5%)",
               "PCA", "All", "PC2(7.2%)",
               "PCA", "All", "PC3(5.4%)",
               "General", "All", "Cell.line"),
             ncol = 3, byrow = TRUE)
ra <- omicsViewer:::.scatter_resolve_axis

ok(ut_cmp_equal(paste(ra("PCA|All|PC2(7.2%)", ts), collapse = "|"),
                "PCA|All|PC2(7.2%)"),
   "resolve - exact match passes through")
ok(ut_cmp_equal(paste(ra("PCA|All|PC1(", ts), collapse = "|"),
                "PCA|All|PC1(10.5%)"),
   "resolve - truncated label falls back to live PC1")
ok(ut_cmp_equal(paste(ra("PCA|All|PC1(21.0%)", ts), collapse = "|"),
                "PCA|All|PC1(10.5%)"),
   "resolve - stale percentage falls back to live PC1")
ok(is.null(ra("Nope|All|PC1", ts)), "resolve - unknown analysis is NULL")
ok(is.null(ra("PCA|Nope|PC1", ts)), "resolve - unknown subset is NULL")
ok(is.null(ra("PCA|All|PC9", ts)), "resolve - unmatched stem is NULL")
ok(is.null(ra("garbage", ts)), "resolve - malformed string is NULL")
ok(is.null(ra(NULL, ts)), "resolve - NULL string is NULL")

## ---- 2. demo.RDS default axes are exact live columns ------------------------
ts_live <- omicsViewer:::trisetter(expr = ex, meta = pd, combine = "pheno")
ok(ut_cmp_equal(attr(dat, "sx"), "PCA|All|PC1(10.5%)"), "demo - sx is the full PC1 label")
ok(ut_cmp_equal(attr(dat, "sy"), "PCA|All|PC2(7.2%)"), "demo - sy is the full PC2 label")
var_of <- function(s) strsplit(s, "\\|")[[1]][3]
ok(1L == sum(ts_live[, 1] == "PCA" & ts_live[, 2] == "All" &
            ts_live[, 3] == var_of(attr(dat, "sx"))),
   "demo - sx matches exactly one triset variable")
ok(1L == sum(ts_live[, 1] == "PCA" & ts_live[, 2] == "All" &
            ts_live[, 3] == var_of(attr(dat, "sy"))),
   "demo - sy matches exactly one triset variable")

## ---- 3. demoSurv.RDS mock ----------------------------------------------------
dsv <- readRDS("inst/extdata/demoSurv.RDS")
pdv <- pData(dsv)
ok(all(c("Surv|all|OS", "Surv|all|DFS") %in% colnames(pdv)),
   "demoSurv - Surv endpoints present")
ok(is.character(pdv[["Surv|all|OS"]]), "demoSurv - OS stored as character (keeps '+')")
ok(ut_cmp_equal(sum(!is.na(pdv[["Surv|all|OS"]])), 47L), "demoSurv - 47 samples with OS")
ok(ut_cmp_equal(sum(grepl("\\+$", pdv[["Surv|all|OS"]][!is.na(pdv[["Surv|all|OS"]])])), 17L),
   "demoSurv - 17 OS events marked '+'")
ok(ut_cmp_equal(sum(grepl("\\+$", pdv[["Surv|all|DFS"]][!is.na(pdv[["Surv|all|DFS"]])])), 28L),
   "demoSurv - 28 DFS events marked '+'")
ok(ut_cmp_equal(attr(dsv, "sx"), "PCA|All|PC1(10.5%)"), "demoSurv - sx default intact")

## ---- 4. testServer probe: seeding resolves the truncated default -------------
## The seed observer depends only on reactive_x/reactive_y + triset() (no
## browser acks needed); a warm axisMode + flushes suffice - the ack machinery
## from test_renderStability is unnecessary here and interacts badly with
## MockShinySession's unprefixed-id update relays.
boot_probe <- function(sx, sy) {
  store <- omicsViewer:::widget_store_new()
  st <- omicsViewer:::widget_store_child(store, "dataspace.sample_space")
  out <- NULL
  shiny::testServer(
    function(input, output, session) omicsViewer:::meta_scatter_module(
      "ds", reactive_meta = reactive(pd), reactive_expr = reactive(ex),
      combine = "pheno", source = "scatter_meta_sample",
      reactive_x = reactive(sx), reactive_y = reactive(sy),
      store = st, selection = omicsViewer:::selection_port(
        omicsViewer:::selection_store_new(), "sample")),
    {
      session$setInputs("ds-axisMode" = "quick")
      for (i in 1:10) session$flushReact()
      ## testServer does not return expr's value - write into the closure
      out <<- list(
        x = paste(unlist(omicsViewer:::store_read(
          st, c("x_analysis", "x_subset", "x_variable"))), collapse = "|"),
        y = paste(unlist(omicsViewer:::store_read(
          st, c("y_analysis", "y_subset", "y_variable"))), collapse = "|"))
    })
  out
}

r <- boot_probe("PCA|All|PC1(", "PCA|All|PC2(")   # the historical demo default
ok(ut_cmp_equal(r$x, "PCA|All|PC1(10.5%)"),
   "probe - truncated default seeds the live PC1 axis")
ok(ut_cmp_equal(r$y, "PCA|All|PC2(7.2%)"),
   "probe - truncated default seeds the live PC2 axis")

r <- boot_probe("PCA|All|PC1(10.5%)", "PCA|All|PC2(7.2%)")  # the fixed demo default
ok(ut_cmp_equal(r$x, "PCA|All|PC1(10.5%)"),
   "probe - exact default seeds unchanged")
ok(ut_cmp_equal(r$y, "PCA|All|PC2(7.2%)"),
   "probe - exact default seeds unchanged (y)")
