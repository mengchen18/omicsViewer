library(omicsViewer)
library(unittest, quietly = TRUE)
library(Matrix)
library(fastmatch)

# =========================== conversion ========================
gs <- cbind(
  s1 = c(1, 1, 1, 1, 0, 0, 0, 0),
  s2 = c(0, 0, 0, 0, 1, 1, 1, 1),
  s3 = c(0, 0, 1, 1, 1, 1, 0, 0))
stats <- 1:8
rownames(gs) <- names(stats) <- paste0("g", 1:8)
gss <- as(gs, "dgCMatrix")
ok(ut_cmp_equal(omicsViewer:::csc2list(gs), omicsViewer:::csc2list(gss)), 
  "convert matrix to data.frame for fgsea")
df <- omicsViewer:::csc2list(gs)
ok(
  ut_cmp_equal(
    gss, omicsViewer:::list2csc(df, dimnames = dimnames(gs)), check.attributes = FALSE
    ), "convert data.frame to matrix"
   )

vectORA <- omicsViewer:::vectORA
r1 <- vectORA( gs, i = 1:4, minOverlap = 1, minSize = 1 )
r2 <- vectORA( gss, i = 1:4, minOverlap = 1, minSize = 1 )
vectORATall <- omicsViewer:::vectORATall
r3 <- vectORATall( df, i=rownames(gs)[1:4], background=8 )

ok(
  ut_cmp_equal(r1$pathway, c("s1", "s3"), check.attributes = FALSE),
  "vectORA 1"
  )
ok(ut_cmp_equal(r1, r2, check.attributes = FALSE), "vectORA 2")
ok(ut_cmp_equal(r1, r3), "vectORATall")



res <- omicsViewer:::fgsea1(
  gs = gss, stats = stats, minSize = 1, maxSize = 500, sampleSize = 3)
ut_cmp_identical(res$pathway, colnames(gss))
ok(ut_cmp_equal(res$pval[1] < 0.2, TRUE), "fgsea significant pathway 1")
ok(ut_cmp_equal(res$pval[2] < 0.2, TRUE), "fgsea significant pathway 2")
ok(ut_cmp_equal(res$pval[3] > 0.6, TRUE), "fgsea insignificant pathway")

totall <- omicsViewer:::totall
tgs <- totall(gs)
ok(ut_cmp_identical(colnames(tgs), c("featureId", "gsId", "weight")), 
  "totall 1")
ok(ut_cmp_equal(nrow(tgs), 12), "totall 2")


# ====================== trisetters ==========================

expr <- matrix(1:6, 3, 2)
rownames(expr) <- c("g1", "g2", "g3")
colnames(expr) <- c("s1", "s2")

m1 <- data.frame(
  "F1|pos2|pos3" = 1:3,
  "F2|pos2|pos3" = 1:3,
  check.names = FALSE
)
rownames(m1) <- rownames(expr)

m2 <- data.frame(
  "var1|pos2|pos3" = 1:2,
  "var2|pos2|pos3" = 1:2,
  check.names = FALSE
)
rownames(m2) <- colnames(expr)

trisetter <- omicsViewer:::trisetter

f1 <- rbind(c("F1", "pos2", "pos3"),
            c("F2", "pos2", "pos3"))
ok(
  ut_cmp_equal(
    trisetter(meta = m1, expr=NULL, combine="none"), f1,
    check.attributes = FALSE
  ),
  "trisetter feature meta"
)

var1 <- rbind(c("var1", "pos2", "pos3"),
              c("var2", "pos2", "pos3"))
ok(
  ut_cmp_equal(
    trisetter(meta = m2, expr=NULL, combine="none"), var1,
    check.attributes = FALSE
  ),
  "trisetter pheno meta"
)

cf1 <- rbind(f1, c("Sample", "Auto", "s1"),
             c("Sample", "Auto", "s2"))
ok(
  ut_cmp_equal(
    trisetter(meta = m1, expr=expr, combine="feature"), cf1,
    check.attributes = FALSE
  ),
  "trisetter feature combined with expr"
)
cvar1 <- rbind(var1, 
               c("Feature", "Auto", "g1"),
               c("Feature", "Auto", "g2"),
               c("Feature", "Auto", "g3"))
ok(
  ut_cmp_equal(
    trisetter(meta = m2, expr=expr, combine="pheno"), cvar1,
    check.attributes = FALSE
  ),
  "trisetter pheno combined with expr"
)

varSelector <- omicsViewer:::varSelector
l1 <- list(analysis = "Feature", subset= "Auto", variable = "g1")
ok(
ut_cmp_equal(
  varSelector(x = l1, expr = expr, meta = m2),
  c(1, 4),
  check.attributes = FALSE
  ),
"select from triselector - feature"
)
  
l1 <- list(analysis = "Sample", subset= "Auto", variable = "s1")
ok(
ut_cmp_equal(
  varSelector(x = l1, expr = expr, meta = m2), 
  1:3, 
  check.attributes = FALSE),
"select from triselector - sample"
)
text2num <- omicsViewer:::text2num
ok(ut_cmp_equal(text2num("-log10(0.01)"), -log10(0.01)), "text2num - 1")
ok(ut_cmp_equal(text2num("0.05"), 0.05), "text2num - 2")

######################################
terms <- data.frame(
  id = c("ID1", "ID2", "ID1", "ID2", "ID8", "ID10"),
 term = c("T1", "T1", "T2", "T2", "T2", "T2"),
  stringsAsFactors = FALSE
)
features <- list(c("ID1", "ID2"), c("ID13"), c("ID4", "ID8", "ID10"))
gsAnnotIdList(idList = features, gsIdMap = terms, minSize = 1, maxSize = 500)

terms <- data.frame(
id = c("ID1", "ID2", "ID1", "ID2", "ID8", "ID10", "ID4", "ID4"),
term = c("T1", "T1", "T2", "T2", "T2", "T2", "T1", "T2"),
stringsAsFactors = FALSE
)
features <- list(F1 = c("ID1", "ID2", "ID4"), F2 = c("ID13"), F3 = c("ID4", "ID8", "ID10"))

res <- data.frame(
  featureId = c(1, 1, 3, 3),
  gsId = c("T1", "T2", "T2", "T1"),
  weight = rep(1, 4),
  stringsAsFactors = FALSE
)
r1 <- gsAnnotIdList(features, gsIdMap = terms, data.frame = TRUE, minSize = 1)
ok(ut_cmp_equal(r1, res, check.attributes = FALSE), 
   "gsAnnotIdList data.frame"
   )

res <- sparseMatrix(i = c(1, 1, 3, 3), j = c(1, 2, 1, 2), x = 1)
colnames(res) <- c("T1", "T2")
r2 <- gsAnnotIdList(features, gsIdMap = terms, data.frame = FALSE, minSize = 1)
ok(ut_cmp_equal(r2, res), "gsAnnotIdList sparseMatrix")

# ==============================================
xq <- rbind(c(4, 2, 4),
            c(20, 40, 10),
            c(11, 234, 10),
            c(200, 1000, 100))

vectORA.core <- omicsViewer:::vectORA.core
r1 <- vectORA.core(xq[1, ], xq[2, ], xq[3, ], xq[4, ])
r2 <- vectORA.core(xq[1, ], xq[2, ], xq[3, ], xq[4, ], unconditional.or = FALSE)

ok(ut_cmp_equal(r1$p.value, r2$p.value), "vecORA.core 1")

# fisher's test
pv <- t(apply(xq, 2, function(x1) {
  m <- rbind(c(x1[1], x1[2]-x1[1]),
             c(x1[3]-x1[1], x1[4] - x1[2] - x1[3] + x1[1]))
  v <- fisher.test(m, alternative = "greater")
  c(p.value = v$p.value, v$estimate)
}))

ok(ut_cmp_equal(r1$p.value, pv[, "p.value"]), "vectORA.core - conditional OR")
ok(ut_cmp_equal(r2$OR, pv[, "odds ratio"]), "vectORA - unconditional OR")


# ================= module-level regression ======================
# ORA must run on the raw feature ids when the collapse triple is unset.
# Cold start: the unified triselector never commits a variable (analysis
# defaults, variable sits at "--select--", no allow_unset) and returns
# NULL - the old req(v1()$variable) suspended the entry observer
# permanently, so rii() was never written and the whole tab rendered
# blank on every gene selection (regression vs the pre-unification
# triselector, which preselected "--select--" and drove the no-collapse
# branch). The fix treats NULL like the historical "--select--".
#
# Observable: the results DT ("stab-table") renders only when OT()
# produced a data.frame, which requires rii() to have been written. On
# the broken build the render suspends (req) and reading the output
# errors. The cold start guarantees the no-collapse branch (v1() NULL).
omv_reg_fd <- local({
  dat <- readRDS(system.file("extdata", "demo.RDS", package = "omicsViewer"))
  fd <- Biobase::fData(dat)
  attr(fd, "GS") <- data.frame(
    featureId = factor(rownames(fd)[1:40]),
    gsId = factor(rep(paste0("gs", 1:4), each = 10)),
    weight = 1)
  fd
})
omv_reg_ids <- head(rownames(omv_reg_fd), 5)

omv_ora_host <- function(sel_ids) function(input, output, session) {
  omicsViewer:::enrichment_analysis_module(
    "ora",
    reactive_featureData = shiny::reactive(omv_reg_fd),
    reactive_i = shiny::reactive(sel_ids))
  shiny::outputOptions(output, "ora-stab-table", suspendWhenHidden = FALSE)
}
omv_reg_res <- NULL
shiny::testServer(omv_ora_host(omv_reg_ids), {
  for (i in 1:10) session$flushReact()
  omv_reg_res <<- tryCatch(
    output[["ora-stab-table"]],
    error = function(e) conditionMessage(e))
})
ok(inherits(omv_reg_res, "json") &&
     grepl("size_backgroung", omv_reg_res, fixed = TRUE),
   "ORA: results table renders without a manual collapse pick (no-collapse branch)")

omv_reg_neg <- NULL
shiny::testServer(omv_ora_host(character(0)), {
  for (i in 1:10) session$flushReact()
  omv_reg_neg <<- tryCatch(
    output[["ora-stab-table"]],
    error = function(e) conditionMessage(e))
})
ok(!inherits(omv_reg_neg, "htmlwidget"),
   "ORA: no selection -> no results table (control)")

# ---------------- R-H2: NA in the collapse column ----------------
# A LEADING NA in the collapsed values used to reach rii() as rii()[1] ==
# NA, and `if (NA == "notest")` in OT() errored the eager oraTab observer
# ("missing value where TRUE/FALSE needed" -> session close). Non-leading
# NAs survived as a pseudo-gene (background + overlap polluted). The fix
# drops NA/"" from val/ck BEFORE collapsing and compares with identical().
omv_na_fd <- local({
  dat <- readRDS(system.file("extdata", "demo.RDS", package = "omicsViewer"))
  fd <- Biobase::fData(dat)[1:40, , drop = FALSE]
  fd
})
attr(omv_na_fd, "GS") <- data.frame(
  featureId = factor(rownames(omv_na_fd)),
  gsId = factor(rep(paste0("gs", 1:4), each = 10)),
  weight = 1)
# f1/f2 collapse to NA (leading NA in the selected set), f3-f6 to A-D
omv_na_fd$`cat|sub|col` <- c(NA, NA, "A", "B", "C", "D", rep("Z", 34))
omv_na_sel <- head(rownames(omv_na_fd), 6)

omv_na_host <- function() function(input, output, session) {
  session$userData$ora_ret <- omicsViewer:::enrichment_analysis_module(
    "ora",
    reactive_featureData = shiny::reactive(omv_na_fd),
    reactive_i = shiny::reactive(omv_na_sel),
    reactive_status = shiny::reactive(list(xax = c("cat", "sub", "col"))))
  shiny::outputOptions(output, "ora-stab-table", suspendWhenHidden = FALSE)
  shiny::outputOptions(output, "ora-errorMsg", suspendWhenHidden = FALSE)
}
omv_na_res <- NULL; omv_na_tri <- NULL
shiny::testServer(omv_na_host(), {
  for (i in 1:5) session$flushReact()
  # drive the collapse cascade directly (testServer relays no
  # updateSelectInput messages; the commit observer reads the inputs)
  session$setInputs(`ora-tris_ora-analysis` = "cat",
                    `ora-tris_ora-subset` = "sub",
                    `ora-tris_ora-variable` = "col")
  for (i in 1:10) tryCatch(session$flushReact(), error = function(e) NULL)
  omv_na_tri <<- tryCatch(session$userData$ora_ret(), error = function(e) NULL)
  omv_na_res <<- tryCatch(
    output[["ora-stab-table"]],
    error = function(e) conditionMessage(e))
})
ok(!is.null(omv_na_tri) && identical(omv_na_tri$xax$variable, "col"),
   "ORA: collapse cascade committed the NA-bearing variable (R-H2 branch evidence)")
ok(inherits(omv_na_res, "json"),
   "ORA: leading-NA collapse column still renders results (R-H2, no session crash)")

# ---------------- R-M6: selection shrinks below testable minimum ----------------
# The old observer only wrote oraTab when rii() was non-NULL, so shrinking
# the selection to <= 1 collapsed entity left the PREVIOUS results table
# on screen. The fix writes an explicit no-test message once something was
# shown (and stays quiet on cold start).
omv_m6_host <- function() function(input, output, session) {
  omicsViewer:::enrichment_analysis_module(
    "ora",
    reactive_featureData = shiny::reactive(omv_reg_fd),
    reactive_i = shiny::reactive(
      if (isTRUE(input$shrink %in% 1)) head(rownames(omv_reg_fd), 1)
      else omv_reg_ids))
  shiny::outputOptions(output, "ora-errorMsg", suspendWhenHidden = FALSE)
}
omv_m6_pre <- NULL; omv_m6_post <- NULL
shiny::testServer(omv_m6_host(), {
  for (i in 1:10) session$flushReact()
  omv_m6_pre <<- tryCatch(output[["ora-errorMsg"]],
                          error = function(e) conditionMessage(e))
  session$setInputs(shrink = 1)
  for (i in 1:10) session$flushReact()
  omv_m6_post <<- tryCatch(output[["ora-errorMsg"]],
                           error = function(e) conditionMessage(e))
})
ok(!is.null(omv_m6_pre) && !grepl("Too few", omv_m6_pre, fixed = TRUE),
   "ORA: no spurious no-test message before the selection shrinks (R-M6 control)")
ok(grepl("Too few feature IDs", omv_m6_post, fixed = TRUE),
   "ORA: shrunk selection shows explicit no-test message (R-M6)")

# ---------------- R-H3: jaccardList sparse reimplementation ----------------
jaccard_ref <- function(x) {
  ax <- unique(unlist(x))
  m <- sapply(x, function(x1) as.integer(ax %in% x1))
  ist <- crossprod(m)
  uni <- apply(m, 2, function(x) colSums(x + m > 0))
  as.dist(1 - ist / uni)
}
jl_sets <- replicate(60, sample(paste0("g", 1:300), sample(20:80, 1)),
                     simplify = FALSE)
ok(ut_cmp_equal(
  as.matrix(omicsViewer:::jaccardList(jl_sets)),
  as.matrix(jaccard_ref(jl_sets))),
  "jaccardList: sparse version equals the dense reference (R-H3)")
jl_sets2 <- replicate(40, sample(paste0("g", 1:50), sample(5:50, 1)),
                      simplify = FALSE)
jl_sets2[[3]] <- c(jl_sets2[[3]], jl_sets2[[3]][1])  # duplicate entry in a set
ok(ut_cmp_equal(
  as.matrix(omicsViewer:::jaccardList(jl_sets2)),
  as.matrix(jaccard_ref(jl_sets2))),
  "jaccardList: duplicate entries within a set are neutralized (R-H3)")
ok(inherits(omicsViewer:::jaccardList(list(a = character(0))), "dist"),
  "jaccardList: degenerate input returns a dist object (R-H3)")
