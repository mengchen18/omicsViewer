library(unittest)
library(drc)

ds <- c(0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4)
y <- omicsViewer:::.modelFormula(ds, b = -1, c = 0.25, d = 1.1, f = 1, e = 2) + rnorm(length(ds), sd = 1e-5)
mod1 <- drm(y ~ ds, fct = LL.4())
f1 <- omicsViewer:::extractParamDC(mod1)
mod2 <- drm(y ~ ds, fct = LL2.4())
f2 <- omicsViewer:::extractParamDC(mod2)

ok(
  all.equal(f1, f2, tolerance = 7.5e-5),
  "functions - .modelFormula, extractParamDC"
) 


testPar <- function( pars, cid = "(Intercept)", b, c, d, e, f, tol = 0.1 ) {
  r <- TRUE
  nr <- nrow(pars)
  id <- paste("ResponseCurve", cid, sep = "|")
  
  iae <- function(x, y, tolerance) isTRUE(all.equal(x, y, tolerance = tolerance))
  if (!missing(b)) {
    rb <- iae(pars[, sprintf("%s|b_hill.slope", id)], rep(b, nr), tolerance = tol)
    r <- r && rb
  }
  
  if (!missing(c)) {
    rc <- iae(pars[, sprintf("%s|c_min.response", id)], rep(c, nr), tolerance = tol)
    r <- r && rc
  }
  
  if (!missing(d)) {
    rd <- iae(pars[, sprintf("%s|d_max.response", id)], rep(d, nr), tolerance = tol)
    r <- r && rd
  }
  
  if (!missing(e)) {
    re <- iae(pars[, sprintf("%s|e_inflection", id)], rep(e, nr), tolerance = tol)
    r <- r && re
  }
  
  if (!missing(f)) {
    rf <- iae(pars[, sprintf("%s|f_asym.factor", id)], rep(f, nr), tolerance = tol)
    r <- r && rf
  }
    
  ok(r)
}

d_log <- -(10:-1)
dd <- 2^d_log

for (eb in c(1, -1)) {
  ec <- 0
  ed <- 1
  ee <- 0.005
  ef <- 1
  
  m1 <- sapply(1:20, function(n) {
    y_sim <- omicsViewer:::.modelFormula(dd, b = eb, c = ec, d = ed, e = ee, f = ef)
    y_sim + rnorm(length(y_sim), sd = 0.0001)
  })
  m1 <- t(m1)
  
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.5()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, d = ed, e = ee, f = ef)
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.5()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, d = ed, e = ee, f = ef)
  
  #
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.4()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, d = ed, e = ee)
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.4()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, d = ed, e = ee)
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.3()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb,  d = ed, e = ee)
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.3()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, d = ed, e = ee)
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.3u()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, e = ee)
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.3u()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, e = ee)
  
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.2()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, e = ee)
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.2()")
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, e = ee)
}

d_log <- rep(-(10:-1), 2)
dd <- 2^d_log

for (eb in c(1, -1)) {
  ec <- 0
  ed <- 1
  ee <- 0.005
  ef <- 1
  
  m1 <- sapply(1:5, function(n) {
    y_sim <- omicsViewer:::.modelFormula(dd, b = eb, c = ec, d = ed, e = ee, f = ef)
    y_sim + rnorm(length(y_sim), sd = 0.0001)
  })
  m1 <- t(m1)
  curve <- rep(c("a", "b"), each = 12)
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.5()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, d = ed, e = ee, f = ef, cid = "a")
  testPar(pp, b = eb, c = ec, d = ed, e = ee, f = ef, cid = "b")
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.5()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, d = ed, e = ee, f = ef, cid = "a")
  testPar(pp, b = eb, c = ec, d = ed, e = ee, f = ef, cid = "b")
  #
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.4()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, d = ed, e = ee, cid = "a")
  testPar(pp, b = eb, c = ec, d = ed, e = ee, cid = "b")
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.4()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, d = ed, e = ee, cid = "a")
  testPar(pp, b = eb, c = ec, d = ed, e = ee, cid = "b")
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.3()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb,  d = ed, e = ee, cid = "a")
  testPar(pp, b = eb,  d = ed, e = ee, cid = "b")
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.3()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, d = ed, e = ee, cid = "a")
  testPar(pp, b = eb, d = ed, e = ee, cid = "b")
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.3u()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, e = ee, cid = "a")
  testPar(pp, b = eb, c = ec, e = ee, cid = "b")
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.3u()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, c = ec, e = ee, cid = "a")
  testPar(pp, b = eb, c = ec, e = ee, cid = "b")
  
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL.2()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, e = ee, cid = "a")
  testPar(pp, b = eb, e = ee, cid = "b")
  
  mod <- omicsViewer:::drmMat(m1, fitvar = dd, fitvar.name = "Dose", fct.name = "LL2.2()", curveid = curve)
  pp <- omicsViewer:::extractParamDCList(mod)
  testPar(pp, b = eb, e = ee, cid = "a")
  testPar(pp, b = eb, e = ee, cid = "b")
}

