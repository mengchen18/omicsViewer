library(unittest)
library(drc)

ds <- c(0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4)
y <- .modelFormula(ds, b = -1, c = 0.25, d = 1.1, f = 1, e = 2) + rnorm(length(ds), sd = 1e-5)
mod1 <- drm(y ~ ds, fct = LL.4())
f1 <- extractParamDC(mod1)
mod2 <- drm(y ~ ds, fct = LL2.4())
f2 <- extractParamDC(mod2)

ok(
  all.equal(f1, f2, tolerance = 7.5e-5),
  "functions - .modelFormula, extractParamDC"
) 