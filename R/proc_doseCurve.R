#' Fitting dose-response models for omics data matrix
#' @description A convenient function to fit dose response models for every row in 
#'   an omics matrix using \code{drm} function in the \code{drc} package. 
#' @param x a numerical matrix where the rows are features and columns are samples.
#' @param fitvar a numerical variable has the same length as \code{ncol(x)} to 
#'   indicate the dose/time/concentration conditions. 
#' @param fitvar.name the name of the \code{fitvar}, a length one character. Will be used as the
#'   label for x-axis when drawing the dose curve. 
#' @param curveid a numeric vector or factor containing the grouping of the columns in x. 
#' @param fct.name the function name, e.g. "LL.4()", "LL.3()", "LL.2()" and "LL.5()", 
#'   which are defined in the \code{drc} package. 
#' @importFrom drc drm LL.2 LL2.2, LL.3 LL.3u, LL2.3, LL2.3u LL.4 LL2.4 LL.5 LL2.5
#' @return a list of \code{drc} object
#' @export

drmMat <- function(
    x, fitvar, fitvar.name = c("Dose", "Time", "Concentration")[1], curveid = NA, 
    fct.name = c("LL.4()", "LL.3()", "LL.2()", "LL.5()")[1]
) {
  
  # fct.name <- match.arg(fct.name, choices = c("LL.4()", "LL.3()", "LL.2()", "LL.5()"))
  
  fitvar.name <- make.names( fitvar.name[[1]] )
  
  if (inherits(x, "data.frame")) {
    rn <- rownames(x)
    x <- apply(x, 2, as.numeric)
    rownames(x) <- rn
  }
  
  if (!inherits(x, "matrix"))
    stop("x needs to be a numeric matrix!")
  
  if (length(fitvar) != ncol(x))
    stop("Length of fitvar should be the same as ncol(x)!")
  
  if (length(curveid) == 1 && is.na(curveid))
    curveid <- rep("Observed", ncol(x))
  
  if (length(curveid) != length(fitvar))
    stop("curveid should be the same length as the fitvar!")
  
  df <- data.frame(
    v1 = fitvar,
    curveid = curveid,
    feature = NA
  )
  colnames(df)[1] <- fitvar.name
  
  form <- as.formula(sprintf("feature ~ %s", fitvar.name))
  fct <- eval(parse(text = fct.name))
  
  l <- apply(x, 1, function(x1, dat) {
    dat$feature <- x1
    rr <- try( drm(form, curveid = curveid, data = dat, fct = fct), silent = TRUE)
    if (inherits(rr, "drc"))
      rr$fct <- fct.name else
        rr <- NA
    rr
  }, dat = df, simplify = FALSE)
  
  if (!is.null(rownames(x)))
    names(l) <- rownames(x)
  l
}

#' Draw dose-response curves
#' @param mod an \code{drc} object
#' @param ylab ylab in plot function
#' @param lty lty in plot function
#' @param pch pch in plot function
#' @param cex cex in plot function
#' @param ... other arguments passed to plot function
#' 
plotDC <- function(mod, ylab = "Abundance", lty = 2, pch = 19, cex = 1, ...) {
  
  if (!inherits(mod, "drc"))
    stop("mod needs to be an object of drc!")
  
  op <- par(no.readonly = TRUE)
  on.exit(op)
  
  if (is.character(mod$fct))
    mod$fct <- eval( parse(text = mod$fct) )
  nl <- length(unique( mod$origData[, mod$curveVarNam] ))
  cc <- nColors( nl )
  plot( mod, type = c("confidence"), col = cc, legend = FALSE, ylab = ylab, lty = lty, ...)
  ya <- par()$yaxp
  abline(h = seq(ya[1], ya[2], length.out = ya[3]+1), lty = 3, col = "gray85")
  plot( mod, type = c("obs"), col = cc, pch = pch, add = TRUE, legend = TRUE, cex = cex)
}

#' convert e (inflection point) to EC50
#' @note
#' Only has an effect when using LL.5 and LL2.5 model
.e2EC50 <- function(b, d, e, f) {
  e * (2^(1/f) - 1) ^ (1/b)
}


#' Extracting parameters from drc models
#' @param mod a drc object
#' @param prefix for column header, the column will be named as prefix|curveid|curveparameter
#' @details
#' modelFormula <- function(x, b, c, d, e, f) {
#'   c + (d - c) / (1 + (x/e)^b)^f
#' }
#' d - Minimum asymptote. In a bioassay where you have a standard curve, this can be thought of as the response value at 0 standard concentration.
#' c - Maximum asymptote. In an bioassay where you have a standard curve, this can be thought of as the response value for infinite standard concentration.
#' e - Inflection point. The inflection point is defined as the point on the curve where the curvature changes direction or signs. e is the concentration of analyte where y=(c-d)/2.
#' b - Hill's slope. The Hill's slope refers to the steepness of the curve. It could either be positive or negative.
#' f - Asymmetry factor. When f=1 we have a symmetrical curve around inflection point and so we have a four-parameters logistic equation.
#' @note
#' when LL2.X is used, e is estimated as log(e), this function will return e in linear scale instead. 

extractParamDC <- function(mod, prefix = "ResponseCurve") {
  
  pm <- mod$parmMat
  rownames(pm) <- unique(mod$parNames[[2]])
  
  # if using LL2 model (e is in the form of log(e) ), transform back
  ee <- pm["e", ]
  if (is.character(mod$fct)) {
    if (grepl("LL2.", mod$fct))
      ee <- exp( pm["e", ] )
  } else {
    if (is.null(mod$fct$scaleFct))
      ee <- exp( pm["e", ])
  }
  pm["e", ] <- ee
  
  mod_parName = mod$parNames[[2]]
  mod_id = mod$parNames[[3]]
  mod_par <- c(t(pm))
  
  uv <- unique(iv <- mod$origData[, mod$curveVarNam])
  m2 <- sapply(uv, function(i) {
    d <- mod$predres[iv == i, ]
    cc <- cor.test(d[, 1], rowSums(d))
    c(Pval = cc$p.value, Pseudo.rsq = cc[["estimate"]][[1]]^2)
  })
  colnames(m2) <- uv
  m2 <- rbind(m2, range = pm["c", ] - pm["d", ])
  
  n <- length(uv)
  mod_parName <- c(mod_parName, rep(c("Pval", "Pseudo.rsq", "range"), each = n))
  mod_id <- c(mod_id, rep(mod_id[1:n], 3))
  mod_par <- c(mod_par, c(t(m2)))

  if (any(mod$parNames[[2]] == "f")) { # calculate EC50
    ec50 <- .e2EC50(pm["b", ], pm["d", ], pm["e", ], pm["f", ])
    mod_parName <- c(mod_parName, rep("ec50", n))
    mod_id <- c(mod_id, mod_id[1:n])
    mod_par <- c(mod_par, ec50)
  }
  dfp <- data.frame(
    parName = mod_parName,
    id = mod_id,
    par = mod_par
  )
  rn <- c(b = "b_hill.slope", c = "c_max.response", d = "d_min.response", e = "e_inflection", 
          f = "f_asym.factor", ec50 = "EC50", Pval = "pval", Pseudo.rsq = "pseudo.rsq", range = "range" )
  dfp$parName <- rn[dfp$parName]
  dfp <- dfp[order(dfp$id), ]
  structure(dfp$par, names = paste(prefix, dfp$id, dfp$parName, sep = "|"))
}

#' Extracting parameter from a list of drc object
#' @description
#' Extracting parameter from a list of drc object and return a data.frame,
#' which can be incorporated into the object visualized by omicsViewer
#' @param x a list of drc object
#' @param prefix for column header
#' @return a data.frame
#' 
extractParamDCList <- function(x, prefix = "ResponseCurve") {
  tx <- lapply(x, function(x) {
    if (!inherits(x, "drc"))
      return(NA)
    extractParamDC(x, prefix = prefix)
  })
  uu <- unique( unlist(sapply(tx, names)) )
  t(sapply(tx, "[", uu))
}