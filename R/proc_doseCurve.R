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
#' @importFrom drc drm LL.2 LL2.2 LL.3 LL.3u LL2.3 LL2.3u LL.4 LL2.4 LL.5 LL2.5
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
  pl <- try(
    plot( mod, type = c("confidence"), col = cc, legend = FALSE, ylab = ylab, lty = lty, ...),
    silent = TRUE
  )
  if (inherits(pl, "try-error"))
    plot( mod, type = "none", col = cc, legend = FALSE, ylab = ylab, lty = lty, ...)
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

#' model fitted by drc
#' @details 
#' func(x) = c + (d - c) / (1 + (x/e)^b)^f
#' @param d Minimum asymptote. In a bioassay where you have a standard curve, this can be thought of as the response value at 0 standard concentration.
#' @param c Maximum asymptote. In an bioassay where you have a standard curve, this can be thought of as the response value for infinite standard concentration.
#' @param e Inflection point. The inflection point is defined as the point on the curve where the curvature changes direction or signs. e is the concentration of analyte where y=(c-d)/2.
#' @param b Hill's slope. The Hill's slope refers to the steepness of the curve. It could either be positive or negative.
#' @param f Asymmetry factor. When f=1 we have a symmetrical curve around inflection point and so we have a four-parameters logistic equation.
.modelFormula <- function(x, b, c = 0, d = 1, e, f = 1) {
  c + (d - c) / (1 + (x/e)^b)^f
}


#' Extracting parameters from drc models
#' @param mod a drc object
#' @param prefix for column header, the column will be named as prefix|curveid|curveparameter
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
  
  uv <- unique(iv <- mod_id)
  m2 <- sapply(uv, function(i) {
    d <- mod$predres[iv == i, ]
    cc <- cor.test(d[, 1], rowSums(d))
    c(Pval = cc$p.value, log.Pval = -log10(cc$p.value), Pseudo.rsq = cc[["estimate"]][[1]]^2)
  })
  colnames(m2) <- uv
  m2 <- rbind(m2, range = pm["c", ] - pm["d", ])
  
  n <- length(uv)
  mod_parName <- c(mod_parName, rep(c("Pval", "log.Pval", "Pseudo.rsq", "range"), each = n))
  mod_id <- c(mod_id, rep(mod_id[1:n], 4))
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
          f = "f_asym.factor", ec50 = "EC50", Pval = "pval", log.Pval = "log.pval", 
          Pseudo.rsq = "pseudo.rsq", range = "range" )
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

#' Draw dose response curve given parameters in the omicsViewer object
#' @description
#' Draw dose response curve given the feature Data/rowData, phenotype data/colData 
#' and expression matrix. The function is usually used in shinyApp. 
#' @param expr expression matrix
#' @param pd phenotype data or colData
#' @param fd feature data or rowData
#' @param featid feature id to be visualized
#' @param dose.var the column header indicating the dose/time/concentration
#' @param curve.var the column header indicating the curve ids
#' @param return.par logical value. If true, no plot generated,
#'   the function only returns the parameters of models.

plotDCMat <- function(expr, pd, fd, featid, dose.var, curve.var, only.par = FALSE) {
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  v <- grep("ResponseCurve", colnames(fd), value = TRUE, ignore.case = TRUE)
  vs <- sub("ResponseCurve\\|", "", v)
  se <- str_split_fixed(vs, pattern = "\\|", n = 2)
  rd <- data.frame(
    se, unlist( fd[featid, v] )
  )
  colnames(rd) <- c("group", "param", "value")
  rdl <- split(rd, rd$group)
  
  if (only.par) {
    ll <- list(
      par = data.frame(
        param = rdl[[1]]$param, 
        sapply(rdl, "[[", "value"),
        check.names = FALSE
      ),
      featInfo = fd[featid, ]
    )
    return(ll)
  }
  
  cid <- pd[, curve.var]
  dose <- pd[, dose.var]
  dd <- data.frame(
    feat = expr[featid, ],
    dose = dose,
    curve = cid
  )
  fdl <- split(dd, dd$curve)
  cc <- nColors( length(fdl) )
  names(cc) <- unique(unique(rd$group))
  ###
  plot(dd$dose, dd$feat, col = cc[dd$curve], pch = 19, 
       xlab = dose.var, ylab = "Response / abundance")
  ya <- par()$yaxp
  abline( h = seq(ya[1], ya[2], length.out = ya[3]+1), lty = 3, col = "gray")
  
  lc <- seq( min(dose), max(dose), length.out = 100 )
  i <- 1
  leg <- c()
  gid <- c()
  
  for (i in seq_along(fdl)) {
    d1 <- fdl[[i]]
    d2 <- rdl[[i]]
    ss <- structure(d2$value, names = d2$param)
    c = 0
    d = 1
    f = 1
    b <- ss[["b_hill.slope"]]
    if (any(grepl("c_max.response", names(ss))))
      c <- ss[["c_max.response"]]
    if (any(grepl("d_min.response", names(ss))))
      d <- ss[["d_min.response"]]
    e <- ss[["e_inflection"]]
    if (any(grepl("f_asym.factor", names(ss))))
      f <- ss[["f_asym.factor"]]
    lg <- paste(
      d2$group[[1]],
      sprintf("p-value=%s", signif(ss[["pval"]], 2)),
      sprintf("pseudo.rsq=%s", signif(ss[["pseudo.rsq"]], 2)),
      sprintf("range=%s", signif(ss[["range"]], 2)),
      sep = ";"
    )
    leg <- c(leg, lg)
    gid <- c(gid, d2$group[[1]])
    yy <- .modelFormula( lc, b = b, c = c, d = d, e = e, f = f )
    lines(lc, yy, type = "l", lty = 2, col = cc[d1$curve[1]], lwd = 2)
  }
  par(xpd = TRUE)
  legend(
    y = par("usr")[4],
    x = par()$xaxp[1], legend = leg, col = cc[gid], pch = 15,
    yjust = 0,
    bty = "n")
}

