#' Drawing ROC and PR curve
#' @param value a numerical vector indicates the predictions
#' @param label true class labels, could be two or more unique values
#' @importFrom ROCR performance prediction
#' @importFrom stats model.matrix
#' @importFrom graphics legend
#' @export
#' @examples 
#' v <- sort(rnorm(100))
#' l <- sample(1:2, size = 100, replace = TRUE)
#' draw_roc_pr(v, l)
#' l <- rep(c("b", "c", "a", "d"), each = 25)
#' draw_roc_pr(v, l)
#' draw_roc_pr(v, sample(l))

draw_roc_pr <- function(value, label) {
  
  pred_one <- function(value, label) {
    cal_auc <- function(pred) {
      perf <- performance(pred, measure = "auc")
      perf@y.values[[1]]
    }
    
    pred <- prediction(value, label)
    auc <- cal_auc(pred)
    if (auc < .5) { # switch label to make AUC > 0.5
      label <- 1 - label
      pred <- prediction(value, label)
      auc <- cal_auc(pred) 
    }
    
    list(
      roc = performance(pred, measure = "tpr", x.measure = "fpr"),
      auc = auc, 
      pr = performance(pred, measure="prec", x.measure="rec")
    )
  }

  label <- as.character(label)
  i <- !is.na(value) & !is.na(label)
  value <- value[i]
  label <- label[i]
  nlab <- length(unique(label))
  if (nlab == 1 || nlab > 12) {
    message("label should have 2 - 12 unique values!")
    plot(0, col = NA, xlab = "", ylab = "", axes = FALSE)
    mtext(text = "label should have 2 - 12 unique values!", side = 3, line = -5)
    return()
  }
  
  dum <- model.matrix( ~ label - 1)
  if (ncol(dum) == 2) {
    dum <- dum[, 1, drop = FALSE]
    colnames(dum) <- "Label"
  } else
    colnames(dum) <- sub("label", "", colnames(dum))
  
    
  curves <- lapply(seq_len(ncol(dum)), function(i) pred_one(value, dum[, i]))
  
  layout(matrix(1:2, 1, 2))
  cols <- nColors(ncol(dum))
  for (i in seq_along(curves)) 
    plot(curves[[i]]$roc, add = i > 1, col = cols[i], main = "ROC Curve", lwd = 2)
  lab <- paste0(colnames(dum), " | ", colSums(dum), " | ", signif(sapply(curves, "[[", "auc"), 2))
  legend("bottomright", col = cols, legend = lab, lty = 1, lwd = 2, bty = "n", title = "label | n | AUC")
  
  for (i in seq_along(curves)) 
    plot(curves[[i]]$pr, add = i > 1, col = cols[i], xlim = c(0, 1), ylim = c(0, 1), main = "PR Curve", lwd = 2)
}



