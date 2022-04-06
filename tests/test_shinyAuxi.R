library(unittest)
library(RColorBrewer)

nColors <- omicsViewer:::nColors
ok(ut_cmp_equal(nColors(k=0, stop = FALSE), NULL), "nColors - return NULL when not stop")
ok(ut_cmp_error(nColors(k=0, stop = TRUE), "k should be an integer between 1 and 60!"), 
   "nColors - stop only when stop")

null2empty <- omicsViewer:::null2empty
ok(ut_cmp_identical(null2empty(NULL), ""), "null2empty")

getSearchCols <- omicsViewer:::getSearchCols
ok(ut_cmp_identical(getSearchCols(NULL), NULL), "getSearchCols - return NULL when NULL given")

getOrderCols <- omicsViewer:::getOrderCols
ok(ut_cmp_identical(getOrderCols(NULL), NULL), "getOrderCols - return NULL when NULL given")

exprsImpute <- omicsViewer:::getExprsImpute
ok(ut_cmp_identical(exprsImpute("5"), NULL), "exprsImpute - return NULL when error")


value2color <- omicsViewer:::value2color
l <- value2color(1:100, n=10)
ok(
   ut_cmp_equal(names(l), c("color", "key")),
   "value2color - return color and key"
   ) 
ok(
   ut_cmp_equal(length(l$color), 100),
   "value2color - return color correct length"
   ) 
ok(
   ut_cmp_equal(c(table(l$color)), rep(10, 10), check.attributes = FALSE),
   "value2color - return key correct length"
   ) 
ok(
   ut_cmp_equal(length(unique(l$key)), 10),
   "value2color - return color correct unique length"
   ) 
