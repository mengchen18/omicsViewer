# omicsViewer

"SummarizedExperiment" and the historical "ExpressionSet" are S4 objects storing high throughput omics data. The core component of the objects is an expression matrix, where the rows are features, such as genes, proteins, and columns are samples. The values in the matrix represent the abundance or presence/absence of features. The meta-information about features (rows) and samples (columns) are stored in *data.frames*-like object called "feature data" (or "row data") and "phenotype data" (or "col data), respectively. More detailed instructions of _ExpressionSet_ and _SummarizeExperiment_ could be found [here](https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) and [here](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html).

omicsViewer visualizes ExpressionSet/SummarizedExperiment in an interactive way to facilitate data exploration. The object to be visualized needs to be prepared in R, usually, it contains all the necessary information for the data interpretation, e.g. functional annotation of features, meta information of samples, and the statistical results (See [Prepare the ExpressionSet]). Using R/Bioconductor in the data preparation stage guarantees maximum flexibility in the statistical analysis. Once the object is prepared, it can be visualized using the omicsViewer interface, implemented using shiny and plotly. At this stage, coding in R is not required, therefore, it can be shared with anyone who does not have experience with any programming knowledge. 

Both features and samples could be selected from (data) tables or graphs (scatter plot/heatmap) in the omicsViewer. Different types of analyses, such as enrichment analysis (using Bioconductor package fgsea or fisher's exact test) and STRING network analysis, will be performed on the fly when a subset of features is selected. The statistical results are visualized simultaneously. When a subset of samples and a phenotype variable is selected, a significance test on mean difference (t-test or ranked based test; when phenotype variable is quantitative) or test of independence (chi-square or fisher’s exact test; when phenotype data is categorical) will be performed to test the association between the phenotype of interest with the selected samples. Therefore, the immediate analyses and result visualization in omicsViewer greatly facilitates data exploration, many different hypotheses can be explored in a short time without the need for knowledge of R.

In addition, the resulted data could be easily shared using a shiny server. Otherwise, a standalone version 
of omicsViewer together with designated omics data could be easily created by integrating it with 
portable R or with docker, which can be shared with collaborators or submitted as supplementary data together with a manuscript.

## Vignette
[Link](https://mengchen18.github.io/omicsViewer/index.html)

## Install package
### Quick installation:
```
devtools::install_github("mengchen18/omicsViewer")
```
If you see an error message, perhaps because you are missing some dependencies, please follow the next section.
### Installation from scratch
First, installing dependencies as follow:
```
install.packages("V8")
install.packages("SparseM")
# cran
s1 <- c(
  "survminer",
  "survival",
  "fastmatch",
  "reshape2",
  "beeswarm",
  "grDevices",
  "shinycssloaders",
  "shinythemes",
  "networkD3",
  "httr",
  "RColorBrewer",
  "psych",
  "stringr",
  "shiny",
  "shinydashboard",
  "shinyWidgets",
  "shinybusy",
  "matrixStats",
  "flatxml",
  "excelR",
  "shinyjs",
  "shinyFiles",
  "DT",
  "plotly",
  "openxlsx",
  "yaml",
  "curl", 
  "sortable",
  "BiocManager",
  "password",
  "ggseqlogo",
  "devtools",
  "RSQLite"
  )

# # BIOC
s2 <- c(
  "Biobase", "fgsea",
  "S4Vectors",
  "SummarizedExperiment"
  )

# 
lapply(s1, function(x) {
  if (x %in% installed.packages()[, 1])
    return()
  install.packages(x)
})
# 
lapply(s2, function(x) {
  if (x %in% installed.packages()[, 1])
    return()
  BiocManager::install(x, update = FALSE)
})
# 
a <- installed.packages()[,1 ]
# 
xs <- c(s1, s2)
missingPkg <- setdiff(xs, a)

if (length(missingPkg) > 0)
  stop(paste("this packages are missing", paste(missingPkg, collapse = " ")))
```

If you see error messages, please solve them first. 
Then install the `omicsViewer` packages
```
devtools::install_github("mengchen18/omicsViewer")
```


## Start the viewer

```
library(omicsViewer)
omicsViewer(system.file("extdata", package = "omicsViewer"))
```


## Docker image


https://hub.docker.com/repository/docker/mengchen18/omicsviewer
