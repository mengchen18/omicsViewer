# ExpressionSetViewer

ExpressionSet is a standard S4 object storing high-throughput omics data in Bioconductor. 
It consists mainly of three components 1) the expression matrix where the rows are features (genes, 
mRNAs, proteins, etc) and columns are the samples (cell lines, patients, etc); 2) phenotype data - 
the annotation information of samples and 3) feature data – the annotation information of features. 
Interpreting biology from omics data relies on the integrative analysis of the data triplet and a 
wide range of statistical methods are often used to answer the different biological questions. 

"ExpressionSetViewer" visualizes ExpressionSet in an interactive way. The ExpressionSetViewer has a 
separate back- and front-end. In the back-end, users need to prepare an ExpressionSet that contains all 
the necessary information for the downstream data interpretation. The pure dependency on R/Bioconductor 
guarantees maximum flexibility in the statistical analysis in the back-end. Once the ExpressionSet is prepared, 
it can be visualized using the front-end, implemented by shiny and plotly. Both features and samples could 
be selected from (data) tables or graphs (scatter plot/heatmap). Different types of analyses, such as 
enrichment analysis (using Bioconductor package fgsea or fisher’s exact test) and STRING network analysis, 
will be performed on-the-fly and the results are visualized simultaneously. When a subset of samples and a 
phenotype variable is selected, a significance test on means (t-test or ranked based test; when phenotype 
variable is quantitative) or test of independence (chi-square or fisher’s exact test; when phenotype data 
is categorical) will be performed to test the association between the phenotype of interest with the selected 
samples. Therefore, ExpressionSetViewer will greatly facilitate data exploration, many different hypotheses 
can be explored in a short time without the need for knowledge of R.

In addition, the resulted data could be easily shared using a shiny server. Otherwise, a standalone version 
of ExpressionSetViewer together with designated omics data could be easily created by integrating it with 
portable R or with docker, which can be shared with collaborators or submitted as supplementary data together 
with a manuscript.

## Vignette
[Link](https://mengchen18.github.io/ExpressionSetViewer/index.html)

## Install package
### Quick installation:
```
devtools::install_github("mengchen18/ExpressionSetViewer/package")
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
  "randomcoloR",
  "plotly",
  "openxlsx",
  "yaml",
  "curl", 
  "sortable",
  "BiocManager",
  "password",
  "ggseqlogo",
  "devtools"
  )

# # BIOC
s2 <- c(
  "Biobase", "fgsea"
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

devtools::install_github('omarwagih/rmotifx')

if (length(missingPkg) > 0)
  stop(paste("this packages are missing", paste(missingPkg, collapse = " ")))
```

If you see error messages, please solve them first. 
Then install the `ExpressionSetViewer` packages
```
devtools::install_github("mengchen18/ExpressionSetViewer/package")
```


## Start the viewer

```
library(ExpressionSetViewer)
ExpressionSetViewer(system.file("extdata", package = "ExpressionSetViewer"))
```




