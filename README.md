# omicsViewer

"SummarizedExperiment" and the historical "ExpressionSet" are S4 objects storing high throughput omics data. The core component of the objects is an expression matrix, where the rows are features, such as genes, proteins, and columns are samples. The values in the matrix represent the abundance or presence/absence of features. The meta-information about features (rows) and samples (columns) are stored in *data.frames*-like object called "feature data" (or "row data") and "phenotype data" (or "col data), respectively. More detailed instructions of _ExpressionSet_ and _SummarizeExperiment_ could be found [here](https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) and [here](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html).

omicsViewer visualizes ExpressionSet/SummarizedExperiment in an interactive way to facilitate data exploration. The object to be visualized needs to be prepared in R, usually, it contains all the necessary information for the data interpretation, e.g. functional annotation of features, meta information of samples, and the statistical results. Using R/Bioconductor in the data preparation stage guarantees maximum flexibility in the statistical analysis. Once the object is prepared, it can be visualized using the omicsViewer interface, implemented using shiny and plotly. At this stage, coding in R is not required, therefore, it can be shared with anyone who does not have experience with any programming knowledge. 

Both features and samples could be selected from (data) tables or graphs (scatter plot/heatmap) in the omicsViewer. Different types of analyses, such as enrichment analysis (using Bioconductor package fgsea or fisher's exact test) and STRING network analysis, will be performed on the fly when a subset of features is selected. The statistical results are visualized simultaneously. When a subset of samples and a phenotype variable is selected, a significance test on mean difference (t-test or ranked based test; when phenotype variable is quantitative) or test of independence (chi-square or fisher’s exact test; when phenotype data is categorical) will be performed to test the association between the phenotype of interest with the selected samples. Therefore, the immediate analyses and result visualization in omicsViewer greatly facilitates data exploration, many different hypotheses can be explored in a short time without the need for knowledge of R.

In addition, the resulted data could be easily shared using a shiny server. Otherwise, a standalone version of omicsViewer together with designated omics data could be easily created by integrating it with portable R or with docker, which can be shared with collaborators or submitted as supplementary data together with a manuscript.

## Links:

[Vignette](https://mengchen18.github.io/omicsViewer/index.html)

[A live example](http://138.246.235.174:3838/sample-apps/omicsViewer/)

[Docker image](https://hub.docker.com/repository/docker/mengchen18/omicsviewer)

[Share using standalone data package on windows](https://zenodo.org/record/6337392)

[Reference]()

## Further references

OmicsViewer collects a set of useful resources created by other researchers for interpreting omics data. Please make sure to cite the corresponding references when these analyses are used:


- __Fast Gene Set Enrichment Analysis (fGSEA):__

Korotkevich, Gennady, Vladimir Sukhov, Nikolay Budin, Boris Shpak, Maxim N. Artyomov, and Alexey Sergushichev. 2016. “Fast Gene Set Enrichment Analysis.” BioRxiv. bioRxiv. https://doi.org/10.1101/060012.

- __Geneshot:__

Lachmann, Alexander, Brian M. Schilder, Megan L. Wojciechowicz, Denis Torre, Maxim V. Kuleshov, Alexandra B. Keenan, and Avi Ma’ayan. 2019. “Geneshot: Search Engine for Ranking Genes from Arbitrary Text Queries.” Nucleic Acids Research 47 (W1): W571–77.

- __STRING Network analysis:__

Szklarczyk, Damian, Annika L. Gable, Katerina C. Nastou, David Lyon, Rebecca Kirsch, Sampo Pyysalo, Nadezhda T. Doncheva, et al. 2021. “The STRING Database in 2021: Customizable Protein-Protein Networks, and Functional Characterization of User-Uploaded Gene/Measurement Sets.” Nucleic Acids Research 49 (D1): D605–12.

- __SeqLogo:__

Wagih, Omar. 2017. “Ggseqlogo: A Versatile R Package for Drawing Sequence Logos.” Bioinformatics  33 (22): 3645–47.


## Usage

### Quick installation:
```
devtools::install_github("mengchen18/omicsViewer")
```

### Start the viewer
```
library(omicsViewer)
omicsViewer(system.file("extdata", package = "omicsViewer"))
```

### Installation from scratch

OmicsViewer has been submitted to Bioconductor so it can be installed via BiocManager in the future. Currently, if you see an error message in the installation, perhaps because you are missing some dependencies. Please follow the next section.

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


