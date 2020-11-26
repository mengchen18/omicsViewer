# user level functions:
# - general information
# - t-test
# - PCA
# - gene set preparation
# - prepare string

# auxi functions:
# ORA
# fgsea

# tutorial structure

# # introduction
#  - What is ExpresssionSetViewer
#  - What is ExpressionSet
#  - Who should use it, what it used for/not used for
# 
# # Data preparation (backend user)
#  - Extention of ExpressionSetViewer
#     - three level column headers
#     - require row/colnames
#  - Reserved and recommanded keywords in headers
#  - Data preparation
#   - General (feature and sample)
#   - GS (feature)
#   - t-test (feature)
#   - PCA (feature and sample)
#   - StringDB (feature)
#   - Surv (sample)
#   - others
#     - elastic net/PLS/clustering results (free keywords)
#  - working pipeline
#   - compare two groups, enable enrichment analysis
#   - multiple group comparison


setwd("/media/share_baybioms/Projects/008_Bioinformatics/B032_ExpressionSetViewer/Git/bookdown/")
library(bookdown)
render_book(input = "./index.Rmd", output_format = "bookdown::gitbook")
