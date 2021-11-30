### =========================================================================
### GSE62944 metadata 
### -------------------------------------------------------------------------
###

meta <- data.frame(
    Title = c(
        "Demo data",
        "Expression matrix",
        "Feature annotation",
        "Feature gene set annotation",
        "Sample drug sensitivity data",
        "Sapmle meta data",
        "Sample survival data"
        ),
    Description = c(
        "Example data to be visualized by ExpressionSetViewer",
        "Expression matrix used to generate Demo data",
        "Feature annotations, no gene set information included",
        "Feature annotations, gene set information",
        "Sample meta data, drug sensitivity data",
        "Sample meta data, general information",
        "Sample meta data, survival information"
        ), 
    BiocVersion = "3.10",
    Genome = NA, 
    SourceType = c("RDS", rep("TSV", 6)), 
    SourceUrl = "",
    SourceVersion = "Sep 9 2021",
    Species = "Homo sapiens",
    TaxonomyId = 9606,
    Coordinate_1_based = NA,
    DataProvider = "BayBioMS",
    Maintainer = "Chen Meng <mengchen18@gmail.com>",
    RDataClass = c("ExpressionSet", rep("data.frame", 6)) ,
    DispatchClass = "EpiMetadata",
    RDataPath = c(
        "demo.RDS",
        "expressionMatrix.tsv",
        "featureGeneral.tsv",
        "featureGS.tsv",
        "sampleDrug.tsv",
        "sampleGeneral.tsv",
        "sampleSurv.tsv"
        ),
    Tags = ""
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
