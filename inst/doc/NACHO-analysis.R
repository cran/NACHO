## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  eval = TRUE,
  collapse = TRUE,
  # results = "asis",
  include = TRUE,
  echo = TRUE,
  warning = TRUE,
  message = TRUE,
  error = TRUE,
  # tidy = FALSE,
  # crop = TRUE,
  # autodep = TRUE,
  fig.align = "center",
  cache = FALSE
)

## ----logo, echo = FALSE, out.width = "150px"----------------------------------
knitr::include_graphics(path = "nacho_hex.png")

## ---- eval = FALSE------------------------------------------------------------
#  # Install NACHO from CRAN:
#  install.packages("NACHO")
#  
#  # Or the the development version from GitHub:
#  # install.packages("remotes")
#  remotes::install_github("mcanouil/NACHO")

## ---- echo = FALSE, results = "asis"------------------------------------------
cat(readLines(system.file("app", "www", "about-nacho.md", package = "NACHO"))[-c(1, 2)], sep = "\n")

## ---- echo = FALSE, results = "asis"------------------------------------------
print(citation("NACHO"), "html")

## ---- echo = FALSE, comment = ""----------------------------------------------
print(citation("NACHO"), "bibtex")

## -----------------------------------------------------------------------------
library(NACHO)
library(GEOquery, quietly = TRUE, warn.conflicts = FALSE)

## -----------------------------------------------------------------------------
data_directory <- file.path(tempdir(), "GSE70970", "Data")

# Download data
gse <- getGEO("GSE70970")
getGEOSuppFiles(GEO = "GSE70970", baseDir = tempdir())
# Unzip data
untar(
  tarfile = file.path(tempdir(), "GSE70970", "GSE70970_RAW.tar"),
  exdir = data_directory
)
# Get phenotypes and add IDs
targets <- pData(phenoData(gse[[1]]))
targets$IDFILE <- list.files(data_directory)

## -----------------------------------------------------------------------------
GSE70970 <- load_rcc(data_directory, targets, id_colname = "IDFILE")

## -----------------------------------------------------------------------------
library(limma)

## -----------------------------------------------------------------------------
selected_pheno <- GSE70970[["nacho"]][
  j = lapply(unique(.SD), function(x) ifelse(x == "NA", NA, x)),
  .SDcols = c("IDFILE", "age:ch1", "gender:ch1", "chemo:ch1", "disease.event:ch1")
]
selected_pheno <- na.exclude(selected_pheno)

## ---- echo = FALSE------------------------------------------------------------
head(selected_pheno)

## -----------------------------------------------------------------------------
expr_counts <- GSE70970[["nacho"]][
  i = grepl("Endogenous", CodeClass),
  j = as.matrix(
    dcast(.SD, Name ~ IDFILE, value.var = "Count_Norm"),
    "Name"
  ),
  .SDcols = c("IDFILE", "Name", "Count_Norm")
]

## ---- echo = FALSE------------------------------------------------------------
expr_counts[1:5, 1:5]

## ---- eval = FALSE------------------------------------------------------------
#  GSE70970[["nacho"]][
#    i = grepl("Endogenous", CodeClass),
#    j = as.matrix(
#      dcast(.SD, Accession ~ IDFILE, value.var = "Count_Norm"),
#      "Accession"
#    ),
#    .SDcols = c("IDFILE", "Accession", "Count_Norm")
#  ]

## -----------------------------------------------------------------------------
samples_kept <- intersect(selected_pheno[["IDFILE"]], colnames(expr_counts))
expr_counts <- expr_counts[, samples_kept]
selected_pheno <- selected_pheno[IDFILE %in% c(samples_kept)]

## -----------------------------------------------------------------------------
design <- model.matrix(~ `disease.event:ch1`, selected_pheno)

## -----------------------------------------------------------------------------
eBayes(lmFit(expr_counts, design))

## -----------------------------------------------------------------------------
GSE70970[["nacho"]][
  i = grepl("Endogenous", CodeClass),
  j = lapply(unique(.SD), function(x) ifelse(x == "NA", NA, x)),
  .SDcols = c(
    "IDFILE", "Name", "Accession", "Count", "Count_Norm",
    "age:ch1", "gender:ch1", "chemo:ch1", "disease.event:ch1"
  )
][
  Name %in% head(unique(Name), 10)
][
  j = as.data.table(
    coef(summary(lm(
      formula = Count_Norm ~ `disease.event:ch1`,
      data = na.exclude(.SD)
    ))),
    "term"
  ),
  by = c("Name", "Accession")
]

