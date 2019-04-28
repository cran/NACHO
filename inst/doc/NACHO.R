## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  eval = TRUE,
  collapse = TRUE,
  results = "asis",
  include = TRUE,
  echo = TRUE,
  warning = TRUE,
  message = TRUE,
  error = TRUE,
  tidy = FALSE,
  crop = TRUE,
  autodep = TRUE,
  fig.align = 'center',
  fig.pos = '!h',
  cache = FALSE
)

## ----logo, echo = FALSE, out.width = "150px"-----------------------------
knitr::include_graphics(path = "nacho_hex.png")

## ---- eval = FALSE-------------------------------------------------------
#  # Install NACHO from CRAN:
#  install.packages("NACHO")
#  
#  # Or the the development version from GitHub:
#  # install.packages("devtools")
#  devtools::install_github("mcanouil/NACHO")

## ----ex1, eval = FALSE---------------------------------------------------
#  library(NACHO)
#  data(GSE74821)
#  visualise(GSE74821)

## ----ex2, results = "hide", message = FALSE, warning = FALSE-------------
library(GEOquery)
# Download data
gse <- getGEO("GSE70970")
# Get phenotypes
targets <- pData(phenoData(gse[[1]]))
getGEOSuppFiles(GEO = "GSE70970", baseDir = tempdir())
# Unzip data
untar(
  tarfile = paste0(tempdir(), "/GSE70970/GSE70970_RAW.tar"), 
  exdir = paste0(tempdir(), "/GSE70970/Data")
)
# Add IDs
targets$IDFILE <- list.files(paste0(tempdir(), "/GSE70970/Data"))

## ----ex3-----------------------------------------------------------------
library(NACHO)
GSE70970_sum <- summarise(
  data_directory = paste0(tempdir(), "/GSE70970/Data"), # Where the data is
  ssheet_csv = targets, # The samplesheet
  id_colname = "IDFILE", # Name of the column that contains the identfiers
  housekeeping_genes = NULL, # Custom list of housekeeping genes
  housekeeping_predict = TRUE, # Predict the housekeeping genes based on the data?
  normalisation_method = "GEO", # Geometric mean or GLM
  n_comp = 5 # Number indicating the number of principal components to compute. 
)

## ---- echo = FALSE, results = "hide"-------------------------------------
unlink(paste0(tempdir(), "/GSE70970"), recursive = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  visualise(GSE70970_sum)

## ----ex5-----------------------------------------------------------------
print(GSE70970_sum[["housekeeping_genes"]])

## ----ex6-----------------------------------------------------------------
my_housekeeping <- GSE70970_sum[["housekeeping_genes"]][-c(1, 2)]
print(my_housekeeping)

## ----ex7-----------------------------------------------------------------
GSE70970_norm <- normalise(
  nacho_object = GSE70970_sum,
  housekeeping_genes = my_housekeeping,
  housekeeping_norm = TRUE,
  normalisation_method = "GEO", 
  remove_outliers = TRUE
)

