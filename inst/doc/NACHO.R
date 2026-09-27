## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  eval = TRUE,
  collapse = TRUE,
  include = TRUE,
  echo = TRUE,
  warning = TRUE,
  message = TRUE,
  error = TRUE,
  fig.align = "center",
  fig.pos = "!h",
  cache = FALSE
)

## -----------------------------------------------------------------------------
knitr::include_graphics(path = "nacho_hex.png")

## -----------------------------------------------------------------------------
# # Install NACHO from CRAN:
# install.packages("NACHO")
# 
# # Or the development version from GitHub:
# # install.packages("pak")
# pak::pak("mcanouil/NACHO")

## -----------------------------------------------------------------------------
# Load NACHO
library(NACHO)

## -----------------------------------------------------------------------------
cat(readLines(system.file("app", "www", "about-nacho.md", package = "NACHO"))[-c(1, 2)], sep = "\n")

## -----------------------------------------------------------------------------
print(citation("NACHO"), "html")

## -----------------------------------------------------------------------------
print(citation("NACHO"), "bibtex")

## -----------------------------------------------------------------------------
# library(NACHO)
# data(GSE74821)
# visualise(GSE74821)

## -----------------------------------------------------------------------------
knitr::include_graphics(path = "README-visualise.png")

## -----------------------------------------------------------------------------
gse <- try(GEOquery::getGEO("GSE70970"), silent = TRUE)
if (inherits(gse, "try-error")) { # when GEOquery is down
  cons <- showConnections(all = TRUE)
  icons <- which(grepl("GSE70970", cons[, "description"])) - 1
  for (icon in icons) close(getConnection(icon))
  cat(
    "Note: `GEOquery` seems to be currently down. Thus, the following code was not executed.\n"
  )
}

## -----------------------------------------------------------------------------
library(GEOquery)
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
rcc_files <- list.files(data_directory, pattern = "\\.RCC(\\.gz)?$", ignore.case = TRUE)
targets$IDFILE <- rcc_files[match(targets$geo_accession, sub("_.*", "", rcc_files))]
targets <- targets[!is.na(targets$IDFILE), ]
# Keep the samples measured with the same CodeSet
codeset <- vapply(
  X = file.path(data_directory, targets$IDFILE),
  FUN = function(file) {
    header <- grep("^GeneRLF,", readLines(file, n = 40), value = TRUE)
    if (length(header) == 0) NA_character_ else header[1]
  },
  FUN.VALUE = character(1)
)
targets <- targets[codeset %in% "GeneRLF,NS_H_miR_1.4", ]

## -----------------------------------------------------------------------------
targets[1:5, unique(c("IDFILE", names(targets)))]

## -----------------------------------------------------------------------------
GSE70970_sum <- load_rcc(
  data_directory = data_directory, # Where the data is
  ssheet_csv = targets, # The samplesheet
  id_colname = "IDFILE", # Name of the column that contains the unique identifiers
  housekeeping_genes = NULL, # Custom list of housekeeping genes
  housekeeping_predict = TRUE, # Whether or not to predict the housekeeping genes
  normalisation_method = "GEO", # Geometric mean or GLM
  n_comp = 5 # Number indicating how many principal components should be computed.
)

## -----------------------------------------------------------------------------
unlink(file.path(tempdir(), "GSE70970"), recursive = TRUE)

## -----------------------------------------------------------------------------
# visualise(GSE70970_sum)

## -----------------------------------------------------------------------------
print(GSE70970_sum[["housekeeping_genes"]])

## -----------------------------------------------------------------------------
cat(
  "Let's say _", GSE70970_sum[["housekeeping_genes"]][1],
  "_ and _", GSE70970_sum[["housekeeping_genes"]][2],
  "_ are not suitable, therefore, you want to exclude these genes from the normalisation process.",
  sep = ""
)

## -----------------------------------------------------------------------------
my_housekeeping <- GSE70970_sum[["housekeeping_genes"]][-c(1, 2)]
print(my_housekeeping)

## -----------------------------------------------------------------------------
GSE70970_norm <- normalise(
  nacho_object = GSE70970_sum,
  housekeeping_genes = my_housekeeping,
  housekeeping_predict = FALSE,
  housekeeping_norm = TRUE,
  normalisation_method = "GEO",
  remove_outliers = TRUE
)

## -----------------------------------------------------------------------------
# autoplot(
#   object = GSE74821,
#   x = "BD",
#   colour = "CartridgeID",
#   size = 0.5,
#   show_legend = TRUE
# )

## -----------------------------------------------------------------------------
metrics <- c(
  "BD" = "Binding Density",
  "FoV" = "Imaging",
  "PCL" = "Positive Control Linearity",
  "LoD" = "Limit of Detection",
  "Positive" = "Positive Controls",
  "Negative" = "Negative Controls",
  "Housekeeping" = "Housekeeping Genes",
  "PN" = "Positive Controls vs. Negative Controls",
  "ACBD" = "Average Counts vs. Binding Density",
  "ACMC" = "Average Counts vs. Median Counts",
  "PCA12" = "Principal Component 1 vs. 2",
  "PCAi" = "Principal Component scree plot",
  "PCA" = "Principal Components planes",
  "PFNF" = "Positive Factor vs. Negative Factor",
  "HF" = "Housekeeping Factor",
  "NORM" = "Normalisation Factor"
)

for (imetric in seq_along(metrics)) {
  cat("\n\n###", metrics[imetric], "\n\n")
  print(autoplot(object = GSE74821, x = names(metrics[imetric])))
  cat("\n")
}

## -----------------------------------------------------------------------------
# deploy(directory = "/srv/shiny-server", app_name = "NACHO")

## -----------------------------------------------------------------------------
# shiny::runApp(system.file("app", package = "NACHO"))

## -----------------------------------------------------------------------------
knitr::include_graphics(path = "README-app.png")

## -----------------------------------------------------------------------------
# render(
#   nacho_object = GSE74821,
#   colour = "CartridgeID",
#   output_file = "NACHO_QC.html",
#   output_dir = ".",
#   size = 0.5,
#   show_legend = TRUE,
#   clean = TRUE
# )

## -----------------------------------------------------------------------------
print(
  x = GSE74821,
  colour = "CartridgeID",
  size = 0.5,
  show_legend = TRUE,
  echo = TRUE,
  title_level = 3
)

