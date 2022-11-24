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
  fig.pos = "!h",
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

## ---- message = FALSE---------------------------------------------------------
# Load NACHO
library(NACHO)

## ---- echo = FALSE, results = "asis"------------------------------------------
cat(readLines(system.file("app", "www", "about-nacho.md", package = "NACHO"))[-c(1, 2)], sep = "\n")

## ---- echo = FALSE, results = "asis"------------------------------------------
print(citation("NACHO"), "html")

## ---- echo = FALSE, comment = ""----------------------------------------------
print(citation("NACHO"), "bibtex")

## ----ex1, eval = FALSE--------------------------------------------------------
#  library(NACHO)
#  data(GSE74821)
#  visualise(GSE74821)

## ----ex1-fig, echo = FALSE, out.width = "650px"-------------------------------
knitr::include_graphics(path = "README-visualise.png")

## ----geo-down, echo = FALSE, warning = FALSE, message = FALSE, error = FALSE----
gse <- try(GEOquery::getGEO("GSE70970"), silent = TRUE)
if (inherits(gse, "try-error")) { # when GEOquery is down
  cons <- showConnections(all = TRUE)
  icons <- which(grepl("GSE70970", cons[, "description"])) - 1
  for (icon in icons) close(getConnection(icon))
  cat(
    "Note: `GEOquery` seems to be currently down. Thus, the following code was not executed.\n"
  )
}

## ----ex2, results = "hide", message = FALSE, warning = FALSE, eval = !inherits(gse, "try-error")----
#  library(GEOquery)
#  # Download data
#  gse <- getGEO("GSE70970")
#  getGEOSuppFiles(GEO = "GSE70970", baseDir = tempdir())
#  # Unzip data
#  untar(
#    tarfile = file.path(tempdir(), "GSE70970", "GSE70970_RAW.tar"),
#    exdir = file.path(tempdir(), "GSE70970", "Data")
#  )
#  # Get phenotypes and add IDs
#  targets <- pData(phenoData(gse[[1]]))
#  targets$IDFILE <- list.files(file.path(tempdir(), "GSE70970", "Data"))

## ---- echo = FALSE, message = FALSE, warning = FALSE, eval = !inherits(gse, "try-error")----
#  targets[1:5, unique(c("IDFILE", names(targets)))]

## ----ex3, eval = !inherits(gse, "try-error")----------------------------------
#  GSE70970_sum <- load_rcc(
#    data_directory = file.path(tempdir(), "GSE70970", "Data"), # Where the data is
#    ssheet_csv = targets, # The samplesheet
#    id_colname = "IDFILE", # Name of the column that contains the unique identfiers
#    housekeeping_genes = NULL, # Custom list of housekeeping genes
#    housekeeping_predict = TRUE, # Whether or not to predict the housekeeping genes
#    normalisation_method = "GEO", # Geometric mean or GLM
#    n_comp = 5 # Number indicating how many principal components should be computed.
#  )

## ---- echo = FALSE, results = "hide", eval = !inherits(gse, "try-error")------
#  unlink(file.path(tempdir(), "GSE70970"), recursive = TRUE)

## -----------------------------------------------------------------------------
#  visualise(GSE70970_sum)

## ----ex5, eval = !inherits(gse, "try-error")----------------------------------
#  print(GSE70970_sum[["housekeeping_genes"]])

## ----intext, eval = !inherits(gse, "try-error"), echo = FALSE, results = "asis"----
#  cat(
#    "Let's say _", GSE70970_sum[["housekeeping_genes"]][1],
#    "_ and _", GSE70970_sum[["housekeeping_genes"]][2],
#    "_ are not suitable, therefore, you want to exclude these genes from the normalisation process.",
#    sep = ""
#  )

## ---- eval = !inherits(gse, "try-error")--------------------------------------
#  my_housekeeping <- GSE70970_sum[["housekeeping_genes"]][-c(1, 2)]
#  print(my_housekeeping)

## ----ex7, eval = !inherits(gse, "try-error")----------------------------------
#  GSE70970_norm <- normalise(
#    nacho_object = GSE70970_sum,
#    housekeeping_genes = my_housekeeping,
#    housekeeping_predict = FALSE,
#    housekeeping_norm = TRUE,
#    normalisation_method = "GEO",
#    remove_outliers = TRUE
#  )

## ---- eval = FALSE------------------------------------------------------------
#  autoplot(
#    object = GSE74821,
#    x = "BD",
#    colour = "CartridgeID",
#    size = 0.5,
#    show_legend = TRUE
#  )

## ---- echo = FALSE, results = "asis"------------------------------------------
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

## ----deploy, eval = FALSE-----------------------------------------------------
#  deploy(directory = "/srv/shiny-server", app_name = "NACHO")

## ----app, eval = FALSE--------------------------------------------------------
#  shiny::runApp(system.file("app", package = "NACHO"))

## ----app-fig, echo = FALSE, out.width = "650px"-------------------------------
knitr::include_graphics(path = "README-app.png")

## ---- eval = FALSE------------------------------------------------------------
#  render(
#    nacho_object = GSE74821,
#    colour = "CartridgeID",
#    output_file = "NACHO_QC.html",
#    output_dir = ".",
#    size = 0.5,
#    show_legend = TRUE,
#    clean = TRUE
#  )

## ----print, results = "asis"--------------------------------------------------
print(
  x = GSE74821,
  colour = "CartridgeID",
  size = 0.5,
  show_legend = TRUE,
  echo = TRUE,
  title_level = 3
)

