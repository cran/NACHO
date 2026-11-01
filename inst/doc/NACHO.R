## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  eval = TRUE,
  collapse = TRUE,
  include = TRUE,
  echo = TRUE,
  warning = TRUE,
  message = TRUE,
  fig.align = "center",
  fig.pos = "!h",
  cache = FALSE
)

## -----------------------------------------------------------------------------
knitr::include_graphics(path = "nacho_hex.png")

## -----------------------------------------------------------------------------
# # Install the release from CRAN:
# install.packages("NACHO")
# 
# # Or install the development version from R-universe:
# install.packages(
#   "NACHO",
#   repos = c("https://mcanouil.r-universe.dev", "https://cloud.r-project.org")
# )
# 
# # Or install the latest commit from GitHub with pak:
# # install.packages("pak")
# pak::pak("mcanouil/NACHO")

## -----------------------------------------------------------------------------
# Load NACHO
library(NACHO)

## -----------------------------------------------------------------------------
cat(readLines(system.file("about", "about-nacho.md", package = "NACHO"))[-c(1, 2)], sep = "\n")

## -----------------------------------------------------------------------------
print(citation("NACHO"), "html")

## -----------------------------------------------------------------------------
print(citation("NACHO"), "bibtex")

## -----------------------------------------------------------------------------
data_directory <- system.file("extdata", "GSE178516", package = "NACHO")
nacho <- load_rcc(
  data_directory = data_directory,
  ssheet_csv = file.path(data_directory, "samplesheet.csv"),
  id_colname = "IDFILE",
  instrument = "max",
  n_comp = 5
)
print(nacho)

## -----------------------------------------------------------------------------
head(nacho_counts(nacho), n = c(5, 3))

## -----------------------------------------------------------------------------
head(nacho_samples(nacho)[, c("IDFILE", "title", "CartridgeID", "BD", "FoV")], n = 3)

## -----------------------------------------------------------------------------
head(nacho_probes(nacho), n = 3)
table(nacho_probes(nacho)[["CodeClass"]])

## -----------------------------------------------------------------------------
head(nacho_qc(nacho)[, c("IDFILE", "BD", "FoV", "PCL", "LoD")], n = 3)

## -----------------------------------------------------------------------------
str(nacho_thresholds(), max.level = 1)

## -----------------------------------------------------------------------------
nacho_qc(nacho)[, c("IDFILE", "status", "n_flags", "reason")]

## -----------------------------------------------------------------------------
legacy <- load_rcc(
  data_directory = data_directory,
  ssheet_csv = file.path(data_directory, "samplesheet.csv"),
  id_colname = "IDFILE",
  instrument = "max",
  n_comp = 5,
  preset = "legacy"
)
table(nacho_qc(legacy)[["status"]])

## -----------------------------------------------------------------------------
thresholds_file <- tempfile(fileext = ".yml")
write_thresholds(nacho_thresholds(preset = "legacy"), thresholds_file)
legacy <- normalise(nacho, outliers_thresholds = read_thresholds(thresholds_file))
table(nacho_qc(legacy)[["status"]])

## -----------------------------------------------------------------------------
geo_norm <- normalise(nacho, normalisation_method = "GEO")
glm_norm <- normalise(nacho, normalisation_method = "GLM")
head(nacho_counts(glm_norm, normalised = TRUE), n = c(3, 3))

## -----------------------------------------------------------------------------
background <- normalise(
  nacho,
  background = "mean_2sd",
  background_mode = "subtract"
)

## -----------------------------------------------------------------------------
cleaned <- exclude_outliers(geo_norm)

## -----------------------------------------------------------------------------
mirna_directory <- system.file("extdata", "GSE270837", package = "NACHO")
mirna <- load_rcc(
  data_directory = mirna_directory,
  ssheet_csv = file.path(mirna_directory, "samplesheet.csv"),
  id_colname = "IDFILE",
  instrument = "sprint",
  n_comp = 5
)
print(mirna)
table(nacho_probes(mirna)[["CodeClass"]])

## -----------------------------------------------------------------------------
# thresholds <- nacho_thresholds(instrument = "sprint")
# thresholds[["Ligation_NEG"]] <- c(-Inf, 100)
# mirna <- normalise(mirna, outliers_thresholds = thresholds)

## -----------------------------------------------------------------------------
autoplot(nacho, type = "BD")

## -----------------------------------------------------------------------------
autoplot(geo_norm, type = "RLE")

## -----------------------------------------------------------------------------
autoplot(geo_norm, type = "PCA12")

## -----------------------------------------------------------------------------
# render(nacho, format = "html", output_dir = tempdir())

## -----------------------------------------------------------------------------
# tuned <- visualise(nacho)

## -----------------------------------------------------------------------------
knitr::include_graphics(path = "README-app.png")

## -----------------------------------------------------------------------------
se <- as_summarized_experiment(nacho)
se

