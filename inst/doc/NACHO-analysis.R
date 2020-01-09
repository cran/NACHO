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
  fig.align = 'center',
  fig.pos = '!h',
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
library(dplyr)
library(tidyr)
library(tibble)
library(NACHO)
library(GEOquery)

## -----------------------------------------------------------------------------
data_directory <- file.path(tempdir(), "GSE70970", "Data")

# Download data
gse <- getGEO("GSE70970")
# Get phenotypes
targets <- pData(phenoData(gse[[1]]))
getGEOSuppFiles(GEO = "GSE70970", baseDir = tempdir())
# Unzip data
untar(
  tarfile = file.path(tempdir(), "GSE70970", "GSE70970_RAW.tar"), 
  exdir = data_directory
)
# Add IDs
targets$IDFILE <- list.files(data_directory)

## -----------------------------------------------------------------------------
GSE70970 <- load_rcc(data_directory, targets, id_colname = "IDFILE")

## -----------------------------------------------------------------------------
library(limma)

## -----------------------------------------------------------------------------
selected_pheno <- GSE70970[["nacho"]] %>% 
  select(IDFILE, `age:ch1`, `gender:ch1`, `chemo:ch1`, `disease.event:ch1`) %>% 
  distinct() %>% 
  mutate_all(~ na_if(.x, "NA")) %>% 
  drop_na()

## ---- echo = FALSE------------------------------------------------------------
head(selected_pheno)

## -----------------------------------------------------------------------------
expr_counts <- GSE70970[["nacho"]] %>% 
  filter(grepl("Endogenous", CodeClass)) %>% 
  select(IDFILE, Name, Count_Norm) %>% 
  pivot_wider(names_from = "Name", values_from = "Count_Norm") %>% 
  column_to_rownames("IDFILE") %>% 
  t()

## ---- echo = FALSE------------------------------------------------------------
expr_counts[1:5, 1:5]

## ---- eval = FALSE------------------------------------------------------------
#  GSE70970[["nacho"]] %>%
#    filter(grepl("Endogenous", CodeClass)) %>%
#    select(IDFILE, Accession, Count_Norm) %>%
#    pivot_wider(names_from = "Accession", values_from = "Count_Norm") %>%
#    column_to_rownames("IDFILE") %>%
#    t()

## -----------------------------------------------------------------------------
samples_kept <- intersect(selected_pheno[["IDFILE"]], colnames(expr_counts))
expr_counts <- expr_counts[, samples_kept]
selected_pheno <- filter(selected_pheno, IDFILE %in% !!samples_kept)

## -----------------------------------------------------------------------------
design <- model.matrix(~ `disease.event:ch1`, selected_pheno)

## -----------------------------------------------------------------------------
eBayes(lmFit(expr_counts, design))

## -----------------------------------------------------------------------------
library(purrr)
res <- GSE70970[["nacho"]] %>% 
  # filter( QC PARAMETER) %>% # possible additional QC filters
  filter(grepl("Endogenous", CodeClass)) %>% 
  select(IDFILE, Name, Accession, Count, Count_Norm, `age:ch1`, `gender:ch1`, `chemo:ch1`, `disease.event:ch1`) %>% 
  mutate_all(~ na_if(.x, "NA")) %>% 
  drop_na() %>% 
  group_by(Name, Accession) %>% 
  nest() %>% 
  ungroup() %>% 
  slice(1:10) %>% # the ten first genes
  mutate(
    lm = map(.x = data, .f = function(idata) {# lm or whatever model you want
      as.data.frame(summary(lm(formula = Count_Norm ~ `disease.event:ch1`, idata))$coef)
    })
  ) %>% 
  select(-data)

unnest(res, lm)

