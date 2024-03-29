#' (re)Normalise a "nacho" object
#'
#' This function creates a list in which your settings, the raw counts and normalised counts are stored,
#' using the result from a call to [`load_rcc()`].
#'
#' @param nacho_object [[list]] A list object of class `"nacho"` obtained
#'   from [`load_rcc()`] or [`normalise()`].
#' @inheritParams load_rcc
#' @param remove_outliers [[logical]] A boolean to indicate if outliers should be excluded.
#' @param outliers_thresholds [[list]] List of thresholds to exclude outliers.
#'
#' @details Outliers definition (`remove_outliers = TRUE`):
#'
#'  * Binding Density (`BD`) < 0.1
#'  * Binding Density (`BD`) > 2.25
#'  * Field of View (`FoV`) < 75
#'  * Positive Control Linearity (`PCL`) < 0.95
#'  * Limit of Detection (`LoD`) < 2
#'  * Positive normalisation factor (`Positive_factor`) < 0.25
#'  * Positive normalisation factor (`Positive_factor`) > 4
#'  * Housekeeping normalisation factor (`house_factor`) < 1/11
#'  * Housekeeping normalisation factor (`house_factor`) > 11
#'
#' @return [[list]] A list containing parameters and data.
#' \describe{
#'   \item{`access`}{[[character]] Value passed to [`load_rcc()`] in `id_colname`.}
#'   \item{`housekeeping_genes`}{[[character]] Value passed to [`load_rcc()`] or [`normalise()`].}
#'   \item{`housekeeping_predict`}{[[logical]] Value passed to [`load_rcc()`].}
#'   \item{`housekeeping_norm`}{[[logical]] Value passed to [`load_rcc()`] or [`normalise()`].}
#'   \item{`normalisation_method`}{[[character]] Value passed to [`load_rcc()`] or [`normalise()`].}
#'   \item{`remove_outliers`}{[[logical]] Value passed to [`normalise()`].}
#'   \item{`n_comp`}{[[numeric]] Value passed to [`load_rcc()`].}
#'   \item{`data_directory`}{[[character]] Value passed to [`load_rcc()`].}
#'   \item{`pc_sum`}{[[data.frame]] A `data.frame` with `n_comp` rows and four columns:
#'     "Standard deviation", "Proportion of Variance", "Cumulative Proportion" and "PC".}
#'   \item{`nacho`}{[[data.frame]] A `data.frame` with all columns from the sample sheet `ssheet_csv`
#'     and all computed columns, *i.e.*, quality-control metrics and counts, with one sample per row.}
#'   \item{`outliers_thresholds`}{[[list]] A `list` of the quality-control thresholds used.}
#'   \item{`raw_counts`}{[[data.frame]] Raw counts with probes as rows and samples as columns.
#'     With `"CodeClass"` (first column), the type of the probes and
#'     `"Name"` (second column), the Name of the probes.}
#'   \item{`normalised_counts`}{[[data.frame]] Normalised counts with probes as rows and samples as columns.
#'     With `"CodeClass"` (first column)), the type of the probes and
#'     `"Name"` (second column), the name of the probes.}
#' }
#'
#' @export
#'
#' @examples
#'
#' data(GSE74821)
#' GSE74821_norm <- normalise(
#'   nacho_object = GSE74821,
#'   housekeeping_norm = TRUE,
#'   normalisation_method = "GEO",
#'   remove_outliers = TRUE
#' )
#'
#' if (interactive()) {
#'   library(GEOquery)
#'   library(NACHO)
#'
#'   # Import data from GEO
#'   gse <- GEOquery::getGEO(GEO = "GSE74821")
#'   targets <- Biobase::pData(Biobase::phenoData(gse[[1]]))
#'   GEOquery::getGEOSuppFiles(GEO = "GSE74821", baseDir = tempdir())
#'   utils::untar(
#'     tarfile = file.path(tempdir(), "GSE74821", "GSE74821_RAW.tar"),
#'     exdir = file.path(tempdir(), "GSE74821")
#'   )
#'   targets$IDFILE <- list.files(
#'     path = file.path(tempdir(), "GSE74821"),
#'     pattern = ".RCC.gz$"
#'   )
#'   targets[] <- lapply(X = targets, FUN = iconv, from = "latin1", to = "ASCII")
#'   utils::write.csv(
#'     x = targets,
#'     file = file.path(tempdir(), "GSE74821", "Samplesheet.csv")
#'   )
#'
#'   # Read RCC files and format
#'   nacho <- load_rcc(
#'     data_directory = file.path(tempdir(), "GSE74821"),
#'     ssheet_csv = file.path(tempdir(), "GSE74821", "Samplesheet.csv"),
#'     id_colname = "IDFILE"
#'   )
#'
#'   # (re)Normalise data by removing outliers
#'   nacho_norm <- normalise(
#'     nacho_object = nacho,
#'     remove_outliers = TRUE
#'   )
#'
#'   # (re)Normalise data with "GLM" method and removing outliers
#'   nacho_norm <- normalise(
#'     nacho_object = nacho,
#'     normalisation_method = "GLM",
#'     remove_outliers = TRUE
#'   )
#' }
#'
normalise <- function(
  nacho_object,
  housekeeping_genes = nacho_object[["housekeeping_genes"]],
  housekeeping_predict = nacho_object[["housekeeping_predict"]],
  housekeeping_norm = nacho_object[["housekeeping_norm"]],
  normalisation_method = nacho_object[["normalisation_method"]],
  n_comp = nacho_object[["n_comp"]],
  remove_outliers = nacho_object[["remove_outliers"]],
  outliers_thresholds = nacho_object[["outliers_thresholds"]]
) {
  if (missing(nacho_object)) {
    stop(
      '[NACHO] "nacho_object" is missing, results from "load_rcc()" and/or "normalise()" is mandatory!'
    )
  }
  if (!attr(nacho_object, "RCC_type") %in% c("n1", "n8")) {
    stop('[NACHO] RCC type must be either "n1" or "n8"!')
  }
  mandatory_fields <- c(
    "access",
    "housekeeping_genes",
    "housekeeping_predict",
    "housekeeping_norm",
    "normalisation_method",
    "remove_outliers",
    "n_comp",
    "data_directory",
    "pc_sum",
    "nacho",
    "outliers_thresholds"
  )
  if (!all(mandatory_fields %in% names(nacho_object))) {
    stop(
      '[NACHO] Mandatory fields are missing in "', substitute(nacho_object), '"!\n',
      '  "load_rcc()" must be called before "normalise()".'
    )
  }

  id_colname <- nacho_object[["access"]]
  type_set <- attr(nacho_object, "RCC_type")

  params_changed <- c(
    "housekeeping_genes" = !isTRUE(all.equal(sort(nacho_object[["housekeeping_genes"]]), sort(housekeeping_genes))),
    "housekeeping_predict" = nacho_object[["housekeeping_predict"]] != housekeeping_predict,
    "housekeeping_norm" = nacho_object[["housekeeping_norm"]] != housekeeping_norm,
    "normalisation_method" = nacho_object[["normalisation_method"]] != normalisation_method,
    "n_comp" = nacho_object[["n_comp"]] != n_comp,
    "remove_outliers" = nacho_object[["remove_outliers"]] != remove_outliers,
    "outliers_thresholds" = !isTRUE(all.equal(nacho_object[["outliers_thresholds"]], outliers_thresholds))
  )

  if (all(!params_changed)) {
    message(
      '[NACHO] Nothing was done. Parameters in "normalise()", were the same as in "', substitute(nacho_object), '".'
    )
    return(nacho_object)
  } else {
    message(
      '[NACHO] Normalising "', substitute(nacho_object), '" with new value for parameters:\n',
      paste(
        paste0("  - ", names(params_changed[which(params_changed)]), " = ", params_changed[which(params_changed)]),
        collapse = "\n"
      )
    )
  }

  if (remove_outliers & !nacho_object[["remove_outliers"]]) {
    nacho_object[["outliers_thresholds"]] <- outliers_thresholds
    nacho_object <- check_outliers(nacho_object)

    if (any(nacho_object[["nacho"]][, "is_outlier"]) | any(params_changed)) {
      nacho_object <- qc_rcc(
        data_directory = nacho_object[["data_directory"]],
        nacho_df = nacho_object[["nacho"]][which(!nacho_object[["nacho"]][, "is_outlier"]), ],
        id_colname = id_colname,
        housekeeping_genes = housekeeping_genes,
        housekeeping_predict = housekeeping_predict,
        housekeeping_norm = housekeeping_norm,
        normalisation_method = normalisation_method,
        n_comp = nacho_object[["n_comp"]]
      )
    }
    nacho_object[["remove_outliers"]] <- remove_outliers
  } else {
    message("[NACHO] Outliers have already been removed!")

    if (any(params_changed)) {
      nacho_object <- qc_rcc(
        data_directory = nacho_object[["data_directory"]],
        nacho_df = nacho_object[["nacho"]],
        id_colname = id_colname,
        housekeeping_genes = housekeeping_genes,
        housekeeping_predict = housekeeping_predict,
        housekeeping_norm = housekeeping_norm,
        normalisation_method = normalisation_method,
        n_comp = nacho_object[["n_comp"]]
      )
    }
  }

  nacho_object[["outliers_thresholds"]] <- outliers_thresholds

  nacho_object[["nacho"]][["Count_Norm"]] <- normalise_counts(
    data = nacho_object[["nacho"]],
    housekeeping_norm = housekeeping_norm
  )

  if (!"RCC_type" %in% names(attributes(nacho_object))) {
    attributes(nacho_object) <- c(attributes(nacho_object), RCC_type = type_set)
  }

  message(paste(
    "[NACHO] Returning a list.",
    "  $ access              : character",
    "  $ housekeeping_genes  : character",
    "  $ housekeeping_predict: logical",
    "  $ housekeeping_norm   : logical",
    "  $ normalisation_method: character",
    "  $ remove_outliers     : logical",
    "  $ n_comp              : numeric",
    "  $ data_directory      : character",
    "  $ pc_sum              : data.frame",
    "  $ nacho               : data.frame",
    "  $ outliers_thresholds : list",
    sep = "\n"
  ))

  class(nacho_object) <- "nacho"

  nacho_object
}


#' @export
#' @rdname normalise
#' @usage NULL
normalize <- normalise
