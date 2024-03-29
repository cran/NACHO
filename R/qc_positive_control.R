#' qc_positive_control
#'
#' @param counts [[data.frame]] A `data.frame` with the count data.
#'
#' @keywords internal
#' @usage NULL
#' @noRd
#'
#' @return [[numeric]]
qc_positive_control <- function(counts) {
  if (any(counts[["Count"]] %in% 0)) {
    measured <- log2(counts[["Count"]] + 1)
  } else {
    measured <- log2(counts[["Count"]])
  }
  known <- log2(as.numeric(sub("^[^(]*\\((.*)\\)$", "\\1", counts[["Name"]]))) # plexset value: "32"
  correlation <- summary(stats::lm(measured ~ known))$r.squared
  unname(round(correlation, 5))
}
