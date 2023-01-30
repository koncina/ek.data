#' Filter out samples of interest
#' based on col data variables (or .sample as the column identifier)
#'
#' @importFrom dplyr filter
#' @import tidyselect
#'
#' @param x A SummarizedExperiment or ExpressionSet object.
#' @param e An expression to filter features based on available variables
#' @param ... Arguments to be passed to methods
#'
#' @export
filter_samples <- function (x,  e, ...) {
  UseMethod("filter_samples", x)
}

#' @importFrom Biobase pData
#'
#' @export
filter_samples.ExpressionSet <- function(x, e, ...) {
  Biobase::pData(x) |>
    as_tibble(rownames = ".sample") |>
    filter({{e}}) |>
    (\(i) x[, i$.sample])()
}

#' @importFrom SummarizedExperiment colData
#'
#' @export
filter_samples.SummarizedExperiment <- function(x, e, ...) {
  SummarizedExperiment::colData(x) |>
    as_tibble(rownames = ".sample") |>
    filter({{e}}) |>
    (\(i) x[, i$.sample])()
}
