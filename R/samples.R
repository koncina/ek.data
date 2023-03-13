#' Filter out samples of interest
#' based on col data variables (or .sample as the column identifier)
#'
#' @importFrom dplyr filter
#' @import tidyselect
#'
#' @param x A SummarizedExperiment or ExpressionSet object.
#' @param ... Arguments passed to filter
#'
#' @export
filter_samples <- function (x, ...) {
  UseMethod("filter_samples", x)
}

#' @importFrom Biobase pData
#'
#' @export
filter_samples.ExpressionSet <- function(x, ...) {
  Biobase::pData(x) |>
    as_tibble(rownames = ".sample") |>
    filter(...) |>
    (\(i) x[, i$.sample])()
}

#' @importFrom SummarizedExperiment colData
#'
#' @export
filter_samples.SummarizedExperiment <- function(x, ...) {
  SummarizedExperiment::colData(x) |>
    as_tibble(rownames = ".sample") |>
    filter(...) |>
    (\(i) x[, i$.sample])()
}


#' Extract available samples as a tibble
#' A wrapper around colData() for SummarizedExperiment
#' and pData() for ExpressionSet objects
#' returning a tibble
#'
#' @param x A SummarizedExperiment or ExpressionSet object.
#' @param ... Arguments to be passed to methods
#'
#' @export
sample_data <- function (x, ...) {
  UseMethod("sample_data", x)
}

#' @importFrom Biobase pData
#' @importFrom tibble as_tibble
#'
#' @export
sample_data.ExpressionSet <- function(x, ...) {
  pData(x) |>
    as_tibble(rownames = ".sample", ...)
}

#' @importFrom SummarizedExperiment colData
#' @importFrom tibble as_tibble
#'
#' @export
sample_data.SummarizedExperiment <- function(x, ...) {
  SummarizedExperiment::colData(x) |>
    as_tibble(rownames = ".sample", ...)
}
