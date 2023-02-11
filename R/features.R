#' Filters out expression data from
#' SummarizedExperiment or ExpressionSet objects
#' matching feature criteria
#'
#' @importFrom dplyr filter
#' @import tidyselect
#'
#' @param x A SummarizedExperiment or ExpressionSet object.
#' @param ... Arguments passed to filter
#'
#' @export
filter_features <- function (x, ...) {
  UseMethod("filter_features", x)
}

#' @importFrom Biobase fData
#'
#' @export
filter_features.ExpressionSet <- function(x, ...) {
  filter(Biobase::fData(x), ...) |>
    (\(f) x[rownames(f),])()
}

#' @importFrom SummarizedExperiment rowData
#'
#' @export
filter_features.SummarizedExperiment <- function(x, ...) {
  SummarizedExperiment::rowData(x) |>
    as_tibble(rownames = ".feature") |>
    filter(...) |>
    (\(f) x[f$.feature,])()
}

#' "Also known as" splits collapsed list of gene identifiers and tests if
#' the expression is part of it (could be replaced by `str_detect()` or similar)
#'
#' @importFrom stringr str_split
#' @importFrom purrr map_lgl
#'
#' @param x Collapsed list of genes (string)
#' @param y Value to search for
#' @param sep the separator to be used (defaults to ` */+ *`)
#'
#' @export
aka <- function(x, y, sep = " */+ *") {
  str_split(x, pattern = sep) |>
    map_lgl(\(x) y %in% x)
}

#' "Also known as" splits collapsed list of gene identifiers and tests if
#' the expression is part of it (could be replaced by `str_detect()` or similar)
#'
#' @param x Collapsed list of genes (string)
#' @param y Value to search for
#'
#' @export
`%aka%` <- function(x, y) {
  aka(x, y, sep =  " */+ *")
}


#' Extract available features as a tibble
#' A wrapper around rowData() for SummarizedExperiment
#' and fData() for ExpressionSet objects
#' returning a tibble
#'
#' @param x A SummarizedExperiment or ExpressionSet object.
#' @param ... Arguments to be passed to methods
#'
#' @export
feature_data <- function (x, ...) {
  UseMethod("feature_data", x)
}


#' @importFrom Biobase fData
#' @importFrom tibble as_tibble
#'
#' @export
feature_data.ExpressionSet <- function(x, ...) {
  fData(x) |>
    as_tibble(rownames = ".feature")
}

#' @importFrom SummarizedExperiment rowData
#' @importFrom tibble as_tibble
#'
#' @export
feature_data.SummarizedExperiment <- function(x, ...) {
  SummarizedExperiment::rowData(x) |>
    as_tibble(rownames = ".feature")
}

