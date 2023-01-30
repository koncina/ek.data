#' Filters out expression data from
#' SummarizedExperiment or ExpressionSet objects
#' matching feature criteria
#'
#' @importFrom dplyr filter
#' @import tidyselect
#'
#' @param x A SummarizedExperiment or ExpressionSet object.
#' @param e An expression to filter features based on available variables
#' @param ... Arguments to be passed to methods
#'
#' @export
filter_features <- function (x,  e, ...) {
  UseMethod("filter_features", x)
}

#' @importFrom Biobase fData
#'
#' @export
filter_features.ExpressionSet <- function(x, e, ...) {
  filter(Biobase::fData(x), {{e}}) |>
    (\(f) x[rownames(f),])()
}

#' @importFrom SummarizedExperiment rowData
#'
#' @export
filter_features.SummarizedExperiment <- function(x, e, ...) {
  SummarizedExperiment::rowData(x) |>
    as_tibble(rownames = ".feature") |>
    filter({{e}}) |>
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
