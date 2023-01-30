globalVariables(c(".feature", ".sample"))

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom Biobase fData pData exprs
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select inner_join
#'
#' @export
tidy.ExpressionSet <- function(x, feature_vars = NULL, sample_vars = NULL, ...) {

  row_data <- fData(x) |>
    as_tibble(rownames = ".feature") |>
    select(c(.feature, {{feature_vars}}))

  col_data <- pData(x) |>
    as_tibble(rownames = ".sample") |>
    select(c(.sample, {{sample_vars}}))

  exprs(x) |>
    as_tibble(rownames = ".feature") |>
    pivot_longer(names_to = ".sample",
                 values_to = "expression",
                 -`.feature`) |>
    inner_join(row_data, by = ".feature") |>
    inner_join(col_data, by = ".sample")

}

#' @importFrom SummarizedExperiment assay assayNames rowData colData
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select inner_join
#' @importFrom purrr map reduce
#'
#' @param assay Assay to use as values. By default all assays are returned as different columns.
#'
#' @export
tidy.SummarizedExperiment <- function(x, feature_vars = NULL, sample_vars = NULL, assay = NULL, ...) {

  row_data <- rowData(x) |>
    as_tibble(rownames = ".feature") |>
    select(c(.feature, {{feature_vars}}))

  col_data <- colData(x) |>
    as_tibble(rownames = ".sample") |>
    select(c(.sample, {{sample_vars}}))

  {{assay}} %n% assayNames(x) %n% 1 |>
    map(\(i)
        assay(x, i) |>
          as_tibble(rownames = ".feature") |>
          pivot_longer(names_to = ".sample",
                       values_to = i,
                       -`.feature`)

    ) |>
    reduce(inner_join, by = c(".feature", ".sample")) |>
    inner_join(row_data, by = ".feature") |>
    inner_join(col_data, by = ".sample")
}
