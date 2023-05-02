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
          # Doesn't harm and needs to be done for HDF5
          as.matrix() |>
          as_tibble(rownames = ".feature") |>
          pivot_longer(names_to = ".sample",
                       values_to = as.character(i),
                       # This is certainly still not the correct way to handle this...
                       names_repair = \(x) str_replace(x, "^(\\d+)$", "expression_\\1"),
                       -`.feature`)

    ) |>
    reduce(\(x, y) inner_join(x, y, by = c(".feature", ".sample"))) |>
    rename_with(\(x) if (length(x) == 1) "expression" else x, starts_with("expression")) |>
    inner_join(row_data, by = ".feature") |>
    inner_join(col_data, by = ".sample")
}

# Experimental:
# Adding / changing methods to handle DESeqDataSet classes
# and provide direct access to normalized counts using the methods
# implemented here

#' @export
setMethod("assayNames", signature(x = "DESeqDataSet"), function(x) {
  c(callNextMethod(), "mrn")
})

#' @importFrom methods callNextMethod
#' @importFrom BiocGenerics counts
#' @export
setMethod("assay", signature(x = "DESeqDataSet", i = "character"), function(x, i, ...) {
  if (i == "mrn") return(counts(x, normalize = TRUE))
  else callNextMethod()
})


