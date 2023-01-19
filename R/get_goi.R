globalVariables(c("feature_id", "gene_symbol", "sample_id"))

#' Function to extract expression data from
#' SummarizedExperiment or ExpressionSet objects
#'
#' @importFrom tibble as_tibble
#' @importFrom tidyr complete separate_rows pivot_longer
#' @importFrom dplyr select filter right_join left_join case_when
#' @import tidyselect
#'
#' @param x A SummarizedExperiment object.
#' @param goi A character vector containing the gene symbols of interest
#' @param ... Arguments to be passed to methods
#'
#' @export
get_goi <- function (x,  goi, ...) {
  UseMethod("get_goi", x)
}

#' @importFrom SummarizedExperiment rowData colData assay assayNames
#' @importFrom dplyr recode
#'
#' @param assay SummarizedExperiment assay to be used
#' @param metadata Optional metadata columns to join (tidyselect).
#'
#' @export
get_goi.SummarizedExperiment <- function(x, goi, metadata = NULL, assay) {

  assay_names <- assayNames(x) %n% "expression"

  value_name <- case_when(
    is.character({{assay}}) ~ as.character({{assay}}),
    is.numeric({{assay}}) ~ assay_names[{{assay}}])

  gene_key <- rowData(x) |>
    gene_symbol_to_id(goi)

  col_data <- colData(x) |>
    # Not convinced that it's the best way to handle an already defined sample_id...
    as_tibble(rownames = "sample_id",
              .name_repair = \(x) recode(x, sample_id = "sample_id.2")) |>
    select(sample_id, {{metadata}})

  assay(x, assay) |>
    as_tibble(rownames = "feature_id") |>
    right_join(gene_key, by = "feature_id") |>
    pivot_longer(names_to = "sample_id",
                 values_to = value_name,
                 -c(feature_id, gene_symbol)) |>
    left_join(col_data, by = "sample_id")
}

#' @importFrom Biobase fData pData exprs
#'
#' @param metadata Optional metadata columns to join (tidyselect).
#'
#' @export
get_goi.ExpressionSet <- function(x, goi, metadata = NULL) {

  gene_key <- fData(x) |>
    gene_symbol_to_id(goi)

  col_data <- pData(x) |>
    # Not convinced that it's the best way to handle an already defined sample_id...
    as_tibble(rownames = "sample_id",
              .name_repair = \(x) recode(x, sample_id = "sample_id.2")) |>
    select(sample_id, {{metadata}})

  exprs(x) |>
    as_tibble(rownames = "feature_id") |>
    right_join(gene_key, by = "feature_id") |>
    pivot_longer(names_to = "sample_id",
                 values_to = "expression",
                 -c(feature_id, gene_symbol)) |>
    left_join(col_data, by = "sample_id")
}

gene_symbol_to_id <- function(x, goi) {
  as_tibble(x, rownames = "feature_id") |>
    select("feature_id", "gene_symbol") |>
    separate_rows("gene_symbol", sep = " /// ") |>
    filter(gene_symbol %in% {{goi}}) |>
    complete(gene_symbol = {{goi}})
}

