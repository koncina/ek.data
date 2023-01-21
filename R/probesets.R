
#' Identify probeset exhibiting the highest expression variance
#' and being the most specific (based on probeset IDs)
#'
#' @importFrom tibble tibble rowid_to_column
#' @importFrom dplyr mutate select group_by summarise
#' @importFrom stats var
#'
#' @param x An ExpressionSet or SummarizedExperiment object
#' @param gene_id_var The gene identifier to group by and associated to
#' multiple probesets (row IDs)
#'
#' @export
select_probeset <- function (x, gene_id_var = gene_symbol, ...) {
  UseMethod("select_probeset", x)
}


#' @importFrom SummarizedExperiment rowData assay
#'
#' @param assay SummarizedExperiment assay to be used
#'
#' @export
select_probeset.SummarizedExperiment <- function(x, gene_id_var, assay) {
  apply(assay(x, i = {{assay}}), 1, var) |>
    enframe(".feature_id", "var") |>
    inner_join(as_tibble(rowData(x),
                         rownames = ".feature_id"),
               by = ".feature_id") |>
    .select_probeset({{gene_id_var}})
}


#' @importFrom Biobase exprs fData
#'
#' @export
select_probeset.ExpressionSet <- function(x, gene_id_var) {
  apply(exprs(x), 1, var) |>
    enframe(".feature_id", "var") |>
    inner_join(as_tibble(fData(x),
                         rownames = ".feature_id"),
               by = ".feature_id") |>
    .select_probeset({{gene_id_var}})
}


#' @importFrom dplyr group_by
#'
#' @export
select_probeset.data.frame <- function(x, gene_id_var, probeset_id_var = feature_id, expression_var) {
  rename(x, ".feature_id" = {{probeset_id_var}}) |>
    group_by({{gene_id_var}}, .feature_id) |>
    summarise(var = var({{expression}}), .groups = "drop") |>
    .select_probeset({{gene_id_var}})
}

.select_probeset <- function(x, gene_id_var) {
  group_by(x, {{gene_id_var}}) |>
    mutate(
      is_unspecific = str_detect(.feature_id, "[xs]_at$"),
      only_unspecific = all(is_unspecific),
      is_candidate = !xor(is_unspecific, only_unspecific),
      selected_probeset = var == max(var[is_candidate])) |>
    ungroup() |>
    select(.feature_id, {{gene_id_var}}, is_unspecific:selected_probeset) |>
    column_to_rownames(".feature_id")

}
