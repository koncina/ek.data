globalVariables(c(".feature", "rowid", "n", "probeset_penalty"))

#' Identify probeset exhibiting the highest expression variance
#' and being the most specific (based on probeset IDs)
#'
#' @importFrom stats var
#'
#' @param x An ExpressionSet or SummarizedExperiment object
#' @param feature_var The gene identifier to group by and associated to
#' multiple probesets (row IDs)
#' @param most_specific remove suboptimal affymetrix probesets (`at > a_at > s_at > x_at`). Defaults to `TRUE`.
#' @param ... Arguments to be passed to methods
#'
#' @export
filter_probeset <- function (x, feature_var, most_specific = TRUE, ...) {
  UseMethod("filter_probeset", x)
}

#' @importFrom Biobase fData exprs
#'
#' @export
filter_probeset.ExpressionSet <- function(x, feature_var,
                                          most_specific = TRUE, ...) {
  .filter_probeset(x, {{feature_var}},
                   most_specific,
                   row_data_f = fData,
                   assay_f = exprs)
}

#' @importFrom SummarizedExperiment rowData assay
#'
#' @export
filter_probeset.SummarizedExperiment <- function(x, feature_var,
                                                 most_specific = TRUE, assay, ...) {
  .filter_probeset(x, {{feature_var}},
                   most_specific,
                   row_data_f = rowData,
                   assay_f = \(x) assay(x, i = {{assay}}))
}

#' @importFrom assertr verify
#' @importFrom tibble as_tibble rowid_to_column enframe
#' @importFrom dplyr select mutate add_count recode across filter if_else group_by ungroup
#' @importFrom stringr str_length str_extract str_trim
#' @importFrom tidyr separate_rows
#'
.filter_probeset <- function(x, id, most_specific = TRUE, row_data_f, assay_f) {

  n_features_step1 <- nrow(x)

  row_data <- row_data_f(x) |>
    as_tibble(rownames = ".feature") |>
    select(.feature, {{id}}) |>
    rowid_to_column() |>
    separate_rows({{id}}, sep =  " */+ *") |>
    mutate(across({{id}}, str_trim)) |>
    add_count(rowid) |>
    verify(!(str_length({{id}}) == 0 & n > 1)) |>
    mutate({{id}} := if_else(str_length({{id}}) == 0,
                                 as.character(rowid), {{id}}))

  if (isTRUE(most_specific)) {
    row_data <-  mutate(row_data,
                        probeset_penalty = str_extract(.feature,
                                                       "(_[asx])?_at$"),
                        across(probeset_penalty, recode,
                               `_at` = 0,
                               `_a_at` = 1,
                               `_s_at` = 2,
                               `_x_at` = 3)) |>
      group_by({{id}}) |>
      filter(probeset_penalty == min(probeset_penalty) |
               is.na(probeset_penalty)) |>
      ungroup() |>
      select(-probeset_penalty)
  }

  x <- x[row_data$.feature,]
  n_features_step2 <- nrow(x)

  if (n_features_step1 != n_features_step2)
    message("Removed ", n_features_step1 - n_features_step2,
            " suboptimal affymetrix probesets")

  row_data <- assay_f(x) |>
    apply(1, var) |>
    enframe(".feature", "var") |>
    inner_join(row_data, by = ".feature") |>
    group_by({{id}}) |>
    filter(var == min(var))

  x <- x[row_data$.feature,]
  n_features_step3 <- nrow(x)

  if (n_features_step3 != n_features_step2)
    message("Removing ", n_features_step2 - n_features_step3,
            " feature(s): keeping features ensuring most variance")


  if (n_features_step1 != n_features_step3)
    message("Keeping ", n_features_step3, " out of ",
            n_features_step1, " initial features")

  x
}
