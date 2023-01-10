NULL

#' Log transform a vector adding and add a constant
#'
#' @param x A numeric vector.
#' @param i Numeric value to add before performing the log-transformation. By default half of the smallest non-zero value is added.
#' @param base A positive or complex number: the base with respect to which logarithms are computed. Defaults to 2.
#'
#' @export
log_transform <- function(x, i = min(x[x > 0]) / 2, base = 2) {
  stopifnot(is.numeric(i) && length(i) == 1)
  log(x + i, base = base)
}
