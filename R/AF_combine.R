#' AF Combination
#'
#' @description This function combines P-values using
#' Adaptive Fisher (AF) method.
#'
#' @param p P-values to be combined. A matrix with dimensions \code{K} by \code{N}.
#' If an object of P-values from perm_score is used,
#' \code{N} is the number of permutations plus 1.
#'
#' @param log Indicator of whether P-values are on the log scale.
#'
#' @param weight Weights given to the P-values. A vector with dimension
#' \code{K}. Constant weights are used if it is not specified.
#'
#' @param n0 Tuning parameter. Discard the first \code{n0-1} P-values of each
#' sample.
#'
#' @return A list object.
#' \describe{
#'   \item{pv}{P-value of AF test.}
#'   \item{stat}{AF statistic.}
#'   \item{indexes}{Indexes of P-values combined into the test
#'   statistic. Indexes are sorted so that P-values are
#'   in ascending order.}
#'   \item{stat_all}{AF statistics for all permuted samples.}
#'   \item{pv_all}{P-values of AF statistics for all permuted samples.}
#' }
#'
#'
#' @export
#'
#' @examples
#' # Combine P-values of normally distributed test statistics
#' U <- matrix(rnorm(10 * 100), ncol=100)
#' p <- 2 * (1 - pnorm(abs(U)))
#' wt <- (1:10)/55
#' test <- AF_combine(p, weight = wt)
#'
AF_combine <- function(p, log = FALSE, weight = NULL, n0 = 1) {

  N <- nrow(p)
  T <- ncol(p)

  if (is.null(weight)) {weight <- rep(1, N)}

  if (log == FALSE) {
    r <- weight * log(p)
  } else {
    r <- weight * p
  }

  if (N == 1) {
    stat <- as.vector(p)
    combine.index <- 1
    p.AF <- rank(stat, ties.method = "max")/T
  } else {
    s <- apply(apply(r, 2, sort), 2, cumsum)
    p.perm <- t(apply(s, 1, rank, ties.method="max")/T)
    stat <- apply(p.perm[n0:N,], 2, min)
    p.AF <- rank(stat, ties.method = "max")/T

    # P-values combined in AF statistic
    obs.order <- order(r[,1])
    combine.index <- obs.order[1:which.min(p.perm[n0:N, 1])]
  }

  result <- list(pv = p.AF[1], stat = stat[1],
                 indexes = combine.index,
                 stat_all = stat, pv_all = p.AF)

  return(result)
}
