#' Adaptive Fisher Test
#'
#' @description This function performs Adaptive Fisher (AF) test on GLM-based
#' Score Statistics.
#'
#' @param Y Response Variable. A vector or length \code{n} (the number of samples).
#' @param X Explanatory Variables to be tested. A matrix with dimensions \code{n} by \code{K}
#' (the number of explanatory variables to be tested).
#' @param binary Indicator of whether \code{Y} is binary.
#' @param cov Covariates. A matrix with dimensions \code{n}
#' by \code{J} (the number of explanatory variables that are not to be tested).
#' @param nperm Number of permutations. Also the starting number of
#' permutations for "step-up" algorithm. Default is \code{1,000}.
#' @param adapt_perm Whether "step-up" algorithm is used for P-value
#' calculation. If FALSE, function permutes \code{nperm} times and stops.
#' If TRUE, \code{nperm} will be increased \code{10} times each round if P-value
#' <= \code{5/nperm}. Algorithm stops if P-value > \code{5/nperm} or <= \code{cutoff}.
#' @param cutoff Cutoff for "step-up" algorithm.
#' @param n0 Tuning parameter. Discard the first \code{n0-1} P-values of each sample.
#' @param seed Specify the seed for permutations.
#'
#' @return An object of "AF" class.
#' \describe{
#'  \item{pv}{P-value of AF test.}
#'  \item{stat}{AF statisitc.}
#'  \item{indexes}{Indexes of P-values combined into the test
#'   statistic. Indexes are sorted so that P-values are
#'   in ascending order.}
#'   \item{stat_all}{AF statistics for all permuted samples.}
#'   \item{pv_all}{P-values of AF statistics for all permuted samples.}
#'   \item{method}{Method used.}
#' }
#'
#' @export
#'
#' @seealso \code{\link{set.seed}}
#'
#' @examples
#' Y <- bs_dense$trait
#' methyl <- bs_dense$methyl
#' pos <- bs_dense$pos
#' test <- AF(Y, methyl, binary = TRUE)
#' summary(test)
#'
AF <- function(Y, X, binary = FALSE, cov = NULL, nperm = 1e3,
               adapt_perm = FALSE, cutoff = 2.5e-6, n0 = 1,
               seed = NULL){

  score <- perm_score(Y, X, binary = binary, cov = cov,
                      nperm = nperm, seed = seed)
  test <- AF_combine(score$pvs, n0 = n0)

  if (adapt_perm) {
    pvs.all <- score$pvs
    while(test$pv <= 5/nperm & test$pv > cutoff) {
      nperm.add <- nperm * 9
      nperm <- nperm + nperm.add
      print(paste("start", nperm, "permutations"))
      score.add <- perm_score(Y, X, binary = binary, cov = cov,
                              nperm = nperm.add)
      pvs.all <- cbind(pvs.all, score.add$pvs[,-1])
      test <- AF_combine(pvs.all, n0 = n0)
    }
  }

  result <- list(pv = test$pv, stat = test$stat,
                 indexes = test$indexes,
                 pv_all = test$pv_all,
                 stat_all = test$stat_all,
                 method = "AF")

  class(result) <- "AF"

  return(result)
}
