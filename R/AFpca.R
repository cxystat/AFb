#' Adaptive Fisher Method with Principal Components (AFpca)
#' for Trait-Methylation Set Association
#'
#' @param Y Y Phenotype data. It can be a continuous trait or a binary trait.
#' A vector of length \code{n} (number of subjects).
#' @param M A matrix of methylation levels with dimensions \code{n} by \code{K}
#' (\code{n} subjects, \code{K} CpG sites).
#' @param binary Indicator of whether \code{Y} is binary.
#' @param cov Covariates. A matrix with dimensions \code{n}
#' by \code{J} (number of covariates).
#' @param varprop Cutoff for proportion of variance to be expalined.
#' The first Kp principal components (PCs) are chosen such that \code{varprop}
#' of the total variance is explained.
#' @param nperm Number of permutations. Also the starting number of permutations
#' for "step-up" algorithm. Default is \code{1,000}.
#' @param n0 Tuning parameter. Discard the first \code{n0-1} P-values of each sample.
#' @param adapt_perm Whether "step-up" algorithm is used for P-value
#' calculation. If FALSE, function permutes \code{nperm} times and stops. If TRUE,
#' \code{nperm} will be increased \code{10} times each round if P-value <= \code{5/nperm}.
#' Algorithm stops if P-value > \code{5/nperm} or <= \code{cutoff}.
#' @param cutoff Cutoff for "step-up" algorithm.
#' @param seed Specify the seed for permutations.
#'
#' @return An object of "AFpca" class.
#' \describe{
#'  \item{pv}{P-value of AFb test.}
#'  \item{stat}{Test statistic of AFb test.}
#'  \item{indexes}{Indexes of PCs combined into the test
#'   statistic. Indexes are sorted so that P-values are
#'   in ascending order.}
#'   \item{stat_all}{AFb statistics for all permuted samples.}
#'   \item{pv_all}{P-values of AFb statistics for all permuted samples.}
#'   \item{method}{Method used.}
#' }
#'
#' @export
#'
#' @examples
#' Y <- bs_dense$trait
#' methyl <- bs_dense$methyl
#' pos <- bs_dense$pos
#' test <- AFpca(Y, methyl, varprop = 0.9,
#'             binary = TRUE, adapt_perm = TRUE)
#' summary(test)
#'
AFpca <- function(Y, M, binary = FALSE, cov = NULL,
                  varprop = 0.95,
                  nperm = 1000, n0 = 1,
                  adapt_perm = FALSE, cutoff = 2.5e-6,
                  seed = NULL) {

  pcmat <- pcaMat(M, varprop = varprop)

  test <- AF(Y, pcmat, binary = binary, cov = cov,
             nperm = nperm, n0 = n0,
             adapt_perm = adapt_perm, cutoff = cutoff,
             seed = seed)
  result <- list(pv = test$pv, stat = test$stat,
                 indexes = test$indexes,
                 pv_all = test$pv_all,
                 stat_all = test$stat_all,
                 method = "AFpca")
  class(result) <- "AF"
  return(result)
}
