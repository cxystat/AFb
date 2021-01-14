#' Smoothed Adaptive Fisher (AFb) Test for Trait-Methylation Set Association
#'
#' @param Y Y Phenotype data. It can be a continuous trait or a binary trait.
#' A vector or length n (number of subjects).
#' @param M A matrix of methylation levels. A matrix with dimensions n by K
#' (n subjects, K methylation sites). Each row for a subject, and each column
#' for a methylated CpG site.
#' @param pos A vector of methylated site locations (in base pairs). Elements
#' should be of the same order as the columns of M.
#' @param nbasis The number of B-spline basis functions.
#' @param binary Indicator of whether Y is binary.
#' @param cov Covariates. A matrix with dimensions n (number of subjects)
#' by J (number of covariates).
#' @param adapt_perm Whether "step-up" algorithm is used for P-value
#' calculation. If FALSE, function permutes nperm times and stops. If TRUE,
#' nperm will be increased 10 times each round if P-value <= 5/nperm.
#' Algorithm stops if P-value > 5/nperm or <= cutoff.
#' @param cutoff Cutoff for "step-up" algorithm.
#' @param nperm Number of permutations. Also the starting number of permutations
#' for "step-up" algorithm. Default is 1,000.
#' @param seed Specify random seed for permutations.
#' @param n0 Tuning parameter. Discard the first n0-1 P-values of each column.
#' @param ... Optional arguments \code{create.bspline.basis}.
#'
#' @return An object of "AFb" class.
#' \describe{
#'  \item{pv}{P-value of AFb test.}
#'  \item{stat}{Test statistic of AFb test.}
#'  \item{loci_combined}{Sites which are combined into the test
#'   statistic. The index of included sites are returned in the ascending
#'   order of their P-values.}
#'   \item{method}{Method used.}
#' }
#'
#' @export
#'
#' @seealso \code{\link[splines]{bs}}, \code{\link[fda]{fourier}},
#' \code{\link{set.seed}}
#'
#' @import wAF
#'
#' @examples
#' Y <- bs_dense$trait
#' methyl <- bs_dense$methyl
#' pos <- bs_dense$pos
#' test <- AFb(Y, methyl, pos, nbasis = 50, binary = TRUE, adapt_perm = TRUE)
#' summary(test)
#'
AFb <- function(Y, M, pos, nbasis,
                binary = FALSE, cov = NULL,
                adapt_perm = FALSE, cutoff = 2.5e-6,
                nperm = 1000, seed = NULL,
                n0 = 1, ...) {

  X <- basisMat(M, pos, nbasis, ...)

  test <- wAF(Y, X, binary = binary, cov = cov, w = "flat",
              adapt_perm = adapt_perm, cutoff = cutoff, nperm = nperm,
              n0 = n0, seed = seed)

  result <- list(pv = test$pv, stat = test$stat,
                 loci_combined = test$loci_combined,
                 stat_all = test$stat_all, pv_all = test$pv_all,
                 method = "AFb")
  class(result)<-"AFb"

  return(result)
}
