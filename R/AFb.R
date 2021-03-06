#' Adaptive Fisher Method with B-spline basis funcitons (AFb)
#' for Trait-Methylation Set Association
#'
#' @param Y Phenotype data. It can be a continuous trait or a binary trait.
#' A vector of length \code{n} (number of subjects).
#' @param M A matrix of methylation levels with dimensions \code{n} by \code{K}
#' (\code{n} subjects, \code{K} CpG sites).
#' @param pos A vector of CpG locations (in base pairs). Elements
#' should be of the same order as the columns of \code{M}.
#' @param nbasis The number of B-spline basis functions.
#' @param start Start location of the region.
#' @param end End location of the region.
#' @param binary Indicator of whether \code{Y} is binary.
#' @param cov Covariates. A matrix with dimensions \code{n}
#' by \code{J} (number of covariates).
#' @param nperm Number of permutations. Also the starting number of permutations
#' for "step-up" algorithm. Default is \code{1,000}.
#' @param adapt_perm Whether "step-up" algorithm is used for P-value
#' calculation. If FALSE, function permutes \code{nperm} times and stops. If TRUE,
#' \code{nperm} will be increased 10 times each round if P-value <= \code{5/nperm}.
#' Algorithm stops if P-value > \code{5/nperm} or <= \code{cutoff}.
#' @param cutoff Cutoff for "step-up" algorithm.
#' @param seed Specify the seed for permutations.
#' @param n0 Tuning parameter. Discard the first \code{n0-1} P-values of each sample.
#' @param ... Optional arguments for \code{create.bspline.basis}.
#'
#' @return An object of "AFb" class.
#' \describe{
#'  \item{pv}{P-value of AFb test.}
#'  \item{stat}{Test statistic of AFb test.}
#'  \item{indexes}{Indexes of basis functions combined into the test
#'   statistic. Indexes are sorted so that P-values are
#'   in ascending order.}
#'   \item{stat_all}{AFb statistics for all permuted samples.}
#'   \item{pv_all}{P-values of AFb statistics for all permuted samples.}
#'   \item{method}{Method used.}
#' }
#'
#' @export
#'
#' @seealso \code{\link[fda]{create.bspline.basis}},
#' \code{\link{set.seed}}
#'
#' @examples
#' Y <- bs_dense$trait
#' methyl <- bs_dense$methyl
#' pos <- bs_dense$pos
#' test <- AFb(Y, methyl, pos, nbasis = 10,
#'             binary = TRUE, adapt_perm = TRUE)
#' summary(test)
#'
AFb <- function(Y, M, pos, nbasis, start = NULL, end = NULL,
                binary = FALSE, cov = NULL, nperm = 1000,
                adapt_perm = FALSE, cutoff = 2.5e-6,
                seed = NULL,
                n0 = 1, ...) {

  X <- basisMat(M, pos, nbasis = nbasis, start = start,
                end = end, ...)

  test <- AF(Y, X, binary = binary, cov = cov, nperm = nperm,
             adapt_perm = adapt_perm, cutoff = cutoff,
             n0 = n0, seed = seed)

  result <- list(pv = test$pv, stat = test$stat,
                 indexes = test$indexes,
                 stat_all = test$stat_all,
                 pv_all = test$pv_all,
                 method = "AFb")
  class(result)<-"AF"

  return(result)
}
