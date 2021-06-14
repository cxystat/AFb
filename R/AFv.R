#' Adaptive Fisher Method for Testing Methylation Variation
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
#' @param bandwidth bandwidth The bandwidth for \code{ksmooth}. Default is \code{5509},
#' 20 \% of the median gene length according to UCSC Genome build hg38.
#' @param kernel kernel The kernel to be used for \code{ksmooth}
#' @param adapt_perm Whether "step-up" algorithm is used for P-value
#' calculation. If FALSE, function permutes \code{nperm} times and stops. If TRUE,
#' \code{nperm} will be increased 10 times each round if P-value <= \code{5/nperm}.
#' Algorithm stops if P-value > \code{5/nperm} or <= \code{cutoff}.
#' @param cutoff Cutoff for "step-up" algorithm.
#' @param seed Specify the seed for permutations.
#' @param n0 Tuning parameter. Discard the first \code{n0-1} P-values of each sample.
#' @param ... Optional arguments \code{create.bspline.basis}.
#'
#' @return An object of "AFv" class.
#' \describe{
#'  \item{pv}{P-value of AFv test.}
#'  \item{stat}{Test statistic of AFv test.}
#'  \item{indexes}{Indexes of basis functions combined into the test
#'   statistic. Indexes are sorted so that P-values are
#'   in ascending order.}
#'   \item{stat_all}{AFv statistics for all permuted samples.}
#'   \item{pv_all}{P-values of AFv statistics for all permuted samples.}
#'   \item{method}{Method used.}
#' }
#'
#'
#' @export
#'
#' @seealso \code{\link[fda]{create.bspline.basis}},
#' \code{\link{set.seed}}, \code{\link[stats]{ksmooth}}
#'
#' @examples
#' Y <- bs_dense$trait
#' methyl <- bs_dense$methyl
#' pos <- bs_dense$pos
#' test <- AFv(Y, methyl, pos, nbasis =10, bandwidth = 500,
#'             kernel = "box", binary = TRUE)
#' summary(test)
#'
AFv <- function(Y, M, pos, nbasis, start = NULL, end = NULL,
                binary = FALSE, cov = NULL,
                nperm = 1000, bandwidth = 5509,
                kernel = c("normal", "box"),
                adapt_perm = FALSE, cutoff = 2.5e-6,
                seed = NULL,
                n0 = 1, ...) {

  V <- dispMat(M, pos = pos, bandwidth = bandwidth,
               kernel = kernel)

  test <- AFb(Y, V, pos = pos, nbasis = nbasis,
              start = start, end = end,
              binary = FALSE, cov = NULL, nperm = 1000,
              adapt_perm = FALSE, cutoff = 2.5e-6,
              seed = NULL, n0 = 1, ...)

  result <- list(pv = test$pv, stat = test$stat,
                 indexes = test$indexes,
                 stat_all = test$stat_all,
                 pv_all = test$pv_all,
                 method = "AFv")
  class(result) <- "AF"

  return(result)

}
