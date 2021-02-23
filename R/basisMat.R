#' Methylation Levels to Basis Function Values
#'
#' @param M A matrix of methylation levels. A matrix with dimensions n by K
#' (n subjects, K methylation sites). Each row for a subject, and each column
#' for a methylated CpG site.
#' @param pos A vector of methylated site locations (in base pairs). Elements
#' should be of the same order as the columns of M.
#' @param nbasis The number of B-spline basis functions (default is cubic
#' splines). Note that nbasis should not be smaller than norder (default is 4).
#' @param start Start location of the region.
#' @param end End location of the region.
#' @param ... Optional arguments for \code{create.bspline.basis}.
#'
#' @return A matrix of basis function values. A matrix with dimensions n by Kb,
#' where Kb is the number of basis functions.
#' @export
#'
#' @seealso \code{\link[fda]{create.bspline.basis}}
#'
#' @import fda
#'
#' @examples
#' M <- ar_dense$methyl
#' pos <- ar_dense$pos
#' l <- length(pos)
#'
#' ## B-spline basis functions
#' X <- basisMat(M, pos, nbasis = 50, start = min(pos),
#'               end = max(pos), norder = 3)
basisMat <- function (M, pos, start, end, nbasis, ...)
{
  dl <- length(dim(M))
  if (dl != 2)
    stop("M must be a matrix")
  if (is.null(pos))
    stop("no positions given")
  if (is.unsorted(pos)) {
    cat("pos is unsorted, sort into ascending order. \n")
    o <- order(pos)
    M <- M[, o]
    pos <- sort(pos)
  }
  K <- length(pos)
  dpos <- c(pos[2] - pos[1], diff(pos, lag = 2)/2,
            pos[K] - pos[K - 1])
  basis <- create.bspline.basis(rangeval = c(start, end),
                                nbasis = nbasis, ...)
  bfMat <- eval.basis(pos, basis)
  X <- M %*% (bfMat * dpos)
  return(X)
}
