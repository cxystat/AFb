#' Methylation Levels to Basis Function Values
#'
#' @param M A matrix of methylation levels with dimensions \code{n} by \code{K}
#' (\code{n} subjects, \code{K} CpG sites).
#' @param pos A vector of CpG locations (in base pairs). Elements
#' should be of the same order as the columns of \code{M}.
#' @param nbasis The number of B-spline basis functions (default is cubic
#' splines). Note that \code{nbasis} should not be smaller than norder (default is 4).
#' @param start Start location of the region.
#' @param end End location of the region.
#' @param ... Optional arguments for \code{create.bspline.basis}.
#'
#' @return A matrix of basis function values with dimensions \code{n} by \code{Kb},
#' where \code{Kb} is the number of basis functions.
#' @export
#'
#' @seealso \code{\link[fda]{create.bspline.basis}}
#'
#' @import fda
#'
#' @examples
#' Y <- bs_dense$trait
#' methyl <- bs_dense$methyl
#' pos <- bs_dense$pos
#' X <- basisMat(methyl, pos, nbasis = 10)
#'
basisMat <- function (M, pos, nbasis, start = NULL, end = NULL, ...) {
  dl <- length(dim(M))
  K <- length(pos)
  if (dl != 2)
    stop("M must be a matrix.\n")
  if (K == 1 & is.null(start) & is.null(end))
    stop("single site, cannot calculate basis functions.\n")
  if (is.null(pos))
    stop("no positions given")
  if (is.unsorted(pos)) {
    cat("pos is unsorted, sort into ascending order. \n")
    o <- order(pos)
    M <- M[, o]
    pos <- pos[o]
  }
  if(is.null(start)) start <- 2 * pos[1] - pos[2]
  if(is.null(end)) end <- 2 * pos[K] - pos[K-1]
  pos.ext <- c(start, pos, end)
  dpos <- diff(pos.ext, lag = 2)/2
  basis <- create.bspline.basis(rangeval = c(start, end),
                                nbasis = nbasis, ...)
  bfMat <- eval.basis(pos, basis)
  X <- M %*% (bfMat * dpos)
  return(X)
}
