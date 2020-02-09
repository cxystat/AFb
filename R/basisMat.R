#' Methylation Levels to Basis Function Values
#'
#' @param M A matrix of methylation levels. A matrix with dimensions n by K
#' (n subjects, K methylation sites). Each row for a subject, and each column
#' for a methylated CpG site.
#' @param pos A vector of methylated site locations (in base pairs). Elements
#' should be of the same order as the columns of M.
#' @param basis Basis functions used, "bs" for B-spline, "fourier" for Fourier.
#' @param ... Optional arguments for constructing basis functions.
#'
#' @return A matrix of basis function values. A matrix with dimensions n by Kb,
#' where Kb is the number of basis functions.
#' @export
#'
#' @seealso \code{\link[splines]{bs}}, \code{\link[fda]{fourier}}
#'
#' @examples
#' M <- ar_dense$methyl
#' pos <- ar_dense$pos
#' l <- length(pos)
#'
#' ## B-spline basis functions
#' knots <-seq(1, pos[l], length.out = 50)[-c(1, 50)]
#' X1 <- basisMat(M, pos, basis = "bs", degree = 2, knots = knots)
#'
#' ## Fourier basis functions
#' # if nbasis is an even number, function fourier takes nbasis = nbasis + 1
#' X2 <- basisMat(M, pos, basis = "fourier", nbasis = 20)
#'
basisMat <- function(M, pos, basis = c("bs", "fourier"), ...) {

  dl <- length(dim(M))
  if (dl != 2)
    stop("M must be a matrix")
  if (is.null(pos))
    stop("no positions given")

  if(is.unsorted(pos)) {
    cat("pos is unsorted, sort into ascending order. \n")
    o <- order(pos)
    M <- M[,o]
    pos <- sort(pos)
  }
  K <- length(pos)
  dpos <- c(pos[2]-pos[1], diff(pos, lag = 2)/2, pos[K]-pos[K-1])

  basis <- match.arg(basis)
  FUN <- match.fun(basis)
  bfMat <- FUN(pos, ...)

  X <- M %*% (bfMat * dpos)

  return(X)
}

