#' Methylation Levels to Basis Function Values
#'
#' @param M A matrix of methylation levels. A matrix with dimensions n by K
#' (n subjects, K methylation sites). Each row for a subject, and each column
#' for a methylated CpG site.
#' @param pos A vector of methylated site locations (in base pairs). Elements
#' should be of the same order as the columns of M.
#' @param nbasis The number of B-spline basis functions (defult is cubic
#' splines).
#' @param ... Optional arguments for \code{create.bspline.basis} to
#' construct basis functions.
#'
#' @return A matrix of basis function values. A matrix with dimensions n by Kb,
#' where Kb is the number of basis functions.
#' @export
#'
#' @seealso \code{\link[fda]{create.bspline.basis}},
#'
#' @examples
#' M <- ar_dense$methyl
#' pos <- ar_dense$pos
#' l <- length(pos)
#'
#' ## B-spline basis functions
#' X1 <- basisMat(M, pos, nbasis = 50, norder = 3)
basisMat <- function(M, pos, nbasis, ...) {

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

  c <- min(pos)
  s <- diff(range(pos))
  pos.std <- as.vector(scale(pos, center = c, scale = s))

  K <- length(pos)
  dpos <- c(pos.std[2]-pos.std[1], diff(pos.std, lag = 2)/2,
            pos.std[K]-pos.std[K-1])

  basis <- create.bspline.basis(nbasis = nbasis, ...)
  bfMat <- eval.basis(pos.std, basis)

  X <- M %*% (bfMat * dpos)

  return(X)
}


