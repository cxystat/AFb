myks <- function(y, x, kernel, bandwidth) {
  if(is.unsorted(x)) {
    o <- order(x)
    x <- x[o]
    y <- y[o]
  }
  ks <- ksmooth(x, y, kernel = kernel,
                bandwidth = bandwidth, x.points = x)
  disp <- abs(y - ks$y)

  return(disp)

}



#' Methylation Levels to Dispersion
#'
#' @param M A matrix of methylation levels with dimensions \code{n} by \code{K} (\code{n} subjects, \code{K} CpG sites).
#' @param pos A vector of CpG locations (in base pairs). Elements should be of the same
#' order as the columns of \code{M}.
#' @param bandwidth The bandwidth for \code{ksmooth}. Default is \code{5509},
#' 20 \% of the median gene length according to UCSC Genome build hg38.
#' @param kernel The kernel to be used for \code{ksmooth}.
#'
#' @return A matrix of methylation dispersion with dimensions \code{n} by \code{K}.
#'
#' @export
#'
#' @import stats
#'
#' @seealso \code{\link[stats]{ksmooth}}
#'
dispMat <- function(M, pos, bandwidth = 5509,
                    kernel = c("normal", "box")) {
  ksres <- apply(M, 1, FUN = myks, x = pos,
                 bandwidth = bandwidth, kernel = kernel)
  Mdisp <- t(ksres)
  return(Mdisp)
}
