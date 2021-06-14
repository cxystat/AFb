#' Methylation Levels to Principal Components
#'
#' @param M A matrix of methylation levels with dimensions \code{n} by \code{K}
#' (\code{n} subjects, \code{K} CpG sites).
#' @param varprop Cutoff for proportion of variance to be expalined.
#'  The first \code{Kp} principal components (PCs) are chosen such that \code{varprop}
#'  of the total variance is explained.
#'
#' @return A matrix of PCs with dimension \code{n} by \code{Kp}.
#' @export
#'
pcaMat <- function(M, varprop = 0.95) {

  K <- ncol(M)

  decomp <- eigen(t(M) %*% M)
  PC <- M %*% decomp$vectors
  varcs <- cumsum(decomp$values)
  pve <- varcs/varcs[K]
  npc <- sum(pve < varprop) + 1
  pcmat <- PC[, 1:npc, drop = FALSE]

  colnames(pcmat) <- paste0("PC", 1:npc)

  return(pcmat)
}
