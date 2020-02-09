#' Summary Function for Objects of "AFb" Class
#'
#' @param x An object of "AFb" class.
#'
#' @return Method used; P-value; basis functions used;
#' CpG sites combined into the test statistic.
#'
#' @export
#'
#' @examples
#'
summary.AFb <- function(x, ...){
  cat("Method:\n")
  cat(paste(x$method, "\n"))
  cat("\n")
  cat("P-value:\n")
  print(x$pv)
  cat("\n")
  cat("Basis: \n")
  cat(paste(x$basis, "\n"))
  cat("\n")
  cat("CpG sites combined into test statistic:\n")
  print(x$loci_combined)
}


#' Print Function for Objects of "AFb" Class
#'
#' @param x An object of "AFb" class.
#'
#' @return Method used; P-value; basis functions used;
#' CpG sites combined into the test statistic.
#'
#' @export
#'
#' @examples
#'
print.AFb <- function(x, ...){
  cat("Method:\n")
  cat(paste(x$method, "\n"))
  cat("\n")
  cat("P-value:\n")
  print(x$pv)
  cat("\n")
  cat("Basis: \n")
  cat(paste(x$basis, "\n"))
  cat("\n")
  cat("CpG sites combined into test statistic:\n")
  print(x$loci_combined)
}
