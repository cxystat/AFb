#' Summary Function for Objects of "AF" Class
#'
#' @param object An object of "AF" class.
#' @param ... Optional arguments for \code{summary}.
#'
#' @return Method used; P-value; basis functions used;
#' CpG sites combined into the test statistic.
#'
#' @export
#'
#' @usage \method{summary}{AF}(object, ...)
#'
summary.AF <- function(object, ...){
  cat("Method:\n")
  cat(paste(object$method, "\n"))
  cat("\n")
  cat("P-value:\n")
  print(object$pv)
  cat("\n")
  cat("CpG sites combined into test statistic:\n")
  print(object$loci_combined)
}


#' Print Function for Objects of "AF" Class
#'
#' @param x An object of "AF" class.
#' @param ... Optional arguments for \code{print}.
#'
#' @return Method used; P-value; basis functions used;
#' CpG sites combined into the test statistic.
#'
#' @export
#'
print.AF <- function(x, ...){
  cat("Method:\n")
  cat(paste(x$method, "\n"))
  cat("\n")
  cat("P-value:\n")
  print(x$pv)
  cat("\n")
  cat("CpG sites combined into test statistic:\n")
  print(x$loci_combined)
}
