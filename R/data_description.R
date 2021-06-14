#' Simulated Methylation Data with B-spline Effects in Dense Scenario
#'
#' A simulated dataset containing 500 subjects and 100 methylation sites.
#' Methylation effects are smoothed by B-spline basis functions.
#' 20\% of the basis functions are simulated to be causal/associated. Effect
#' strengths are randomly sampled from U(-0.05, 0.05) distribution.
#'
#' @docType data
#'
#' @usage data(bs_dense)
#'
#' @format A list object.
#' \describe{
#'   \item{trait}{A vector of length 500 containing disease labels
#'   (250 cases and 250 controls).}
#'   \item{methyl}{A 500 by 100 matrix containing methylation levels. Each
#'   row represents the methylation profile for one subject. Each
#'   element of the matrix denotes BS-seq methylation level at a CpG site.}
#'   \item{pos}{Locations for CpG sites (in base pairs).}
#' }
"bs_dense"

#' Simulated Methylation Data with Fourier Effects in Sparse Scenario
#'
#' A simulated dataset containing 500 subjects and 93 methylation sites.
#' Methylation effects are smoothed by Fourier basis functions.
#' 1\% of the basis functions are simulated to be causal/associated. Effect
#' strengths are randomly sampled from U(-5, 5) distribution.
#'
#' @docType data
#'
#' @usage data(fourier_sparse)
#'
#' @format A list object.
#' \describe{
#'   \item{trait}{A vector of length 500 containing disease labels
#'   (254 cases and 246 controls).}
#'   \item{methyl}{A 500 by 93 matrix containing methylation levels. Each
#'   row represents the methylation profile for one subject. Each
#'   element of the matrix denotes BS-seq methylation level at a CpG site.}
#'   \item{pos}{Locations for CpG sites (in base pairs).}
#' }
"fourier_sparse"

#' Simulated Methylation Data with Autoregressive Effects in Dense Scenario
#'
#' A simulated dataset containing 500 subjects and 81 methylation sites.
#' Methylation effects are constructed by first-order autoregressive model.
#'
#' @docType data
#'
#' @usage data(ar_dense)
#'
#' @format A list object.
#' \describe{
#'   \item{trait}{A vector of length 500 containing disease labels
#'   (251 cases and 249 controls).}
#'   \item{methyl}{A 500 by 81 matrix containing methylation levels. Each
#'   row represents the methylation profile for one subject. Each
#'   element of the matrix denotes BS-seq methylation level at a CpG site.}
#'   \item{pos}{Locations for CpG sites (in base pairs).}
#' }
"ar_dense"
