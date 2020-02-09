#' Simulated Methylation Data with B-spline Effects in Dense Scenario
#'
#' A simulated dataset containing 1,000 subjects and 527 methylation sites.
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
#'   \item{trait}{A vector of length 1,000 containing disease labels
#'   (492 cases and 508 controls).}
#'   \item{methyl}{A 1,000 by 527 matrix containing methylation levels. Each
#'   row represents the methylaiton profile for one subject. Each
#'   element of the matrix denotes BS-seq methylation level at a CpG site.}
#'   \item{pos}{Locations for CpG sites (in base pairs).}
#' }
"bs_dense"

#' Simulated Methylation Data with Fourier Effects in Sparse Scenario
#'
#' A simulated dataset containing 1,000 subjects and 490 methylation sites.
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
#'   \item{trait}{A vector of length 1,000 containing disease labels
#'   (523 cases and 477 controls).}
#'   \item{methyl}{A 1,000 by 490 matrix containing methylation levels. Each
#'   row represents the methylaiton profile for one subject. Each
#'   element of the matrix denotes BS-seq methylation level at a CpG site.}
#'   \item{pos}{Locations for CpG sites (in base pairs).}
#' }
"fourier_sparse"

#' Simulated Methylation Data with Autoregressive Effects in Dense Scenario
#'
#' A simulated dataset containing 1,000 subjects and 477 methylation sites.
#' Methylation effects are constructed by first-order autoregressive model.
#'
#' @docType data
#'
#' @usage data(ar_dense)
#'
#' @format A list object.
#' \describe{
#'   \item{trait}{A vector of length 1,000 containing disease labels
#'   (510 cases and 490 controls).}
#'   \item{methyl}{A 1,000 by 477 matrix containing methylation levels. Each
#'   row represents the methylaiton profile for one subject. Each
#'   element of the matrix denotes BS-seq methylation level at a CpG site.}
#'   \item{pos}{Locations for CpG sites (in base pairs).}
#' }
"ar_dense"
