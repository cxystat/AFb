#' Permutation of GLM-based Score Statistics
#'
#' @description This function calculates and permutes GLM-based
#' score statistics.
#'
#' @param Y Response Variable. A vector or length \code{n} (the number of samples).
#' @param X Explanatory Variables to be tested. A matrix with dimensions \code{n} by \code{K}
#' (the number of explanatory variables to be tested).
#' @param binary Indicator of whether \code{Y} is binary.
#' @param cov Covariates. A matrix with dimensions \code{n} (number of subjects)
#' by \code{J} (the number of explanatory variables that are not to be tested).
#' @param nperm Number of permutations. Default is \code{1,000}.
#' @param seed Specify the seed for permutations.
#'
#' @return A list object.
#' \describe{
#'   \item{Us}{Score statistics. The first column contains the original
#'   data.}
#'   \item{pvs}{P-values of two-tailed score tests;
#'   can be further used by AF_combine function.}
#'   \item{pvs_left}{Left-side P-values of one-tailed score tests;
#'   can be further used by AF_combine function.}
#'}
#'
#' @export
#'
#' @seealso {\code{\link{set.seed}}}
#'
#' @import stats statmod
#'
#' @examples
#' Y <- ar_dense$trait
#' X <- ar_dense$methyl
#' result <- perm_score(Y, X, binary = TRUE, nperm = 100)
#' names(result)
#'
perm_score <- function(Y, X, binary = FALSE, cov = NULL,
                       nperm = 1e3, seed = NULL) {
  K <- ncol(X)
  n <- nrow(X)

  if (length(Y) != n)
    stop("length of Y should be the same as row numbers of X")

  s <- apply(X, 2, sd)
  if (any (s == 0)) {
    ind <- which(s == 0)
    print("column indexes of X with no variation")
    print(ind)
    X <- X[, -ind]
  }

  if(is.null(cov)) cov <- rep(1, length(Y))

  if (all(s == 0)) {
    print("all columns of X have no variation; set all P-values as 1.")
    Up <- matrix(Inf, K, nperm + 1)
    p <- p1 <- matrix(1, K, nperm + 1)
  } else {
    ## fit the null model
    model <- ifelse(binary, "binomial", "gaussian")
    data1 <- data.frame(trait = Y, cov)
    null.model <- glm(trait ~ ., family = model, data = data1)

    cov <- data.frame(cov)
    X.null <- lm(X ~ ., data = cov)
    Xres <- resid(X.null)
    if(K == 1) {Xres <- matrix(Xres, ncol = K)}

    ## calculate score test statistics U
    U <- glm.scoretest(null.model, Xres)

    ## permute X residuals and calculate permutation Us
    if (!is.null(seed)) {
      set.seed(seed)
    }
    order.perm <- replicate(nperm, sample(1:n))
    U.perm <- apply(order.perm, 2,
                    function(x) glm.scoretest(null.model, Xres[x, ]))

    ## Calculate p-values
    if (K == 1) {
      Up <- matrix(c(U, U.perm), nrow = K)
    } else {
      Up <- cbind(U, U.perm)
    }
    p <- 2 * (1 - pnorm(abs(Up)))
    p1 <- pnorm(Up)
  }

  result <- list(Us = Up, pvs = p, pvs_left = p1)
  return(result)
}
