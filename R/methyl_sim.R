## Data Simulation
##
## @param n Number of subjects.
## @param tot Total length (in base pairs) of the region being tested.
## @param il Length (in base pairs) of CpG island.
## @param prop Proportion of true causal/associated basis.
## @param strg Effect size. Effects are selected from U(-strg, strg)
## distribution.
## @param beta Model used to construct methylation effects, "bs" for B-spline
## basis functions, "fourier" for Fourier basis functions,
## "ar" for AR(1) model. Default is "bs".
## @param weight Weight of baseline methylation levels. Default is 0.9.
## @param p.CG A vector of length three containing CpG site proportions on CpG
## island, shore and desert.
## @param knots A vector of knots for B-spline basis funcitons. Used only when
## argument beta = "bs".
## @param n.basis Number of Fourier basis functions. Used only when argument
## beta = "fourier".
## @param trans.prob Transition matrix for Markov chain. Used only when
## argument beta = "ar".
## @param seed Random seed.
## @param ... Optional argument for "bs" and "fourier".
##
## @return A list object.
## \describe{
##  \item{trait}{A vector of values for a binary trait for n subjects.}
##  \item{methyl}{A matrix with methylation levels for n subjects.}
##  \item{pos}{Locations of methylation sites.}
## }
## @export
##
## @examples
## K <- 329 # Q1 of CpG island length
## t <- 39139 # Q1 of CpG + gap + gene
## n <- 1000
## p.sig <- 0.2
## delta <- 0.05
##
## bs_dense <- methyl_sim(n, t, K, p.sig, delta,
##                        beta = "bs", seed = 7)
methyl_sim <- function(n, tot, il, prop, strg,
                       beta = c("bs", "fourier", "ar"),
                       weight = 0.9,
                       p.CG = c(0.09443, 0.019, 0.012462),
                       knots = NULL,
                       n.basis = NULL,
                       trans.prob = c(0.8, 0.95),
                       seed = NULL, ...) {
  source("R/methyl_sim_components.R")

  temp <- CG(tot, il, p.CG, seed = seed)
  pos <- temp$pos

  designMat <- methylprop(n, tot, pos, weight = weight)

  beta <- match.arg(beta)
  FUN <- match.fun(paste0("beta_", beta))
  if (beta == "ar") {
    coefs <- FUN(tot, pos, trans.prob)
  } else {
    if (beta == "fourier") {
      if (is.null(n.basis)) {
        coefs <- FUN(pos, prop, strg, ...)
      } else {
        coefs <- FUN(pos, prop, strg, n.basis, ...)
      }
    } else {
      coefs <- FUN(tot, pos, prop, strg, knots = NULL, ...)
    }
  }

  methyl.std <- scale(designMat$methyl, scale = FALSE)
  n.cg <- length(pos)
  hd <- c(pos[2] - pos[1] ,diff(pos, lag = 2)/2,
          pos[n.cg]-pos[n.cg-1])

  w.coef <- coefs * hd
  d.prob <- inv.logit(methyl.std %*% w.coef)
  Y <- rbinom(n, 1, d.prob)

  result <- list(trait = Y, methyl = designMat$methylobs, pos = pos)

  return(result)

}
