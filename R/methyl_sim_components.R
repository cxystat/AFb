# CG <- function(tot, il, p.CG = c(0.09443, 0.019, 0.012462),
#                seed = NULL) {
#   set.seed(seed)
#   cg <- rep(0, tot)
#   sh <- 2000
#   j <- 1
#   p <- p.CG[1]
#   while (j <= tot) {
#     cg[j] <- rbinom(1, 1, p)
#     j <- j + cg[j] + 1
#     if (j > il & j <= il + sh) {
#       p <- p.CG[2]
#     }
#     if (j > il + sh) {
#       p <- p.CG[3]
#     }
#   }
#   pos <- which(cg == 1)
#   result <- list(cg = cg, pos = pos)
#   return(result)
# }
#
# ts_to_rate <- function(x) {
#   t <- atan(x)
#   m <- apply(t, 1, min)
#   s <- apply(apply(t, 1, range), 2, diff)
#   result <- (t - m)/s
#   return(result)
# }
#
# methylprop <- function(n, tot, pos, weight = 0.9) {
#
#   n.cg <- length(pos) # number of CpG sites
#
#   # overall trend (time series)
#   ar.all <- 0.99^(1/10)
#   burnin <- 1000
#   ts0.all <- arima.sim(list(ar = ar.all), n = tot + burnin, sd = 0.15)- 2
#   ts.all <- ts0.all[pos + burnin]
#
#   # individual trend (AR covariance)
#   rho <- 0.9
#   corr.mat <- rho ^ abs(outer(pos, pos, "-"))
#   each <- mvrnorm(n = n, mu = rep(0, n.cg), Sigma = corr.mat)
#
#   # design matrix
#   ts <- t(weight * ts.all + (1-weight) * t(each))
#   ts.noise <- ts + matrix(rnorm(n * n.cg, sd = 0.05), n, n.cg)
#
#   methyl <- ts_to_rate(ts)
#   methylobs <- ts_to_rate(ts.noise)
#
#   result <- list(methyl = methyl, methylobs = methylobs)
#
#   return(result)
# }
#
# beta_bs <- function(tot, pos, prop, strg, knots = NULL, ...){
#   nk <- ceiling(length(pos) * 0.2) - 1
#   if (is.null(knots)) {
#     knots <- seq(1, tot, length.out = nk)[-c(1, nk)]
#   }
#   bs.beta <- bs(pos, knots = knots)
#   n.basis <- ncol(bs.beta)
#   nsig <- ceiling(n.basis * prop)
#   nz <- sample(1:n.basis, nsig)
#   gamma <- rep(0, n.basis)
#   gamma[nz] <- runif(nsig, -strg, strg)
#   beta <- bs.beta %*% gamma
#
#   return(beta)
# }
#
# beta_fourier <- function(pos, prop, strg,
#                          n.basis = ceiling(length(pos) * 0.2), ...) {
#   bs.beta <- fourier(pos, nbasis = n.basis, ...)
#   n.basis <- ncol(bs.beta)
#   nsig <- ceiling(n.basis * prop)
#   nz <- sample(1:n.basis, nsig)
#   gamma <- rep(0, n.basis)
#   gamma[nz] <- runif(nsig, -strg, strg)
#   beta <- bs.beta %*% gamma
#   return(beta)
# }
#
# beta_ar <- function(tot, pos, trans.prob = c(0.8, 0.95)) {
#   ar.beta <- 0.98^(1/5)
#   sd.beta <- 0.005
#   burnin <- 1000
#   ts.beta <- arima.sim(list(ar=ar.beta), n = tot + burnin, sd = sd.beta)
#   beta.all <- ts.beta[pos + burnin]
#   state <- c("y", "n")
#   trans.mat <- matrix(c(trans.prob[1], 1-trans.prob[1],
#                         1-trans.prob[2], trans.prob[2]),
#                       byrow = TRUE, nrow = 2)
#   mc <- new("markovchain", states = state,
#             transitionMatrix = trans.mat,
#             name = "beta")
#   n.cg <- length(pos)
#   beta.mc <- rmarkovchain(n.cg, mc, t0 = sample(state, 1))
#   nz.ind <- which(beta.mc == "y") # signal locations
#   beta <- rep(0, n.cg)
#   beta[nz.ind] <- beta.all[nz.ind]
#
#   return(beta)
# }
