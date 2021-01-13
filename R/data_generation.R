# library(AFb)
# library(devtools)
# source("R/methyl_sim.R")
#
# ## B-spline dense
# # K <- 329 # Q1 of CpG island length
# # t <- 39139 # Q1 of CpG + gap + gene
# # n <- 1000
#
# K <- 300
# t <- 5000
# n <- 500
# p.sig <- 0.2
# delta <- 0.05
#
# bs_dense <- methyl_sim(n, t, K, p.sig, delta,
#                        beta = "bs", seed = 7)
# use_data(bs_dense)
#
#
# ## Fourier sparse
# K <- 300 # Q1 of CpG island length
# t <- 5000 # Q1 of CpG + gap + gene
# n <- 500
# p.sig <- 0.01
# delta <- 5
#
# fourier_sparse <- methyl_sim(n, t, K, p.sig, delta,
#                             beta = "fourier", seed = 28)
#
# use_data(fourier_sparse)
#
# ## AR(1)
# K <- 300 # Q1 of CpG island length
# t <- 5000 # Q1 of CpG + gap + gene
# n <- 500
# prob <- c(0.8, 0.98)
#
# ar_dense <- methyl_sim(n, t, K, trans.prob = prob,
#                        beta = "ar", seed = 56)
#
# use_data(ar_dense)
