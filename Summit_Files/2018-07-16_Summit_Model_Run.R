
load("2018-07-16_occConnC_Summit_Data.RData")

library(raster)
library(doParallel)
library(gdistance)
library(dplyr)

# run sampler
occConn_mout <- occuConnC(
  x = data$coords,
  y = data$y_mat,
  j = data$j_mat,
  site_covs = data$sitecovs,
  obs_covs = data$season_vec,
  r_covs = data$res_covs,
  disp_dist = 50000,
  iters = 50000,
  report = 100,
  monitor.z = TRUE,
  plot.z = TRUE,
  # tune order (alpha[1], alpha[2], alpha[3], b0.gam, b.gam[1], b.gam[2], b.gam[3], b.gam[4], sigma, b0.psi1, b.psi1[1], b.psi1[2], b.psi1[3],
  # b0.eps, b.eps[1], b.eps[2], b.eps[3], b.eps[4], b.eps[5], b.eps[6], b.eps[7], b.eps[8], a0, season[2], season[3], season[4]
  tune = c(0.3, 0.3, 0.3, # alpha coefficients
           0.3, 0.3, 0.3, 0.3, 0.3, # gamma coefficients
           0.3, # sigma
           0.3, 0.3, 0.3, 0.3, # psi1 coefficients
           0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, # epsilon coefficients
           0.2, 0.2, 0.2, 0.2 #a0 and season parameters (detection)
  ),
  param_mon = c("alpha[1]","alpha[2]", "alpha[3]", "sigma", "b0.gam", 
                "b.gam[1]", "b.gam[2]", "b.gam[3]", "b.gam[4]", "b0.psi1", 
                "b.psi1[1]", "b.psi1[3]", "b.psi1[3]", "b0.eps", "b.eps[1]", 
                "b.eps[2]", "b.eps[3]", "b.eps[4]", "b.eps[5]", "b.eps[6]", 
                "b.eps[7]", "b.eps[8]", "a0", "season[2]","season[3]", 
                "season[4]", "zk", "deviance")
)
