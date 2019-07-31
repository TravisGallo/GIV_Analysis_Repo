
files <- list.files("./ARGO_Return/run1_output/")

samples_list <- vector("list", length(files))
for(i in 1:length(files)){
  samples_list[[i]] <- readRDS(paste0("./ARGO_Return/run1_output/",files[i]))$samples
}

library(coda)
library(mcmcplots)
mcmc_list <- lapply(samples_list, function(x){ as.mcmc(x) })

mcmc_chains <- mcmc.list(mcmc_list)


plot(mcmc_chains)
