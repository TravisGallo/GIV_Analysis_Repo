#setwd("~/Documents/GitHub Repos/GIV_Analysis_Repo")
## MCMC algorithm

# load libraries
library(sp, lib.loc = "/scratch/hgallo/Rlibs/")
library(raster, lib.loc = "/scratch/hgallo/Rlibs/")
library(foreach, lib.loc = "/scratch/hgallo/Rlibs/")
library(iterators, lib.loc = "/scratch/hgallo/Rlibs/")
library(doParallel, lib.loc = "/scratch/hgallo/Rlibs/")
library(igraph, lib.loc = "/scratch/hgallo/Rlibs/")
library(gdistance, lib.loc = "/scratch/hgallo/Rlibs/")
library(crayon, lib.loc = "/scratch/hgallo/Rlibs/")
library(dplyr, lib.loc = "/scratch/hgallo/Rlibs/")
library(data.table, lib.loc = "/scratch/hgallo/Rlibs/")

load("/home/hgallo/2018-08-27_occConnC_Data.RData")

## species specific settings
# dispersal distance of species
disp_dist <- 50000

## some early settings for inspection
report <- 100 # report every x amount of iterations
monitor.z <- FALSE # track z
plot.z <- FALSE # plot z with each report

## parameters for mcmc sample
# mcmc iterations
iters <- 3
# parameters to monitor
param_mon <- c("alpha[1]","alpha[2]", "alpha[3]", "sigma", "b0.gam", 
              "b.gam[1]", "b.gam[2]", "b.gam[3]", "b.gam[4]", "b0.psi1", 
              "b.psi1[1]", "b.psi1[3]", "b.psi1[3]", "b0.eps", "b.eps[1]", 
              "b.eps[2]", "b.eps[3]", "b.eps[4]", "b.eps[5]", "b.eps[6]", 
              "b.eps[7]", "b.eps[8]", "a0", "season[2]","season[3]", 
              "season[4]", "zk", "deviance")
# tuning parameters
tune <- c(0.3, 0.3, 0.3, # alpha coefficients
         0.3, 0.3, 0.3, 0.3, 0.3, # gamma coefficients
         0.3, # sigma
         0.3, 0.3, 0.3, 0.3, # psi1 coefficients
         0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, # epsilon coefficients
         0.2, 0.2, 0.2, 0.2 #a0 and season parameters (detection)
)

## Dimensions
nsite <- nrow(data$coords) # Number of possible sites instead of only the sites sampled
nseason <- ncol(data$y_mat)
nsampled <- nrow(data$y_mat)

## indicating real detections from data
anyDetections <- matrix(FALSE, nsite, nseason)
anyDetections[1:nsampled,] <- as.numeric(data$y_mat > 0)

## initial values for sampler
# resistence
alpha <- rgamma(3, 1, 4)
# initial occupancy
b0.psi1 <- rnorm(1, 0, 0.5)
b.psi1 <- rnorm(3, 0, 0.5)
# colonization
b0.gam <- rnorm(1,0,0.5) # intercept
b.gam <- rnorm(4,0,0.5) # 4 covariates
gamma0 <- plogis(b0.gam + b.gam[1]*data$sitecovs[,"size"] + 
                   b.gam[2]*data$sitecovs[,"park"] + 
                   b.gam[3]*data$sitecovs[,"cem"] + 
                   b.gam[4]*data$sitecovs[,"golf"])
sigma <- runif(1,1,10)
# extinction
b0.eps <- rnorm(1, 0, 0.5) # intercept
b.eps <- rnorm(8, 0, 0.5) # 8 covariates
epsilon <- plogis(b0.eps + b.eps[1]*data$sitecovs[ ,"tree"] + 
                           b.eps[2]*data$sitecovs[ ,"total_veg"] + 
                           b.eps[3]*data$sitecovs[ ,"size"] + 
                           b.eps[4]*data$sitecovs[ ,"pop10"] +
                           b.eps[5]*data$sitecovs[ ,"water"] +
                           b.eps[6]*data$sitecovs[ ,"park"] +
                           b.eps[7]*data$sitecovs[ ,"cem"] +
                           b.eps[8]*data$sitecovs[ ,"golf"] )
# detection
p <- rep(0, nseason)
a0 <- rnorm(1, 0, 0.5)
season <- rnorm(4,0 , 0.5)
season[1] <- 0
p <- plogis(a0 + season[data$season_vec])
# need to run through the model to generate starting values for z[,k-1] based off z[,1] starting values,
# psi, gamma, and likelihoods for z and y
gamma <- matrix(NA, nsite, nseason-1)
psi1 <- plogis(b0.psi1 + b.psi1[1]*data$sitecovs[,"park"] +
                         b.psi1[2]*data$sitecovs[,"cem"]  +
                         b.psi1[3]*data$sitecovs[,"golf"])

# psi matrix holds all occupancy probabilities in t = 1 and then the 
# associated colonization / extinction probabilities at t > 1.
psi <- matrix(NA, nsite, nseason)

# occupancy prob season 1 from initial values
psi[,1] <- psi1

# log likelihood matrices for latent and observed states
ll.z <- matrix(0, nsite, nseason)
ll.y <- matrix(0, nsampled, nseason)

# create resistance surface
cost <- exp(alpha[1]*data$res_covs[[1]] + 
            alpha[2]*data$res_covs[[2]] + 
            alpha[3]*data$res_covs[[3]])

# create the transition function and cell lists to be used in calculating the
# shortest paths
adj_graph <- getAllEdgeWeights(cost, directions=8, coords = data$coords)

# have to make variables their own stand alone variable to export to cluster
adjacencyGraph <- adj_graph$adjacencyGraph
cell_list <- adj_graph$cell_list
toCells <- adj_graph$toCells
fromCells <- adj_graph$fromCells
rm(adj_graph) # remove function object to save memory

# create wrapper for shortest.paths function to go into lApply
shortestPaths <- function(x){ shortest.paths(adjacencyGraph,
                                             v = x[[1]],
                                             to = x[[2]],
                                             mode = "out",
                                             algorithm = "dijkstra") }
cl <- makeCluster(detectCores())
# export function and variables to clusters
clusterExport(cl=cl, list("adjacencyGraph", "shortestPaths", 
                          "shortest.paths"))
# send libraries to cluster
clusterEvalQ(cl,library(sp, lib.loc = "/scratch/hgallo/Rlibs/"))
clusterEvalQ(cl,library(raster, lib.loc = "/scratch/hgallo/Rlibs/"))
clusterEvalQ(cl,library(iterators, lib.loc = "/scratch/hgallo/Rlibs/"))
clusterEvalQ(cl,library(igraph, lib.loc = "/scratch/hgallo/Rlibs/"))
clusterEvalQ(cl,library(gdistance, lib.loc = "/scratch/hgallo/Rlibs/"))
clusterEvalQ(cl,library(crayon, lib.loc = "/scratch/hgallo/Rlibs/"))
clusterEvalQ(cl,library(dplyr, lib.loc = "/scratch/hgallo/Rlibs/"))
clusterEvalQ(cl,library(data.table, lib.loc = "/scratch/hgallo/Rlibs/"))

# use shortestPaths function in parallel
costDist <- parLapply(cl=cl,cell_list, shortestPaths)

# to relate this back to the full distance matrix 
# we have to find the cells that were removed as duplicates and join the two tables
# create our own funciton that uses the left join in data.table package
l_join <- function(x){data.table(to=toCells, val=x[1,])[data.table(to=fromCells), on="to"]$val}
# export variables and functions needed
clusterExport(cl, c("costDist", "toCells", "fromCells", "l_join", "data.table"))
# left join across each element of the costDist list
# this returns a list of vectors that is length nrow(coords)
D <- parLapply(cl, costDist, l_join)
# unlist and rbind so that it matches our full distance matrix
D <- do.call("rbind", D)

# divide by 1000 to scale to kilometers from meters
D <- D/1000
# D goes into our gamma calculation
G <- gamma0*exp(-D^2/(2*sigma^2))
G[is.na(G)] <- 0

# incorporate spatially-explicit gamma into occupancy model
z <- matrix(0, nsite, nseason)

# initial z at t = 1
z[,1] <- rbinom(nsite, 1, psi[,1]) 

# species is there if we detected it
z[which(anyDetections[,1] == 1),1] <- 1 

# likelihood of first season
ll.z[,1] <- dbinom(z[,1], 1, psi[,1], log = TRUE)
ll.y[,1] <- dbinom(data$y_mat[,1], data$j_mat[,1], z[1:nsampled,1]*p[1], log=TRUE)

# generate z and y for season t > 1 and get log likelihood
for(k in 2:nseason) {
  zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
  PrNotColonizedByNeighbor <- 1 - G*zkt
  PrNotColonizedAtAll <- parApply(cl,PrNotColonizedByNeighbor, 1, prod, na.rm=TRUE)
  gamma[,k-1] <- 1 - PrNotColonizedAtAll
  psi[,k] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) + (1-z[,k-1])*gamma[,k-1] #Rescue effect
  z[,k] <- rbinom(nsite, 1, psi[,k])
  z[which(anyDetections[,k] == 1),k] <- 1
  ll.z[,k] <- dbinom(z[,k], 1, psi[,k], log=TRUE)
  # observation model
  ll.y[,k] <- dbinom(data$y_mat[,k], data$j_mat[,k], z[1:nsampled,k]*p[k], log=TRUE)
}

# Initial candidate values based off of initials
ll.z.cand <- ll.z
ll.z.sum <- sum(ll.z, na.rm=TRUE)
ll.y.cand <- ll.y
ll.y.sum <- sum(ll.y, na.rm=TRUE)
gamma.cand <- gamma
psi.cand <- psi

# Used to compute expected occupancy at each site
nz1 <- z

##### STARTING UPDATING PROCESS ####

# The number of parameters to estimate and zk for 13 seasons
npar <- length(param_mon) + nseason-1

# Objects to hold posterior samples
# Holds posterior samples
samples <- matrix(NA, iters, npar)
colnames(samples) <- c(paste0("alpha",1:length(alpha)),
                              "sigma",
                              "b0.gam", paste0("b.gam",1:length(b.gam)),
                              "b0.psi1", paste0("b.psi1.",1:length(b.psi1)), 
                              "b0.eps", paste0("b.eps",1:length(b.eps)),
                              "a0", paste0("season",2:4), 
                              paste0("zk",1:nseason),
                              "deviance")

# monitor z estimates for each patch
zk <- matrix(NA, nsite, iters)

# indicator that report it parameter is true
reportit <- report > 0

# create array that will hold z values
zA <- NULL
if(monitor.z){
    zA <- array(NA, c(nsite, nseason, iters))
}
# report the first iteration of starting values
if(reportit){
  cat("iter 1\n")
  cat("params =", round(c(alpha, b0.psi1, b.psi1, b0.gam, b.gam, sigma, 
                          b0.eps, b.eps, a0, season[2:4]), 5), "\n")
  cat("z[k] =", round(colSums(z), 2), "\n")
  cat("sum_ll.z =", round(sum(ll.z), 2), "\n")
  cat("deviance =", round(-2*ll.y.sum, 2), "\n")
  cat("time =", format(Sys.time()), "\n")
  if(plot.z) {
      library(lattice)
      zd <- data.frame(z=as.integer(z), season=rep(1:nseason, each=nsite),
                       x=as.numeric(data$coords[,1])/1000, y=as.numeric(data$coords[,2])/1000)
      print(xyplot(y ~ x | season, zd, groups=z, aspect="iso", pch=c(1,16), 
                   as.table=TRUE))
  }
}

#### START SAMPLING FROM POSTERIOR #####

for(s in 1:iters) {
  
  # update ll.z.sum from the previous iteration
  ll.z.sum <- sum(ll.z, na.rm=TRUE) ## This is important!
  
  # report information from previous iteration
  if(reportit) {
    if(s %in% c(2:100) | s %% 10 == 0){
      cat("iter", s, "\n")
      cat("theta =", round(samples[s-1,], 5), "\n")
      cat("z[k] =", zk, "\n")
      cat("accepted", round((zkup/nsite)*100, 1), "percent of z[k] proposals \n")
      cat("sum(ll.z) =", ll.z.sum, "\n")
      cat("deviance =", round(samples[s-1,38], 2), "\n")
      cat("time =", format(Sys.time()), "\n")
      if(plot.z) {
        library(lattice)
        zd$z <- as.integer(z)
        print(xyplot(y ~ x | season, zd, groups=z, aspect="iso", pch=c(1,16), 
                     as.table=TRUE))
      }
    }
  }
  #
  # Metropolis updates. We need to randomly subsample the parameters
  #  in a different order each step in the MCMC chain. To do so we
  #  generate a sequence from 1 to length(tune) and randomly order it. 
  #  The object 'tune' is the number of tuning parameters set for the
  #  metropolis algorithms
  #  We've numbered each of algorithms below from 1 to length(tune) and they
  #  will be sampled at step `subiter` of the vector sampling_order.
  sampling_order <- sample(seq(1, length(tune), 1), length(tune))
  for(subiter in 1:length(tune)){
    ## Metropolis update for alpha[1]
    if(sampling_order[subiter] == 1){
      alpha1.cand <- rnorm(1, alpha[1], tune[1])
      # create resistance surface
      cost <- exp(alpha1.cand*data$res_covs[[1]] + # 1-NDVI
                  alpha[2]*data$res_covs[[2]] +    # population density
                  alpha[3]*data$res_covs[[3]])     # patch indicator
      
      # calculate the ecological distance matrix in parallel
      # divide by 1000 to scale to km
      # note that transition function is 1/mean(x)

      # create the transition function and cell lists to be used in calculating the
      # shortest paths
      adj_graph <- getAllEdgeWeights(cost, directions=8, coords = data$coords)
      
      # have to make variables their own stand alone variable to export to cluster
      adjacencyGraph <- adj_graph$adjacencyGraph
      cell_list <- adj_graph$cell_list
      toCells <- adj_graph$toCells
      fromCells <- adj_graph$fromCells
      rm(adj_graph) # remove function object to save memory
      
      # create wrapper for shortest.paths function to go into lApply
      shortestPaths <- function(x){ shortest.paths(adjacencyGraph,
                                                   v = x[[1]],
                                                   to = x[[2]],
                                                   mode = "out",
                                                   algorithm = "dijkstra") }
      
      # export function and variables to clusters
      clusterExport(cl=cl, list("adjacencyGraph", "shortestPaths", 
                                "shortest.paths"))
      
      # use shortestPaths function in parallel
      costDist <- parLapply(cl=cl,cell_list, shortestPaths)
      
      # to relate this back to the full distance matrix 
      # export new variables and functions needed
      clusterExport(cl, c("costDist", "toCells", "fromCells", "l_join", "data.table"))
      # left join across each element of the costDist list
      # this returns a list of vectors that is length nrow(coords)
      D.cand <- parLapply(cl, costDist, l_join)
      # unlist and rbind so that it matches our full distance matrix
      D.cand <- do.call("rbind", D.cand)
      
      # divide by 1000 to scale to kilometers from meters
      D.cand <- D.cand/1000
      
      # incorporte into gamma calculation
      G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2)) # NA's = distances too far to colonize
      G.cand[is.na(G.cand)] <- 0 # change NA's to 0
      # model
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
      for(k in 2:nseason) {
        zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) +
          (1-z[,k-1])*gamma.cand[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      
      # priors
      prior.alpha1.cand <- dnorm(alpha1.cand, 0, 2, log=TRUE)
      prior.alpha1 <- dnorm(alpha[1], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.alpha1.cand) -
                          (ll.z.sum + prior.alpha1))) {
          alpha[1] <- alpha1.cand
          D <- D.cand
          G <- G.cand
          gamma <- gamma.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close alpha[1] sampler
    
    ## Metropolis update for alpha[2]
    if(sampling_order[subiter] == 2){
      alpha2.cand <- rnorm(1, alpha[2], tune[2])
      # create resistance surface
      cost <- exp(alpha[1]*data$res_covs[[1]] + 
                    alpha2.cand*data$res_covs[[2]] + 
                    alpha[3]*data$res_covs[[3]]) 
      # calculate the ecological distance matrix in parallel
      # divide by 1000 to scale to km
      # note that transition function is 1/mean(x)
      # create the transition function and cell lists to be used in calculating the
      # shortest paths
      adj_graph <- getAllEdgeWeights(cost, directions=8, coords = data$coords)
      
      # have to make variables their own stand alone variable to export to cluster
      adjacencyGraph <- adj_graph$adjacencyGraph
      cell_list <- adj_graph$cell_list
      toCells <- adj_graph$toCells
      fromCells <- adj_graph$fromCells
      rm(adj_graph) # remove function object to save memory
      
      # create wrapper for shortest.paths function to go into lApply
      shortestPaths <- function(x){ shortest.paths(adjacencyGraph,
                                                   v = x[[1]],
                                                   to = x[[2]],
                                                   mode = "out",
                                                   algorithm = "dijkstra") }
      
      # export function and variables to clusters
      clusterExport(cl=cl, list("adjacencyGraph", "shortestPaths", 
                                "shortest.paths"))
      
      # use shortestPaths function in parallel
      costDist <- parLapply(cl=cl,cell_list, shortestPaths)
      
      # to relate this back to the full distance matrix 
      # export new variables and functions needed
      clusterExport(cl, c("costDist", "toCells", "fromCells", "l_join", "data.table"))
      # left join across each element of the costDist list
      # this returns a list of vectors that is length nrow(coords)
      D.cand <- parLapply(cl, costDist, l_join)
      # unlist and rbind so that it matches our full distance matrix
      D.cand <- do.call("rbind", D.cand)
      # divide by 1000 to scale to kilometers from meters
      D.cand <- D.cand/1000
      G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2)) # NA's are distances that were too far to colonize
      G.cand[is.na(G.cand)] <- 0 # change NA's to 0
      # model
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
      for(k in 2:nseason) {
        zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) +
          (1-z[,k-1])*gamma.cand[,k-1]
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      
      # priors
      prior.alpha2.cand <- dnorm(alpha2.cand, 0, 2, log=TRUE)
      prior.alpha2 <- dnorm(alpha[2], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.alpha2.cand) -
                          (ll.z.sum + prior.alpha2))) {
          alpha[2] <- alpha2.cand
          D <- D.cand
          G <- G.cand
          gamma <- gamma.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler for alpha[2]
  
    ## Metropolis update for alpha[3]
    if(sampling_order[subiter] == 3){
      alpha3.cand <- rnorm(1, alpha[3], tune[3])
      # create resistance surface
      cost <- exp(alpha[1]*data$res_covs[[1]] + 
                    alpha[2]*data$res_covs[[2]] + 
                    alpha3.cand*data$res_covs[[3]])
      # calculate the ecological distance matrix in parallel
      # divide by 1000 to scale to km
      # note that transition function is 1/mean(x)
      # create the transition function and cell lists to be used in calculating the
      # shortest paths
      adj_graph <- getAllEdgeWeights(cost, directions=8, coords = data$coords)
      
      # have to make variables their own stand alone variable to export to cluster
      adjacencyGraph <- adj_graph$adjacencyGraph
      cell_list <- adj_graph$cell_list
      toCells <- adj_graph$toCells
      fromCells <- adj_graph$fromCells
      rm(adj_graph) # remove function object to save memory
      
      # create wrapper for shortest.paths function to go into lApply
      shortestPaths <- function(x){ shortest.paths(adjacencyGraph,
                                                   v = x[[1]],
                                                   to = x[[2]],
                                                   mode = "out",
                                                   algorithm = "dijkstra") }
      
      # export function and variables to clusters
      clusterExport(cl=cl, list("adjacencyGraph", "shortestPaths", 
                                "shortest.paths"))
      
      # use shortestPaths function in parallel
      costDist <- parLapply(cl=cl,cell_list, shortestPaths)
      
      # to relate this back to the full distance matrix 
      # export new variables and functions needed
      clusterExport(cl, c("costDist", "toCells", "fromCells", "l_join", "data.table"))
      # left join across each element of the costDist list
      # this returns a list of vectors that is length nrow(coords)
      D.cand <- parLapply(cl, costDist, l_join)
      # unlist and rbind so that it matches our full distance matrix
      D.cand <- do.call("rbind", D.cand)
      # divide by 1000 to scale to kilometers from meters
      D.cand <- D.cand/1000
      
      G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2)) # NA's are distances that were too far to colonize
      G.cand[is.na(G.cand)] <- 0 # change NA's to 0
      # model
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
      for(k in 2:nseason) {
        zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) +
          (1-z[,k-1])*gamma.cand[,k-1] 
       ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      
      # priors
      prior.alpha3.cand <- dnorm(alpha3.cand, 0, 2, log=TRUE)
      prior.alpha3 <- dnorm(alpha[3], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.alpha3.cand) -
                          (ll.z.sum + prior.alpha3))) {
          alpha[3] <- alpha3.cand
          D <- D.cand
          G <- G.cand
          gamma <- gamma.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler for alpha[3]
  
    ## Metropolis update for b0.gam - part of gamma0 linear model
    if(sampling_order[subiter] == 4){
      b0.gam.cand <- rnorm(1, b0.gam, tune[4])
      gamma0.cand <- plogis(b0.gam.cand + b.gam[1]*data$sitecovs[,"size"] +
                            b.gam[2]*data$sitecovs[,"park"] +
                            b.gam[3]*data$sitecovs[,"cem"] +
                            b.gam[4]*data$sitecovs[,"golf"])
      G.cand <- gamma0.cand*exp(-D^2/(2*sigma^2))
      G.cand[is.na(G.cand)] <- 0 # change NA's to 0
      # model
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
      for(k in 2:nseason) {
        zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + # w/ rescue effect 
          (1-z[,k-1])*gamma.cand[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b0.gam.cand <- dnorm(b0.gam.cand, 0, 2, log=TRUE)
      prior.b0.gam <- dnorm(b0.gam, 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b0.gam.cand) -
                          (ll.z.sum + prior.b0.gam))) {
          b0.gam <- b0.gam.cand
          gamma0 <- gamma0.cand
          G <- G.cand
          gamma <- gamma.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b0.gam
  
    ## Metropolis update for b.gam[1] - part of gamma0 linear model
    if(sampling_order[subiter] == 5){
      b1.gam.cand <- rnorm(1, b.gam[1], tune[5])
      gamma0.cand <- plogis(b0.gam + b1.gam.cand*data$sitecovs[,"size"] +
                              b.gam[2]*data$sitecovs[,"park"] + 
                              b.gam[3]*data$sitecovs[,"cem"] + 
                              b.gam[4]*data$sitecovs[,"golf"])
      G.cand <- gamma0.cand*exp(-D^2/(2*sigma^2))
      G.cand[is.na(G.cand)] <- 0 # change NA's to 0
      # model
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
      for(k in 2:nseason) {
        zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) +
          (1-z[,k-1])*gamma.cand[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b1.gam.cand <- dnorm(b1.gam.cand, 0, 2, log=TRUE)
      prior.b1.gam <- dnorm(b.gam[1], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b1.gam.cand) -
                          (ll.z.sum + prior.b1.gam))) {
          b.gam[1] <- b1.gam.cand
          gamma0 <- gamma0.cand
          G <- G.cand
          gamma <- gamma.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.gam[1]
  
    ## Metropolis update for b.gam[2] - part of gamma0 linear model
    if(sampling_order[subiter] == 6){
      b2.gam.cand <- rnorm(1, b.gam[2], tune[6])
      gamma0.cand <- plogis(b0.gam + b.gam[1]*data$sitecovs[,"size"] + 
                                    b2.gam.cand*data$sitecovs[,"park"] +
                                    b.gam[3]*data$sitecovs[,"cem"] +
                                    b.gam[4]*data$sitecovs[,"golf"])
      G.cand <- gamma0.cand*exp(-D^2/(2*sigma^2))
      G.cand[is.na(G.cand)] <- 0 # change NA's to 0
      # model
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
      for(k in 2:nseason) {
        zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + # w/ rescue effect
          (1-z[,k-1])*gamma.cand[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b2.gam.cand <- dnorm(b2.gam.cand, 0, 2, log=TRUE)
      prior.b2.gam <- dnorm(b.gam[2], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b2.gam.cand) -
                          (ll.z.sum + prior.b2.gam))) {
          b.gam[2] <- b2.gam.cand
          gamma0 <- gamma0.cand
          G <- G.cand
          gamma <- gamma.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.gam[2]
  
    ## Metropolis update for b.gam[3] - part of gamma0 linear model
    if(sampling_order[subiter] == 7){
      b3.gam.cand <- rnorm(1, b.gam[3], tune[7])
      gamma0.cand <- plogis(b0.gam + b.gam[1]*data$sitecovs[,"size"] + 
                                    b.gam[2]*data$sitecovs[,"park"] + 
                                    b3.gam.cand*data$sitecovs[,"cem"] +
                                    b.gam[4]*data$sitecovs[,"golf"])
      G.cand <- gamma0.cand*exp(-D^2/(2*sigma^2))
      G.cand[is.na(G.cand)] <- 0 # change NA's to 0
      # model
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
      for(k in 2:nseason) {
        zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + # w/ rescue effect 
          (1-z[,k-1])*gamma.cand[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b3.gam.cand <- dnorm(b3.gam.cand, 0, 2, log=TRUE)
      prior.b3.gam <- dnorm(b.gam[3], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
     
        if(runif(1) < exp((ll.z.sum.cand + prior.b3.gam.cand) -
                          (ll.z.sum + prior.b3.gam))) {
          b.gam[3] <- b3.gam.cand
          gamma0 <- gamma0.cand
          G <- G.cand
          gamma <- gamma.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.gam[3]
  
    ## Metropolis update for b.gam[4] - part of gamma0 linear model
    if(sampling_order[subiter] == 8){
      b4.gam.cand <- rnorm(1, b.gam[4], tune[8])
      gamma0.cand <- plogis(b0.gam + b.gam[1]*data$sitecovs[,"size"] + 
                                    b.gam[2]*data$sitecovs[,"park"] +
                                    b.gam[3]*data$sitecovs[,"cem"] + 
                                    b4.gam.cand*data$sitecovs[,"golf"])
      G.cand <- gamma0.cand*exp(-D^2/(2*sigma^2))
      G.cand[is.na(G.cand)] <- 0 # change NA's to 0
      # model
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
      for(k in 2:nseason) {
        zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + # w/ rescue effect 
          (1-z[,k-1])*gamma.cand[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b4.gam.cand <- dnorm(b4.gam.cand, 0, 2, log=TRUE)
      prior.b4.gam <- dnorm(b.gam[4], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b4.gam.cand) -
                          (ll.z.sum + prior.b4.gam))) {
          b.gam[4] <- b4.gam.cand
          gamma0 <- gamma0.cand
          G <- G.cand
          gamma <- gamma.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler for b.gam[4]
  
  ## Metropolis update for sigma
    if(sampling_order[subiter] == 9){
      sigma.cand <- rnorm(1, sigma, tune[9])
      if(sigma.cand > 0){
        G.cand <- gamma0*exp(-D^2/(2*sigma.cand^2))
        G.cand[is.na(G.cand)] <- 0 # change NA's to 0
        # model
        ll.z.cand[,1] <- dbinom(z[,1], 1, psi[,1], log=TRUE)
        for(k in 2:nseason) {
          zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
          gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
          psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) +
            (1-z[,k-1])*gamma.cand[,k-1] 
          ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
        }
        # priors
        prior.sigma.cand <- dgamma(sigma.cand, 0.001, 0.001)
        prior.sigma <- dgamma(sigma, 0.001, 0.001)
        # update
        ll.z.sum.cand <- sum(ll.z.cand)
        
          if(runif(1) < exp((ll.z.sum.cand + prior.sigma.cand) -
                            (ll.z.sum + prior.sigma))) {
            sigma <- sigma.cand
            G <- G.cand
            gamma <- gamma.cand
            psi <- psi.cand
            ll.z <- ll.z.cand
            ll.z.sum <- ll.z.sum.cand
          }
        
      }
    } # close sampler for sigma
  
  ## Metropolis update for b0.psi1 - part of linear predictor for initial occupancy
    if(sampling_order[subiter] == 10){
      b0.psi1.cand <- rnorm(1, b0.psi1, tune[10])
      # model
      psi1.cand <- plogis(b0.psi1.cand + b.psi1[1]*data$sitecovs[,"park"] +
                            b.psi1[2]*data$sitecovs[,"cem"] + 
                            b.psi1[3]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1.cand
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) + (1-z[,k-1])*gamma[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b0.psi1.cand <- dnorm(b0.psi1.cand, 0, 2, log=TRUE)
      prior.b0.psi1 <- dnorm(b0.psi1, 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b0.psi1.cand) -
                          (ll.z.sum + prior.b0.psi1))) {
          b0.psi1 <- b0.psi1.cand
          psi1 <- psi1.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # Close sampler for b0.psi1
  
  ## Metropolis update for b.psi1[1] - part of linear predictor for initial occupancy
    if(sampling_order[subiter] == 11){
      b1.psi1.cand <- rnorm(1, b.psi1[1], tune[11])
      # model
      psi1.cand <- plogis(b0.psi1 + b1.psi1.cand*data$sitecovs[,"park"] +
                            b.psi1[2]*data$sitecovs[,"cem"] + 
                            b.psi1[3]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1.cand
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) +
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b1.psi1.cand <- dnorm(b1.psi1.cand, 0, 2, log=TRUE)
      prior.b1.psi1 <- dnorm(b.psi1[1], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b1.psi1.cand) -
                          (ll.z.sum + prior.b1.psi1))) {
          b.psi1[1] <- b1.psi1.cand
          psi1 <- psi1.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # Close sampler for b.psi[1]
  
  ## Metropolis update for b.psi1[2] - part of linear predictor for initial occupancy
    if(sampling_order[subiter] == 12){
      b2.psi1.cand <- rnorm(1, b.psi1[2], tune[12])
      # model
      psi1.cand <- plogis(b0.psi1 + b.psi1[1]*data$sitecovs[,"park"] + 
                            b2.psi1.cand*data$sitecovs[,"cem"] +
                            b.psi1[3]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1.cand
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) + # w/ rescue effect
          (1-z[,k-1])*gamma[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b2.psi1.cand <- dnorm(b2.psi1.cand, 0, 2, log=TRUE)
      prior.b2.psi1 <- dnorm(b.psi1[2], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b2.psi1.cand) -
                          (ll.z.sum + prior.b2.psi1))) {
          b.psi1[2] <- b2.psi1.cand
          psi1 <- psi1.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler for b.psi1[2]
  
  ## Metropolis update for b.psi1[3] - part of linear predictor for initial occupancy
    if(sampling_order[subiter] == 13){
      b3.psi1.cand <- rnorm(1, b.psi1[3], tune[13])
      # model
      psi1.cand <- plogis(b0.psi1 + b.psi1[1]*data$sitecovs[,"park"] + 
                            b.psi1[2]*data$sitecovs[,"cem"] +
                            b3.psi1.cand*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1.cand
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) + # w/ rescue effect
          (1-z[,k-1])*gamma[,k-1] 
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b3.psi1.cand <- dnorm(b3.psi1.cand, 0, 2, log=TRUE)
      prior.b3.psi1 <- dnorm(b.psi1[3], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
     
        if(runif(1) < exp((ll.z.sum.cand + prior.b3.psi1.cand) -
                          (ll.z.sum + prior.b3.psi1))) {
          b.psi1[3] <- b3.psi1.cand
          psi1 <- psi1.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler for b.psi1[3]
  
  ## Metropolis update for b0.eps - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 14){
      b0.eps.cand <- rnorm(1, b0.eps, tune[14])
      epsilon.cand <- plogis(b0.eps.cand + b.eps[1]*data$sitecovs[,"tree"] + 
                                          b.eps[2]*data$sitecovs[,"total_veg"] +
                                          b.eps[3]*data$sitecovs[,"size"] +
                                          b.eps[4]*data$sitecovs[,"pop10"] +
                                          b.eps[5]*data$sitecovs[,"water"] +
                                          b.eps[6]*data$sitecovs[,"park"] +
                                          b.eps[7]*data$sitecovs[,"cem"] +
                                          b.eps[8]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon.cand*(1-gamma[,k-1])) + 
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b0.eps.cand <- dnorm(b0.eps.cand, 0, 2, log=TRUE)
      prior.b0.eps <- dnorm(b0.eps, 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
     
        if(runif(1) < exp((ll.z.sum.cand + prior.b0.eps.cand) -
                          (ll.z.sum + prior.b0.eps))) {
          b0.eps <- b0.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b0.eps
 
  ## Metropolis update for b.eps[1] - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 15){
      b1.eps.cand <- rnorm(1, b.eps[1], tune[15])
      epsilon.cand <- plogis(b0.eps + b1.eps.cand*data$sitecovs[,"tree"] +
                                      b.eps[2]*data$sitecovs[,"total_veg"] +
                                      b.eps[3]*data$sitecovs[,"size"] +
                                      b.eps[4]*data$sitecovs[,"pop10"] +
                                      b.eps[5]*data$sitecovs[,"water"] +
                                      b.eps[6]*data$sitecovs[,"park"] +
                                      b.eps[7]*data$sitecovs[,"cem"] +
                                      b.eps[8]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon.cand*(1-gamma[,k-1])) +
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b1.eps.cand <- dnorm(b1.eps.cand, 0, 2, log=TRUE)
      prior.b1.eps <- dnorm(b.eps[1], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b1.eps.cand) -
                          (ll.z.sum + prior.b1.eps))) {
          b.eps[1] <- b1.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.eps[1]
  
  ## Metropolis update for b.eps[2] - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 16){
      b2.eps.cand <- rnorm(1, b.eps[2], tune[16])
      epsilon.cand <- plogis(b0.eps + b.eps[1]*data$sitecovs[,"tree"] +
                                      b2.eps.cand*data$sitecovs[,"total_veg"] +
                                      b.eps[3]*data$sitecovs[,"size"] +
                                      b.eps[4]*data$sitecovs[,"pop10"] +
                                      b.eps[5]*data$sitecovs[,"water"] +
                                      b.eps[6]*data$sitecovs[,"park"] +
                                      b.eps[7]*data$sitecovs[,"cem"] +
                                      b.eps[8]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon.cand*(1-gamma[,k-1])) +
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b2.eps.cand <- dnorm(b2.eps.cand, 0, 2, log=TRUE)
      prior.b2.eps <- dnorm(b.eps[2], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b2.eps.cand) -
                          (ll.z.sum + prior.b2.eps))) {
          b.eps[2] <- b2.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.eps[2]
  
  ## Metropolis update for b.eps[3] - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 17){
      b3.eps.cand <- rnorm(1, b.eps[3], tune[17])
      epsilon.cand <- plogis(b0.eps + b.eps[1]*data$sitecovs[,"tree"] +
                                      b.eps[2]*data$sitecovs[,"total_veg"] +
                                      b3.eps.cand*data$sitecovs[,"size"] +
                                      b.eps[4]*data$sitecovs[,"pop10"] +
                                      b.eps[5]*data$sitecovs[,"water"] +
                                      b.eps[6]*data$sitecovs[,"park"] +
                                      b.eps[7]*data$sitecovs[,"cem"] +
                                      b.eps[8]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon.cand*(1-gamma[,k-1])) + 
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b3.eps.cand <- dnorm(b3.eps.cand, 0, 2, log=TRUE)
      prior.b3.eps <- dnorm(b.eps[3], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b3.eps.cand) -
                          (ll.z.sum + prior.b3.eps))) {
          b.eps[3] <- b3.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.eps[3]
  
  ## Metropolis update for b.eps[4] - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 18){
      b4.eps.cand <- rnorm(1, b.eps[4], tune[18])
      epsilon.cand <- plogis(b0.eps + b.eps[1]*data$sitecovs[,"tree"] +
                                      b.eps[2]*data$sitecovs[,"total_veg"] +
                                      b.eps[3]*data$sitecovs[,"size"] +
                                      b4.eps.cand*data$sitecovs[,"pop10"] +
                                      b.eps[5]*data$sitecovs[,"water"] +
                                      b.eps[6]*data$sitecovs[,"park"] +
                                      b.eps[7]*data$sitecovs[,"cem"] +
                                      b.eps[8]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon.cand*(1-gamma[,k-1])) +
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b4.eps.cand <- dnorm(b4.eps.cand, 0, 2, log=TRUE)
      prior.b4.eps <- dnorm(b.eps[4], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
     
        if(runif(1) < exp((ll.z.sum.cand + prior.b4.eps.cand) -
                          (ll.z.sum + prior.b4.eps))) {
          b.eps[4] <- b4.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.eps[4]
    
  ## Metropolis update for b.eps[5] - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 19){
      b5.eps.cand <- rnorm(1, b.eps[5], tune[19])
      epsilon.cand <- plogis(b0.eps + b.eps[1]*data$sitecovs[,"tree"] +
                                      b.eps[2]*data$sitecovs[,"total_veg"] +
                                      b.eps[3]*data$sitecovs[,"size"] +
                                      b.eps[4]*data$sitecovs[,"pop10"] +
                                      b5.eps.cand*data$sitecovs[,"water"] +
                                      b.eps[6]*data$sitecovs[,"park"] +
                                      b.eps[7]*data$sitecovs[,"cem"] +
                                      b.eps[8]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon.cand*(1-gamma[,k-1])) +
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b5.eps.cand <- dnorm(b5.eps.cand, 0, 2, log=TRUE)
      prior.b5.eps <- dnorm(b.eps[5], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b5.eps.cand) -
                          (ll.z.sum + prior.b5.eps))) {
          b.eps[5] <- b5.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.eps[5]

  ## Metropolis update for b.eps[6] - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 20){
      b6.eps.cand <- rnorm(1, b.eps[6], tune[20])
      epsilon.cand <- plogis(b0.eps + b.eps[1]*data$sitecovs[,"tree"] +
                                      b.eps[2]*data$sitecovs[,"total_veg"] +
                                      b.eps[3]*data$sitecovs[,"size"] +
                                      b.eps[4]*data$sitecovs[,"pop10"] +
                                      b.eps[5]*data$sitecovs[,"water"] +
                                      b6.eps.cand*data$sitecovs[,"park"] +
                                      b.eps[7]*data$sitecovs[,"cem"] +
                                      b.eps[8]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon.cand*(1-gamma[,k-1])) +
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b6.eps.cand <- dnorm(b6.eps.cand, 0, 2, log=TRUE)
      prior.b6.eps <- dnorm(b.eps[6], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b6.eps.cand) -
                          (ll.z.sum + prior.b6.eps))) {
          b.eps[6] <- b6.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.eps[6]
  
  ## Metropolis update for b.eps[7] - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 21){
      b7.eps.cand <- rnorm(1, b.eps[7], tune[21])
      epsilon.cand <- plogis(b0.eps + b.eps[1]*data$sitecovs[,"tree"] +
                                      b.eps[2]*data$sitecovs[,"total_veg"] +
                                      b.eps[3]*data$sitecovs[,"size"] +
                                      b.eps[4]*data$sitecovs[,"pop10"] +
                                      b.eps[5]*data$sitecovs[,"water"] +
                                      b.eps[6]*data$sitecovs[,"park"] +
                                      b7.eps.cand*data$sitecovs[,"cem"] +
                                      b.eps[8]*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon.cand*(1-gamma[,k-1])) +
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b7.eps.cand <- dnorm(b7.eps.cand, 0, 2, log=TRUE)
      prior.b7.eps <- dnorm(b.eps[7], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
     
        if(runif(1) < exp((ll.z.sum.cand + prior.b7.eps.cand) -
                          (ll.z.sum + prior.b7.eps))) {
          b.eps[7] <- b7.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.eps[7]
  
  ## Metropolis update for b.eps[8] - part of the linear predictor for epsilon
    if(sampling_order[subiter] == 22){
      b8.eps.cand <- rnorm(1, b.eps[8], tune[22])
      epsilon.cand <- plogis(b0.eps + b.eps[1]*data$sitecovs[,"tree"] +
                                      b.eps[2]*data$sitecovs[,"total_veg"] +
                                      b.eps[3]*data$sitecovs[,"size"] +
                                      b.eps[4]*data$sitecovs[,"pop10"] +
                                      b.eps[5]*data$sitecovs[,"water"] +
                                      b.eps[6]*data$sitecovs[,"park"] +
                                      b.eps[7]*data$sitecovs[,"cem"] +
                                      b8.eps.cand*data$sitecovs[,"golf"])
      psi.cand[,1] <- psi1
      ll.z.cand[,1] <- dbinom(z[,1], 1, psi.cand[,1], log=TRUE)
      for(k in 2:nseason) {
        psi.cand[,k] <- z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) +
          (1-z[,k-1])*gamma[,k-1] # w/ rescue effect
        ll.z.cand[,k] <- dbinom(z[,k], 1, psi.cand[,k], log=TRUE)
      }
      # priors
      prior.b8.eps.cand <- dnorm(b8.eps.cand, 0, 2, log=TRUE)
      prior.b8.eps <- dnorm(b.eps[8], 0, 2, log=TRUE)
      # update
      ll.z.sum.cand <- sum(ll.z.cand)
      
        if(runif(1) < exp((ll.z.sum.cand + prior.b8.eps.cand) -
                          (ll.z.sum + prior.b8.eps))) {
          b.eps[8] <- b8.eps.cand
          epsilon <- epsilon.cand
          psi <- psi.cand
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
        }
      
    } # close sampler b.eps[8]
  
  ## Metropolis update for z
    if(sampling_order[subiter] == 23){
      zkup <- rep(0, nseason)
      zknown1 <- anyDetections[,1]==1
      zknown1[is.na(zknown1)] <- FALSE
  # update for seasons k+1
      for(k in 1:nseason) {
        zknown <- anyDetections[,k]==1
        zknown[is.na(zknown)] <- FALSE
        for(i in 1:nsite) {
          if(zknown[i]) next
          # z candidate values
          zk.cand <- z[,k] # grab a vector of z candidates from initial values for each season
          zk.cand[i] <- 1-z[i,k] # give it the oppposite, does that decrease or increase likelihood?
          # ll.y candidate
          ll.y.tmp <- ll.y.cand.tmp <- 0
          if(i <= nsampled) { # only use sampled sites
            ll.y.cand.tmp <- dbinom(data$y_mat[i,k], data$j_mat[i,k], zk.cand[i]*p[k], log=TRUE)
            ll.y.tmp <- sum(ll.y[i,k], na.rm=TRUE) # trick to turn NA's to 0
          }
          # prior for z
          # note from Howell et al (2018): Prior must be calculated for time k and k+1 b/c change in z affects both
          ll.z.cand[i,k] <- dbinom(zk.cand[i], 1, psi[i,k], log=TRUE)
          ll.z2 <- ll.z2.cand <- 0
          if(k < nseason) {
            zkt.cand <- matrix(zk.cand, nsite, nsite, byrow=TRUE)
            gamma.cand[,k] <- 1 - exp(rowSums(log(1-G*zkt.cand)))
            psi.cand[,k+1] <- (zk.cand*(1-epsilon*(1-gamma.cand[,k])) +
                               (1-zk.cand)*gamma.cand[,k])
            ll.z.cand[,k+1] <- dbinom(z[,k+1], 1, psi.cand[,k+1], log=TRUE)
            ll.z2 <- sum(ll.z[,k+1])
            ll.z2.cand <- sum(ll.z.cand[,k+1])
          }
          
            if(runif(1) < exp((sum(ll.y.cand.tmp, na.rm=TRUE) +
                               ll.z.cand[i,k] + ll.z2.cand) -
                              (ll.y.tmp + ll.z[i,k] + ll.z2))) {
              z[,k] <- zk.cand
              ll.z[i,k] <- ll.z.cand[i,k]
              if(k < nseason) {
                gamma[,k] <- gamma.cand[,k]
                psi[,k] <- psi.cand[,k]
                ll.z[,k] <- ll.z.cand[,k]
              }
              if(i <= nsampled) {
                ll.y[i,k] <- ll.y.cand.tmp
              }
              zkup[k] <- zkup[k] + 1
            }
        }
      }
      
      nz1 <- nz1+z
    } # close sampler z
  
  ## Metropolis update for a0 - part of detection probability
    if(sampling_order[subiter] == 23){
      a0.cand<-rnorm(1, a0, tune[23])
      p.cand <- plogis(a0.cand + season[data$season_vec])
      p.mat <- matrix(p, nsampled, nseason, byrow=TRUE)
      p.cand.mat <- matrix(p.cand, nsampled, nseason, byrow=TRUE)
  
      ll.y <- dbinom(data$y_mat, data$j_mat, z[1:nsampled,]*p.mat, log=TRUE)
      ll.y.cand <- dbinom(data$y_mat, data$j_mat, z[1:nsampled,]*p.cand.mat, log=TRUE)
      # priors
      prior.a0.cand <- dnorm(a0.cand, 0, 2, log=TRUE) 
      prior.a0 <- dnorm(a0, 0, 2, log=TRUE)
  
      ll.y.sum <- sum(ll.y, na.rm=TRUE)
      ll.y.sum.cand <- sum(ll.y.cand, na.rm=TRUE)
      
        if(runif(1) < exp((ll.y.sum.cand + prior.a0.cand) -
                          (ll.y.sum + prior.a0))) {
            a0 <- a0.cand
            p <- p.cand
            ll.y <- ll.y.cand
            ll.y.sum <- ll.y.sum.cand
        }
      
    } # close sampler a0

  #Update for season 2
    if(sampling_order[subiter] == 24){
      season2.cand.vec <- season
      season2.cand.vec[2] <- rnorm(1, season[2], tune[24])
      p.cand <- plogis(a0 + season2.cand.vec[data$season_vec])
      p.mat <- matrix(p, nsampled, nseason, byrow=TRUE)
      p.cand.mat <- matrix(p.cand, nsampled, nseason, byrow=TRUE)
      
      ll.y <- dbinom(data$y_mat, data$j_mat, z[1:nsampled,]*p.mat, log=TRUE)
      ll.y.cand <- dbinom(data$y_mat, data$j_mat, z[1:nsampled,]*p.cand.mat, log=TRUE)
      # priors
      prior.season2.cand <- dnorm(season2.cand.vec[2], 0, 2, log=TRUE) 
      prior.season2 <- dnorm(season[2], 0, 2, log=TRUE)
      
      ll.y.sum <- sum(ll.y, na.rm=TRUE)
      ll.y.sum.cand <- sum(ll.y.cand, na.rm=TRUE)
      
        if(runif(1) < exp((ll.y.sum.cand + prior.season2.cand) -
                          (ll.y.sum + prior.season2))) {
          season <- season2.cand.vec
          p <- p.cand
          ll.y <- ll.y.cand
          ll.y.sum <- ll.y.sum.cand
        }
      
    } # close sampler season2

  # Update for season 3
    if(sampling_order[subiter] == 25){
      season3.cand.vec <- season
      season3.cand.vec[3] <- rnorm(1, season[3], tune[25])
      p.cand <- plogis(a0 + season3.cand.vec[data$season_vec])
      p.mat <- matrix(p, nsampled, nseason, byrow=TRUE)
      p.cand.mat <- matrix(p.cand, nsampled, nseason, byrow=TRUE)
      
      ll.y <- dbinom(data$y_mat, data$j_mat, z[1:nsampled,]*p.mat, log=TRUE)
      ll.y.cand <- dbinom(data$y_mat, data$j_mat, z[1:nsampled,]*p.cand.mat, log=TRUE)
      # priors
      prior.season3.cand <- dnorm(season3.cand.vec[3], 0, 2, log=TRUE) 
      prior.season3 <- dnorm(season[3], 0, 2, log=TRUE)
      
      ll.y.sum <- sum(ll.y, na.rm=TRUE)
      ll.y.sum.cand <- sum(ll.y.cand, na.rm=TRUE)
      
        if(runif(1) < exp((ll.y.sum.cand + prior.season3.cand) -
                          (ll.y.sum + prior.season3))) {
          season <- season3.cand.vec
          p <- p.cand
          ll.y <- ll.y.cand
          ll.y.sum <- ll.y.sum.cand
        }
      
    } # close sampler season3
  
  # Update for season 4
    if(sampling_order[subiter] == 26){
      season4.cand.vec <- season
      season4.cand.vec[4] <- rnorm(1, season[4], tune[26])
      p.cand <- plogis(a0 + season4.cand.vec[data$season_vec])
      p.mat <- matrix(p, nsampled, nseason, byrow=TRUE)
      p.cand.mat <- matrix(p.cand, nsampled, nseason, byrow=TRUE)
      
      ll.y <- dbinom(data$y_mat, data$j_mat, z[1:nsampled,]*p.mat, log=TRUE)
      ll.y.cand <- dbinom(data$y_mat, data$j_mat, z[1:nsampled,]*p.cand.mat, log=TRUE)
      # priors
      prior.season4.cand <- dnorm(season4.cand.vec[4], 0, 2, log=TRUE) 
      prior.season4 <- dnorm(season[4], 0, 2, log=TRUE)
      
      ll.y.sum <- sum(ll.y, na.rm=TRUE)
      ll.y.sum.cand <- sum(ll.y.cand, na.rm=TRUE)
     
        if(runif(1) < exp((ll.y.sum.cand + prior.season4.cand) -
                          (ll.y.sum + prior.season4))) {
          season <- season4.cand.vec
          p <- p.cand
          ll.y <- ll.y.cand
          ll.y.sum <- ll.y.sum.cand
        }
      
    } # close sampler season4
  } # close subiter loop

  zk <- colSums(z)
  
  samples[s,] <- c(alpha, sigma, b0.gam, b.gam, b0.psi1, 
                   b.psi1, b0.eps, b.eps, a0, season[2:4], zk=zk, 
                   deviance=-2*ll.y.sum)
  
  #zK[,s] <- z[,nseason]
  if(monitor.z){
      zA[,,s] <- z
  }
} # closes sampler

final_state <- list(z=z, D=D, samples=samples[s,])

# remember to put zK back if we keep it
return_list <- list(samples=samples, final_state=final_state,
            zA=zA, Ez=nz1/iters, seed=.Random.seed)

saveRDS(return_list,paste0(Sys.time(),"_argo_run.rds"))
