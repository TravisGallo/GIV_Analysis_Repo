
## MCMC algorithm
occuConn <- function(y,            # nsampled x nseason matrix of detection data*
                     j,            # nsampled x nseason matrix of sampling occasions
                     x,            # nSites x 2 matrix of site coordinates. Note that nSampled will usually be <nSites*
                     disp_dist,    # dispersal distance
                     r_covs,       # list of resistance covariate
                     site_covs,    # site-level covariate
                     obs_covs,     # observation-level covariate
                     iters=10,     # MCMC iterations
                     tune,         # Tuning order: sigma,gamma0.i,gamma0.s,gamma0.p,eps.i,eps.s,eps.p,
                                   #               beta0,beta1,beta2,alpha[1],alpha[2] (12 in total)
                     param_mon,    # character vector of parameters to monitor
                     monitor.z=FALSE, # store each iteration of the z matrix?
                     report=0,      # Only report progress if >0
                     plot.z=FALSE,  # Plot the latent presence-absence state (if report>0)
                     n.cores)       # number of cores for parallel 
{
  ## Dimensions
  nsite <- nrow(x) # Number of possible sites instead of only the sites sampled
  nseason <- ncol(y)
  nsampled <- nrow(y)
  
  ## indicating real detections from data
  anyDetections <- matrix(FALSE, nsite, nseason)
  anyDetections[1:nsampled,] <- as.numeric(y > 0)
  
  ## initial values for sampler
  # resistence
  alpha <-   rgamma(3, 1, 4)
  # initial occupancy
  b0.psi1 <- rnorm(1, 0, 0.5)
  b.psi1 <-  rnorm(3, 0, 0.5)
  # colonization
  b0.gam <- rnorm(1,0,0.5) # intercept
  b.gam <- rnorm(4,0,0.5) # 4 covariates
  gamma0 <- plogis(b0.gam + b.gam[1]*site_covs[,"size"] + 
                     b.gam[2]*site_covs[,"park"] + 
                     b.gam[3]*site_covs[,"cem"] + 
                     b.gam[4]*site_covs[,"golf"])
  sigma <- runif(1,1,10)
  # extinction
  b0.eps <-  rnorm(1, 0, 0.5) # intercept
  b.eps <-   rnorm(8, 0, 0.5) # 8 covariates
  epsilon <- plogis(b0.eps + b.eps[1]*site_covs[ ,"tree"] + 
                             b.eps[2]*site_covs[ ,"total_veg"] + 
                             b.eps[3]*site_covs[ ,"size"] + 
                             b.eps[4]*site_covs[ ,"pop10"] +
                             b.eps[5]*site_covs[ ,"water"] +
                             b.eps[6]*site_covs[ ,"park"] +
                             b.eps[7]*site_covs[ ,"cem"] +
                             b.eps[8]*site_covs[ ,"golf"] )
  # detection
  p <- rep(0, nseason)
  a0 <- rnorm(1, 0, 0.5)
  season <- rnorm(4,0 , 0.5)
  season[1] <- 0
  p <- plogis(a0 + season[obs_covs])
  # need to run through the model to generate starting values for z[,k-1] based off z[,1] starting values,
  # psi, gamma, and likelihoods for z and y
  gamma <- matrix(NA, nsite, nseason-1)
  psi1 <- plogis(b0.psi1 + b.psi1[1]*site_covs[,"park"] +
                           b.psi1[2]*site_covs[,"cem"]  +
                           b.psi1[3]*site_covs[,"golf"])
  
  # psi matrix holds all occupancy probabilities in t = 1 and then the 
  #  associated colonization / extinction probabilities at t > 1.
  psi <- matrix(NA, nsite, nseason)
  
  #  occupancy prob season 1 from initial values
  psi[,1] <- psi1
  
  # log likelihood matrices for latent and observed states
  ll.z <- matrix(0, nsite, nseason)
  ll.y <- matrix(0, nsampled, nseason)
  
  # create resistance surface
  cost <- exp(alpha[1]*r_covs[[1]] + 
              alpha[2]*r_covs[[2]] + 
              alpha[3]*r_covs[[3]])
  
  cl <- makeCluster(detectCores())
  # use modified function that creates transition matrix, corrects for diagnols and
    # calculates the ecological distance matrix in parallel
  D <- distanceMatrix(cost, directions=16, symm = TRUE,
                      fromCoords=test_data$coords, toCoords=test_data$coords,
                      dist.cutoff=disp_dist, cluster=cl)
  
  # divide by 1000 to scale to kilometers from meters
  D <- D/1000
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
  ll.y[,1] <- dbinom(y[,1], j[,1], z[1:nsampled,1]*p[1], log=TRUE)
  
  # generate z and y for season t > 1 and get log likelihood
  for(k in 2:nseason) {
    zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
    PrNotColonizedByNeighbor <- 1 - G*zkt
    PrNotColonizedAtAll <- apply(PrNotColonizedByNeighbor, 1, prod, na.rm=TRUE)
    gamma[,k-1] <- 1 - PrNotColonizedAtAll
    psi[,k] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) + (1-z[,k-1])*gamma[,k-1] #Rescue effect
    z[,k] <- rbinom(nsite, 1, psi[,k])
    z[which(anyDetections[,k] == 1),k] <- 1
    ll.z[,k] <- dbinom(z[,k], 1, psi[,k], log=TRUE)
    # observation model
    ll.y[,k] <- dbinom(y[,k], j[,k], z[1:nsampled,k]*p[k], log=TRUE)
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
  
  # STARTING UPDATING PROCESS
  # Objects to hold posterior samples
  
  # The number of parameters to estimate and zk for 13 seasons
  npar <- length(param_mon) + nseason-1

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
  if(reportit) {
      cat("iter 1\n")
      cat("    params =", round(c(alpha, b0.psi1, b.psi1, b0.gam, b.gam, 
                                 sigma, b0.eps, b.eps, a0, season[2:4]), 5),
          "\n")
      cat("    z[k] =", round(colSums(z), 2), "\n")
      cat("    sum_ll.z =", round(sum(ll.z), 2), "\n")
      cat("    deviance =", round(-2*ll.y.sum, 2), "\n")
      cat("    time =", format(Sys.time()), "\n")
      if(plot.z) {
          library(lattice)
          zd <- data.frame(z=as.integer(z), season=rep(1:nseason, each=nsite),
                           x=as.numeric(x[,1])/1000, y=as.numeric(x[,2])/1000)
          print(xyplot(y ~ x | season, zd, groups=z, aspect="iso", pch=c(1,16), 
                       as.table=TRUE))
      }
  }

  
  # START SAMPLING FROM POSTERIOR
  for(s in 1:iters) {
    
    # update ll.z.sum from the previous iteration
    ll.z.sum <- sum(ll.z, na.rm=TRUE) ## This is important!
    
    # report information from previous iteration
    if(reportit) {
      if(s %in% c(2:100) || s %% report == 0) {
        cat("iter", s, "\n")
        cat("    theta =", round(samples[s-1,], 5), "\n")
        cat("    z[k] =", zk, "\n")
        cat("    accepted", round((zkup/nsite)*100, 1), 
            "percent of z[k] proposals \n")
        cat("    sum(ll.z) =", ll.z.sum, "\n")
        cat("    deviance =", round(samples[s-1,38], 2), "\n")
        cat("    time =", format(Sys.time()), "\n")
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
        cost <- exp(alpha1.cand*r_covs[[1]] + # 1-NDVI
                    alpha[2]*r_covs[[2]] +    # population density
                    alpha[3]*r_covs[[3]])     # patch indicator
        
        # calculate the ecological distance matrix in parallel
        # divide by 1000 to scale to km
        # note that transition function is 1/mean(x)
        D.cand <- distanceMatrix(cost, directions=16, symm = TRUE,
                                 fromCoords=test_data$coords, toCoords=test_data$coords,
                                 dist.cutoff=disp_dist, cluster=cl)/1000
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
        cost <- exp(alpha[1]*r_covs[[1]] + alpha2.cand*r_covs[[2]] + alpha[3]*r_covs[[3]]) 
        # calculate the ecological distance matrix in parallel
        # divide by 1000 to scale to km
        # note that transition function is 1/mean(x)
        D.cand <- distanceMatrix(cost, directions=16, symm = TRUE,
                                 fromCoords=test_data$coords, toCoords=test_data$coords,
                                 dist.cutoff=disp_dist, cluster=cl)/1000
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
        cost <- exp(alpha[1]*r_covs[[1]] + alpha[2]*r_covs[[2]] + alpha3.cand*r_covs[[3]])
        # calculate the ecological distance matrix in parallel
        # divide by 1000 to scale to km
        # note that transition function is 1/mean(x)
        D.cand <- distanceMatrix(cost, directions=16, symm = TRUE,
                                 fromCoords=test_data$coords, toCoords=test_data$coords,
                                 dist.cutoff=disp_dist, cluster=cl)/1000
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
        gamma0.cand <- plogis(b0.gam.cand + b.gam[1]*site_covs[,"size"] +
                              b.gam[2]*site_covs[,"park"] +
                              b.gam[3]*site_covs[,"cem"] +
                              b.gam[4]*site_covs[,"golf"])
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
        gamma0.cand <- plogis(b0.gam + b1.gam.cand*site_covs[,"size"] +
                                b.gam[2]*site_covs[,"park"] + 
                                b.gam[3]*site_covs[,"cem"] + 
                                b.gam[4]*site_covs[,"golf"])
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
        gamma0.cand <- plogis(b0.gam + b.gam[1]*site_covs[,"size"] + 
                                      b2.gam.cand*site_covs[,"park"] +
                                      b.gam[3]*site_covs[,"cem"] +
                                      b.gam[4]*site_covs[,"golf"])
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
        gamma0.cand <- plogis(b0.gam + b.gam[1]*site_covs[,"size"] + 
                                      b.gam[2]*site_covs[,"park"] + 
                                      b3.gam.cand*site_covs[,"cem"] +
                                      b.gam[4]*site_covs[,"golf"])
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
        gamma0.cand <- plogis(b0.gam + b.gam[1]*site_covs[,"size"] + 
                                      b.gam[2]*site_covs[,"park"] +
                                      b.gam[3]*site_covs[,"cem"] + 
                                      b4.gam.cand*site_covs[,"golf"])
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
        psi1.cand <- plogis(b0.psi1.cand + b.psi1[1]*site_covs[,"park"] +
                              b.psi1[2]*site_covs[,"cem"] + 
                              b.psi1[3]*site_covs[,"golf"])
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
        psi1.cand <- plogis(b0.psi1 + b1.psi1.cand*site_covs[,"park"] +
                              b.psi1[2]*site_covs[,"cem"] + 
                              b.psi1[3]*site_covs[,"golf"])
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
        psi1.cand <- plogis(b0.psi1 + b.psi1[1]*site_covs[,"park"] + 
                              b2.psi1.cand*site_covs[,"cem"] +
                              b.psi1[3]*site_covs[,"golf"])
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
        psi1.cand <- plogis(b0.psi1 + b.psi1[1]*site_covs[,"park"] + 
                              b.psi1[2]*site_covs[,"cem"] +
                              b3.psi1.cand*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps.cand + b.eps[1]*site_covs[,"tree"] + 
                                            b.eps[2]*site_covs[,"total_veg"] +
                                            b.eps[3]*site_covs[,"size"] +
                                            b.eps[4]*site_covs[,"pop10"] +
                                            b.eps[5]*site_covs[,"water"] +
                                            b.eps[6]*site_covs[,"park"] +
                                            b.eps[7]*site_covs[,"cem"] +
                                            b.eps[8]*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps + b1.eps.cand*site_covs[,"tree"] +
                                        b.eps[2]*site_covs[,"total_veg"] +
                                        b.eps[3]*site_covs[,"size"] +
                                        b.eps[4]*site_covs[,"pop10"] +
                                        b.eps[5]*site_covs[,"water"] +
                                        b.eps[6]*site_covs[,"park"] +
                                        b.eps[7]*site_covs[,"cem"] +
                                        b.eps[8]*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps + b.eps[1]*site_covs[,"tree"] +
                                        b2.eps.cand*site_covs[,"total_veg"] +
                                        b.eps[3]*site_covs[,"size"] +
                                        b.eps[4]*site_covs[,"pop10"] +
                                        b.eps[5]*site_covs[,"water"] +
                                        b.eps[6]*site_covs[,"park"] +
                                        b.eps[7]*site_covs[,"cem"] +
                                        b.eps[8]*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps + b.eps[1]*site_covs[,"tree"] +
                                        b.eps[2]*site_covs[,"total_veg"] +
                                        b3.eps.cand*site_covs[,"size"] +
                                        b.eps[4]*site_covs[,"pop10"] +
                                        b.eps[5]*site_covs[,"water"] +
                                        b.eps[6]*site_covs[,"park"] +
                                        b.eps[7]*site_covs[,"cem"] +
                                        b.eps[8]*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps + b.eps[1]*site_covs[,"tree"] +
                                        b.eps[2]*site_covs[,"total_veg"] +
                                        b.eps[3]*site_covs[,"size"] +
                                        b4.eps.cand*site_covs[,"pop10"] +
                                        b.eps[5]*site_covs[,"water"] +
                                        b.eps[6]*site_covs[,"park"] +
                                        b.eps[7]*site_covs[,"cem"] +
                                        b.eps[8]*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps + b.eps[1]*site_covs[,"tree"] +
                                        b.eps[2]*site_covs[,"total_veg"] +
                                        b.eps[3]*site_covs[,"size"] +
                                        b.eps[4]*site_covs[,"pop10"] +
                                        b5.eps.cand*site_covs[,"water"] +
                                        b.eps[6]*site_covs[,"park"] +
                                        b.eps[7]*site_covs[,"cem"] +
                                        b.eps[8]*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps + b.eps[1]*site_covs[,"tree"] +
                                        b.eps[2]*site_covs[,"total_veg"] +
                                        b.eps[3]*site_covs[,"size"] +
                                        b.eps[4]*site_covs[,"pop10"] +
                                        b.eps[5]*site_covs[,"water"] +
                                        b6.eps.cand*site_covs[,"park"] +
                                        b.eps[7]*site_covs[,"cem"] +
                                        b.eps[8]*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps + b.eps[1]*site_covs[,"tree"] +
                                        b.eps[2]*site_covs[,"total_veg"] +
                                        b.eps[3]*site_covs[,"size"] +
                                        b.eps[4]*site_covs[,"pop10"] +
                                        b.eps[5]*site_covs[,"water"] +
                                        b.eps[6]*site_covs[,"park"] +
                                        b7.eps.cand*site_covs[,"cem"] +
                                        b.eps[8]*site_covs[,"golf"])
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
        epsilon.cand <- plogis(b0.eps + b.eps[1]*site_covs[,"tree"] +
                                        b.eps[2]*site_covs[,"total_veg"] +
                                        b.eps[3]*site_covs[,"size"] +
                                        b.eps[4]*site_covs[,"pop10"] +
                                        b.eps[5]*site_covs[,"water"] +
                                        b.eps[6]*site_covs[,"park"] +
                                        b.eps[7]*site_covs[,"cem"] +
                                        b8.eps.cand*site_covs[,"golf"])
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
      # # update for seasons 1
      # for(i in 1:nsite){
      #   if(zknown1[i]) next # only estimate sites where we did not detect the species
      #   z1.cand <- z[,1]
    #   z1.cand[i] <- 1-z[i,1]
    #   ll.y1.tmp <- ll.y1.cand.tmp <- 0
    #   if(i <= nsampled) { # only use sampled sites
    #     ll.y1.cand.tmp <- dbinom(y[i,1], j[i,1], z1.cand[i]*p[1], log=TRUE)
    #     ll.y1.tmp <- sum(ll.y[i,1], na.rm=TRUE) # trick to turn all NA's to 0
    #   }
    #   ll.z1.cand <- dbinom(z1.cand[i], 1, psi[i,1], log=TRUE)
    #   ll.z1.tmp <- sum(ll.z[i,1], na.rm=TRUE)
    #   # proposal
    #   if(runif(1) < exp((sum(ll.y1.cand.tmp, na.rm=TRUE) + ll.z1.cand) -
    #                     (ll.y1.tmp + ll.z1.tmp))) {
    #     z[,1] <- z1.cand
    #     ll.z[i,1] <- ll.z1.cand
    #     if(i <= nsampled) {
    #       ll.y[i,1] <- ll.y1.cand.tmp
    #     }
    #     zkup[1] <- zkup[1] + 1
    #   }
    # }
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
              ll.y.cand.tmp <- dbinom(y[i,k], j[i,k], zk.cand[i]*p[k], log=TRUE)
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
        p.cand <- plogis(a0.cand + season[obs_covs])
        p.mat <- matrix(p, nsampled, nseason, byrow=TRUE)
        p.cand.mat <- matrix(p.cand, nsampled, nseason, byrow=TRUE)
    
        ll.y <- dbinom(y, j, z[1:nsampled,]*p.mat, log=TRUE)
        ll.y.cand <- dbinom(y, j, z[1:nsampled,]*p.cand.mat, log=TRUE)
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
        p.cand <- plogis(a0 + season2.cand.vec[obs_covs])
        p.mat <- matrix(p, nsampled, nseason, byrow=TRUE)
        p.cand.mat <- matrix(p.cand, nsampled, nseason, byrow=TRUE)
        
        ll.y <- dbinom(y, j, z[1:nsampled,]*p.mat, log=TRUE)
        ll.y.cand <- dbinom(y, j, z[1:nsampled,]*p.cand.mat, log=TRUE)
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
        p.cand <- plogis(a0 + season3.cand.vec[obs_covs])
        p.mat <- matrix(p, nsampled, nseason, byrow=TRUE)
        p.cand.mat <- matrix(p.cand, nsampled, nseason, byrow=TRUE)
        
        ll.y <- dbinom(y, j, z[1:nsampled,]*p.mat, log=TRUE)
        ll.y.cand <- dbinom(y, j, z[1:nsampled,]*p.cand.mat, log=TRUE)
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
        p.cand <- plogis(a0 + season4.cand.vec[obs_covs])
        p.mat <- matrix(p, nsampled, nseason, byrow=TRUE)
        p.cand.mat <- matrix(p.cand, nsampled, nseason, byrow=TRUE)
        
        ll.y <- dbinom(y, j, z[1:nsampled,]*p.mat, log=TRUE)
        ll.y.cand <- dbinom(y, j, z[1:nsampled,]*p.cand.mat, log=TRUE)
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
  }
  
  # stop cluster
  stopCluster(cl)
  
  final_state <- list(z=z, D=D, samples=samples[s,])
  
  # remember to put zK back if we keep it
  return(list(samples=samples, final_state=final_state,
              zA=zA, Ez=nz1/iters, seed=.Random.seed))
}

occuConnC <- cmpfun(occuConn)
