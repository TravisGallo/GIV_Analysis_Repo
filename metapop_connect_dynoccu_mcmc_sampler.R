library(compiler)

## MCMC algorithm
dynroccH <- function(y,            # nsampled x nseason matrix of detection data*
                     j,            # nsampled x nseason matrix of sampling occasions
                     x,            # nSites x 2 matrix of site coordinates. Note that nSampled will usually be <nSites*
                     disp_cutoff,  # dispersal distance cutoff
                     r_covs,       # list of resistance covariate
                     site_covs,    # site-level covariate
                     obs_covs,     # observation-level covariate
                     iters=10,     # MCMC iterations
                     tune,         # Tuning order: sigma,gamma0.i,gamma0.s,gamma0.p,eps.i,eps.s,eps.p,
                                   #               beta0,beta1,beta2,alpha[1],alpha[2] (12 in total)
                     param_mon,    # character vector of parameters to monitor
                    
                    inits=NULL,    # until you run algorithm, inits are based on what is given.
                    zProp=c("ind","vec"), # Update z matrix by either proposing z(i,k) or z(,k), respectively
                    zProbs=NULL,   # matrix of proposal probs use if zProp="vec"
                    monitor.z=FALSE, # store each iteration of the z matrix?
                    report=0,      # Only report progress if >0
                    plot.z=FALSE,  # Plot the latent presence-absence state (if report>0)
                    tol=0)      # This will reject a proposal of z(i,k)=1 if mu(i,k-1)<tol
{

  zProp <- match.arg(zProp)

  ## Dimensions
  nsite <- nrow(x) # Number of possible sites instead of only the sites sampled
  nseason <- ncol(y)
  nsampled <- nrow(y)
  
  ## Using this to avoid likelihood calculations for sites not sampled
  dataYears <- apply(!is.na(y), 3, any) #********************* need to see if this stays or is adapted ****************************
  
  ## indicating real detections from data
  anyDetections <- matrix(FALSE, nsite, nseason)
  anyDetections[1:nsampled,] <- as.numeric(y > 0)

  ## initial values
  # resistence
  alpha <- rnorm(3)
  # initial occupancy
  b0 <- rnorm(1)
  b <- rnorm(3)
  psi1 <- plogis(b0 + b[1]*site_covs[,"park"] + b[2]*site_covs[,"cem"] + b[3]*site_covs[,"golf"])
  z <- matrix(0, nsite, nseason)
  z[,1] <- rbinom(nsite, 1, psi1) # first season
  z[1:nsampled,1] <- as.numeric(y_mat[,1] > 0) # set actual detections
  # colonization
  g0 <- rnorm(1) # intercept
  g <- rnorm(4) # 4 covariates
  gamma0 <- plogis(g0 + g[1]*site_covs[,"size"] + g[2]*site_covs[,"park"] + g[3]*site_covs[,"cem"] + g[4]*site_covs[,"golf"])
  sigma <- runif(1,0,8)
  # extinction
  e0 <- rnorm(1) # intercept
  e <- rnorm(8) # 8 covariates
  epsilon <- plogis(e0 + e[1]*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                      e[5]*site_covs[,"water"] + e[6]*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
  # detection
  p <- rep(0, nseason)
  a0 <- rnorm(1)
  season <- rnorm(nseason)
  for(k in 1:nseason){
    p[k] <- plogis(a0 + obs_covs[season_vec[k]])
  }
  # need to run through the model to generate starting values for z[,k-1] based off z[,1] starting values,
  # psi, gamma, and likelihoods for z and y
  gamma <- psi <- matrix(NA, nsite, nseason-1)
  ll.z <- matrix(0, nsite, nseason)
  ll.y <- matrix(0, nsampled, nseason)
  # create resistance surface
  cost <- exp(alpha[1]*r_covs[[1]] + alpha[2]*r_covs[[2]] + alpha[3]*r_covs[[3]])*r_covs[[4]]
  # calculate conductances among neighbors
  tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16) 
  #adjust diag.conductances
  tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE, scl=FALSE) 
  #calculate the ecological distance matrix
  D <- costDistance(tr1CorrC,x,x)/1000
  G <- gamma0*exp(-D^2/(2*sigma^2))
  # incorporate spatially-explicit gamma into occupancy model
  for(k in 2:nseason) { 
    PrNotColonizedByNeighbor <- 1 - gamma0*exp(-D^2/(2*sigma^2)) * t(z[,rep(k-1, nsite)])
    PrNotColonizedAtAll <- apply(PrNotColonizedByNeighbor, 1, prod, na.rm=TRUE)
    gamma[,k-1] <- 1 - PrNotColonizedAtAll
    psi[,k-1] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) + (1-z[,k-1])*gamma[,k-1] #Rescue effect
    z[,k] <- rbinom(nsite, 1, psi[,1])
    z[which(anyDetections[,k] == 1),k] <- 1
    ll.z[,k-1] <- dbinom(z[,k], 1, psi[,k-1], log=TRUE)
    # observation model
    ll.y[,k-1] <- dbinom(y[,k-1], j[,k-1], z[1:nsampled,k-1]*p[k-1], log=TRUE)
  }
  ll.z.cand <- ll.z
  ll.z.sum <- sum(ll.z, na.rm=TRUE)
  ll.y.cand <- ll.y
  ll.y.sum <- sum(ll.y, na.rm=TRUE)
  gamma.cand <- gamma
  psi.cand <- psi
  
  # used to compute expected occupancy at each site
  nz1 <- z
  
  # used to monitor accpetance rate
  zkup <- rep(0, nseason-1)
  
  # STARTING UPDATING PROCESS
  # objects to hold posterior samples
  npar <- length(param_mon) + nseason # The number of parameters to estimate for each season
  samples <- array(NA, c(iters, npar))
  colnames(samples) <- param_mon
  # monitor z estimates for each patch
  zK <- matrix(NA, nsite, iters)
  
  # indicator that report it parameter is true
  reportit <- report > 0
  
  
  ###############################################
  ###############################################
 
  zA <- NULL
  if(monitor.z){
      zA <- array(NA_integer_, c(nsite, nseason, iters))
  }
  # report the first iteration of starting values
  if(reportit) {
      cat("iter 1\n")
      cat("    theta =", round(c(alpha, ps1, b0, b, gamma0, g0, g, sigma, epsilon, e0, e, p, a0, season), 5), "\n")
      cat("    z[k] =", round(colSums(z), 2), "\n")
      cat("    ll.z =", round(sum(ll.z), 2), "\n")
      cat("    deviance =", round(-2*ll.y.sum, 2), "\n")
      cat("    time =", format(Sys.time()), "\n")
      if(plot.z) { ###############*********************######################
          library(lattice)
          zd <- data.frame(z=as.integer(z), year=factor(rep(2003:2017, each=nSites)),
                           x=as.numeric(x[,1])/1000, y=as.numeric(x[,2])/1000)
          print(xyplot(y ~ x | year, zd, groups=z, aspect="iso", pch=c(1,16), as.table=TRUE))
      }
  }

  # START SAMPLING FROM POSTERIOR
  for(s in 1:nIter) {
    
    # update ll.z.sum from the previous iteration
    ll.z.sum <- sum(ll.z) ## This is important!
    
    # report information from previous iteration
    if(reportit) {
      if(s %in% c(2:100) || s %% report == 0) {
        cat("iter", s, "\n")
        cat("    theta =", round(samples[s-1,1:12], 5), "\n")
        cat("    z[k] =", zk, "\n")
        cat("    accepted", round(zkup/(nSites)*100, 1), "percent of z[k] proposals \n")
        cat("    sum(ll.z) =", ll.z.sum, "\n")
        cat("    deviance =", round(samples[s-1,"deviance"], 2), "\n")
        cat("    time =", format(Sys.time()), "\n")
        if(plot.z) {
          library(lattice)
          zd$z <- as.integer(z)
          print(xyplot(y ~ x | year, zd, groups=z, aspect="iso", pch=c(1,16), as.table=TRUE))
        }
        }
      }

    #Metropolis update for alpha[1]
    library(gdistance)
    alpha1.cand <- rnorm(1, alpha[1], tune[11])
    # create resistance surface
    cost <- exp(alpha.cand*r_covs[[1]] + alpha[2]*r_covs[[2]] + alpha[3]*r_covs[[3]])*r_covs[[4]]
    # calculate conductances among neighbors
    tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16) 
    # adjust diag. conductances
    tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE, scl=FALSE) 
    # calculate least cost distance between all pairs of sites.
    D.cand <- costDistance(tr1CorrC,x,x)/1000 #calculate the ecological distance matrix
    G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2))
    # model
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.alpha1.cand <- dnorm(alpha1.cand, 0, 10, log=TRUE)
    prior.alpha1 <- dnorm(alpha[1], 0, 10, log=TRUE)
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

    # Metropolis update for alpha[2]
    alpha2.cand <- rnorm(1, alpha[2], tune[12])
    # create resistance surface
    cost <- exp(alpha[1]*r_covs[[1]] + alpha2.cand*r_covs[[2]] + alpha[3]*r_covs[[3]])*r_covs[[4]]    
    # calculate conductances among neighbors
    tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16) 
    tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE,scl=FALSE) #adjust diag.conductances
    # calculate least cost distance between all pairs of sites.
    D.cand <- costDistance(tr1CorrC,x,x)/1000 #calculate the ecological distance matrix
    G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2))
    # model
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.alpha2.cand <- dnorm(alpha2.cand, 0, 10, log=TRUE)
    prior.alpha2 <- dnorm(alpha[2], 0, 10, log=TRUE)
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
    
    # Metropolis update for alpha[3]
    alpha3.cand <- rnorm(1, alpha[2], tune[12])
    # create resistance surface
    cost <- exp(alpha[1]*r_covs[[1]] + alpha[2]*r_covs[[2]] + alpha3.cand*r_covs[[3]])*r_covs[[4]]    ## calculate conductances among neighbors
    tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16) 
    tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE,scl=FALSE) #adjust diag.conductances
    # calculate least cost distance between all pairs of sites.
    D.cand <- costDistance(tr1CorrC,x,x)/1000 #calculate the ecological distance matrix
    G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2))
    # model
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.alpha3.cand <- dnorm(alpha3.cand, 0, 10, log=TRUE)
    prior.alpha3 <- dnorm(alpha[3], 0, 10, log=TRUE)
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

    # Metropolis update for g0 - part of gamma0 linear model
    g0.cand <- rnorm(1, b0, tune[12])
    gamma0.cand <- plogis(g0.cand + g[1]*site_covs[,"size"] + g[2]*site_covs[,"park"] + g[3]*site_covs[,"cem"] + g[4]*site_covs[,"golf"])
    G.cand <- gamma0.cand*exp(-D.cand^2/(2*sigma^2))
    # model
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.g0.cand <- dnorm(g0.cand, 0, 10, log=TRUE)
    prior.g0 <- dnorm(g0, 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.g0.cand) -
                      (ll.z.sum + prior.g0))) {
      g0 <- g0.cand
      D <- D.cand
      gamma0 <- gamma0.cand
      G <- G.cand
      gamma <- gamma.cand
      psi <- psi.cand
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
    }
    
    # Metropolis update for g[1] - part of gamma0 linear model
    g1.cand <- rnorm(1, g[1], tune[12])
    gamma0.cand <- plogis(g0 + g1.cand*site_covs[,"size"] + g[2]*site_covs[,"park"] + g[3]*site_covs[,"cem"] + g[4]*site_covs[,"golf"])
    G.cand <- gamma0.cand*exp(-D.cand^2/(2*sigma^2))
    # model
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.g1.cand <- dnorm(g1.cand, 0, 10, log=TRUE)
    prior.g1 <- dnorm(g[1], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.g1.cand) -
                      (ll.z.sum + prior.g1))) {
      g[1] <- g1.cand
      D <- D.cand
      gamma0 <- gamma0.cand
      G <- G.cand
      gamma <- gamma.cand
      psi <- psi.cand
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
    }
    
    # Metropolis update for g[2] - part of gamma0 linear model
    g2.cand <- rnorm(1, g[2], tune[12])
    gamma0.cand <- plogis(g0 + g[1]*site_covs[,"size"] + g2.cand*site_covs[,"park"] + g[3]*site_covs[,"cem"] + g[4]*site_covs[,"golf"])
    G.cand <- gamma0.cand*exp(-D.cand^2/(2*sigma^2))
    # model
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.g2.cand <- dnorm(g2.cand, 0, 10, log=TRUE)
    prior.g2 <- dnorm(g[2], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.g2.cand) -
                      (ll.z.sum + prior.g2))) {
      g[2] <- g2.cand
      D <- D.cand
      gamma0 <- gamma0.cand
      G <- G.cand
      gamma <- gamma.cand
      psi <- psi.cand
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
    }
    
    # Metropolis update for g[3] - part of gamma0 linear model
    g3.cand <- rnorm(1, g[3], tune[12])
    gamma0.cand <- plogis(g0 + g[1]*site_covs[,"size"] + g[2]*site_covs[,"park"] + g3.cand*site_covs[,"cem"] + g[4]*site_covs[,"golf"])
    G.cand <- gamma0.cand*exp(-D.cand^2/(2*sigma^2))
    # model
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.g3.cand <- dnorm(g3.cand, 0, 10, log=TRUE)
    prior.g3 <- dnorm(g[3], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.g3.cand) -
                      (ll.z.sum + prior.g3))) {
      g[3] <- g3.cand
      D <- D.cand
      gamma0 <- gamma0.cand
      G <- G.cand
      gamma <- gamma.cand
      psi <- psi.cand
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
    }
    
    # Metropolis update for g[4] - part of gamma0 linear model
    g4.cand <- rnorm(1, g[4], tune[12])
    gamma0.cand <- plogis(g0 + g[1]*site_covs[,"size"] + g[2]*site_covs[,"park"] + g[3]*site_covs[,"cem"] + g4.cand*site_covs[,"golf"])
    G.cand <- gamma0.cand*exp(-D.cand^2/(2*sigma^2))
    # model
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.g4.cand <- dnorm(g4.cand, 0, 10, log=TRUE)
    prior.g4 <- dnorm(g[4], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.g4.cand) -
                      (ll.z.sum + prior.g4))) {
      g[4] <- g4.cand
      D <- D.cand
      gamma0 <- gamma0.cand
      G <- G.cand
      gamma <- gamma.cand
      psi <- psi.cand
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
    }
    
    # Metropolis update for sigma
    sigma.cand <- abs(rnorm(1, sigma, tune[1]))
    G.cand <- gamma0.cand*exp(-D^2/(2*sigma.cand^2))
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.sigma.cand <- dgamma(sigma.cand, 0.001, 0.001)
    prior.sigma <- dgamma(sigma, 0.001, 0.001)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.sigma.cand) -
                      (ll.z.sum + prior.sigma))) {
      sigma <- sigma.cand
      gamma0 <- gamma0.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
      }
    
    # Metropolis update for e0 - part of the linear predictor for epsilon
    e0.cand <- rnorm(1, e0, tune[12])
    epsilon.cand <- plogis(e0.cand + e[1]*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                        e[5]*site_covs[,"water"] + e[6]*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e0.cand <- dnorm(e0.cand, 0, 10, log=TRUE)
    prior.e0 <- dnorm(e0, 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e0.cand) -
                      (ll.z.sum + prior.e0))) {
      e0 <- e0.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }
   
    # Metropolis update for e1 - part of the linear predictor for epsilon
    e1.cand <- rnorm(1, e[1], tune[12])
    epsilon.cand <- plogis(e0 + e1.cand*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                             e[5]*site_covs[,"water"] + e[6]*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e1.cand <- dnorm(e1.cand, 0, 10, log=TRUE)
    prior.e1 <- dnorm(e[1], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e1.cand) -
                      (ll.z.sum + prior.e1))) {
      e[1] <- e1.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }
    
    # Metropolis update for e2 - part of the linear predictor for epsilon
    e2.cand <- rnorm(1, e[2], tune[12])
    epsilon.cand <- plogis(e0 + e[1]*site_covs[,"tree"] + e2.cand*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                             e[5]*site_covs[,"water"] + e[6]*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e2.cand <- dnorm(e2.cand, 0, 10, log=TRUE)
    prior.e2 <- dnorm(e[2], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e2.cand) -
                      (ll.z.sum + prior.e2))) {
      e[2] <- e2.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }
    
    # Metropolis update for e3 - part of the linear predictor for epsilon
    e3.cand <- rnorm(1, e[3], tune[12])
    epsilon.cand <- plogis(e0 + e[1]*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e3.cand*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                             e[5]*site_covs[,"water"] + e[6]*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e3.cand <- dnorm(e3.cand, 0, 10, log=TRUE)
    prior.e3 <- dnorm(e[3], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e3.cand) -
                      (ll.z.sum + prior.e3))) {
      e[3] <- e3.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }
    
    # Metropolis update for e4 - part of the linear predictor for epsilon
    e4.cand <- rnorm(1, e[4], tune[12])
    epsilon.cand <- plogis(e0 + e[1]*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e4.cand*site_covs[,"pop10"] +
                             e[5]*site_covs[,"water"] + e[6]*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e4.cand <- dnorm(e4.cand, 0, 10, log=TRUE)
    prior.e4 <- dnorm(e[4], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e4.cand) -
                      (ll.z.sum + prior.e4))) {
      e[4] <- e4.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }
    
    # Metropolis update for e5 - part of the linear predictor for epsilon
    e5.cand <- rnorm(1, e[5], tune[12])
    epsilon.cand <- plogis(e0 + e[1]*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                             e5.cand*site_covs[,"water"] + e[6]*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e5.cand <- dnorm(e5.cand, 0, 10, log=TRUE)
    prior.e5 <- dnorm(e[5], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e5.cand) -
                      (ll.z.sum + prior.e5))) {
      e[5] <- e5.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }

    # Metropolis update for e6 - part of the linear predictor for epsilon
    e6.cand <- rnorm(1, e[6], tune[12])
    epsilon.cand <- plogis(e0 + e[1]*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                             e[5]*site_covs[,"water"] + e6.cand*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e6.cand <- dnorm(e6.cand, 0, 10, log=TRUE)
    prior.e6 <- dnorm(e[6], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e6.cand) -
                      (ll.z.sum + prior.e6))) {
      e[6] <- e6.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }
    
    # Metropolis update for e7 - part of the linear predictor for epsilon
    e7.cand <- rnorm(1, e[7], tune[12])
    epsilon.cand <- plogis(e0 + e[1]*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                             e[5]*site_covs[,"water"] + e[6]*site_covs[,"park"] + e7.cand*site_covs[,"cem"] + e[8]*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e7.cand <- dnorm(e7.cand, 0, 10, log=TRUE)
    prior.e7 <- dnorm(e[7], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e7.cand) -
                      (ll.z.sum + prior.e7))) {
      e[7] <- e7.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }
    
    # Metropolis update for e8 - part of the linear predictor for epsilon
    e8.cand <- rnorm(1, e[8], tune[12])
    epsilon.cand <- plogis(e0 + e[1]*site_covs[,"tree"] + e[2]*site_covs[,"total_veg"] + e[3]*site_covs[,"size"] + e[4]*site_covs[,"pop10"] +
                             e[5]*site_covs[,"water"] + e[6]*site_covs[,"park"] + e[7]*site_covs[,"cem"] + e8.cand*site_covs[,"golf"])
    for(k in 2:nseason) {
      zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon.cand*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    # priors
    prior.e8.cand <- dnorm(e8.cand, 0, 10, log=TRUE)
    prior.e8 <- dnorm(e[8], 0, 10, log=TRUE)
    # update
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.e8.cand) -
                      (ll.z.sum + prior.e8))) {
      e[8] <- e8.cand
      epsilon <- epsilon.cand
      gamma <- gamma.cand 
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      muz <- muz.cand
      G <- G.cand
    }
    
    
    ## update z
    ## We can update each z(i,t) individually, and it results in better mixing than updating a vector of z's
    zkup <- rep(0, nYears-1)
   for(k in 2:nseason) {
       anyDet <- anyDetections[,k]==1 

       prop.back <- prop.cand <- 0
       for(i in 1:nsite) {
           if(zknown[i])
               next
           ## Reject highly unlikely proposals (before proposing them)
           ## This speed trick shouldn't affect anything but
           ## can double check by changing toleranc (tol)
           if(z[i,k]<1 & psi[i,k-1]<tol)
               next
           zk.wide <- matrix(z[,k], nSites, nReps)
           zk.cand <- z[,k]
           zk.cand[i] <- 1-z[i,k]
           zk.cand.wide <- matrix(zk.cand, nSites, nReps)

           ll.y.tmp <- 0
           ll.y.cand.tmp <- 0
           if((k > 4) & (i <= nSampled)) { ## Ignore first 4 years without data
               ll.y.cand.tmp <- dbinom(y[i,k], 1, zk.cand[i]*p[i,,k], log=TRUE)
               ll.y.tmp <- sum(ll.y[i,k], na.rm=TRUE)
           }
           ## RC: Prior must be calculated for time k and k+1 b/c change in z affects both
           ll.z.cand[i,k-1] <- dbinom(zk.cand[i], 1, psi[i,k-1], log=TRUE)
           ll.z2 <- ll.z2.cand <- 0
           if(k < nseason) {
               zkt.cand <- matrix(zk.cand, nsite, nsite, byrow=TRUE)
               gamma.cand[,k] <- 1 - exp(rowSums(log(1-G*zkt.cand)))
               psi.cand[,k] <- (zk.cand*(1-epsilon*(1-gamma.cand[,k])) + (1-zk.cand)*gamma.cand[,k])
               ll.z.cand[,k] <- dbinom(z[,k+1], 1, psi.cand[,k], log=TRUE)
               ll.z2 <- sum(ll.z[,k])
               ll.z2.cand <- sum(ll.z.cand[,k])
           }
           if(runif(1) < exp((sum(ll.y.cand.tmp, na.rm=TRUE) + ll.z.cand[i,k-1] +
                              ll.z2.cand + prop.back) -
                             (ll.y.tmp + ll.z[i,k-1] +
                              ll.z2 + prop.cand))) {
               z[,k] <- zk.cand
               ll.z[i,k-1] <- ll.z.cand[i,k-1]
               if(k < nYears) {
                   gamma[,k] <- gamma.cand[,k]
                   psi[,k] <- psi.cand[,k]
                   ll.z[,k] <- ll.z.cand[,k]
               }
               if((i <= nSampled) & (k>4)) {
                   ll.y[i,,k] <- ll.y.cand.tmp
               }
               zkup[k-1] <- zkup[k-1] + 1
           }
       }
   }

    nz1 <- nz1+z

    #Update for beta0
    beta0.cand<-rnorm(1, beta0, tune[8])
    p.cand <- plogis(beta0.cand + beta1*p.cov1 + beta2*p.cov2)
    z.wide <- z[,rep(1:nYears, each=nReps)]
    z.a <- array(z.wide, c(nSites, nReps, nYears))

    ll.y[,,dataYears] <- dbinom(y[,,dataYears], 1, z.a[1:nSampled,,dataYears]*p[,,dataYears], log=TRUE)
    ll.y.cand[,,dataYears] <- dbinom(y[,,dataYears], 1, z.a[1:nSampled,,dataYears]*p.cand[,,dataYears], log=TRUE)
    prior.beta0.cand <- dnorm(beta0.cand, 0, 10, log=TRUE) 
    prior.beta0 <- dnorm(beta0, 0, 10, log=TRUE)

    ll.y.sum <- sum(ll.y, na.rm=TRUE)
    ll.y.sum.cand <- sum(ll.y.cand, na.rm=TRUE)
    if(runif(1) < exp((ll.y.sum.cand + prior.beta0.cand) -
                      (ll.y.sum + prior.beta0))) {
        beta0 <- beta0.cand
        p <- p.cand
        ll.y <- ll.y.cand
        ll.y.sum <- ll.y.sum.cand
    }

    #Update for beta1
    beta1.cand<-rnorm(1, beta1, tune[9])
    p.cand <- plogis(beta0 + beta1.cand*p.cov1 + beta2*p.cov2)
    z.wide <- z[,rep(1:nYears, each=nReps)]
    z.a <- array(z.wide, c(nSites, nReps, nYears))

    ll.y.cand[,,dataYears] <- dbinom(y[,,dataYears], 1, z.a[1:nSampled,,dataYears]*p.cand[,,dataYears], log=TRUE)
    prior.beta1.cand <- dnorm(beta1.cand, 0, 10, log=TRUE) 
    prior.beta1 <- dnorm(beta1, 0, 10, log=TRUE)

    ll.y.sum.cand <- sum(ll.y.cand, na.rm=TRUE)
    if(runif(1) < exp((ll.y.sum.cand + prior.beta0.cand) -
                      (ll.y.sum + prior.beta0))) {
        beta1 <- beta1.cand
        p <- p.cand
        ll.y <- ll.y.cand
        ll.y.sum <- ll.y.sum.cand
    }

    #Update for beta2
    beta2.cand<-rnorm(1, beta2, tune[10])
    p.cand <- plogis(beta0 + beta1*p.cov1 + beta2.cand*p.cov2)
    z.wide <- z[,rep(1:nYears, each=nReps)]
    z.a <- array(z.wide, c(nSites, nReps, nYears))

    ll.y.cand[,,dataYears] <- dbinom(y[,,dataYears], 1, z.a[1:nSampled,,dataYears]*p.cand[,,dataYears], log=TRUE)
    prior.beta2.cand <- dnorm(beta2.cand, 0, 10, log=TRUE) 
    prior.beta2 <- dnorm(beta2, 0, 10, log=TRUE)

    ll.y.sum.cand <- sum(ll.y.cand, na.rm=TRUE)
    if(runif(1) < exp((ll.y.sum.cand + prior.beta0.cand) -
                      (ll.y.sum + prior.beta0))) {
        beta2 <- beta2.cand
        p <- p.cand
        ll.y <- ll.y.cand
        ll.y.sum <- ll.y.sum.cand
    }

    zk <- colSums(z)

    samples[s,] <- c(sigma, gamma0.i, gamma0.s, gamma0.p,
                     epsilon.i, epsilon.s, epsilon.p,
                     beta0, beta1, beta2, alpha, zk=zk, deviance=-2*ll.y.sum)
    zK[,s] <- z[,nYears]
    if(monitor.z)
        zA[,,s] <- z
  }

  final.state <- list(z=z, D=D, samples=samples[s,])
  library(coda)
  return(list(samples=samples, final.state=final.state,
              zK=zK, zA=zA, Ez=nz1/nIter,
              seed=.Random.seed))
}

dynroccHC <- cmpfun(dynroccH)
