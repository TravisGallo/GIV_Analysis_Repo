library(compiler)

## MCMC algorithm
dynroccH <- function(y,            # nSampled x nSeason matrix of detection data*
                    x,             # nSites x 2 matrix of site coordinates. Note that nSampled will usually be <nSites*
                    patch_ind,     # indicate cells that are considered habitat patches
                    disp_cutoff,   # dispersal distance cutoff
                    r_cov1,        # resistance covariate
                    r_cov2,        # resistance covariate
                    r_cov3,        # resistance covariate
                    e_cov1,        # extinction covariate
                    e_cov2,        # extinction covariate
                    e_cov3,        # extinction covariate
                    e_cov4,        # extinction covariate
                    e_cov5,        # extinction covariate
                    e_cov6,        # extinction covariate
                    p_cov1,        # detection covariate
                    nIter=10,      # MCMC iterations
                    tune,          # Tuning order: sigma,gamma0.i,gamma0.s,gamma0.p,eps.i,eps.s,eps.p,
                                   #               beta0,beta1,beta2,alpha[1],alpha[2] (12 in total)
                    estAlpha=TRUE, # Estimate the resistance coefficient? ###IS THIS NEEDED####
                    inits=NULL,    # until you run algorithm, inits are based on what is given.
                    zProp=c("ind","vec"), # Update z matrix by either proposing z(i,k) or z(,k), respectively
                    zProbs=NULL,   # matrix of proposal probs use if zProp="vec"
                    monitor.z=FALSE, # store each iteration of the z matrix?
                    report=0,      # Only report progress if >0
                    plot.z=FALSE,  # Plot the latent presence-absence state (if report>0)
                    tol=0)      # This will reject a proposal of z(i,k)=1 if mu(i,k-1)<tol
{

  zProp <- match.arg(zProp)

  ## Dimensions*
  nsite <- nrow(x) # Number of possible sites instead of only the sites sampled
  nseason <- ncol(y)
  
  ## Using this to avoid likelihood calculations for sites not sampled
  nsampled <- nrow(y)
  dataYears <- apply(!is.na(y), 3, any) #*

  anyDetections <- matrix(FALSE, nsite, nseason)
  anyDetections[1:nsampled,] <- as.numeric(y > 0)

  #y.wide <- matrix(y, nSampled)

  #isInter <- e.cov=="Intermittent"
  #isSemi <- e.cov=="Semi-permanent"
  #isPerm <- e.cov=="Permanent"
  epsilon <- rep(NA, nsite)  
  gamma0 <- rep(NA, nsite) 

  ## initial values
  gamma <- muz <- matrix(NA, nSites, nYears-1)
  if(is.null(inits)) {
      epsilon.p<-epsilon.s<-epsilon.i<-runif(1) 
      sigma <- runif(1,3,4)
      gamma0.i <- runif(1, 0.01, 0.3)
      gamma0.s <- runif(1, 0.01, 0.3)
      gamma0.p <- runif(1, 0.01, 0.3)
      gamma0[isInter] <- gamma0.i
      gamma0[isSemi] <- gamma0.s
      gamma0[isPerm] <- gamma0.p
      beta0<-runif(1, 0.1, 0.6)
      epsilon[isInter] <- epsilon.i
      epsilon[isSemi] <- epsilon.s
      epsilon[isPerm] <- epsilon.p
      beta0 <- rnorm(1)
      beta1 <- rnorm(1)
      beta2 <- rnorm(1)
      alpha <- c(0, 0) 
      p <- plogis(beta0 + beta1*p.cov1 + beta2*p.cov2) 
      z <- matrix(0, nSites, nYears)
      
      #Which sites were the reintroduction sites
      z[c(15, 33, 274),1] <- 1 

      ## NOTE: For some organisms and systems the maximum dispersal distance within a single time step 
        #(e.g., annual) may be known and it could be considered very unlikely that colonists would arrive from patches exceeding this distance. 
        #To speed up computation, a neighborhood matrix for each patch could be supplied to eliminate consideration of colonization from neighboring patches 
        #that exceed a reasonable distance from the focal patch.     
      
      ## create resistance surface
      cost <- exp(alpha[1]*r_cov1 + alpha[2]*r_cov2 + alpha[3]*r_cov3)*patch_ind
      ## calculate conductances among neighbors
      tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16) 
      #adjust diag.conductances
      tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE,scl=FALSE) 
      ## calculate least cost distance between all pairs of sites.
      if(!estAlpha)
          alpha <- c(0,0) ## Force alpha to be 0 if you aren't estimating it. Results in appox Euclidean dist
      #calculate the ecological distance matrix - assuming going back and forth is equal*
      D <- costDistance(tr1CorrC,x)/1000 
      G <- gamma0*exp(-D^2/(2*sigma^2))

      for(k in 2:nSeasons) { 
          PrNotColonizedByNeighbor <- 1 - gamma0*exp(-D^2/(2*sigma^2)) * t(z[,rep(k-1, nSites)])
          PrNotColonizedAtAll <- apply(PrNotColonizedByNeighbor, 1, prod)
          gamma[,k-1] <- 1 - PrNotColonizedAtAll
          
          psi[,k-1] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) + (1-z[,k-1])*gamma[,k-1] #Rescue effect
          
          z[,k] <- rbinom(nSites, 1, psi[,1]) # why hard coded 1 and not k
          z[known0[,k],k] <- 0 
          z[which(anyDetections[,k]),k] <- 1 
      }
  } else {
      gamma0.i <- inits$samples["gamma0.i"]
      gamma0.s <- inits$samples["gamma0.s"]
      gamma0.p <- inits$samples["gamma0.p"]
      gamma0[isInter] <- gamma0.i
      gamma0[isSemi] <- gamma0.s
      gamma0[isPerm] <- gamma0.p
      sigma <- inits$samples["sigma"]
      epsilon.i <- inits$samples["epsilon.i"]
      epsilon.s <- inits$samples["epsilon.s"]
      epsilon.p <- inits$samples["epsilon.p"]
      epsilon <- rep(NA, nSites)
      epsilon[isInter] <- epsilon.i
      epsilon[isSemi] <- epsilon.s
      epsilon[isPerm] <- epsilon.p
      alpha<-c(inits$samples["alpha1"],inits$samples["alpha2"])
      D <- inits$D
      beta0 <- inits$samples["beta0"]
      beta1 <- inits$samples["beta1"]
      beta2 <- inits$samples["beta2"]
      p <- plogis(beta0 + beta1*p.cov1 + beta2*p.cov2)
      z <- inits$z
      .Random.seed <- inits$seed ## use same random seed as before
  }

  ll.z <- matrix(0, nSites, nYears)
  ll.y <- array(0, c(nSampled, nReps, nYears))
  for(k in 2:nseason) {
      PrNotColonizedByNeighbor <- 1 - gamma0*exp(-D^2/(2*sigma^2  ))*t(z[,rep(k-1, nSites)])
      PrNotColonizedAtAll <- apply(PrNotColonizedByNeighbor, 1, prod)
      gamma[,k-1] <- 1 - PrNotColonizedAtAll
      psi[,k-1] <- (1-z[,k-1])*gamma[,k-1] + z[,k-1]*(1-epsilon*(1-gamma[,k-1])) #PH rescue effect
      ll.z[,k-1] <- dbinom(z[,k], 1, psi[,k-1], log=TRUE)
      # observation model
      ll.y[,k] <- dbinom(y[,k], 1, z[1:nsampled,k]*p[,k], log=TRUE)
      
  }
  ll.z.cand <- ll.z
  ll.z.sum <- sum(ll.z)
  ll.y.cand <- ll.y
  ll.y.sum <- sum(ll.y, na.rm=TRUE)
  gamma.cand <- gamma
  psi.cand <- psi

    nz1 <- z ## Used to compute expected occupancy at each site

    zkup <- rep(0, nseason-1)

  ## posterior samples
  nPar <- 13+nseason # 13 is the number of parameters including intercept to eastimate for each season
  samples <- array(NA, c(nIter, nPar))
  zK <- matrix(NA, nsite, nIter)
  colnames(samples) <- c("sigma", "gamma0.i", "gamma0.s", "gamma0.p",
                         "epsilon.i", "epsilon.s","epsilon.p",
                         "beta0", "beta1", "beta2", "alpha1", "alpha2",
                         paste("zk", 1:nseason, sep=""), "deviance")

  reportit <- report>0
    nzup <- rep(0, nYears-1)
  zA <- NULL
  if(monitor.z){
      zA <- array(NA_integer_, c(nsite, nseason, nIter))
  }

  if(reportit) {
      cat("iter 1\n")
      cat("    theta =", round(c(sigma,gamma0.i,gamma0.s,gamma0.p,epsilon.i,epsilon.s,epsilon.p,beta0,beta1,beta2,alpha), 5), "\n")
      cat("    z[k] =", round(colSums(z), 2), "\n")
      cat("    ll.z =", round(sum(ll.z), 2), "\n")
      cat("    deviance =", round(-2*ll.y.sum, 2), "\n")
      cat("    time =", format(Sys.time()), "\n")
      if(plot.z) {
          library(lattice)
          zd <- data.frame(z=as.integer(z), year=factor(rep(2003:2017, each=nSites)),
                           x=as.numeric(x[,1])/1000, y=as.numeric(x[,2])/1000)
          print(xyplot(y ~ x | year, zd, groups=z, aspect="iso", pch=c(1,16), as.table=TRUE))
      }
  }

  ## Sample from posteriors
  for(s in 1:nIter) {

    ll.z.sum <- sum(ll.z) ## This is important!

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

    if(estAlpha) {
      
    library(gdistance) 
    #Metropolis update for alpha
    alpha1.cand <- rnorm(1, alpha[1], tune[11])
    #create resistance surface
    cost <- exp(alpha1.cand*r.cov1 + alpha[2]*r_cov2 + alpha[3]*r_cov3)*patch_ind
    ## calculate conductances among neighbors
    tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16)
    tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE,scl=FALSE) #adjust diag.conductances

    ## calculate least cost distance between all pairs of sites.
    D.cand <- costDistance(tr1CorrC,x,x)/1000 #calculate the ecological distance matrix
    G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2))

    for(k in 2:nYears) {
      zkt <- matrix(z[,k-1], nSites, nSites, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    prior.alpha.cand <- dnorm(alpha1.cand, 0, 10, log=TRUE)
    prior.alpha <- dnorm(alpha[1], 0, 10, log=TRUE)

    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.alpha.cand) -
                        (ll.z.sum + prior.alpha))) {
      alpha[1] <- alpha1.cand
      D <- D.cand
      G <- G.cand
      gamma <- gamma.cand
      psi <- psi.cand
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      }

    #Metropolis update for alpha[2]
    alpha2.cand <- rnorm(1, alpha[2], tune[12])
    #create resistance surface
    cost <- exp(alpha[1]*r.cov1 + alpha2.cand*r.cov2 + alpha[3]*r_cov3)*patch_ind
    ## calculate conductances among neighbors
    tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16) 
    tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE,scl=FALSE) #adjust diag.conductances
    ## calculate least cost distance between all pairs of sites.
    D.cand <- costDistance(tr1CorrC,x,x)/1000 #calculate the ecological distance matrix
    G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2  ))

    for(k in 2:nYears) {
      zkt <- matrix(z[,k-1], nSites, nSites, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    prior.alpha.cand <- dnorm(alpha2.cand, 0, 10, log=TRUE)
    prior.alpha <- dnorm(alpha[2], 0, 10, log=TRUE)

    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.alpha.cand) -
                        (ll.z.sum + prior.alpha))) {
      alpha[2] <- alpha2.cand
      D <- D.cand
      G <- G.cand
      gamma <- gamma.cand
      muz <- muz.cand
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
    }
    
    #Metropolis update for alpha[3]
    alpha3.cand <- rnorm(1, alpha[2], tune[12])
    #create resistance surface
    cost <- exp(alpha[1]*r.cov1 + alpha[2]*r.cov2 + alpha3.cand*r_cov3)*patch_ind
    ## calculate conductances among neighbors
    tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16) 
    tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE,scl=FALSE) #adjust diag.conductances
    ## calculate least cost distance between all pairs of sites.
    D.cand <- costDistance(tr1CorrC,x,x)/1000 #calculate the ecological distance matrix
    G.cand <- gamma0*exp(-D.cand^2/(2*sigma^2  ))
    
    for(k in 2:nYears) {
      zkt <- matrix(z[,k-1], nSites, nSites, byrow=TRUE)
      gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
      psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
      ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
    }
    prior.alpha.cand <- dnorm(alpha3.cand, 0, 10, log=TRUE)
    prior.alpha <- dnorm(alpha[3], 0, 10, log=TRUE)
    
    ll.z.sum.cand <- sum(ll.z.cand)
    if(runif(1) < exp((ll.z.sum.cand + prior.alpha.cand) -
                      (ll.z.sum + prior.alpha))) {
      alpha[3] <- alpha3.cand
      D <- D.cand
      G <- G.cand
      gamma <- gamma.cand
      muz <- muz.cand
      ll.z <- ll.z.cand
      ll.z.sum <- ll.z.sum.cand
      }
    }

    ## Metropolis update for sigma
    sigma.cand <- rnorm(1, sigma, tune[1])
    if(sigma.cand > 0) {
        G.cand <- gamma0*exp(-D^2/(2*sigma.cand^2))
      for(k in 2:nYears) {
        zkt <- matrix(z[,k-1], nSites, nSites, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        psi.cand[,k-1] <- (z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])
        ll.z.cand[,k-1] <- dbinom(z[,k], 1, psi.cand[,k-1], log=TRUE)
      }
      prior.sigma.cand <- dgamma(sigma.cand, 0.001, 0.001)
      prior.sigma <- dgamma(sigma, 0.001, 0.001)
      ll.z.sum.cand <- sum(ll.z.cand)
      if(runif(1) < exp((ll.z.sum.cand + prior.sigma.cand) -
                        (ll.z.sum + prior.sigma))) {
          sigma <- sigma.cand
          gamma <- gamma.cand 
          ll.z <- ll.z.cand
          ll.z.sum <- ll.z.sum.cand
          muz <- muz.cand
          G <- G.cand
      }
    }

    #Metropolis update for gamma0 at intermittent sites (part of the gammaDist calculation)
    prior.gamma0.cand <- prior.gamma0 <- 0
    gamma0.i.cand <- rnorm(1, gamma0.i, tune[2])
    if(gamma0.i.cand > 0 & gamma0.i.cand < 1) {
        gamma0.cand <- gamma0
        gamma0.cand[isInter] <- gamma0.i.cand
        G.cand <- gamma0.cand*exp(-D^2/(2*sigma^2  ))
      for(k in 2:nYears) { #nYears
        zkt <- matrix(z[,k-1], nSites, nSites, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        muz.cand[,k-1] <-(z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])*notFailed[,k]
        ll.z.cand[,k-1] <- dbinom(z[,k], 1, muz.cand[,k-1], log=TRUE)
    }
      ll.z.sum.cand <- sum(ll.z.cand)
      if(runif(1) < exp((ll.z.sum.cand + prior.gamma0.cand) -
                        (ll.z.sum + prior.gamma0))) {
        gamma0.i <- gamma0.i.cand
        gamma0 <- gamma0.cand
        gamma <- gamma.cand
        muz <- muz.cand
        ll.z <- ll.z.cand
        ll.z.sum <- ll.z.sum.cand
        G <- G.cand
    }
  }

    #Metropolis update for gamma0 at semi-permanent sites (part of the gammaDist calculation)
    gamma0.s.cand <- rnorm(1, gamma0.s, tune[3])
    if(gamma0.s.cand > 0 & gamma0.s.cand < 1) {
        gamma0.cand <- gamma0
        gamma0.cand[isSemi] <- gamma0.s.cand
        G.cand <- gamma0.cand*exp(-D^2/(2*sigma^2  ))
      for(k in 2:nYears) { #nYears
        zkt <- matrix(z[,k-1], nSites, nSites, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        muz.cand[,k-1] <-(z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])*notFailed[,k]
        ll.z.cand[,k-1] <- dbinom(z[,k], 1, muz.cand[,k-1], log=TRUE)
    }
      ll.z.sum.cand <- sum(ll.z.cand)
      if(runif(1) < exp((ll.z.sum.cand + prior.gamma0.cand) -
                        (ll.z.sum + prior.gamma0))) {
        gamma0.s <- gamma0.s.cand
        gamma0 <- gamma0.cand
        gamma <- gamma.cand
        muz <- muz.cand
        ll.z <- ll.z.cand
        ll.z.sum <- ll.z.sum.cand
        G <- G.cand
    }
  }


    #Metropolis update for gamma0 at permanent sites (part of the gammaDist calculation)
    gamma0.p.cand <- rnorm(1, gamma0.p, tune[4])
    if(gamma0.p.cand > 0 & gamma0.p.cand < 1) {
        gamma0.cand <- gamma0
        gamma0.cand[isPerm] <- gamma0.p.cand
        G.cand <- gamma0.cand*exp(-D^2/(2*sigma^2  ))
      for(k in 2:nYears) { #nYears
        zkt <- matrix(z[,k-1], nSites, nSites, byrow=TRUE)
        gamma.cand[,k-1] <- 1 - exp(rowSums(log(1-G.cand*zkt)))
        muz.cand[,k-1] <-(z[,k-1]*(1-epsilon*(1-gamma.cand[,k-1])) + (1-z[,k-1])*gamma.cand[,k-1])*notFailed[,k]
        ll.z.cand[,k-1] <- dbinom(z[,k], 1, muz.cand[,k-1], log=TRUE)
    }
      ll.z.sum.cand <- sum(ll.z.cand)
      if(runif(1) < exp((ll.z.sum.cand + prior.gamma0.cand) -
                        (ll.z.sum + prior.gamma0))) {
        gamma0.p <- gamma0.p.cand
        gamma0 <- gamma0.cand
        gamma <- gamma.cand
        muz <- muz.cand
        ll.z <- ll.z.cand
        ll.z.sum <- ll.z.sum.cand
        G <- G.cand
    }
  }


    ## Metropolis update for epsilon (intermittent sites)
    epsilon.i.cand <- rnorm(1, epsilon.i, tune[5])
    if(epsilon.i.cand > 0 & epsilon.i.cand < 1) {
        for(k in 2:nYears) {
            muz.cand[isInter,k-1] <- z[isInter,k-1]*(1-epsilon.i.cand*(1-gamma[isInter,k-1])) + (1-z[isInter,k-1])*gamma[isInter,k-1] 
            muz.cand[isInter,k-1] <- muz.cand[isInter,k-1]*notFailed[isInter,k]
            ll.z.cand[isInter,k-1] <- dbinom(z[isInter,k], 1, muz.cand[isInter,k-1], log=TRUE)
        }
        prior.epsilon.i.cand <- dbeta(epsilon.i.cand, 1, 1, log=TRUE)
        prior.epsilon.i <- dbeta(epsilon.i, 1, 1, log=TRUE)
        if(runif(1) < exp((sum(ll.z.cand[isInter,]) + prior.epsilon.i.cand) -
                          (sum(ll.z[isInter,]) + prior.epsilon.i))) {
            epsilon.i <- epsilon.i.cand
            epsilon[isInter] <- epsilon.i.cand
            muz[isInter,] <- muz.cand[isInter,]
            ll.z[isInter,] <- ll.z.cand[isInter,]
        }
    }


    ## Metropolis update for epsilon (semi-permanent sites)
    epsilon.s.cand <- rnorm(1, epsilon.s, tune[6])
    if(epsilon.s.cand > 0 & epsilon.s.cand < 1) {
        for(k in 2:nYears) {
            muz.cand[isSemi,k-1] <- z[isSemi,k-1]*(1-epsilon.s.cand*(1-gamma[isSemi,k-1])) + (1-z[isSemi,k-1])*gamma[isSemi,k-1] 
            muz.cand[isSemi,k-1] <- muz.cand[isSemi,k-1]*notFailed[isSemi,k]
            ll.z.cand[isSemi,k-1] <- dbinom(z[isSemi,k], 1, muz.cand[isSemi,k-1], log=TRUE)
        }
        prior.epsilon.s.cand <- dbeta(epsilon.s.cand, 1, 1, log=TRUE)
        prior.epsilon.s <- dbeta(epsilon.s, 1, 1, log=TRUE)
        if(runif(1) < exp((sum(ll.z.cand[isSemi,]) + prior.epsilon.s.cand) -
                          (sum(ll.z[isSemi,]) + prior.epsilon.s))) {
            epsilon.s <- epsilon.s.cand
            epsilon[isSemi] <- epsilon.s.cand
            muz[isSemi,] <- muz.cand[isSemi,]
            ll.z[isSemi,] <- ll.z.cand[isSemi,]
        }
    }


    ## Metropolis update for epsilon (permanent sites)
    epsilon.p.cand <- rnorm(1, epsilon.p, tune[7])
    if(epsilon.p.cand > 0 & epsilon.p.cand < 1) {
        for(k in 2:nYears) {
            muz.cand[isPerm,k-1] <- z[isPerm,k-1]*(1-epsilon.p.cand*(1-gamma[isPerm,k-1])) + (1-z[isPerm,k-1])*gamma[isPerm,k-1] 
            muz.cand[isPerm,k-1] <- muz.cand[isPerm,k-1]*notFailed[isPerm,k]
            ll.z.cand[isPerm,k-1] <- dbinom(z[isPerm,k], 1, muz.cand[isPerm,k-1], log=TRUE)
        }
        prior.epsilon.p.cand <- dbeta(epsilon.p.cand, 1, 1, log=TRUE)
        prior.epsilon.p <- dbeta(epsilon.p, 1, 1, log=TRUE)
        if(runif(1) < exp((sum(ll.z.cand[isPerm,]) + prior.epsilon.p.cand) -
                          (sum(ll.z[isPerm,]) + prior.epsilon.p))) {
            epsilon.p <- epsilon.p.cand
            epsilon[isPerm] <- epsilon.p.cand
            muz[isPerm,] <- muz.cand[isPerm,]
            ll.z[isPerm,] <- ll.z.cand[isPerm,]
        }
    }



    ## update z
    ## We can update each z(i,t) individually, and it results in better mixing than updating a vector of z's
    zkup <- rep(0, nYears-1)
   for(k in 2:nYears) {
       anyDet <- anyDetections[,k]==1 
       zknown <- anyDet | !notFailed[,k] 

       prop.back <- prop.cand <- 0
       for(i in 1:nSites) {
           if(zknown[i])
               next
           ## Reject highly unlikely proposals (before proposing them)
           ## This speed trick shouldn't affect anything but
           ## can double check by changing toleranc (tol)
           if(z[i,k]<1 & muz[i,k-1]<tol)
               next
           zk.wide <- matrix(z[,k], nSites, nReps)
           zk.cand <- z[,k]
           zk.cand[i] <- 1-z[i,k]
           zk.cand.wide <- matrix(zk.cand, nSites, nReps)

           ll.y.tmp <- 0
           ll.y.cand.tmp <- 0
           if((k > 4) & (i <= nSampled)) { ## Ignore first 4 years without data
               ll.y.cand.tmp <- dbinom(y[i,,k], 1, zk.cand[i]*p[i,,k], log=TRUE)
               ll.y.tmp <- sum(ll.y[i,,k], na.rm=TRUE)
           }
           ## RC: Prior must be calculated for time k and k+1 b/c change in z affects both
           ll.z.cand[i,k-1] <- dbinom(zk.cand[i], 1, muz[i,k-1], log=TRUE)
           ll.z2 <- ll.z2.cand <- 0
           if(k < nYears) {
               zkt.cand <- matrix(zk.cand, nSites, nSites, byrow=TRUE)
               gamma.cand[,k] <- 1 - exp(rowSums(log(1-G*zkt.cand)))
               muz.cand[,k] <- (zk.cand*(1-epsilon*(1-gamma.cand[,k])) + (1-zk.cand)*gamma.cand[,k])*notFailed[,k+1]
               ll.z.cand[,k] <- dbinom(z[,k+1], 1, muz.cand[,k], log=TRUE)
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
                   muz[,k] <- muz.cand[,k]
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


#This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.