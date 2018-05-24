### Generate Some Data for Simulations ###

y_gen <- matrix(0, 98, 11)
## Dimensions
nsite <- nrow(x) # Number of possible sites instead of only the sites sampled
nseason <- ncol(y_gen)
nsampled <- nrow(y_gen)

## hard code dummy parameters
# resistence
alpha1 <- c(0.95, 0.95, 0.001)
# initial occupancy
b0.psi1 <- 0.4
b.psi1 <- c(-0.75, 0.25, 0.25)
# colonization
b0.gam <- 0.25
b.gam <- c(0.75, -0.5, 0.25, 0.25)
gamma0 <- plogis(b0.gam + b.gam[1]*site_covs[,"size"] + b.gam[2]*site_covs[,"park"] + b.gam[3]*site_covs[,"cem"] + b.gam[4]*site_covs[,"golf"])
sigma <- 2
# extinction
b0.eps <- -0.75
b.eps <- c(-1, -0.2, -1.5, 1.25, -0.5, 0.25, -0.25, -0.25)
epsilon <- plogis(b0.eps + b.eps[1]*site_covs[,"tree"] + b.eps[2]*site_covs[,"total_veg"] + b.eps[3]*site_covs[,"size"] + b.eps[4]*site_covs[,"pop10"] +
                    b.eps[5]*site_covs[,"water"] + b.eps[6]*site_covs[,"park"] + b.eps[7]*site_covs[,"cem"] + b.eps[8]*site_covs[,"golf"])
# detection - made to be low
a0 <- -1.25
season <- c(0.025, -0.025, 0.025, 0.025)
season[1] <- 0
p <- plogis(a0 + season[season_vec])
# initial occupancy
psi1 <- plogis(b0.psi1 + b.psi1[1]*site_covs[,"park"] + b.psi1[2]*site_covs[,"cem"] + b.psi1[3]*site_covs[,"golf"])
psi <- matrix(NA, nsite, nseason)
psi[,1] <- psi1
# create resistance surface
cost <- exp(alpha[1]*r_covs[[1]] + alpha[2]*r_covs[[2]] + alpha[3]*r_covs[[3]])
# calculate conductances among neighbors
# create transition matrix - here we convert our cost to conductance by doing 1/resistence
tr1 <- transition(cost, transitionFunction=function(x) 1/mean(x), directions=16) 
# adjust diag. conductances
tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE, scl=FALSE) 
# calculate the ecological distance matrix in parallel
D <- costDistance_mod(tr1CorrC, fromCoords=x, toCoords=x, dist.cutoff=disp_dist, n.cores)/1000
G <- gamma0*exp(-D^2/(2*sigma^2))
G[is.na(G)] <- 0
# incorporate spatially-explicit gamma into occupancy model
z <- matrix(NA, nsite, nseason)
gamma <- matrix(NA, nsite, nseason-1)
z[,1] <- rbinom(nsite, 1, psi[,1])
y_gen[,1] <- rbinom(nsampled, j[,1], z[1:nsampled,1]*p[1])
for(k in 2:nseason) {
  zkt <- matrix(z[,k-1], nsite, nsite, byrow=TRUE)
  PrNotColonizedByNeighbor <- 1 - G*zkt
  PrNotColonizedAtAll <- apply(PrNotColonizedByNeighbor, 1, prod, na.rm=TRUE)
  gamma[,k-1] <- 1 - PrNotColonizedAtAll
  psi[,k] <- z[,k-1]*(1-epsilon*(1-gamma[,k-1])) + (1-z[,k-1])*gamma[,k-1] #Rescue effect
  z[,k] <- rbinom(nsite, 1, psi[,k])
  # observation model
  y_gen[,k] <- rbinom(nsampled, j[,k], z[1:nsampled,k]*p[k])
}