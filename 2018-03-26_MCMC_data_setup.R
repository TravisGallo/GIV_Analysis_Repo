#######################################################################################
#Spatial Occupancy Model for COYOTE. Written by Travis Gallo and Mason Fidino
#Modified from Paige Howell, Richard Chandler
#Updated 3/27/2018
#######################################################################################

#load mcmc sampler
#source("mcmcfile.R")
source("GIV_utility_functions.R")

#load covariate data
covs <- readRDS("./Data/patch_covariates.RDS")

#site coordinates
coords <- as.matrix(covs[, c("x","y")])

##########################################
##########################################

#site hydroperiod
watercov1<-read.csv("water.dryad.csv")
watercov<-as.character(watercov1[,2])

#presence-absence data
y1<-read.csv("y.wide.dryad.csv")
y.pre1<-data.matrix(y1[,2:45], 47)
y.pre<-array(y.pre1, dim=c(47, 3, 15))

#occasion covariates - wind, temp
wind1<-read.csv("wind.wide.dryad.csv")
wind.pre1<-data.matrix(wind1[,2:45], 47)
wind.pre<-array(wind.pre1, dim=c(47, 3, 15))

temp1<-read.csv("temp.wide.dryad.csv")
temp.pre1<-data.matrix(temp1[,2:45], 47)
temp.pre<-array(temp.pre1, dim=c(47, 3, 15))

#elevation
dem.900mag <- raster("dem")
dem.900mag.scale <- scale(dem.900mag) 

#Take a look at the raw data

sum(y.pre, na.rm=TRUE)                            # Total dets
apply(y.pre, 3, sum, na.rm=TRUE)                  # Dets per year
colSums(apply(y.pre, c(1,3), sum, na.rm=TRUE)>0)  # Sites occupied

nSampled <- nrow(y.pre)

devAskNewPage(TRUE)
for(t in 1:dim(y.pre)[3]) {
    dets.t <- which(rowSums(y.pre[,,t], na.rm=TRUE)>0)
    plot(dem.900mag.scale, main=paste("Year", (2003:2017)[t], "with", length(dets.t), "sites occupied"))
    points(coords, pch=3)
    points(coords[1:nSampled,], pch=16)
    points(coords[dets.t,], pch=16, col=rgb(0,0,1,0.9))
}


#Code to run sampler with data for 10 MCMC iterations
## Tune order: sigma, gamma0.i, gamma0.s, gamma0.p, eps.i, eps.s, eps.p, z, beta0, beta1, beta2, alpha1, alpha2
out1se <- dynroccH(y=y.pre, 
                x=coords, 
                r.cov1=dem.900mag.scale,
                e.cov=watercov, 
                p.cov1=wind.pre, 
                p.cov2=temp.pre,
                nIter=10,
                tune=c(0.3, 0.01, 0.1, 0.1, 
                       0.2, 0.2, 0.2, 
                       0.9, 0.9, 0.9, 
                       0.2, 0.2),
                estAlpha=FALSE,
##               zProp="vec", zProbs=Ez,            ## This results in slooow mixing
               report=10, plot.z=TRUE, tol=1e-3)


save(out1se, file="out1se.gzip")
mc1 <- mcmc(out1se$samples)
rejectionRate(window(mc1, start=1))