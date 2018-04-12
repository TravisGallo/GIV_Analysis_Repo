#######################################################################################
#Spatial Occupancy Model for COYOTE. Written by Travis Gallo and Mason Fidino
#Modified from Paige Howell, Richard Chandler
#Updated 3/27/2018
#######################################################################################

#load mcmc sampler
#source("mcmcfile.R")
source("GIV_utility_functions.R")

#load covariate data
covs <- readRDS("./Data/2018-03-27_patch_covariates.RDS")

#site coordinates
coords <- as.matrix(covs[, c("x","y")])

# site-level covariates
# site type: 1 = park; 2 = golf course; 3 = cemetery' 4 = natural area/conservation
type <- covs[,"site_type"]
# tree cover (scaled)
tree <- (covs[,"patch_tree"] - mean(covs[,"patch_tree"]))/sd(covs[,"patch_tree"])
# vegetation cover (scaled)
total_veg <- (covs[,"patch_total_veg"] - mean(covs[,"patch_total_veg"]))/sd(covs[,"patch_total_veg"])
# water present
water <- covs[,"water_present"]
# patch size (scaled)
size <- (covs[,"patch_size"] - mean(covs[,"patch_size"]))/sd(covs[,"patch_size"])
# population density 2010 (scaled)
pop10 <- (covs[,"CMAP_pop10_dens"] - mean(covs[,"CMAP_pop10_dens"]))/sd(covs[,"CMAP_pop10_dens"])
# population density 2040 (scaled)
pop40 <- (covs[,"CMAP_pop40_dens"] - mean(covs[,"CMAP_pop40_dens"]))/sd(covs[,"CMAP_pop40_dens"])

# observation covariate - season
# vector indicating which calander season the observation was obtained
season <- c(3,4,1,2,3,4,1,2,3,4,1)

# resistance covariates
# NDVI resistence layer: 1-NDVI = resistance (scaled)
ndvi_res <- scale(raster("./Data/NDVI_to_Resistence.tif"))
# 2010 population density raster (scaled)
pop10_res <- scale(raster("./Data/CMAP_PHH10.tif"))
# 2040 population density raster (scaled)
pop40_res <- scale(raster("./Data/CMAP_PHH40.tif"))
# interstate raster
interstate_res <- raster("./Data/Interstate_Resistance.tif")
# patch indicator indicating that habitat patches have 0 resistance
patch_indicator <- raster("./Data/2018-03-20_patch_indicator_raster.tif")

###################################################
###################################################

#presence-absence data
y1<-read.csv("y.wide.dryad.csv")
y.pre1<-data.matrix(y1[,2:45], 47)
y.pre<-array(y.pre1, dim=c(47, 3, 15))

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