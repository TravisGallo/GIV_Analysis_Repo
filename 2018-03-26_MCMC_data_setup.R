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
# season <- c(2,3,4,1,2,3,4,1,2,3,4,1,2)

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

# observation and sampling data
full_site_names <- read.table("./Data/sites_used_in_sp10_sp13_analysis_6_1_17.txt", header = TRUE)
full_site_names[!full_site_names[,1] %in% covs$Station.ID,] # see the mismatch
# randomly chose 1 site for patches that had multiple sites. Remove unused sites
sites_names <- as.character(full_site_names[-c(23,70,69,95,96),])

# load species specific data - remove sites that we do not use in this analysis
z_array <- df_2_array(read.table("./Data/z_matrix_sp10_sp13_6_1_17.txt", header = TRUE, sep = "\t"))[1,-c(23,70,69,95,96),-c(1,2)]
# build y-array and j-matrix
y_array <- df_2_array(read.table("./Data/y_matrix_sp10_sp13_6_1_17.txt", header = TRUE, sep = "\t"))[1,-c(23,70,69,95,96),-c(1,2)]
j_mat <- read.table("./Data/j_matrix_sp10_sp13_6_1_17.txt", header = TRUE, sep = "\t")[-c(23,70,69,95,96),-c(1,2)]

#Take a look at the raw data
sum(y_array, na.rm=TRUE) # Total dets
colSums(y_array, na.rm=TRUE) # Dets per season
sum(rowSums(y_array, na.rm = TRUE) > 0)/nrow(y_array) # Proportion of sites occupied


###################################################
###################################################



nSampled <- nrow(y.pre)


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