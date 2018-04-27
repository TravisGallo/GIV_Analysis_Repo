#######################################################################################
#Spatial Occupancy Model for COYOTE. Written by Travis Gallo and Mason Fidino
#Modified from Paige Howell, Richard Chandler
#######################################################################################

#load mcmc sampler
#source("mcmcfile.R")
source("GIV_utility_functions.R")

#load patch-level covariate data
covs <- readRDS("./Data/2018-03-27_patch_covariates.RDS")
# order covariates so the sampled sites are in the same order as detection data frames
# rows 1-98 are sampled sites and data is in the same order as y, z, and j matrices (even though Station.ID does not match perfectly)
covs <- covs[order(covs$Station.ID),]

# observation and sampling data
full_site_names <- read.table("./Data/sites_used_in_sp10_sp13_analysis_6_1_17.txt", header = TRUE)[,1]
full_site_names[!full_site_names %in% covs$Station.ID] # see the mismatch
# randomly chose 1 site for patches that had multiple sites. Remove unused sites
sites_names <- as.character(full_site_names[-c(20,23,70,94,96)])

# load species specific data - remove sites that we do not use in this analysis
z_mat <- df_2_array(read.table("./Data/z_matrix_sp10_sp13_6_1_17.txt", header = TRUE, sep = "\t"))[1,-c(20,23,70,94,96),-c(1,2)]
# build y-array and j-matrix
y_mat <- df_2_array(read.table("./Data/y_matrix_sp10_sp13_6_1_17.txt", header = TRUE, sep = "\t"))[1,-c(20,23,70,94,96),-c(1,2)]
j_mat <- read.table("./Data/j_matrix_sp10_sp13_6_1_17.txt", header = TRUE, sep = "\t")[-c(20,23,70,94,96),-c(1,2)]

# take a look at the raw data
sum(y_array, na.rm=TRUE) # Total dets
colSums(y_array, na.rm=TRUE) # Dets per season
sum(rowSums(y_array, na.rm = TRUE) > 0)/nrow(y_array) # Proportion of sites occupied

# site coordinates for all sites
coords <- as.matrix(covs[, c("x","y")])

# site-level covariates
# site type: 1 = park; 2 = golf course; 3 = cemetery' 4 = natural area/conservation
type <- covs[,"site_type"]
park <- golf <- cem <- rep(0, length(type))
park[which(type == 1)] <- 1
golf[which(type == 2)] <- 1
cem[which(type == 3)] <- 1
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
# data matrix
sitecovs <- cbind(tree, total_veg, water, size, pop10, park, golf, cem)

# observation covariate - season
# vector indicating which calander season the observation was obtained
season_vec <- c(4,1,2,3,4,1,2,3,4,1,2)

# resistance covariates
# NDVI resistence layer: 1-NDVI = resistance (scaled)
ndvi_res <- scale(raster("./Data/NDVI_to_Resistence.tif"))
# 2010 population density raster (scaled)
pop10_res <- scale(raster("./Data/CMAP_PHH10.tif"))
# 2040 population density raster (scaled) - for future projections, not initial model run
pop40_res <- scale(raster("./Data/CMAP_PHH40.tif"))
# interstate raster
interstate_res <- raster("./Data/Interstate_Resistance.tif")
# give resistence a really high value but a non-zero probability of crossing
to_replace <- values(interstate_res) == 1
values(interstate_res)[to_replace] <- max(values(ndvi_res), na.rm=TRUE) + max(values(pop10_res), na.rm=TRUE)
# extend raster to match all othe rasters
interstate_res_extend <- crop(extend(ndvi_res, interstate_res), ndvi_res)  
# patch indicator indicating that habitat patches have 0 resistance
patch_indicator <- raster("./Data/2018-03-20_patch_indicator_raster.tif")

# create a list of resistence covariates
res_covs <- as.list(c(ndvi=ndvi_res, pop=pop10_res, interstate=interstate_res, patch=patch_indicator))

###################################################
###################################################

# Test data

# creat small space to test
test_extent <- extent(sites_sampled) + 500
# crop raster and patches shape for practice
ndvi_crop <- crop(ndvi_res, test_extent)
patch_crop <- crop(patch_indicator, test_extent)
pop_crop <- crop(pop10_res, test_extent)
interstate_crop <- crop(interstate_res_extend, test_extent)
global_patches_crop <- crop(global_patches, test_extent)
# list of resistence covariates
res_covs <- as.list(c(ndvi_crop, pop_crop, interstate_crop, patch_crop))

covs2 <- covs[which(covs$patch %in% global_patches_crop@data$Patch),]
covs2 <- covs2[order(covs2$Station.ID),]
covs2 <- covs2[c(1:nsampled, sample(99:12077, 50)),]
# site coordinates for all sites
coords <- as.matrix(covs2[, c("x","y")])
# site-level covariates
# site type: 1 = park; 2 = golf course; 3 = cemetery' 4 = natural area/conservation
type <- covs2[,"site_type"]
park <- golf <- cem <- rep(0, length(type))
park[which(type == 1)] <- 1
golf[which(type == 2)] <- 1
cem[which(type == 3)] <- 1
# tree cover (scaled)
tree <- (covs2[,"patch_tree"] - mean(covs2[,"patch_tree"]))/sd(covs2[,"patch_tree"])
# vegetation cover (scaled)
total_veg <- (covs2[,"patch_total_veg"] - mean(covs2[,"patch_total_veg"]))/sd(covs2[,"patch_total_veg"])
# water present
water <- covs2[,"water_present"]
# patch size (scaled)
size <- (covs2[,"patch_size"] - mean(covs2[,"patch_size"]))/sd(covs2[,"patch_size"])
# population density 2010 (scaled)
pop10 <- (covs2[,"CMAP_pop10_dens"] - mean(covs2[,"CMAP_pop10_dens"]))/sd(covs2[,"CMAP_pop10_dens"])
# population density 2040 (scaled)
pop40 <- (covs2[,"CMAP_pop40_dens"] - mean(covs2[,"CMAP_pop40_dens"]))/sd(covs2[,"CMAP_pop40_dens"])
# data matrix
sitecovs2 <- cbind(tree, total_veg, water, size, pop10, park, golf, cem)

x <- coords
y <- y_mat
j <- j_mat
site_covs <- sitecovs2
obs_covs <- season_vec
r_covs <- res_covs

plot(coords)

###################################################

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