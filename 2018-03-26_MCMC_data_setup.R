#######################################################################################
# Spatial Occupancy Model for COYOTE. Written by Travis Gallo and Mason Fidino
#  Modified from Paige Howell, Richard Chandler
#######################################################################################

# load internal functions. Assumes this script is in your working directory.
source("GIV_utility_functions.R")

# load patch-level covariate data
covs <- readRDS("./Data/2018-03-27_patch_covariates.RDS")

# Order covariates so the sampled sites are in the same order as detection data
#  frames. Rows 1-98 are sampled sites and data is in the same order as y, z, 
#  and j matrices (even though Station.ID does not match perfectly)
covs <- covs[order(covs$Station.ID),]

# observation and sampling data
full_site_names <- 
  read.table("./Data/sites_used_in_sp10_sp13_analysis_6_1_17.txt", 
             header = TRUE)[,1]

full_site_names[!full_site_names %in% covs$Station.ID] # see the mismatch

# randomly chose 1 site for patches that had multiple sites. Remove unused sites
sites_names <- as.character(full_site_names[-c(20,23,70,94,96)])

# load species specific data - remove sites that we do not use in this analysis
z_mat <- df_2_array(read.table("./Data/z_matrix_sp10_sp13_6_1_17.txt", header = TRUE, sep = "\t"))[1,-c(20,23,70,94,96),-c(1,2)]

# build y-array (detection) and j-matrix (n days sampled)
y_mat <- df_2_array(read.table("./Data/y_matrix_sp10_sp13_6_1_17.txt", 
                    header = TRUE, sep = "\t"))[1,-c(20,23,70,94,96),-c(1,2)]

j_mat <- as.matrix(read.table("./Data/j_matrix_sp10_sp13_6_1_17.txt", 
                   header = TRUE, sep = "\t"))[-c(20,23,70,94,96),-c(1,2)]

# take a look at the raw data
sum(y_mat, na.rm=TRUE) # Total dets
colSums(y_mat, na.rm=TRUE) # Dets per season
sum(rowSums(y_mat, na.rm = TRUE) > 0)/nrow(y_mat) # Proportion of sites occupied

# site coordinates for all sites
coords <- as.matrix(covs[, c("x","y")])

# Site-level covariates
#  site types include: 
#  1 = park 
#  2 = golf course
#  3 = cemetery
#  4 = natural area/conservation
type <- covs[,"site_type"]
park <- golf <- cem <- rep(0, length(type))
park[which(type == 1)] <- 1
golf[which(type == 2)] <- 1
cem[which(type == 3)]  <- 1

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

## resistance covariates
# study extent to crop larger extent rasters
study_extent <- extent(global_patches) + 5000

# NDVI resistence layer: 1-NDVI = resistance (scaled)
ndvi_res <- crop(raster("./Data/NDVI_to_Resistence.tif"), study_extent)
ndvi_coarse <- aggregate(ndvi_res, fact=4, fun=mean) # scale up raster to 120x120
ndvi_scale <- scale(ndvi_coarse)

# 2010 population density raster (scaled)
pop10_res <- crop(raster("./Data/CMAP_PHH10.tif"), study_extent)
pop10_coarse <- aggregate(pop10_res, fact=4, fun=mean)
pop10_scale <- scale(pop10_coarse)

# 2040 population density raster (scaled) - for future projections, not initial model run
pop40_res <- crop(raster("./Data/CMAP_PHH40.tif"), study_extent)
pop40_coarse <- aggregate(pop40_res, fact=4, fun=mean)
pop40_scale <- scale(pop40_coarse)

# interstate raster - ULTIMATELY DID NOT USE INTERSTATES AS A POTENTIAL BARRIER
# interstate_res <- raster("./Data/Interstate_Resistance.tif")
# give resistence a really high value but a non-zero probability of crossing
#to_replace <- values(interstate_coarse) == 1
#values(interstate_coarse)[to_replace] <- max(values(ndvi_scale), na.rm=TRUE) + max(values(pop10_scale), na.rm=TRUE)
# extend raster to match all othe rasters
#interstate_extend <- crop(extend(ndvi_coarse, interstate_coarse), ndvi_coarse)  

# patch indicator indicating that habitat patches have 0 resistance
# 0 resistence when converted to conductance (1/resistence) gives us infinity, which is what we want but not computationally possible
# therefore, we give it basically 0 resistence (0.0001) so that it converts to basically infinity (some really high number) when divided by 1
patch_indicator <- fasterize(st_as_sf(global_patches), ndvi_coarse, field = "AREA", fun="max")
val_to_replace0 <- !is.na(values(patch_indicator))
values(patch_indicator)[val_to_replace0] <- 1
val_to_replace1 <- is.na(values(patch_indicator))
values(patch_indicator)[val_to_replace1] <- 0

# cells from the ndvi and population raster become 0 to estimate transition cost within patches seperately from the matrix
to_0 <- values(patch_indicator) == 1
values(ndvi_scale)[to_0] <- 0
values(pop10_scale)[to_0] <- 0

# create a list of resistence covariates
res_covs <- as.list(c(ndvi=ndvi_scale, pop=pop10_scale, patch=patch_indicator))

###################################################
###################################################

# Test data

# creat small space to test
test_extent <- extent(sites_sampled) + 2000

# crop raster and patches shape for practice
ndvi_crop <- crop(ndvi_scale, test_extent)
patch_crop <- crop(patch_indicator, test_extent)
pop_crop <- crop(pop10_scale, test_extent)
#interstate_crop <- crop(interstate_res_extend, test_extent)
global_patches_crop <- crop(global_patches, test_extent)

# list of resistence covariates
res_covs <- as.list(c(ndvi_crop, pop_crop, patch_crop))

covs2 <- covs[which(covs$patch %in% global_patches_crop@data$Patch),]
covs2 <- covs2[order(covs2$Station.ID),]
covs2 <- covs2[c(1:98, sample(99:12077, 50)),]
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
y <- y.gen
j <- j_mat
site_covs <- sitecovs2
obs_covs <- season_vec
r_covs <- res_covs
disp_dist <- 100000
n.cores <- 4
iters <- 1000
report <- 100
monitor.z <- TRUE
plot.z <- TRUE
# tune order (alpha[1], alpha[2], alpha[3], b0.gam, b.gam[1], b.gam[2], b.gam[3], b.gam[4], sigma, b0.psi1, b.psi1[1], b.psi1[2], b.psi1[3],
# b0.eps, b.eps[1], b.eps[2], b.eps[3], b.eps[4], b.eps[5], b.eps[6], b.eps[7], b.eps[8], a0, season[2], season[3], season[4]
tune <- c(0.3, 0.3, 0.3, # alpha coefficients
          1, 1, 1, 1, 1, # gamma coefficients
          0.3, # sigma
          1, 1, 1, 1, # psi1 coefficients
          1, 1, 1, 1, 1, 1, 1, 1, 1, # epsilon coefficients
          1, 1, 1, 1 #a0 and season parameters (detection)
          )

length(tune)


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
                iters=10,
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