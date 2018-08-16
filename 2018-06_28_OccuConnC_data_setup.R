#######################################################################################
# Data set up and function run for Spatial Occupancy Model for COYOTE. 
# Written by Travis Gallo and Mason Fidino
# Modified from Paige Howell, Richard Chandler
#######################################################################################

# load internal functions. Assumes this script is in your working directory.
source("GIV_utility_functions.R")

# observation and sampling data
full_site_names <- 
  read.table("./Data/sites_used_in_sp10_sp13_analysis_6_1_17.txt", 
             header = TRUE)[,1]

# randomly chose 1 site for patches that had multiple sites. Remove unused sites
sites_names <- as.character(full_site_names[-c(20,23,70,94,96)])

# load species specific data - remove sites that we do not use in this analysis
z_mat <- df_2_array(read.table("./Data/z_matrix_sp10_sp13_6_1_17.txt", header = TRUE, 
                               sep = "\t"))[1,-c(20,23,70,94,96),-c(1,2)]

# build y-array (detection) and j-matrix (n days sampled)
y_mat <- df_2_array(read.table("./Data/y_matrix_sp10_sp13_6_1_17.txt", 
                    header = TRUE, sep = "\t"))[1,-c(20,23,70,94,96),-c(1,2)]

j_mat <- as.matrix(read.table("./Data/j_matrix_sp10_sp13_6_1_17.txt", 
                   header = TRUE, sep = "\t"))[-c(20,23,70,94,96),-c(1,2)]

# take a look at the raw data
sum(y_mat, na.rm=TRUE) # Total dets
colSums(y_mat, na.rm=TRUE) # Dets per season
sum(rowSums(y_mat, na.rm = TRUE) > 0)/nrow(y_mat) # Proportion of sites occupied

# build a study extent that is a square around each transect 
  # with a diamter as long as the transect + 10km
# we then dissolve the three squares to be one full study area that we will work within
# load end points of each transect
end_points <-read.csv("./Data/transect_endpoints.csv")
# transect names
transect_names <- as.character(unique(end_points$Transect))
# create transects
transects <- vector("list", length(transect_names))
square <- s <- vector("list", length(transect_names))
theta <- dx <- dy <- rep(0, length(transect_names))
x0_1 <- y0_1 <- x1_1 <- y1_1 <- rep(0, length(transect_names))
x0_2<- y0_2 <- x1_2 <- y1_2 <- rep(0, length(transect_names))
# loop through and make a list of transects for each city
for (i in 1:length(transect_names)) {
  transects[[i]] <- Lines(list(Line(rbind(as.matrix(end_points[which(end_points$Transect == transect_names[i]),2:3]),
                                          as.matrix(end_points[which(end_points$Transect == transect_names[i]),4:5])))),
                          as.character(transect_names[i]))
  # do some trig
  theta[i] <- atan2(transects[[i]]@Lines[[1]]@coords[1,2] - transects[[i]]@Lines[[1]]@coords[2,2],
                    transects[[i]]@Lines[[1]]@coords[1,1] - transects[[i]]@Lines[[1]]@coords[2,1])
  thetaT <- theta+pi/2
  # calculate the distance for both the x and y points
  dx <- ((LineLength(transects[[i]]@Lines[[1]]@coords, longlat = FALSE, sum = TRUE) 
          + 10000)/2) * cos(thetaT)
  dy <- ((LineLength(transects[[i]]@Lines[[1]]@coords, longlat = FALSE, sum = TRUE) 
          + 10000)/2) * sin(thetaT)
  # add those distances to the transect
  # end 1
  x0_1[i] <- transects[[i]]@Lines[[1]]@coords[1,1] + dx[i]
  y0_1[i] <- transects[[i]]@Lines[[1]]@coords[1,2] + dy[i]
  x1_1[i] <- transects[[i]]@Lines[[1]]@coords[1,1] - dx[i]
  y1_1[i] <- transects[[i]]@Lines[[1]]@coords[1,2] - dy[i]
  # end 2
  x0_2[i] <- transects[[i]]@Lines[[1]]@coords[2,1] + dx[i]
  y0_2[i] <- transects[[i]]@Lines[[1]]@coords[2,2] + dy[i]
  x1_2[i] <- transects[[i]]@Lines[[1]]@coords[2,1] - dx[i]
  y1_2[i] <- transects[[i]]@Lines[[1]]@coords[2,2] - dy[i]
  
  # calculate polygon coordinates
  square[[i]] <- matrix(c(x0_1[i], y0_1[i], # corner
                          x0_2[i], y0_2[i],      # corner
                          x1_2[i], y1_2[i],      # corner
                          x1_1[i], y1_1[i],      # corner
                          x0_1[i], y0_1[i]),     # first corner again - closes the polygon
                        ncol = 2, byrow = TRUE)
  # make them a list of polygons
  s[[i]]<-Polygons(list(Polygon(square[[i]])), i)
  
}

# make the polygons spatial
square_poly <- SpatialPolygons(s, proj4string = CRS("+init=epsg:26916"))
# make transects spatial just for plot checks
transects <- SpatialLines(transects, proj4string = CRS("+init=epsg:26916"))

# dissolve the squares
ID <- rep(1,length(transect_names))
study_extent <- unionSpatialPolygons(square_poly, ID)
# we need to add a bit of distance to the outside of the study extent poly
  # we will crop the rasters at the larger extent and later the patches at the smaller extent
  # this will insure that no coordinates land outside the resistence rasters
  # when the above happens we get errors in the shortest.path function
study_extent_ballooned <- buffer(study_extent, width=1000)

## resistance covariates
# NDVI resistence layer: 1-NDVI = resistance (scaled)
ndvi_coarse <- aggregate(raster("./Data/NDVI_to_Resistence.tif"), fact=10, fun=mean) # scale up raster to 300x300
ndvi_study_extent <- mask(crop(ndvi_coarse, study_extent_ballooned), study_extent_ballooned)
ndvi_scale <- scale(ndvi_study_extent)

# 2010 population density raster (scaled)
pop10_coarse <- aggregate(raster("./Data/CMAP_PHH10.tif"), fact=10, fun=mean)
pop10_study_extent <- mask(crop(pop10_coarse, study_extent_ballooned), study_extent_ballooned)
# add 0.001 to all population values to remove 0's
values(pop10_study_extent) <- values(pop10_study_extent) + 0.001
# log transform to remove large right tale and scale mean 0
pop10_scale <- scale(log(pop10_study_extent))

# 2040 population density raster (scaled) - for future projections, not initial model run
pop40_coarse <- aggregate(raster("./Data/CMAP_PHH40.tif"), fact=10, fun=mean)
pop40_study_extent <- mask(crop(pop40_coarse, study_extent_ballooned), study_extent_ballooned)
# add 0.001 to all population values to remove 0's
values(pop40_study_extent) <- values(pop40_study_extent) + 0.001
# log transform to remove large right tale and scale mean 0
pop40_scale <- scale(log(pop40_study_extent))

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
load("2018-04-12_GIV_Workspace.RData")
patch_indicator <- fasterize(st_as_sf(global_patches), ndvi_coarse, field = "AREA", fun="max")
patch_indicator <- mask(crop(patch_indicator, study_extent_ballooned), study_extent_ballooned)
val_to_replace0 <- !is.na(values(patch_indicator))
values(patch_indicator)[val_to_replace0] <- 1
val_to_replace1 <- is.na(values(patch_indicator))
values(patch_indicator)[val_to_replace1] <- 0

# cells from the ndvi and population raster become 0 to estimate transition cost within patches seperately from the matrix
to_0 <- values(patch_indicator) == 1
values(ndvi_scale)[to_0] <- 0
values(pop10_scale)[to_0] <- 0

## create a list of resistence covariates
res_covs <- as.list(c(ndvi=ndvi_scale, pop=pop10_scale, patch=patch_indicator))

# crop all the global patches to the smaller study extent
global_patches_study_extent <- crop(global_patches, study_extent, study_extent)

# load patch-level covariate data
covs <- readRDS("./Data/2018-03-27_patch_covariates.RDS")

# Order covariates so the sampled sites are in the same order as detection data
#  frames. Rows 1-98 are sampled sites and data is in the same order as y, z, 
#  and j matrices (even though Station.ID does not match perfectly)
covs <- covs[order(covs$Station.ID),]

# grab covariates from global data frame that fall in study extent
sitecovs <- covs[which(covs$patch %in% global_patches_study_extent@data$Patch),]
sitecovs <- sitecovs[order(sitecovs$Station.ID),]
#covs2 <- covs2[c(1:98, sample(99:12077, 50)),]
# site coordinates for all sites
coords <- as.matrix(sitecovs[, c("x","y")])
# site-level covariates
# site type: 1 = park; 2 = golf course; 3 = cemetery' 4 = natural area/conservation
type <- sitecovs[,"site_type"]
park <- golf <- cem <- rep(0, length(type))
park[which(type == 1)] <- 1
golf[which(type == 2)] <- 1
cem[which(type == 3)] <- 1
# tree cover (scaled)
tree <- (sitecovs[,"patch_tree"] - mean(sitecovs[,"patch_tree"]))/sd(sitecovs[,"patch_tree"])
# vegetation cover (scaled)
total_veg <- (sitecovs[,"patch_total_veg"] - mean(sitecovs[,"patch_total_veg"]))/sd(sitecovs[,"patch_total_veg"])
# water present
water <- sitecovs[,"water_present"]
# patch size (scaled)
size_km <- log(sitecovs[,"patch_size"]/1000000) # log transform to remove long right tail
size <- (size_km - mean(size_km))/sd(size_km)
# population density 2010 (scaled)
pop_10_log <- log(sitecovs[,"CMAP_pop10_dens"]) # log transform to remove long right tail
pop10 <- (pop_10_log - mean(pop_10_log))/sd(pop_10_log)
# population density 2040 (scaled)
pop40 <- (sitecovs[,"CMAP_pop40_dens"] - mean(sitecovs[,"CMAP_pop40_dens"]))/sd(sitecovs[,"CMAP_pop40_dens"])
# data matrix
sitecovs2 <- cbind(tree, total_veg, water, size, pop10, park, golf, cem)

# observation covariate - season
# vector indicating which calander season the observation was obtained
season_vec <- c(4,1,2,3,4,1,2,3,4,1,2)


data_list <- list(coords=coords,
                  y_mat=y_mat,
                  j_mat=j_mat,
                  sitecovs = sitecovs2,
                  season_vec = season_vec,
                  res_covs=res_covs)

# save so that we can load test data into a smaller environment
saveRDS(data_list, "2018-07-02_occusampler_data.RDS")

# read back in test data and load needed functions
data <- readRDS("2018-07-02_occusampler_data.RDS")

source("GIV_utility_functions.R")
source("metapop_connect_dynoccu_mcmc_sampler.R")

rm(packs,df_2_array,extractLandcover,getIntersect,occuConn,package_load)

save.image("./Summit_Files/2018-07-16_occConnC_Summit_Data.RData")
