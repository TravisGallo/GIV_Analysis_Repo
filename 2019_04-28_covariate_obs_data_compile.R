# load required packages and functions that will be used
source("GIV_utility_functions.R")

# load in master patch layer
global_patches <- st_read(dsn="./Data/GIS",layer="2018-03-15_Global_Patches")
#as.character(unique(global_patches$Station.ID)) # check we have all sampled sites

# for computational reasons we need to find ways to reduce patches that make sense
# here we remove linear righ-o-ways from rail and utility
# could be habitat patches, but less likely than other landuses and we have to make decisions
global_patches2 <- global_patches[-which(global_patches$LANDUSE == "1511" |
                                           global_patches$LANDUSE == "1561"),]

## We need to create a study area and crop the global patches to this study area
# build a study extent that is 25km from origin
# we then buffer out from that 10km

# load end points of each transect
end_points <-read.csv("./Data/transect_endpoints.csv")

# shorten our transect lines by 25km
new_point <- matrix(NA, nrow = nrow(end_points), ncol = 2)
for(i in 1:nrow(end_points)){
  # distance of full line
  d <- sqrt((end_points[i,4] - end_points[i,2])^2 + (end_points[i,5] - end_points[i,3])^2)
  # the proportion our new distance is of the full distance
  t <- 25000/d
  # calculate the new end points
  new_point[i,1] <- (1-t)*end_points[i,2] + t*end_points[i,4]
  new_point[i,2] <- (1-t)*end_points[i,3] + t*end_points[i,5]
}


# create a list of points (3 end points and an origin)
# note the origin needs to be repeated at the end to create closure
# we also created points to the east to allow for our study extent to grab lakefront
# note southeast point is truncated to the southern extent of patch layer
points <- list(rbind(as.matrix(end_points[1,2:3]),
                     c(end_points[1,2],new_point[1,2]),
                     new_point,
                     c(end_points[1,2],end_points[1,3] - 9000),
                     as.matrix(end_points[1,2:3])))
# turn this into a polygon
end_points_poly <- st_polygon(x = points)

# create buffer around polygon
full_extent <- st_buffer(end_points_poly, dist = 3000)

# crop the patches from within the study area boundary
patches_full_extent <- global_patches2[full_extent,]


# To keep our sampled patches as independent patches, we will seperate them out before union
sampled_patches <- patches_full_extent[!is.na(patches_full_extent$Station.ID),"Station.ID"]
# seperate out the non-sampled patches
M_patches <- patches_full_extent[is.na(patches_full_extent$Station.ID),]

# union the non-sampled patches so any that are within 25m of each other become one patch
M_patches_union <- M_patches %>%
  # break apart all multi-part polys
  st_cast("POLYGON") %>%
  # close < 25 m gaps by creating a buffer around the each patch
  st_buffer(dist = 50) %>%
  # combine any that touch or intersect
  st_union(.) %>% 
  # turn them back into individual polygons
  st_cast("POLYGON") %>%
  # remove the 25 meter buffer we created from each
  st_buffer(dist = -50) %>%
  # to cast them all to Polygons, you first have to match them all as Multipolygons
  st_cast("MULTIPOLYGON") %>%
  # now cast as polygons - this will break up the few multipolygons into individual polygons
  st_cast("POLYGON") %>%
  # put data back into sf
  st_sf(.)

# add data frame to new sf object with Station ID column
M_patches_union$Station.ID <- rep(NA, length(M_patches_union))

# bring back our sampled patches
patches_combined <- rbind(M_patches_union, sampled_patches)

# calculate area of each patch
patches_combined$area <- as.numeric(st_area(patches_combined))/1e6
# our smallest sampled patch will be the area cutoff
min_area <- min(patches_combined$area[which(!is.na(patches_combined$Station.ID))])

# remove M_patches smaller than our smallest sampled patch
patches <- patches_combined[which(patches_combined$area >= min_area),]

# order patches so sampled sites are always first
patches<- patches[order(patches$Station.ID),]

# give each polygon a unique name for organization
patches$patch_name <- paste0("p_", seq(1,nrow(patches), 1))

# create starting point in the center of the patch
# our sampled patches are multipolygons and this will put it in the center of the multipolygon
# we will replace the points for sampled patches with our actual camera location
# so these points will be replaced with real locations
patch_centroids <- st_centroid(patches)[,c("Station.ID", "patch_name")]

# load camera locations
camera_points <- st_read(dsn="./Data/GIS", layer="2018_02-15_Sites_Used")
# keep only the ones we use for this analysis
sampled_points <- camera_points[which(camera_points$Station.ID %in% 
                                        sampled_patches$Station.ID),"Station.ID"]
sampled_points <- sampled_points[order(sampled_points$Station.ID),]

# give the points a patch name to match the patch dataset
# since they are ordered it will match p_1:p_97 in the patch covariate dataset
sampled_points[,"patch_name"] <- paste0("p_", seq(1,nrow(sampled_points), 1))

# combine the sampled points with the remaining unsampled points
all_points <- rbind(sampled_points, 
                    patch_centroids[(nrow(sampled_points)+1):nrow(patch_centroids),])
all_points <- all_points[order(all_points$Station.ID),]

# get coordinates to use in model
coords <- st_coordinates(all_points)

###############################################################
# Calculate within patch covariates that describes each patch #
###############################################################

# calculate the proportion of each landcover variable within each patch
rm <- raster("./Data/GIS/2010_HRLC_cropped_compressed.tif")

# create a list of patches so that you can work with the sf object using lapply
patch_list <- split(patches, seq(nrow(patches)))

cl <- makeCluster(detectCores()-6)
clusterEvalQ(cl,library("sf"))
clusterEvalQ(cl, library("raster"))
clusterExport(cl, c("rm", "patch_list"))

# extract from each patch
raster_extract <- parLapply(cl, patch_list, function(x) {
  prop.table(table(raster::extract(rm, x)))
})

stopCluster(cl)

prop_table <- do.call(bind_rows, raster_extract)

prop_table[is.na(prop_table)] <- 0
 
# combine 1 & 2 to get total vegetation cover
total_veg <- rowSums(prop_table[,1:2])

# water present (binary)
water_present <- rep(0, nrow(prop_table))
water_present[prop_table[,4] > 0] <- 1

# Save info
patch_covs <- data.frame(patch = patches$patch_name, 
                         patch_tree = data.frame(prop_table[,1]), 
                         patch_total_veg = total_veg, 
                         water_present = water_present,
                         patch_size = patches$area)

colnames(patch_covs)[2] <- "patch_tree"


# save since this takes a couple days
write.csv(patch_covs, paste0("./Data/",Sys.Date(),"_raster_extract.csv"))


###########################################################
# Calculate population density around the starting point #
##########################################################

# buffer points to take into account urban coyote home range
# see Roland Kays coyote project protocol for citations
points_buffered <- st_buffer(all_points, dist= 2000)

# now intersect population data with buffers to get mean pop density

# create a list of patches so that you can work with the sf object using lapply
point_list <- split(points_buffered, seq(nrow(points_buffered)))

# current and projected population data
# for 2010 we are still going to use CMAP GoTo2040 dataset
# using 2010 data for our base model will match temporally with other site level covariates
CMAP_pop <- st_read(dsn="./Data/GIS", layer="CMAP_2040_Forecast_26916")

cl <- makeCluster(detectCores()-6)
clusterEvalQ(cl,library("sf"))
clusterExport(cl, c("CMAP_pop", "point_list"))

# extract from each patch
pop_extract <- parLapply(cl, point_list, function(x) {
  pop_intersect <- st_intersection(x, CMAP_pop[,"PHH10"])
  pop10_dens <- sum(as.numeric(as.character(pop_intersect$PHH10)))/(pi*(2^2))
  return(pop10_dens)
})

stopCluster(cl)

# turn into a vector
pop10_dens <- do.call(c, pop_extract)
# take the log of 
pop10_log <- as.vector(log(pop10_dens))



# CMAP population data -  CMAP GoTo2050 - to get future projections
# CMAP projected forecast at subzone level (~800m resolution)
CMAP_ONTO2050 <- st_read(dsn="./Data/GIS", layer="CMAP_ONTO2050_ForecastByLAZ") 
# reproject to match all other data
CMAP_ONTO2050_reproject <- st_transform(CMAP_ONTO2050, st_crs(all_points))


cl <- makeCluster(detectCores()-6)
clusterEvalQ(cl,library("sf"))
clusterExport(cl, c("CMAP_ONTO2050_reproject", "point_list"))

# extract from each patch
pop_extract2 <- parLapply(cl, point_list, function(x) {
  pop_intersect <- st_intersection(x, CMAP_ONTO2050_reproject[,c("PHH_2020","PHH_2030",
                                                                 "PHH_2040","PHH_2050")])
  pop20_dens <- sum(as.numeric(as.character(pop_intersect$PHH_2020)))/(pi*(2^2))
  pop30_dens <- sum(as.numeric(as.character(pop_intersect$PHH_2030)))/(pi*(2^2))
  pop40_dens <- sum(as.numeric(as.character(pop_intersect$PHH_2040)))/(pi*(2^2))
  pop50_dens <- sum(as.numeric(as.character(pop_intersect$PHH_2050)))/(pi*(2^2))
  
  return(c(pop20_dens,pop30_dens,pop40_dens,pop50_dens))
})

stopCluster(cl)

# collapse into a matrix
pop_future_dens <- do.call(rbind, pop_extract2)
# take the log
pop_future_log <- log(pop_future_dens)
colnames(pop_future_log) <- c("POP20","POP30","POP40","POP50")
pop_future_log <- data.frame(pop_future_log)


##############################################
# Combine all the site-level data we created #
##############################################

site_covs <- data.frame(patch_covs,
                        POP10 = pop10_log,
                        pop_future_log)


# look at correlation between variables
C <- cor(site_covs[,-1])
corrplot.mixed(C)

# scale all covariates
site_covs$tree_scaled <- with(site_covs, scale(site_covs$patch_tree))
site_covs$total_veg_scaled <- with(site_covs, scale(site_covs$patch_total_veg))
site_covs$size_scaled <- with(site_covs, scale(site_covs$patch_size))
site_covs$POP10_scaled <- with(site_covs, scale(site_covs$POP10))
site_covs$POP20_scaled <- with(site_covs, scale(site_covs$POP20))
site_covs$POP30_scaled <- with(site_covs, scale(site_covs$POP30))
site_covs$POP40_scaled <- with(site_covs, scale(site_covs$POP40))
site_covs$POP50_scaled <- with(site_covs, scale(site_covs$POP50))


# save
saveRDS(site_covs, paste0("./Data/",Sys.Date(),"_patch_covariates.rds"))


#############################
# Create Resistance Layers #
############################

# create a buffered study extent to crop rasters
study_extent_buffered <- st_buffer(full_extent, dist = 5000)
# convert to spatial so we can use it to raster::crop and raster::mask
study_extent_sp <- as(study_extent_buffered, "Spatial")

## NDVI ##
# load NDVI raster and increase to 300x300 resolution (check actual line of code)
ndvi_coarse <- aggregate(raster("./Data/GIS/NDVI_to_Resistence.tif"), fact=10, fun=mean)
# crop to study extent
ndvi_study <- mask(crop(ndvi_coarse, study_extent_sp), study_extent_sp)
# scale raster for use in model
ndvi_scaled <- scale(ndvi_study)

## Population ##
#rasterize 2010 population data
pop10_raster <- fasterize(CMAP_pop, ndvi_coarse, field = "PHH10", fun="max")
# crop to study extent
pop10_study_extent <- mask(crop(pop10_raster, study_extent_sp), study_extent_sp)
# take log values of population with addition of constant 1
pop10_log <- log(pop10_study_extent + 1)
# scale raster to use in model
pop10_scaled <- scale(pop10_log)

# rasterize future population data

# create a list to hold the population rasters
pop_raster_list <- vector("list", 4)
# vector of fields to loop through and make a raster of each
field <- c("PHH_2020","PHH_2030","PHH_2040","PHH_2050")

for(i in 1:length(pop_raster_list)){
  # convert population data to raster using the NDVI raster as a template
  Raster <- fasterize(CMAP_ONTO2050_reproject, ndvi_coarse, field = field[i], fun="max")
  # crop to study extent
  pop_raster_list[[i]] <- mask(crop(Raster, study_extent_sp), study_extent_sp)
}

# take the log of each raster with addition of constant 1
pop_future_log <- lapply(pop_raster_list, function(x) { log(x + 1) })
# scale the rasters to use in later models
pop_future_list_scaled <- lapply(pop_future_log, function(x) { scale(x) })
names(pop_future_list_scaled) <- c("pop20", "pop30", "pop40", "pop50")


# create a raster list with the data we created
raster_list <- list(ndvi = ndvi_scaled, 
                    pop10 = pop10_scaled, 
                    future_pop = pop_future_list_scaled)


########################
# Process Camera Data #
#######################

# observation and sampling data
full_site_names <- as.character(read.table("./Data/sites_used_in_sp10_sp13_analysis_6_1_17.txt", 
                              header = TRUE)[,1])

# randomly chose 1 site for patches that had multiple sites. Remove unused sites
sampled_sites <- as.character(sampled_points$Station.ID)

# book keeping: change the name of full sites to match sampled sites
full_site_names <- gsub('0$', 1, full_site_names)

species_names <- read.table("./Data/species_used_in_sp10_sp13_analysis.txt", header=TRUE,
                            stringsAsFactors = FALSE)[,1]

# load species specific data - remove sites that we do not use in this analysis
zarray <- df_2_array(read.table("./Data/z_matrix_sp10_sp13_6_1_17.txt", header = TRUE, 
                               sep = "\t"))[,,-c(1,2)]
# name the dimensions for book keeping
dimnames(zarray) <- list(species_names,full_site_names,seq(1,11,1))

# isolate a single species to run model - here we use raccoon
zmat <- zarray["Raccoon",,]

# keep only the sites that we used for this analysis
z_mat <- zmat[which(rownames(zmat) %in% sampled_sites),]

# build y-array (detection)
yarray <- df_2_array(read.table("./Data/y_matrix_sp10_sp13_6_1_17.txt", 
                               header = TRUE, sep = "\t"))[,,-c(1,2)]
# name the dimensions for book keeping
dimnames(yarray) <- list(species_names,full_site_names,seq(1,11,1))

# isolate a single species to run model - here we use raccoon
ymat <- yarray["Raccoon",,]

# keep only the sites that we used for this analysis
y_mat <- ymat[which(rownames(ymat) %in% sampled_sites),]

#j-matrix (n days sampled)
jmat <- as.matrix(read.table("./Data/j_matrix_sp10_sp13_6_1_17.txt", 
                              header = TRUE, sep = "\t"))[,-c(1,2)]
rownames(jmat) <- full_site_names
# keep only the sites that we used for this analysis
j_mat <- jmat[which(rownames(jmat) %in% sampled_sites),]


# take a look at the raw data
sum(y_mat, na.rm=TRUE) # Total dets
colSums(y_mat, na.rm=TRUE) # Dets per season
sum(rowSums(y_mat, na.rm = TRUE) > 0)/nrow(y_mat) # Proportion of sites occupied
rowSums(y_mat, na.rm=TRUE) # Det per site


###########################################################################
# Combine all of the data we have created into a list to be used in model #
###########################################################################
site_covs <- readRDS("./Data/2019-04-15_patch_covariates.rds")

data_list <- list(coords=coords,
                  z_mat=z_mat,
                  y_mat=y_mat,
                  j_mat=j_mat,
                  sitecovs = site_covs,
                  res_covs=raster_list)

saveRDS(data_list, paste0(Sys.Date(),"_DynOccu_Connectivity_DataList.rds"))


## Prep data to load to ARGO
# this will be the bare minimum we need to send to ARGO cluster
# remove everything
rm(list=ls())
# load functions
source("GIV_utility_functions.R")
# removoe unneeded functions
rm(list = c("packs","df_2_array","extractLandcover","getIntersect","package_load"))
# load data
data <- readRDS("2019-04-29_DynOccu_Connectivity_DataList.rds")
# save workspace to load into ARGO
save.image(paste0(Sys.Date(),"_Connectivity_Workspace_ARGO.RData"))

