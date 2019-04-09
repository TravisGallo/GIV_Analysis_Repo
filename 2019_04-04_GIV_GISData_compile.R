# load required packages and functions that will be used
source("GIV_utility_functions.R")

# load in master patch layer
global_patches <- st_read(dsn="./Data/GIS",layer="2018-03-15_Global_Patches")
#as.character(unique(global_patches$Station.ID)) # check we have all sampled sites

## We need to create a study area and crop the global patches to this study area

# build a study extent that is a square around each transect 
# with a diameter as long as the transect + 10km
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
  
}

# squares into polygons
square_polys <- lapply(square, function(x){ 
  st_sfc(st_polygon(list(x)), crs = st_crs(global_patches)) 
  })

# combine them into a multi feature polygon
polys_combined <- do.call(c, square_polys)

# create one single polygon of the study area boundaries
study_boundary <- st_union(polys_combined)

# crop the patches from within the study area boundary
study_patches <- global_patches[study_boundary,]

# To keep our sampled patches as independent patches, we will seperate them out before union
sampled_patches <- study_patches[!is.na(study_patches$Station.ID),"Station.ID"]
# seperate out the non-sampled patches
M_patches <- study_patches[is.na(study_patches$Station.ID),]

# union the non-sampled patches so any that are within 25m of each other become one patch
M_patches_union <- M_patches %>%
  # break apart all multi-part polys
  st_cast("POLYGON") %>%
  # close < 25 m gaps by creating a buffer around the each patch
  st_buffer(dist = 25) %>%
  # combine any that touch or intersect
  st_union(.) %>% 
  # turn them back into individual polygons
  st_cast("POLYGON") %>%
  # remove the 25 meter buffer we created from each
  st_buffer(dist = -25) %>%
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
camera_points <- st_read(dsn="./Data/GIS",layer="2018_02-15_Sites_Used")
# keep only the ones we use for this analysis
sampled_points <- camera_points[which(camera_points$Station.ID %in% sampled_patches$Station.ID),"Station.ID"]
sampled_points <- sampled_points[order(sampled_points$Station.ID),]

# give the points a patch name to match the patch dataset
# since they are ordered it will match p_1:p_97 in the patch covariate dataset
sampled_points[,"patch_name"] <- paste0("p_", seq(1,nrow(sampled_points), 1))

# combine the sampled points with the remaining unsampled points
all_points <- rbind(sampled_points, patch_centroids[98:nrow(patch_centroids),])

# get coordinates to use in model
coords <- st_coordinates(all_points)

###############################################################
# Calculate within patch covariates that describes each patch #
###############################################################

# calculate the proportion of each landcover variable within each patch
rm <- raster("./Data/GIS/2010_HRLC_cropped_compressed.tif")

# create a list of patches so that you can work with the sf object using lapply
patch_list <- split(patches, seq(nrow(patches)))

cl <- makeCluster(detectCores()-2)
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


write.csv(patch_covs, "2019-04-04_patch_covs.csv")


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

cl <- makeCluster(detectCores()-2)
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
pop10_df <- data.frame(POP10 = pop10_dens)

# CMAP population data -  CMAP GoTo2050 - to get future projections
# CMAP projected forecast at subzone level (~800m resolution)
CMAP_ONTO2050 <- st_read(dsn="./Data/GIS", layer="CMAP_ONTO2050_ForecastByLAZ") 
# reproject to match all other data
CMAP_ONTO2050_reproject <- st_transform(CMAP_ONTO2050, st_crs(all_points))


cl <- makeCluster(detectCores()-2)
clusterEvalQ(cl,library("sf"))
clusterExport(cl, c("CMAP_ONTO2050_reproject", "point_list"))

# extract from each patch
pop_extract2 <- parLapply(cl, point_list, function(x) {
  pop_intersect <- st_intersection(point_list[[1]], CMAP_ONTO2050_reproject[,c("PHH_2020","PHH_2030",
                                                                 "PHH_2040","PHH_2050")])
  pop20_dens <- sum(as.numeric(as.character(pop_intersect$PHH_2020)))/(pi*(2^2))
  pop30_dens <- sum(as.numeric(as.character(pop_intersect$PHH_2030)))/(pi*(2^2))
  pop40_dens <- sum(as.numeric(as.character(pop_intersect$PHH_2040)))/(pi*(2^2))
  pop50_dens <- sum(as.numeric(as.character(pop_intersect$PHH_2050)))/(pi*(2^2))
  
  return(c(pop20_dens,pop30_dens,pop40_dens,pop50_dens))
})

stopCluster(cl)

pop_future_dens <- do.call(rbind, pop_extract2)
colnames(pop_future_dens) <- c("POP20","POP30","POP40","POP50")


##############################################
# Combine all the site-level data we created #
##############################################

site_covs <- data.frame(patch_covs,
                        POP10 = pop10_dens,
                        pop_future_dens)


# look at correlation between variables
C <- cor(site_covs[,c(2:4,4:5)])
corrplot.mixed(C)

# save to source into model script
saveRDS(site_covs, "./Data/2019-04-04_patch_covariates.rds")


#############################
# Create resistance layers #
############################

# create a buffered study extent to crop rasters
study_extent_buffered <- st_buffer(study_boundary, dist = 5000)
# convert to spatial so we can use it to raster::crop and raster::mask
study_extent_sp <- as(study_extent_buffered, "Spatial")

## NDVI ##
# load NDVI raster and increase to 900x900 resolution
ndvi_coarse <- aggregate(raster("./Data/GIS/NDVI_to_Resistence.tif"), fact=30, fun=mean)
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

## Patch Indicator
# patch indicator indicating that habitat patches have 0 resistance
# 0 resistence when converted to conductance (1/resistence) gives us infinity, 
# which is what we want but not computationally possible
# therefore, we give it basically 0 resistence (0.0001) 
#so that it converts to basically infinity (some really high number) when divided by 1

# all geometries much match to use fasterize
patches_polygons <- st_cast(patches, "MULTIPOLYGON") %>% st_cast("POLYGON")
# fasterize patches
patch_indicator <- mask(crop(fasterize(patches_polygons, ndvi_coarse, field = "area", fun="max"),
                        study_extent_sp), study_extent_sp)
# find NA cells
val_to_replace0 <- !is.na(values(patch_indicator))
values(patch_indicator)[val_to_replace0] <- 1
val_to_replace1 <- is.na(values(patch_indicator))
values(patch_indicator)[val_to_replace1] <- 0

# cells from the ndvi and population raster become 0 to estimate transition cost within patches seperately from the matrix
to_0 <- values(patch_indicator) == 1
values(ndvi_scale)[to_0] <- 0
values(pop10_scale)[to_0] <- 0

