# load required packages and functions that will be used
source("GIV_utility_functions.R")

# load in master patch layer
global_patches <- readOGR(dsn="./Data",layer="2018-03-15_Global_Patches")
as.character(unique(global_patches@data[,"Station.ID"])) # check we have all sampled sites

# name patches
prefix <- "p_"
num <- seq(1,length(global_patches), 1)
patch_name <- paste(prefix,num, sep="")
global_patches$Patch <- patch_name
global_patches@data

# calculate the area of each patch
global_patches$AREA <- area(global_patches)
global_patches$AREA <- as.numeric(as.character(global_patches$AREA))

# remove patches less than 600m^2: Approx the size of our smallest sampled patch
global_patches_reduced <- global_patches[global_patches$AREA > 600,]

###############################################################
# Calculate within patch covariates that describes each patch #
###############################################################

# calculate the proportion of each landcover variable within each patch
rm <- raster("~/Documents/GIS/LandCover2010ChicagoRegion/2010_landcover_NCS_merged2.tif")
  
#setup parallel backend to use many processors
cl <- makeCluster(7) #set number of cores
registerDoParallel(cl) # register backend

# loop through each point
dataLong <- foreach(i=1:length(global_patches),.combine=rbind, .packages=c("raster", "rgdal", "stringr")) %dopar% {
  # extract land cover data for each point, given buffer size
  Landcover <- prop.table(table(extract(rm, global_patches[i,])))
  if(length(Landcover)==0) {
    Landcover <- NA
    names(Landcover) <- "BLANK"
  }
  # summarize each site's data by proportion of each cover type
  # convert to data frame
  data.frame(id = global_patches[i,]$Patch,
              cover = names(Landcover),
              percent = as.numeric(Landcover)
  )
}

#stop cluster
stopCluster(cl)
  
# reshape data
mydf_reshape <- reshape(dataLong,idvar="id",timevar="cover", direction="wide")
mydf_reshape$id <- as.character(mydf_reshape$id)
  
# remove dataLong to save memory
#rm(dataLong)
  
# NA's to 0 
mydf_reshape[is.na(mydf_reshape)] <- 0
if(sum(grep("BLANK",colnames(mydf_reshape))) > 0){
  mydf_reshape <- mydf_reshape[,-grep("BLANK",colnames(mydf_reshape))]
}
  
# create cover name column
# order by site and order columns 1-7
df <- mydf_reshape[order((as.numeric(mydf_reshape$id))),order(names(mydf_reshape))]
  
#setup parallel backend to use many processors
cl <- makeCluster(7) #set number of cores
  
# combine 1 & 2 to get total vegetation cover
df$total_veg <- 0
df$total_veg <- parApply(cl,df[,2:3],1,sum)
  
#stop cluster
stopCluster(cl)

# create vector that indicates if open water is present in the patch
water <- df$percent.4
water[water > 0] <- 1

# calculate patch area
global_patches$AREA

# Save info
patch_covs <- data.frame(patch=df[,1], patch_tree=df[,2], patch_total_veg=df[,"total_veg"], water_present=water,
                         patch_size=as.numeric(as.character(global_patches@data[,"AREA"]))
                         )

write.csv(patch_covs, "2018-03-21_patch_covs.csv")

#################################################################################################################
# Calculate covariates around the random starting point to characterize the urban matrix around each coordinate #
#################################################################################################################

# loop through each polygon and lay a random point down inside
for (i in 1:length(slot(global_patches, "polygons"))) {
  pt = spsample(global_patches[i,], n=1, type="random", iter=100)
  if (i == 1){
    pts = pt
  } else{
    pts = rbind(pts, pt)
  }
}

# combine is.na() with over() to do check that all points are within a polygon; note that we
# need to "demote" patches to a SpatialPolygons object first
# Thanks NCEAS: https://www.nceas.ucsb.edu/scicomp/usecases/point-in-polygon
inside_patch <- !is.na(over(pts,as(global_patches, "SpatialPolygons"))) # flip the objects back and forth to thoroughly check
mean(inside_patch) # mean should be 1 in both cases - all points in patches and all patches catch a point

# double check they match: change around the numbers
plot(global_patches[20000,])
plot(pts[20000,], add=TRUE)

# create columns for xy coordinates
pts$X <-  coordinates(pts)[,1]
pts$Y <-  coordinates(pts)[,2]

# use 'over' again, this time with patches as a SpatialPolygonsDataFrame
# object, to combine the XY coordinates with the patch dataframe
global_patches$X <- over(global_patches, pts)$X
global_patches$Y <- over(global_patches, pts)$Y

# seperate sampled sites
# sampled sites
sites_sampled <- readOGR(dsn="./Data",layer="2018_02-15_Sites_Used")
plot(sites_sampled)


# use spatial overlay to find the sites that fall within a patch
# some sites fall just outside the boundry of the spatially defined patch, so we associated that site to the patch
# here we want to use the actual camera location
global_patches$Station_X <- over(global_patches, sites_sampled)$Easting
global_patches$Station_Y <- over(global_patches, sites_sampled)$Northing

# seperate the sampled patch
# NA's represent where point does not fall within spatially defined boundries of global patches layer
sampled <- global_patches[which(global_patches$Station.ID != "<NA>"),c("Patch","Station.ID", "X", "Y", "Station_X", "Station_Y")]
dim(sampled@data) # check we have all 98 sampled points
# change coords to numeric
sampled@data[,3:6] <- apply(sampled@data[,3:6], 2, as.numeric)
sites_sampled@data[,4:5] <- apply(sites_sampled@data[,4:5], 2, as.numeric)

# find missed sites and replace them with actual coordinates of our cameras
missed_sites <- sampled[is.na(sampled$Station_X),"Station.ID"]
missed_sites <- missed_sites$Station.ID
missed_sites <- as.character(missed_sites)

for (i in 1:length(missed_sites)){
  sampled@data[which(sampled@data$Station.ID==missed_sites[i]), "Station_X"] <- sites_sampled@data[which(sites_sampled@data$Station.ID==missed_sites[i]),"Easting"]
  sampled@data[which(sampled@data$Station.ID==missed_sites[i]), "Station_Y"] <- sites_sampled@data[which(sites_sampled@data$Station.ID==missed_sites[i]),"Northing"]
}

plot(sampled)

# break apart meta population of patches
m_patches <- global_patches[which(!global_patches$Patch %in% sampled$Patch),]

# turn polygons to points based on the random point or the camera location
m_points <- SpatialPointsDataFrame(m_patches@data[,c("X","Y")], data=m_patches@data, proj4string = CRS(proj4string(m_patches)))
m_points <- m_points[,-c(3:6,8:12)] # clean up
s_points <- SpatialPointsDataFrame(sampled@data[,c("Station_X","Station_Y")], data=sampled@data, proj4string = CRS(proj4string(m_patches)))
s_points <- s_points[,1:2] # clean up

# create buffer around each point
# 5km for now: average home range of coyote in Cook County
m_buffer <- gBuffer(m_points, byid=TRUE, width=5000)
s_buffer <- gBuffer(s_points, byid=TRUE, width=5000)

# write as shapefiles so that we can see them in QGIS for visual inspection
writeOGR(obj=m_buffer, dsn=".", layer="m_buffers", driver="ESRI Shapefile")
writeOGR(obj=s_buffer, dsn=".", layer="s_buffers", driver="ESRI Shapefile")


# extract data using getIntersect function

# housing units
# IL population data - SILVIS Lab
IL_pop <- readOGR(dsn="./Data", layer="2010_IL_Housing_GIV_Extent_26916") # to save space, load shapefile here
m_house_500 <- getIntersect(IL_pop, m_buffer, "HU10", n.cores=6) # housing units
s_house_500 <- getIntersect(IL_pop, s_buffer, "HU10", n.cores=6) # housing units
m_pop_500 <- getIntersect(IL_pop, m_buffer, "POP10", n.cores=6) # population
s_pop_500 <- getIntersect(IL_pop, s_buffer, "POP10", n.cores=6) # population
# clip the m_patch patch buffers to the data layer of IL pop: had some invalid geometries in the IL_pop layer so here is a workaround
extent_IL_pop <- extent(IL_pop)
extent_IL_pop@xmax <- extent_IL_pop@xmax - 4700
m_buffer_clip <- raster::crop(m_buffer, extent_IL_pop)
#writeOGR(obj=m_buffer_clip, dsn=".", layer="clip_test2", driver="ESRI Shapefile") # write as shapefile to visually inspect
# calculate area of each buffer in the clipped layer
m_buffer_clip$AREA <- raster::area(m_buffer_clip)/1000000 # convert to sq km's
# houseing density
m_house_500$density <- m_house_500$dat/m_buffer_clip$AREA
s_house_500$density <- s_house_500$dat/max(m_buffer_clip$AREA)
m_pop_500$density <- m_pop_500$dat/m_buffer_clip$AREA
s_pop_500$density <- s_pop_500$dat/max(m_buffer_clip$AREA)
# clean up environment to save memory: remove shapefile when done
#rm(IL_pop, m_buffer_clip, extent_IL_pop)

# projected housing units
# CMAP population data -  CMAP GoTo2040
CMAP_pop <- readOGR(dsn=path.expand("./Data"), layer="CMAP_2040_Forecast_26916") # CMAP projected forecast at subzone level (~800m resolution)
m_PHH10_500 <- getIntersect(CMAP_pop, m_buffer, "PHH10", n.cores=6) # 2010 population in households
m_PHH40_500 <- getIntersect(CMAP_pop, m_buffer, "PHH40", n.cores=6) # projected population in households
s_PHH10_500 <- getIntersect(CMAP_pop, s_buffer, "PHH10", n.cores=6) # 2010 population in households
s_PHH40_500 <- getIntersect(CMAP_pop, s_buffer, "PHH40", n.cores=6) # projected population in households
# clip the m_patches to the data layer
m_buffer_cmap_clip <- raster::crop(m_buffer, CMAP_pop)
#writeOGR(obj=m_buffer_cmap_clip, dsn=".", layer="clip_test_cmap", driver="ESRI Shapefile") # write as shapefile to visually inspect
# calculate area of each buffer in the clipped layer
m_buffer_cmap_clip$AREA <- raster::area(m_buffer_cmap_clip)/1000000 # convert to sq km's
# 2010 and projected 2040 housing and population densities
m_PHH10_500$density <- m_PHH10_500$dat/m_buffer_cmap_clip$AREA
s_PHH10_500$density <- s_PHH10_500$dat/max(m_buffer_cmap_clip$AREA)
m_PHH40_500$density <- m_PHH40_500$dat/m_buffer_cmap_clip$AREA
s_PHH40_500$density <- s_PHH40_500$dat/max(m_buffer_cmap_clip$AREA)
# remove shapefiles when done to save memory
#rm(CMAP_pop, m_buffer_cmap_clip)

# vacant housing units from 2010 census data
census_dat <- readOGR(dsn="./Data", layer="2010_Census_Block_GIV_Extent_26916") # IL Tiger Data Products census data
m_vacant_500 <- getIntersect(census_dat, m_buffer, "DP0180003", n.cores=6)
s_vacant_500 <- getIntersect(census_dat, s_buffer, "DP0180003", n.cores=6)
# all buffers fit inside the census layer so just calculate density
m_vacant_500$density <- m_vacant_500$dat/77.25425 # hard coded the max area from previous layers
s_vacant_500$density <- s_vacant_500$dat/77.25425
# remove shapefile to same memory
#rm(census_dat)

# highway density
# definition of a highway -
# http://www.ilga.gov/legislation/ilcs/ilcs4.asp?DocName=060500050HArt%2E+2+Div%2E+2&ActID=1745&ChapterID=45&SeqStart=1100000&SeqEnd=3200000
roads_dat <- readOGR(dsn="./Data", layer="2012_Roads_GIV_Extent_26916") # IL roads - IDOT
m_road_length_500 <- getIntersect(roads_dat, m_buffer, "LENGTH", getLength=TRUE, n.cores=6) # total length of roads in buffer
s_road_length_500 <- getIntersect(roads_dat, s_buffer, "LENGTH", getLength=TRUE, n.cores=6) # total length of roads in buffer
# clip the m_patch buffers to the data layer of road layer
m_buffer_roads_clip <- raster::crop(m_buffer, roads_dat)
#writeOGR(obj=m_buffer_roads_clip, dsn=".", layer="clip_test_roads", driver="ESRI Shapefile") # write as shapefile to visually inspect
m_buffer_roads_clip$AREA <- raster::area(m_buffer_roads_clip)/1000000 # convert to sq km's
m_road_length_500$Density <- m_road_length_500$dat/m_buffer_roads_clip$AREA
s_road_length_500$Density <- s_road_length_500$dat/max(m_buffer_roads_clip$AREA)
# remove shapefile to same memory
#rm(roads_dat, m_buffer_roads_clip)

# extract dominate landcover using extractLandcover function from the 2011 NLCD
m_dom_lancov <- extractLandcover(m_buffer, n.cores=6)
s_dom_lancov <- extractLandcover(s_buffer, n.cores=6)

# check for duplicate/equal cover at some sites
checkDuplicates <- function(dom){
  dups <- dom[which(duplicated(m_dom_lancov$id)),] # 37 sites
  dups2 <- dom[which(dom$id %in% dups$id),] # duplicates
  
  return(list(dups, dups2))
}

m_dups <- checkDuplicates(m_dom_lancov) # no duplicates
s_dups <- checkDuplicates(s_dom_lancov) # no duplicates

table(m_dom_lancov$cover) # seven types of landcover
table(s_dom_lancov$cover) # only two types of landcover
# ultimately will not use this covariate

# categorize the meta patches layer into 4 land use types following our study design

# new field
m_patches$TYPE <- 0

# park = 1
1 -> m_patches[which(m_patches$CATEGORY == "City Park" |
                  m_patches$CATEGORY == "Open Space, Primarily Recreation" |
                  m_patches$CATEGORY == "Parks and Recreation" |
                  m_patches$CATEGORY == "County Park" |
                  m_patches$LANDUSE == 3100), "TYPE"]
# golf course = 2
2 -> m_patches[which(m_patches$CATEGORY == "Golf Course" | m_patches$LANDUSE == 3200), "TYPE"]

# cemeteries = 3

3 -> m_patches[which(m_patches$CATEGORY == "Cemetery" | m_patches$LANDUSE == 1360), "TYPE"]

# natural areas/open spaces = 4
4 -> m_patches[-which(m_patches$TYPE == 1 | m_patches$TYPE == 2 | m_patches$TYPE ==3), "TYPE"]

table(m_patches$TYPE) # check the distribution

# categorize the meta patches layer into 4 land use types following our study design
# seperate out the patches that we sampled
s_patches <- global_patches[which(global_patches$Patch %in% sampled$Patch),]

# our sampling desing doesn't fall cleanly with land use layer categories, so I hard coded the site_type
s_patches$TYPE <- c(4,4,4,3,4,4,4,4,4,4,4,4,4,4,2,4,4,4,4,4,4,4,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,4,4,4,4,
                    4,4,1,1,1,3,4,4,1,4,4,1,4,4,1,4,4,1,4,1,2,1,4,1,
                    1,1)

table(s_patches$TYPE) # check the distribution

# add to dataset

# covariate data frames
m_covs <- data.frame(patch=m_buffer$Patch,
                     house_dens=m_house_500$density, pop_dens=m_pop_500$density,
                     CMAP_pop10_dens=m_PHH10_500$density, CMAP_pop40_dens=m_PHH40_500$density,
                     vacant_house_dens=m_vacant_500$density,
                     road_dens=m_road_length_500$Density,
                     site_type = as.character(m_patches$TYPE),
                     Station.ID=NA,
                     x=as.numeric(coordinates(m_points)[,1]), y=as.numeric(coordinates(m_points)[,2]))

s_covs <- data.frame(patch=s_buffer$Patch,
                     house_dens=s_house_500$density, pop_dens=s_pop_500$density,
                     CMAP_pop10_dens=s_PHH10_500$density, CMAP_pop40_dens=s_PHH40_500$density,
                     vacant_house_dens=s_vacant_500$density,
                     road_dens=s_road_length_500$Density,
                     site_type = as.character(s_patches$TYPE), 
                     Station.ID=s_patches$Station.ID,
                     x=as.numeric(coordinates(s_points)[,1]), y=as.numeric(coordinates(s_points)[,2]))

# combine the datasets
covs1 <- rbind(m_covs, s_covs)
covs1 <- covs[order(as.numeric(as.character(substr(covs$patch,3,10)))),] # order by patch id
covs <-data.frame(patch_covs, covs1[,-1])

# look at correlation between variables
C <- cor(covs[,c(2:3,5:11)])
corrplot.mixed(C)

# save to source into model script
saveRDS(covs, "./Data/2018-03-27_patch_covariates.RDS")

# create raster layers out of population data to use for resistance

# convert popluation data to sf class
CMAP_pop_sf <- st_as_sf(CMAP_pop)

# convert to raster using the NDVI raster layer as a template
pop10_raster <- fasterize(CMAP_pop_sf, ndvi_res, field = "PHH10", fun="min", background = mean(CMAP_pop_sf$PHH10))
pop40_raster <- fasterize(CMAP_pop_sf, ndvi_res, field = "PHH40", fun="min")

# plot
plot(pop10_raster)

# write to explore in QGIS
writeRaster(pop10_raster, "./Data/CMAP_PHH10", format = "GTiff")
writeRaster(pop40_raster, "CMAP_PHH40", format = "GTiff")

save.image("2018-04-12_GIV_Workspace.RData")

######## HARD HAT AREA - PROCEED WITH CAUTION ###############

# Euclidian Distance
x.loc <- vector("list", nrow(x))
x.list <- vector("list", nrow(x))
for(j in 1:nrow(x)){
  x.loc[[j]] <- which(sqrt((x[j,1]-x[,1])^2 + (x[j,2]-x[,2])^2) < 10000)
  x.list[[j]] <- x[x.loc[[j]],]
}


sapply(x.list, nrow)
