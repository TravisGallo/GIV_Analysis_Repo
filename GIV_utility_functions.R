# load required packages

package_load<-function(packages = NULL, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# the packages needed for this analysis
packs <- c("dplyr", "reshape2", "runjags", "mcmcplots",
           "runjags", 'parallel','raster','rgdal','stringr',
           'foreach','doParallel', 'sf', 'spatialEco', 'rgeos',
           'maptools', 'cleangeo', 'gdistance', 'corrplot',
           'fasterize')

package_load(packs)

# function to intersect background data layer (vector) with buffer layers (vector) and calculate needed data
# layer 1 should be the input layer with data to extract, layer 2 is the intersecting layer
getIntersect <- function(layer1, layer2, columnName, getMean=FALSE, getLength=FALSE, getArea=FALSE, n.cores){
  
  # set and register cores
  cl <- makeCluster(n.cores) #set number of cores
  registerDoParallel(cl) # register backend
  
  # use a foreach loop to parallel and loop through each feature of the buffer layer
  dat <- foreach(i=1:length(layer2), .combine=c, .packages=c("raster", "dplyr")) %dopar% {
    # intersect two vector layers - layer1 is the input layer
    buffer_intersect <- raster::intersect(layer1,layer2[i,])
    # loop through each feature and extract the needed data
    if(is.null(buffer_intersect)){
      buffer_intersect <- 0
    } else {
      if(getLength){
        sum(SpatialLinesLengths(buffer_intersect))
      } else {
        if(getArea) {
          data_intersect$area <- raster::area(data_intersect)
          # calculate the areas of each landuse category
          data_table <- as.data.frame(data_intersect@data %>% 
                                        group_by_at(vars(one_of(columnName))) %>%
                                        summarise(total=sum(as.numeric(as.character(area)))))
          # extract the land use category that has the greatest area
          as.character(data_table[which(data_table$total == max(data_table$total, na.rm=TRUE)),1])
        } else {
          if(getMean){
            mean(as.numeric(as.character((buffer_intersect@data[,columnName]))))
          } else {
            sum(as.numeric(as.character((buffer_intersect@data[,columnName]))))
          }
        }
      }
    }
  }
  
  #stop core cluster
  stopCluster(cl)
  
  output <- data.frame(Patch=layer2@data$Patch, dat)
  
  return(output) # returns vector of values in order of patch_ID
}

# extract the proportion of land cover categories for each polygon
extractLandcover <- function(layer, n.cores){
  
  # load nlcd raster
  nlcd <- raster("./Data/nlcd_2011_GIV_Extent_26916.tif")
  
  # setup parallel backend to use many processors
  cl <- makeCluster(n.cores) #set number of cores
  registerDoParallel(cl) # register backend
  
  # loop through each patch
  dom_landcover <- foreach(i=1:length(layer),.combine=rbind, .packages=c("raster", "rgdal", "stringr")) %dopar% {
    # extract land cover data for each point, given buffer size
    Landcover <- prop.table(table(extract(nlcd, layer[i,])))
    if(length(Landcover)==0) {
      Landcover <- NA
      names(Landcover) <- "BLANK"
    }
    # summarize each site's data by proportion of each cover type
    # convert to data frame
    data.frame(id = layer[i,]$Patch,
               cover = names(Landcover),
               percent = as.numeric(Landcover)
    )
  }
  
  #stop cluster
  stopCluster(cl)
  
  # combine categories
  dom_landcover[which(dom_landcover$cover == 43), "cover"] <- "41" # combine mixed forest with decidious forest (forest)
  dom_landcover[which(dom_landcover$cover == 82), "cover"] <- "81" # combine pasture and cultivated crops (agriculture)
  dom_landcover[which(dom_landcover$cover == 95), "cover"] <- "90" # combine emergent wetland with woody wetland (wetland)
  dom_landcover[which(dom_landcover$cover == 71), "cover"] <- "52" # combine grassland with shrub/scrub (grassland/shrubland)
  dom_landcover[which(dom_landcover$cover == 31), "cover"] <- "24" # combine barren land (quarries and industrial) with high developed
  dom_landcover[which(dom_landcover$cover == 23), "cover"] <- "22" # combine medium and low urban (med/low developed)
  
  # sum across those that were combined
  dom_comb <- as.data.frame(dom_landcover) %>% group_by(id,cover) %>%
    summarise(percent=sum(percent)) %>% as.data.frame
  
  # filter to get the lancover class that has the highest proportion
  dom <- as.data.frame(dom_comb %>% group_by(id) %>%
                         filter(percent == max(percent)))
  
  return(dom)
}