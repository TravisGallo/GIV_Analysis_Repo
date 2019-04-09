data_list <- list(patches = patches, raster = rm)

saveRDS(data_list, "2019-03-23_SpatialData_ARGO.rds")


# For ARGO
library(sf)
library(raster)
library(dplyr)

data_list <- readRDS("2019-03-23_SpatialData_ARGO.rds")
patches <- data_list[[1]]
rm <- data_list[[2]]

#setup parallel backend to use many processors
cl <- makeCluster(detectCores()) #set number of cores
registerDoParallel(cl) # register backend

# loop through each point
dataLong <- foreach(i=1:nrow(patches),.combine=rbind, .packages=c("raster", "rgdal", "stringr")) %dopar% {
  
  for(i in 1:nrow(patches)){
    # extract land cover data for each point, given buffer size
    Landcover <- prop.table(table(extract(rm, patches[i,])))
    if(length(Landcover)==0) {
      Landcover <- NA
      names(Landcover) <- "BLANK"
    }
    # summarize each site's data by proportion of each cover type
    # convert to data frame
    data.frame(id = patches[i,]$patch_name,
               cover = names(Landcover),
               percent = as.numeric(Landcover)
    )
  }
  
  #stop cluster
  stopCluster(cl)
  
  # reshape data
  mydf_reshape <- reshape(dataLong, idvar="id", timevar="cover", direction="wide")
  mydf_reshape$id <- as.character(mydf_reshape$id)
  
  # NA's to 0 
  mydf_reshape[is.na(mydf_reshape)] <- 0
  if(sum(grep("BLANK",colnames(mydf_reshape))) > 0){
    mydf_reshape <- mydf_reshape[,-grep("BLANK",colnames(mydf_reshape))]
  }
  
  # create cover name column
  # order by site and order columns 1-7
  df <- mydf_reshape[order((as.numeric(mydf_reshape$id))),order(names(mydf_reshape))]
  
  #setup parallel backend to use many processors
  cl <- makeCluster(detectCores()) #set number of cores
  
  # combine 1 & 2 to get total vegetation cover
  df$total_veg <- 0
  df$total_veg <- parApply(cl,df[,2:3],1,sum)
  
  #stop cluster
  stopCluster(cl)
  
  # create vector that indicates if open water is present in the patch
  water <- df$percent.4
  water[water > 0] <- 1
  
  # Save info
  patch_covs <- data.frame(patch=df[,1], patch_tree=df[,2], patch_total_veg=df[,"total_veg"], water_present=water,
                           patch_size=as.numeric(as.character(global_patches@data[,"AREA"]))
  )