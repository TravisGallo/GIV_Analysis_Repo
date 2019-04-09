library(raster)
library(sf)
library(parallel)

# calculate the proportion of each landcover variable within each patch
data_list <- readRDS("2019-03-25_raster_patches_ARGO.rds")

patches <- data_list[[1]]
rm <- data_list[[2]]

patch_list <- split(patches, seq(nrow(patches)))

cl <- makeCluster(detecCores())
clusterEvalQ(cl,library("sf"))
clusterEvalQ(cl, library("raster"))
clusterExport(cl, c("rm", "patch_list"))


raster_extract <- parLapply(cl, patch_list, function(x) {
  prop_vec <- prop.table(table(raster::extract(rm, x)))
  })

writeRDS(raster_extract, paste0(Sys.Date(),"_raster_extract.rds"))

