### load required packages ###

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
           'fasterize', 'compiler')

package_load(packs)

### change dataframes back to an array ###

df_2_array <- function(my_df = NULL){
  
  require(reshape2)
  my_array <- acast(my_df, species~site~season, value.var = "count")
  dimnames(my_array) <- NULL
  return(my_array)
}

#### function to intersect background data layer (vector) with buffer layers (vector) and calculate needed data ###
# layer 1 should be the input layer with data to extract, layer 2 is the intersecting layer
getIntersect <- function(layer1, layer2, columnName, getMean=FALSE, 
                         getLength=FALSE, getArea=FALSE, n.cores){
  
  # set and register cores
  cl <- makeCluster(n.cores) # set number of cores
  registerDoParallel(cl)     # register backend
  
  # use a foreach loop to parallel and loop through each feature of 
  #  the buffer layer
  dat <- foreach(i=1:length(layer2), 
                 .combine=c, 
                 .packages=c("raster", "dplyr")) %dopar% {
    # intersect two vector layers - layer1 is the input layer
    buffer_intersect <- raster::intersect(layer1,layer2[i,])
    # loop through each feature and extract the needed data
    if(is.null(buffer_intersect)){
      buffer_intersect <- 0
    } else {
      if(getLength){
        buffer_intersect %>% 
          SpatialLinesLengths(.) %>% 
          sum(.) %>% 
          return(.)
      } else {
        if(getArea) {
          data_intersect$area <- raster::area(data_intersect)
          # calculate the areas of each landuse category
          data_table <- data_intersect@data %>% 
                        group_by_at(vars(one_of(columnName))) %>%
                        summarise(total=sum(as.numeric(as.character(area)))) %>% 
                        as.data.frame(.)
          # extract the land use category that has the greatest area
          data_table[which(data_table$total == 
                           max(data_table$total, na.rm=TRUE)),1] %>% 
                     as.character(.) %>% 
                     return(.)
        } else {
          if(getMean){
            buffer_intersect@data[,columnName] %>% 
              as.character(.) %>% 
              as.numeric(.) %>% 
              mean(.) %>% 
              return(.)
          } else {
            buffer_intersect@data[,columnName] %>% 
              as.character(.) %>% 
              as.numeric(.) %>% 
              sum(.) %>% 
              return(.)
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

### extract the proportion of land cover categories for each polygon ###
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


# function that creates a transition layer, corrects diagnols, and creates a cost distance matrix
getEdgeWeights <- function(cost_raster, transitionFunction, directions, fromCoords, toCoords, dist.cutoff)
{
  ### modified transition function from gdistance package to format the cost raster
  # calculate conductances among neighbors
  # create empty transition matrix
  tr <- new("TransitionLayer",
            nrows=as.integer(nrow(cost_raster)),
            ncols=as.integer(ncol(cost_raster)),
            extent=extent(cost_raster),
            crs=projection(cost_raster, asText=FALSE),
            transitionMatrix = Matrix(0,ncell(cost_raster),ncell(cost_raster)),
            transitionCells = 1:ncell(cost_raster))
  # weird way to make a function, but from the gdistance package
  transitionMatr <- transitionMatrix(tr)
  # find cells that actually have values in them
  Cells <- which(!is.na(getValues(cost_raster)))
  # identify adjacent cells and create adjacency matrix
  adj <- adjacent(cost_raster, cells=Cells, pairs=TRUE, target=Cells, directions=directions)
  # make it a symmetrical matrix
  adj <- adj[adj[,1] < adj[,2],]
  # get cell values from cost raster
  dataVals <- cbind(getValues(cost_raster)[adj[,1]],getValues(cost_raster)[adj[,2]])
  # calcualte the transition (conductance) values using 1/mean(x)
  transition.values <- 1/rowMeans(dataVals)
  # there cant be negative values
  if(!all(transition.values>=0)){warning("transition function gives negative values")}
  # turn the adjacency values into a transition matrix
  transitionMatr[adj] <- as.vector(transition.values)
  transitionMatr <- forceSymmetric(transitionMatr)
  transitionMatrix(tr) <- transitionMatr
  matrixValues(tr) <- "conductance"
  ### geocorrect the transition matrix
  # adjust diag. conductances
  trCorr <- geoCorrection(tr, type="c", multpl = FALSE, scl = FALSE)
  
  ### modified costdistance function from gdistance package to format the to and from coordinates
  ## create a list of sites only within dispersal distance (Euclidian) from each site
  toCoords.loc <- vector("list", nrow(toCoords)) # location of each retained coordinate in the distance matrix
  toCoords.list <- vector("list", nrow(toCoords)) # list of coordinates within dispersal distance of each individual point (list in order of points)
  for(j in 1:nrow(toCoords)){
    # calcualte Euclidian distance from each point to all the rest
    toCoords.loc[[j]] <- which(sqrt((fromCoords[j,1]-fromCoords[,1])^2 + 
                                      (fromCoords[j,2]-fromCoords[,2])^2) < dist.cutoff)
    # save the location of the used coordinates from the entire list of coordinates
    toCoords.list[[j]] <- toCoords[toCoords.loc[[j]],]
  }
  
  # align the coordinates with the cells from the transition matrix
  fromCells <- cellFromXY(trCorr, fromCoords)
  
  # remove duplicates with unique function
  # have to set up toCells a bit different
  toCells <- toCells_gdist <- vector("list", nrow(fromCoords)) # list to hold the toCells for each individual fromCell
  # indicate the toCells, but keep them within the list to hold true to the dispersal cutoff
  for(i in 1:length(toCoords.list)){
    toCells[[i]] <- cellFromXY(trCorr, toCoords.list[[i]])
    toCells_gdist[[i]] <- unique(toCells[[i]])
  }
  
  # make a list object that can be run through parLapply
  cell_list <- vector("list", length(fromCells))
  for(i in 1:length(cell_list)){
    cell_list[[i]] <- list(fromCells[i],toCells_gdist[[i]])
  }
  
  # make the transition matrix an object that igraph can read
  y <- gdistance::transitionMatrix(trCorr)
  # create an adjacencyGraph to get edge weights
  if(isSymmetric(y)) {m <- "undirected"} else{m <- "directed"}
  adjacencyGraph <- graph.adjacency(y, mode=m, weighted=TRUE)
  # reclassify edge weights as cost
  E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
  
  return(list(adjacencyGraph = adjacencyGraph,
              cell_list = cell_list,
              fromCoords = fromCoords,
              toCoords.loc = toCoords.loc,
              toCells = toCells,
              toCells_gdist = toCells_gdist))
}
         
# function that creates a transition layer, corrects diagnols, and creates a cost distance matrix
getAllEdgeWeights <- function(cost_raster, directions, coords){
  ### modified transition function from gdistance package to format the cost raster
  # calculate conductances among neighbors
  # create empty transition matrix
  tr <- new("TransitionLayer",
            nrows=as.integer(nrow(cost_raster)),
            ncols=as.integer(ncol(cost_raster)),
            extent=extent(cost_raster),
            crs=projection(cost_raster, asText=FALSE),
            transitionMatrix = Matrix(0,ncell(cost_raster),ncell(cost_raster)),
            transitionCells = 1:ncell(cost_raster))
  # weird way to make a function, but from the gdistance package
  transitionMatr <- transitionMatrix(tr)
  # find cells that actually have values in them
  Cells <- which(!is.na(getValues(cost_raster)))
  # identify adjacent cells and create adjacency matrix
  adj <- adjacent(cost_raster, cells=Cells, pairs=TRUE, target=Cells, directions=directions)
  # make it a symmetrical matrix
  adj <- adj[adj[,1] < adj[,2],]
  # get cell values from cost raster
  dataVals <- cbind(getValues(cost_raster)[adj[,1]],getValues(cost_raster)[adj[,2]])
  # calcualte the transition (conductance) values using 1/mean(x)
  transition.values <- 1/rowMeans(dataVals)
  # there cant be negative values
  if(!all(transition.values>=0)){warning("transition function gives negative values")}
  # turn the adjacency values into a transition matrix
  transitionMatr[adj] <- as.vector(transition.values)
  transitionMatr <- forceSymmetric(transitionMatr)
  transitionMatrix(tr) <- transitionMatr
  matrixValues(tr) <- "conductance"
  ### geocorrect the transition matrix
  # adjust diag. conductances
  trCorr <- geoCorrection(tr, type="c", multpl = FALSE, scl = FALSE)
  # make the transition matrix an object that igraph can read
  y <- gdistance::transitionMatrix(trCorr)
  # create an adjacencyGraph to get edge weights
  if(isSymmetric(y)) {m <- "undirected"} else{m <- "directed"}
  adjacencyGraph <- graph.adjacency(y, mode=m, weighted=TRUE)
  # reclassify edge weights as cost
  E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
  
  # align the coordinates with the cells from the transition matrix
  fromCells <- cellFromXY(trCorr, coords)
  # there can not be duplications in the toCells
  toCells <- unique(fromCells)
  
  # make a list object that can be run through parLapply
  cell_list <- vector("list", length(fromCells))
  for(i in 1:length(cell_list)){
    cell_list[[i]] <- list(fromCells[i],toCells)
  }
  
  return(list(adjacencyGraph = adjacencyGraph,
              cell_list = cell_list,
              fromCells = fromCells,
              toCells = toCells))
}

  
